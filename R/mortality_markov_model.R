#' Build a Markov probability model from multi-age prediction models of mortality
#'
#' @description All patients at a given age will be binned according to their model score (using quantiles).
#' Each bin is assigned a state, and we are computing the probability for traversing from
#' each state to the next model state. \cr
#' Patients with missing score are also included for this model to reflect actual population numbers.
#'
#' @param models a list of prediction models (output of \code{mldp_mortality_multi_age_predictors}).
#' @param survival_years The number of survival years the model aimed to predict for the oldest age group (see \code{mldp_mortality_multi_age_predictors}).
#' @param qbins a vector of quantile bins for the prediction score from which the markov model will define a state. Default is to use 10 equally sized quantiles.
#' @param required_conditions a string with an expression for any filter to apply to the patients to filter out from model computation. Can be used to filter out patients with missing data or limiting the time window. See \code{\link[dplyr]{filter}} for more details.
#'
#' @return a list of with the following members:
#' \itemize{
#' \item{prob: }{a data frame containing the probability to survive ("alive") or die ("death") before the oldest age group + survival_years for each age, gender and score bin (quantile of score).}
#' \item{models: }{a list of matrices containing the probability for each quantile(score) bin to reach
#' each of the quantile(score) bins of the next model by age.}
#' }
#'
#' @examples
#'
#' library(ggplot2)
#' library(dplyr)
#'
#' mortality <- load_mortality_example_data(N = 100, num_age_groups = 3)
#' predictors <- mldp_mortality_multi_age_predictors(
#'     mortality@patients,
#'     mortality@features,
#'     survival_years = 5,
#'     nfolds = 2,
#'     q_thresh = 0.2,
#'     nthread = 2 # CRAN allows only 2 cores
#' )
#'
#' markov <- mldp_mortality_markov(predictors, 5, qbins = seq(0, 1, by = 0.1))
#'
#' head(markov$prob)
#'
#' # plot the probability of survival for each age group
#' markov$prob %>%
#'     mutate(sex = factor(c("male", "female")[sex])) %>%
#'     ggplot(aes(x = sbin, y = alive, colour = sex, group = sex)) +
#'     geom_point() +
#'     geom_line() +
#'     ylab("p(survival)") +
#'     xlab("score bin") +
#'     scale_color_manual(name = "", values = c("male" = "blue", "female" = "red")) +
#'     facet_wrap(~age, nrow = 1) +
#'     theme_bw() +
#'     theme(axis.text.x = element_text(angle = 90, hjust = 1))
#'
#' @export
mldp_mortality_markov <- function(models, survival_years, qbins = seq(0, 1, by = 0.05), required_conditions = NULL) {
    markov_models <- list()
    steps <- abs(diff(purrr::map_dbl(models, "age")))

    # first model is the oldest model, used to compute the actual risk
    markov_models[[1]] <- mortality_markov_model_for_outcome_model(models[[1]], survival_years, qbins, required_conditions)

    i <- 2
    while (i <= length(models)) {
        markov_models[[i]] <- mortality_markov_model_for_stitch_model(markov_models[[i - 1]], models[[i]], steps[i - 1], qbins, required_conditions)
        i <- i + 1
    }
    names(markov_models) <- names(models)

    probs <- purrr::imap_dfr(markov_models, ~ {
        bind_rows(
            as.data.frame(.x$model[[1]]) %>%
                tibble::rownames_to_column("sbin") %>%
                mutate(sex = 1),
            as.data.frame(.x$model[[2]]) %>%
                tibble::rownames_to_column("sbin") %>%
                mutate(sex = 2)
        ) %>%
            mutate(age = .y)
    }) %>%
        select(age, sbin, sex, everything()) %>%
        mutate(sbin = factor(sbin, levels = c(1:(length(qbins - 1)), "death", "no_score"))) %>%
        mutate(quantile = levels(cut(qbins, qbins, right = FALSE, include.lowest = TRUE))[sbin])

    models <- map(markov_models, ~ purrr::map(.x$local_model, as.matrix))

    return(list(prob = probs, models = models))
}


#' Build an Markov probability model from multi-age prediction models
#'
#' @param markov - the markov model computed for the next age (older). the states in this age will be
#' mapped to the states in this markov layer
#' @param model - prediction model (output of build_cross_validation_time_stitch_classification_models)
#' @param step - time between prediction models
#' @param qbins - quantile bin size of prediction score for which the markov model will define a state
#' @param required_conditions - any filter to apply to the patients to filter out from model computation,
#' for example limiting the time window
#' @param min_obs_for_estimate - minimum of observations required to compute probability per sex/sbin.
#' If minimum is not available, probability will be compuated using all data.
#' @return a markov model, a list with the following members:
#' - model - matrix containing the probability for each quantile(score) bin to reach each of
#' the target_classes provided in the oldest model.
#' - local.model - data.frame containing the probability for each quantile(score) bin to reach
#' each of the quantile(score) bins of the next model by age.
#' - qbins -  bins
#' - target - data frame containing the target bin for this age model (to be used as outcome for
#' the younger age model)
#' @noRd
mortality_markov_model_for_stitch_model <- function(markov, model, step, qbins, required_conditions, min_obs_for_estimate = 10) {
    m <- mortality_set_sbin(model$test, qbins)
    pm <- m %>%
        left_join(markov$target %>% select(id, target_sbin), by = "id")

    if (!is.null(required_conditions)) {
        pm <- pm %>%
            filter(eval(rlang::parse_expr(required_conditions)))
    }

    pm <- pm %>%
        filter(!is.na(target_sbin) | step_outcome != "alive") %>%
        mutate(outcome = ifelse(is.na(target_sbin), step_outcome, as.character(target_sbin)))

    km_model <- mortality_km_sex_sbin(pm, step, min_obs_for_estimate)

    local_model <- km_model %>%
        mutate(outcome = factor(outcome, levels = rownames(markov$local_model[[1]]))) %>%
        arrange(outcome, sbin) %>%
        tidyr::pivot_wider(id_cols = c(sbin, sex), names_from = outcome, values_from = est) %>%
        arrange(sex, sbin) %>%
        replace(is.na(.), 0)

    local_model_by_sex <- mortality_local_model_by_sex(local_model, qbins)
    # adding missing possible outcomes
    local_model_by_sex <- purrr::map2(local_model_by_sex, markov$model, function(a, b) {
        missing_outcome <- setdiff(rownames(b), colnames(a))
        a[, missing_outcome] <- 0
        return(a[, rownames(b)])
    })

    model_by_sex <- purrr::map2(local_model_by_sex, markov$model, ~ as.matrix(.x) %*% .y)

    return(list(
        model = model_by_sex,
        local_model = local_model_by_sex,
        qbins = qbins,
        target = m %>% mutate(target_sbin = sbin)
    ))
}

mortality_markov_model_for_outcome_model <- function(model, step, qbins, required_conditions, min_obs_for_estimate = 10) {
    m <- mortality_set_sbin(model$test, qbins)

    if (!is.null(required_conditions)) {
        pm <- m %>%
            filter(eval(rlang::parse_expr(required_conditions)))
    }

    pm <- m %>%
        filter(!is.na(step_outcome)) %>%
        rename(outcome = step_outcome)

    # compute competing risk models for each sex indeptendently according to source bin(sbin)
    km_model <- mortality_km_sex_sbin(pm, step, min_obs_for_estimate)
    # adding estimates for censored data (patient is alive)
    censored_estimate <- km_model %>%
        group_by(sbin, sex) %>%
        summarize(est = 1 - sum(est), .groups = "drop") %>%
        mutate(outcome = "alive")
    # combining all outcome states
    local_model <- km_model %>%
        bind_rows(censored_estimate) %>% # change target to factor?
        tidyr::pivot_wider(id_cols = c(sbin, sex), names_from = outcome, values_from = est) %>%
        arrange(sex, sbin)

    local_model_by_sex <- mortality_local_model_by_sex(local_model, qbins)
    model_by_sex <- purrr::map(local_model_by_sex, ~ as.matrix(.x))

    return(list(
        model = model_by_sex,
        local_model = local_model_by_sex,
        qbins = qbins,
        target = m %>% mutate(target_sbin = sbin)
    ))
}

mortality_set_sbin <- function(test_score, qbins) {
    return(test_score %>% # contains only patients with score (some were dropped)
        mutate(sbin = as.numeric(cut(qpredict, qbins, right = FALSE, include.lowest = TRUE))) %>%
        mutate(sbin = ifelse(is.na(sbin), "no_score", sbin)) %>%
        mutate(sbin = factor(sbin, levels = c(1:(length(qbins) - 1), "no_score"))))
}

mortality_km_sex_sbin <- function(sbin_outcome_obsT, step, min_obs_for_estimate) {
    km_model <- plyr::adply(expand.grid(age = unique(sbin_outcome_obsT$age), sbin = unique(sbin_outcome_obsT$sbin), sex = c(1, 2)), 1, function(x) {
        data <- sbin_outcome_obsT %>%
            filter(sbin == x$sbin, sex == x$sex)

        if (nrow(data) < min_obs_for_estimate) {
            # bin is too small to compute probabilities, will use entire pop prob to fill in
            cli::cli_alert("Insufficient stats (n={.val {nrow(data)}}) for age: {.val {x$age}}, sbin: {.val {x$sbin}}, sex: {.val {x$sex}}, using entire population")
            data <- sbin_outcome_obsT %>% filter(sex == data$sex[1])
        }
        if (nrow(data %>% filter(outcome != "alive")) == 0) {
            return(data.frame(time = step, est = 0, var = 0, group = 1, outcome = "death"))
        }
        fit <- cmprsk::cuminc(data$obsT, data$outcome, cencode = "alive")
        ret <- purrr::map2_df(fit, names(fit), ~ as.data.frame(.x[1:3]) %>% mutate(name = .y)) %>%
            tidyr::separate(name, into = c("group", "outcome"), sep = " (?=[^ ]*$)") %>%
            filter(time <= step) %>%
            arrange(desc(time), desc(est)) %>% # this is the estimate at latest measured time before outcome years.
            distinct(outcome, .keep_all = TRUE)
        return(ret)
    }) %>%
        select(sbin, sex, est, outcome) %>%
        arrange(sex, sbin, outcome)
    return(km_model)
}

mortality_local_model_by_sex <- function(local_model, qbins) {
    local_model_by_sex <- plyr::dlply(local_model, plyr::.(sex), function(s) {
        s %>%
            select(-sex) %>%
            bind_rows(data.frame(sbin = c("death"), death = c(1))) %>%
            mutate(sbin = factor(sbin, levels = c(1:(length(qbins) - 1), "no_score", "death"))) %>%
            replace(is.na(.), 0) %>%
            arrange(sbin) %>%
            tibble::column_to_rownames("sbin")
    })
    return(local_model_by_sex)
}
