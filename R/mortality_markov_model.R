#' build a Markov probability model from multi-age prediction models
#' To implement this, all patients at a given age will be binned according to their model score (using quantiles).
#' Each bin is assigned a state, and we are computing the probability for traversing from each state to the next model state
#' Patients with missing score are also included for this model to reflect actual population numbers
#' @param models - list of prediction models (output of mldpEHR.cv_train_stitch_outcome)
#' @param outcome - time from oldest model (first) to target outcome
#' @param step - time between prediction models
#' @param qbins - quantile bin size of prediction score for which the markov model will define a state
#' @return a list of markov models (per age), each is a list with the following members:
#' - model - matrix containing the probability for each quantile(score) bin to reach each of the target_classes provided in the oldest model.
#' - local.model - data.frame containing the probability for each quantile(score) bin to reach each of the quantile(score) bins of the next model by age.
#' - qbins -  bins
#' - target - data frame containing the target bin for this age model (to be used as outcome for the younger age model)

#' @examples
#'
#' library(dplyr)
#' library(ggplot2)
#' # build base predictor
#' outcome <- data.frame(id = 1:1000, sex = rep(0:1, 500), age=rep(80, 1000), death = rep(c(NA, 82), each = 500), followup=rep(5, 1000))
#' patients <- c(list(outcome), purrr::map(1:4, ~ data.frame(
#'      id = 1:1000, 
#'      sex = rep(0:1, 500),
#'      age = 80 - .x * 5,
#'      death = rep(c(NA, 82), each = 500),
#'      followup = .x*5+5))) %>% setNames(seq(80, by=-5, length.out=5))
#' features <- purrr::map(1:5, ~ data.frame(id = 1:1000, a = c(rnorm(500), rnorm(500, mean=2, sd=1)), b = c(rnorm(500), rnorm(500, mean=02, sd=1)))) %>%
#'     setNames(seq(80, by=-5, length.out=5))
#' predictors <- mldpEHR.mortality_multi_age_predictors(patients, features, 3, q_thresh = 0.5)
#' markov <- mldpEHR.mortality_markov(predictors, 5, 5, qbins = seq(0, 1, by = 0.25))
#' prob <- purrr::map2_df(markov, names(markov), ~
#'     as_tibble(.x$model[[1]], rownames = "sbin") %>%
#'         mutate(sex = 0, model = .y) %>%
#'         bind_rows(as_tibble(.x$model[[2]], rownames = "sbin") %>% mutate(sex = 1, model = .y))) %>%
#'     mutate(sbin = factor(sbin, levels = c(1:4, "death")))
#'
#' ggplot(prob, aes(x = sbin, y = death, colour = factor(sex), group = factor(sex))) +
#'     geom_point() +
#'     geom_line() +
#'     facet_wrap(~model, nrow = 1) +
#'     theme_bw()
#'
#' @export


mldpEHR.mortality_markov <- function(models, outcome, step, qbins = seq(0, 1, by = 0.05), required_conditions="id==id") {
    markov_models <- list()

    # first model is the oldest model, used to compute the actual risk
    markov_models[[1]] <- .mortality_markov_model_for_outcome_model(models[[1]], outcome, qbins, required_conditions)
    for (i in 2:length(models)) {
        markov_models[[i]] <- .mortality_markov_model_for_stitch_model(markov_models[[i - 1]], models[[i]], step, qbins, required_conditions)
    }
    names(markov_models) <- names(models)
    return(markov_models)
}


#' build an Markov probability model from multi-age prediction models
#' @param markov - the markov model computed for the next age (older). the states in this age will be mapped to the states in this markov layer
#' @param model - prediction model (output of build_cross_validation_time_stitch_classification_models)
#' @param follow -  data.frames  defining for each patient how much follow-up time they had. Dataframe must include the following columns:
#' - id
#' - time_in_system
#' - target_class - the classification if known within step (i.e. patient died). NA for unknown step outcome.
#' @param step - time between prediction models
#' @param qbins - quantile bin size of prediction score for which the markov model will define a state
#' @return a markov model, a list with the following members:
#' - model - matrix containing the probability for each quantile(score) bin to reach each of the target_classes provided in the oldest model.
#' - local.model - data.frame containing the probability for each quantile(score) bin to reach each of the quantile(score) bins of the next model by age.
#' - qbins -  bins
#' - target - data frame containing the target bin for this age model (to be used as outcome for the younger age model)

.mortality_markov_model_for_stitch_model <- function(markov, model, step, qbins, required_conditions, min_obs_for_estimate = 10) {
    m <- model$test %>% # contains only patients with score (some were dropped)
        mutate(sbin = as.numeric(cut(qpredict, qbins, right = FALSE, include.lowest = TRUE))) %>%
        mutate(sbin = replace_na(sbin, "no_score")) %>%
        mutate(sbin = factor(sbin, levels = c(1:(length(qbins) - 1), "no_score")))
    pm <- m %>% 
        left_join(markov$target %>% select(id, target_sbin), by = "id") %>% 
        filter(!is.na(target_sbin) | step_outcome !="alive", eval(rlang::parse_expr(required_conditions))) %>% 
        mutate(outcome=ifelse(is.na(target_sbin), step_outcome, as.character(target_sbin)))

    km_model <- plyr::ddply(pm, plyr::.(sbin, sex), function(data) {
        if (nrow(data) < min_obs_for_estimate) {
            # bin is too small to compute probabilities, will use entire pop prob to fill in
            data <- pm %>% filter(sex == data$sex[1])
            warning("Insufficient stats for age:", data$age[1], " for sbin:", data$sbin[1], " for sex:", c("male", "female")[data$sex[1]], " ,using pop")
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

    local_model <- km_model %>%
        mutate(outcome = factor(outcome, levels = rownames(markov$local_model[[1]]))) %>%
        arrange(outcome, sbin) %>%
        pivot_wider(id_cols = c(sbin, sex), names_from = outcome, values_from = est) %>%
        arrange(sex, sbin) %>%
        replace(is.na(.), 0)

    local_model_by_sex <- plyr::dlply(local_model, plyr::.(sex), function(s) {
        s %>%
            select(-sex) %>%
            bind_rows(data.frame(sbin = c("death"), death = c(1))) %>%
            mutate(sbin = factor(sbin, levels = c(1:(length(qbins) - 1), "no_score", "death"))) %>%
            replace(is.na(.), 0) %>%
            arrange(sbin) %>%
            column_to_rownames("sbin")
    })

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

.mortality_markov_model_for_outcome_model <- function(model, step, qbins, required_conditions, min_obs_for_estimate = 10) {
    m <- model$test %>% # contains only patients with score (some were dropped)
        mutate(sbin = as.numeric(cut(qpredict, qbins, right = FALSE, include.lowest = TRUE))) %>%
        mutate(sbin = replace_na(sbin, "no_score")) %>%
        mutate(sbin = factor(sbin, levels = c(1:(length(qbins) - 1), "no_score")))
    pm <- m %>% filter(!is.na(step_outcome), eval(rlang::parse_expr(required_conditions))) %>% 
        rename(outcome=step_outcome) 

    # compute competing risk models for each sex indeptendently according to source bin(sbin)
    km_model <- plyr::ddply(pm, plyr::.(sbin, sex), function(data) {
        #message(data$sbin[1], " :: " , data$sex[1])
        if (nrow(data) < min_obs_for_estimate) {
            # bin is too small to compute probabilities, will use entire pop prob to fill in
            data <- pm %>% filter(sex == data$sex[1])
            warning("Insufficient stats for age:", data$age[1], " for sbin:", data$sbin[1], " for sex:", c("male", "female")[data$sex[1]], " ,using pop")
        }
        if (nrow(data %>% filter(outcome == "death")) == 0) {
            return(data.frame(time=step, est=0, var=0, group=1, outcome="death"))
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
    # adding estimates for censored data (patient is alive)
    censored_estimate <- km_model %>%
        group_by(sbin, sex) %>%
        summarize(est = 1 - sum(est), .groups = "drop") %>%
        mutate(outcome = "alive")
    # combining all outcome states
    local_model <- km_model %>%
        bind_rows(censored_estimate) %>% # change target to factor?
        pivot_wider(id_cols = c(sbin, sex), names_from = outcome, values_from = est) %>%
        arrange(sex, sbin)

    local_model_by_sex <- plyr::dlply(local_model, plyr::.(sex), function(s) {
        s %>%
            select(-sex) %>%
            bind_rows(data.frame(sbin = c("death"), death = c(1))) %>%
            mutate(sbin = factor(sbin, levels = c(1:(length(qbins) - 1), "no_score", "death"))) %>%
            replace(is.na(.), 0) %>%
            arrange(sbin) %>%
            column_to_rownames("sbin")
    })
    model_by_sex <- purrr::map(local_model_by_sex, ~ as.matrix(.x))

    return(list(
        model = model_by_sex,
        local_model = local_model_by_sex,
        qbins = qbins,
        target = m %>% mutate(target_sbin = sbin)
    ))
}
