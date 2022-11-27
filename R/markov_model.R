#' build a Markov probability model from multi-age prediction models
#' To implement this, all patients at a given age will be binned according to their model score (using quantiles).
#' Each bin is assigned a state, and we are computing the probability for traversing from each state to the next model state
#' Patients with missing score are also included for this model to reflect actual population numbers
#' @param models - list of prediction models (output of build_cross_validation_time_stitch_classification_models)
#' @param follow_time - list of data.frames per age, defining for each patient how much follow-up time they had. Dataframe must include the following columns:
#' - id
#' - time_in_system
#' - target_class - the classification if known for outcome years or within step (i.e. patient died). NA for unknown step outcome.
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
#' target <- data.frame(id = 1:1000, target_class = rep(c(0, 1), each = 500), sex = rep(0:1, 500))
#' features <- data.frame(id = 1:500, a = rnorm(500), b = rnorm(500)) %>%
#'     bind_rows(
#'         data.frame(id = 501:1000, a = rnorm(500, mean = 2, sd = 2), b = rnorm(500, mean = -2, sd = 2))
#'     )
#' predictor <- mldpEHR.cv_train_outcome(target, features, folds = 3)
#' target_list <- purrr::map(1:3, ~ data.frame(id = 1:1000, target_class = NA, sex = rep(0:1, 500))) %>% setNames(1:3)
#' feature_list <- purrr::map(1:3, ~ data.frame(id = 1:500, a = rnorm(500), b = rnorm(500)) %>%
#'     bind_rows(
#'         data.frame(id = 501:1000, a = rnorm(500, mean = 2, sd = 1), b = rnorm(500, mean = -2, sd = 1))
#'     )) %>% setNames(1:3)
#' models <- mldpEHR.cv_train_stitch_outcome(0, predictor, target_list, feature_list, q_thresh = 0.5)
#' follow_time <- c(list(target %>% mutate(time_in_system = 5)), purrr::map(target_list, ~ .x %>%
#'     select(id, sex, target_class) %>%
#'     mutate(time_in_system = 5)))
#' markov <- mldpEHR.mortality_markov(models, follow_time, 5, 5, qbins = seq(0, 1, by = 0.25))
#' prob <- purrr::map2_df(markov, names(markov), ~
#'     as_tibble(.x$model[[1]], rownames = "sbin") %>%
#'         mutate(sex = 0, model = .y) %>%
#'         bind_rows(as_tibble(.x$model[[2]], rownames = "sbin") %>% mutate(sex = 1, model = .y))) %>%
#'     mutate(sbin = factor(sbin, levels = c(1:4, "dead")))
#'
#' ggplot(prob, aes(x = sbin, y = dead, colour = factor(sex), group = factor(sex))) +
#'     geom_point() +
#'     geom_line() +
#'     facet_wrap(~model, nrow = 1) +
#'     theme_bw()
#'
#' @export

mldpEHR.mortality_markov <- function(models, follow_time, outcome, step, qbins = seq(0, 1, by = 0.05)) {
    markov_models <- list()

    # first model is the oldest model, used to compute the actual risk
    markov_models[[1]] <- .mortality_markov_model_for_outcome_model(models[[1]], follow_time[[1]], outcome, qbins)
    for (i in 2:length(models)) {
        markov_models[[i]] <- .mortality_markov_model_for_stitch_model(markov_models[[i - 1]], models[[i]], follow_time[[i]], step, qbins)
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

.mortality_markov_model_for_stitch_model <- function(markov, model, follow, step, qbins, min_obs_for_estimate = 10) {
    m <- model$test %>% # contains only patients with score (some were dropped)
        mutate(sbin = as.numeric(cut(qpredict, qbins, right = FALSE, include.lowest = TRUE))) %>%
        mutate(sbin = factor(sbin, levels = c(1:(length(qbins) - 1), "no_score")))
    pm <- follow %>% # contains all patients
        #        filter(!is.na(survival) | time_in_system > step) %>% #assumes that people that leave the system are not more/less prone to die
        left_join(m %>% select(-target_class), by = c("id", "sex")) %>%
        mutate(sbin = replace_na(sbin, "no_score")) %>%
        left_join(markov$target %>% select(id, target_sbin), by = "id")

    # filtering out patients with missing target. The assumption is that patients who are censored (due to database ending)
    # are not health state dependent and therefore will not affect the estimations
    pm <- pm %>% filter(!is.na(target_sbin) | !is.na(target_class) | time_in_system >= step)

    # setting the observed time to event for patients that were classified in the next step
    pm$time_in_system[!is.na(pm$target_sbin)] <- step

    # taking care of the patients without target in next age
    # setting the target for patients that die
    levels(pm$target_sbin) <- c(levels(pm$target_sbin), "dead")
    pm$target_sbin[!is.na(pm$target_class) & pm$target_class == 1] <- "dead"

    # setting the target for patients that are alive but without score in next age
    pm$target_sbin[is.na(pm$target_class) & is.na(pm$target_sbin)] <- "no_score"

    km_model <- plyr::ddply(pm, plyr::.(sbin, sex), function(data) {
        if (nrow(data) < min_obs_for_estimate) {
            # bin is too small to compute probabilities, will use entire pop prob to fill in
            data <- pm %>% filter(sex == data$sex[1])
            warning("Insufficient stats for age:", data$age[1], " for sbin:", data$sbin[1], " for sex:", c("male", "female")[data$sex[1]], " ,using pop")
        }
        fit <- cmprsk::cuminc(data$time_in_system, data$target_sbin, cencode = "alive")
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
            bind_rows(data.frame(sbin = c("dead"), dead = c(1))) %>%
            mutate(sbin = factor(sbin, levels = c(1:(length(qbins) - 1), "no_score", "dead"))) %>%
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

.mortality_markov_model_for_outcome_model <- function(model, follow, outcome, qbins, min_obs_for_estimate = 10) {
    # in cases of mortality target_class is either 0 - patient is alive, or 1 - patient died
    target <- model$target %>%
        filter(!is.na(target_class)) %>%
        distinct(target_class) %>%
        mutate(outcome = factor(c("alive", "dead")[target_class + 1], levels = c("alive", "dead")))
    # note for mortality, target_class == 1 --> "death"

    m <- model$test %>% # contains only patients with score (some were dropped)
        mutate(sbin = as.numeric(cut(qpredict, qbins, right = FALSE, include.lowest = TRUE))) %>%
        mutate(sbin = factor(sbin, levels = c(1:(length(qbins) - 1), "no_score")))
    pm <- follow %>% # contains all patients
        left_join(m %>% select(-target_class), by = c("id", "sex")) %>%
        mutate(sbin = replace_na(sbin, "no_score")) %>% # adding missing data patients
        left_join(target, by = "target_class") %>%
        filter(!is.na(outcome))

    # compute competing risk models for each sex indeptendently according to source bin(sbin)
    km_model <- plyr::ddply(pm, plyr::.(sbin, sex), function(data) {
        if (nrow(data) < min_obs_for_estimate) {
            # bin is too small to compute probabilities, will use entire pop prob to fill in
            data <- pm %>% filter(sex == data$sex[1])
            warning("Insufficient stats for age:", data$age[1], " for sbin:", data$sbin[1], " for sex:", c("male", "female")[data$sex[1]], " ,using pop")
        }
        fit <- cmprsk::cuminc(data$time_in_system, data$outcome, cencode = "alive")
        ret <- purrr::map2_df(fit, names(fit), ~ as.data.frame(.x[1:3]) %>% mutate(name = .y)) %>%
            tidyr::separate(name, into = c("group", "outcome"), sep = " (?=[^ ]*$)") %>%
            filter(time <= outcome) %>%
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
            bind_rows(data.frame(sbin = c("dead"), dead = c(1))) %>%
            mutate(sbin = factor(sbin, levels = c(1:(length(qbins) - 1), "no_score", "dead"))) %>%
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
