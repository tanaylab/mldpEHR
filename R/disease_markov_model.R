
#' Build a Markov probability model from multi-age prediction models of a disease
#'
#' @description All patients at a given age will be binned according to their model score (using quantiles).
#' Each bin is assigned a state, and we are computing the probability for traversing from
#' each state to the next model state. \cr
#' Patients with missing score are also included for this model to reflect actual population numbers.
#'
#' @param models a list of prediction models (output of \code{mldp_disease_multi_age_predictors}).
#' @param target_years the number of years to predict the disease for (e.g. 5 years).
#'
#' @return a list of with the following members:
#' \itemize{
#' \item{prob: }{a data frame containing the life-long probability to get the disease ("disease"), die ("death") die with the disease ("death_disease"), or not get the disease and live ("healthy") for each age, gender and score bin (quantile of score).}#'
#' \item{models: }{a list of matrices containing the probability for each quantile(score) bin to reach
#' each of the quantile(score) bins of the next model by age.}
#' }
#'
#'
#' @examples
#'
#'
#' # Load a small example data
#' disease_data <- load_disease_example_data(N = 100, num_age_groups = 3)
#'
#' # Build predictors
#' predictors <- mldp_disease_multi_age_predictors(
#'     disease_data@patients,
#'     disease_data@features,
#'     target_years = 5,
#'     nfolds = 2
#' )
#'
#' markov <- mldp_disease_markov(predictors, 5, qbins = seq(0, 1, by = 0.1))
#'
#' # plot the life-long probability for the disease for each age and gender
#'
#' library(ggplot2)
#' library(dplyr)
#'
#' markov$prob %>%
#'     mutate(sex = factor(c("male", "female")[sex])) %>%
#'     ggplot(aes(x = sbin, y = disease + disease_death, colour = sex, group = sex)) +
#'     geom_point() +
#'     geom_line() +
#'     ylab("p(disease)") +
#'     xlab("score bin") +
#'     scale_color_manual(name = "", values = c("male" = "blue", "female" = "red")) +
#'     facet_wrap(~age, nrow = 1) +
#'     theme_bw() +
#'     theme(axis.text.x = element_text(angle = 90, hjust = 1))
#'
#' @inheritParams mldp_mortality_markov
#' @export
mldp_disease_markov <- function(models, target_years, qbins = seq(0, 1, by = 0.05), required_conditions = NULL) {
    markov_models <- list()
    steps <- abs(diff(purrr::map_dbl(models, "age")))

    # first model is the oldest model, used to compute the actual risk
    markov_models[[1]] <- disease_markov_model_for_outcome_model(models[[1]], target_years, qbins, required_conditions)
    i <- 2
    while (i <= length(models)) {
        markov_models[[i]] <- disease_markov_model_for_stitch_model(markov_models[[i - 1]], models[[i]], steps[i - 1], qbins, required_conditions)
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
        mutate(sbin = factor(sbin, levels = c(1:(length(qbins - 1)), "disease", "disease_death", "death", "no_score"))) %>%
        mutate(quantile = levels(cut(qbins, qbins, right = FALSE, include.lowest = TRUE))[sbin])

    models <- purrr::map(markov_models, ~ purrr::map(.x$local_model, as.matrix))

    return(list(prob = probs, models = models))
}


#' build an Markov probability model from multi-age prediction models
#' @param markov - the markov model computed for the next age (older). the states in this age will be mapped to the states in this markov layer
#' @param model - prediction model (output of build_cross_validation_time_stitch_classification_models)
#' @param step - time between prediction models
#' @param qbins - quantile bin size of prediction score for which the markov model will define a state
#' @param required_conditions - any filter to apply to the patients to filter out from model computation,
#' for example limiting the time window
#' @param min_obs_for_estimate - minimum of observations required to compute probability per sex/sbin.
#' If minimum is not available, probability will be compuated using all data.
#' @return a markov model, a list with the following members:
#' - model - matrix containing the probability for each quantile(score) bin to reach each of the target_classes provided in the oldest model.
#' - local.model - data.frame containing the probability for each quantile(score) bin to reach each of the quantile(score) bins of the next model by age.
#' - qbins -  bins
#' - target - data frame containing the target bin for this age model (to be used as outcome for the younger age model)
#' @noRd
disease_markov_model_for_stitch_model <- function(markov, model, step, qbins, required_conditions, min_obs_for_estimate = 10) {
    m <- disease_set_sbin(model$test, qbins)

    pm <- m %>%
        left_join(markov$target %>% select(id, target_sbin), by = "id")

    if (!is.null(required_conditions)) {
        pm <- pm %>% filter(eval(rlang::parse_expr(required_conditions)))
    }

    pm <- pm %>%
        filter(!is.na(target_sbin) | step_outcome != "healthy") %>%
        mutate(outcome = ifelse(is.na(target_sbin), step_outcome, as.character(target_sbin)))

    km_model <- disease_km_sex_sbin(pm, step, min_obs_for_estimate)
    # adding estimates for censored data (patient remains sick)
    censored_estimate <- km_model %>%
        filter(sbin == "disease") %>%
        group_by(sbin, sex) %>%
        summarize(est = 1 - sum(est), .groups = "drop") %>%
        mutate(outcome = "disease")
    # combining all outcome states

    local_model <- km_model %>%
        bind_rows(censored_estimate) %>% # change target to factor?
        # correcting for numerical errors
        group_by(sbin, sex) %>%
        mutate(s = sum(est)) %>%
        ungroup() %>%
        mutate(est = est / s) %>%
        select(-s) %>%
        mutate(outcome = factor(outcome, levels = c(levels(markov$target$target_sbin), "disease_death", "death"))) %>%
        arrange(outcome) %>%
        tidyr::pivot_wider(id_cols = c(sbin, sex), names_from = outcome, values_from = est) %>%
        arrange(sex, sbin)

    local_model_by_sex <- disease_local_model_by_sex(local_model, qbins)

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

disease_markov_model_for_outcome_model <- function(model, step, qbins, required_conditions, min_obs_for_estimate = 10) {
    m <- disease_set_sbin(model$test, qbins)

    pm <- m %>% filter(!is.na(step_outcome))

    if (!is.null(required_conditions)) {
        pm <- pm %>% filter(eval(rlang::parse_expr(required_conditions)))
    }

    pm <- pm %>%
        rename(outcome = step_outcome) %>%
        filter(outcome != "healthy" | obsT >= step)

    # compute competing risk models for each sex indeptendently according to source bin(sbin)
    km_model <- disease_km_sex_sbin(pm, step, min_obs_for_estimate)
    # adding estimates for censored data (patient is healthy)
    censored_estimate <- km_model %>%
        group_by(sbin, sex) %>%
        summarize(est = 1 - sum(est), .groups = "drop") %>%
        mutate(outcome = ifelse(sbin == "disease", "disease", "healthy"))
    # combining all outcome states

    local_model <- km_model %>%
        bind_rows(censored_estimate) %>% # change target to factor?
        # correcting for numerical errors
        group_by(sbin, sex) %>%
        mutate(s = sum(est)) %>%
        ungroup() %>%
        mutate(est = est / s) %>%
        select(-s) %>%
        mutate(outcome = factor(outcome, levels = c("disease", "disease_death", "death", "healthy"))) %>%
        arrange(outcome) %>%
        tidyr::pivot_wider(id_cols = c(sbin, sex), names_from = outcome, values_from = est) %>%
        arrange(sex, sbin)

    local_model_by_sex <- disease_local_model_by_sex(local_model, qbins)
    model_by_sex <- purrr::map(local_model_by_sex, ~ as.matrix(.x))

    return(list(
        model = model_by_sex,
        local_model = local_model_by_sex,
        qbins = qbins,
        target = m %>% mutate(target_sbin = sbin)
    ))
}

disease_set_sbin <- function(test_score, qbins) {
    test_score %>% # contains only patients with score (some were dropped)
        mutate(sbin = as.numeric(cut(qpredict, qbins, right = FALSE, include.lowest = TRUE))) %>%
        mutate(sbin = ifelse(!is.na(disease) & disease < age, "disease", sbin)) %>%
        mutate(sbin = ifelse(is.na(sbin), "no_score", sbin)) %>%
        mutate(sbin = factor(sbin, levels = c(1:(length(qbins) - 1), "no_score", "disease")))
}

disease_km_sex_sbin <- function(sbin_outcome_obsT, step, min_obs_for_estimate) {
    km_model <- plyr::adply(levels(sbin_outcome_obsT$sbin), 1, function(cur_sbin) {
        return(plyr::adply(unique(sbin_outcome_obsT$sex), 1, function(cur_sex) {
            # message(cur_sbin, " :: ", cur_sex)
            data <- sbin_outcome_obsT %>% filter(sbin == cur_sbin, sex == cur_sex)
            cencode <- "healthy"
            # message(data$sbin[1], " :: " , data$sex[1])
            if (nrow(data) < min_obs_for_estimate & cur_sbin != "disease") {
                # bin is too small to compute probabilities, will use entire pop prob to fill in
                cli::cli_alert("Insufficient stats (n={.val {nrow(data)}}) for age: {.val {sbin_outcome_obsT$age[1]}}, sbin: {.val {cur_sbin}}, sex: {.val {c('male', 'female')[cur_sex]}}, using entire population")
                data <- sbin_outcome_obsT %>% filter(sex == cur_sex)
            }
            if (cur_sbin == "disease") {
                # setting censor outcome to disease as possible outcomes are disease or death
                cencode <- "disease"
            }
            # checking if no events occur, just censoring (this will cause an error)
            if (sum(data$outcome == cencode) == nrow(data) | nrow(data) == 0) {
                return(data.frame(sbin = cur_sbin, sex = cur_sex, time = step, est = 0, var = 0, group = 1, outcome = "death"))
            }

            fit <- cmprsk::cuminc(data$obsT, data$outcome, cencode = cencode)
            ret <- purrr::map2_df(fit, names(fit), ~ as.data.frame(.x[1:3]) %>% mutate(name = .y)) %>%
                tidyr::separate(name, into = c("group", "outcome"), sep = " (?=[^ ]*$)") %>%
                filter(time <= step) %>%
                arrange(desc(time), desc(est)) %>% # this is the estimate at latest measured time before outcome years.
                distinct(outcome, .keep_all = TRUE)
            return(ret %>% mutate(sbin = cur_sbin, sex = cur_sex))
        }))
    }) %>%
        select(sbin, sex, est, outcome) %>%
        arrange(sex, sbin, outcome)
    return(km_model)
}


disease_local_model_by_sex <- function(local_model, qbins) {
    local_model_by_sex <- plyr::dlply(local_model, plyr::.(sex), function(s) {
        s %>%
            select(-sex) %>%
            bind_rows(data.frame(sbin = c("disease_death", "death"), disease_death = c(1, 0), death = c(0, 1))) %>%
            mutate(sbin = factor(sbin, levels = c(1:(length(qbins) - 1), "no_score", "disease", "disease_death", "death"))) %>%
            replace(is.na(.), 0) %>%
            arrange(sbin) %>%
            tibble::column_to_rownames("sbin")
    })
    return(local_model_by_sex)
}

#' Compute empirical probabilities for disease
#'
#' @description This method computes the probability for future disease according to the availablitiy of data
#' by basically computing a stitched markov model
#' as score availability is not unbiased (patients with less measurements tend to be healthier),
#' the probaility will be computed seperately for patients with / without score
#' we will define the following modes:
#' - 0 has required conditions at young age and known score / outcome at older age
#' - 1 has required conditions at young age and no restrictions at older age
#' - 2 no restrictions at young/older age
#' @param population - list of data.frames of all the patients in the system going back in time. For example the first data.frame represents age 80, next is 75 and so forth.
#' Each patient data.frame contains the following columns:
#' - patient id
#' - sex
#' - target_class (does this patient have a known outcome)
#' - disease age of disease
#' - death age of death
#' - followup
#' - any other columns that can be used in required conditions
#' @param steps - time between populations
#' @param required_conditions - conditions that will be applied on patients that will be have a prediction score
#' @noRd
mldp_disease_empirical_prob_for_disease <- function(population, steps, required_conditions = NULL) {
    km <- list()
    if (is.null(required_conditions)) {
        required_conditions <- "id==id"
    }
    for (i in 1:length(population)) {
        k <- plyr::adply(0:2, 1, function(state) {
            if (state == 2) {
                young_pop <- population[[i]] %>% filter(is.na(disease) | disease > age) # ?do we need ot filter out sick?
            } else {
                young_pop <- population[[i]] %>% filter(is.na(disease) | disease > age, eval(rlang::parse_expr(required_conditions)))
                if (state == 0) {
                    if (i == 1) {
                        old_pop <- population[[i]] %>%
                            filter(!is.na(target_class)) %>%
                            pull(id)
                    } else {
                        old_pop <- c(
                            population[[i]] %>% filter(!is.na(target_class)) %>% pull(id),
                            population[[i - 1]] %>% filter(eval(rlang::parse_expr(required_conditions))) %>% pull(id)
                        )
                    }
                    young_pop <- young_pop %>% filter(id %in% old_pop)
                }
            }
            young_pop <- young_pop %>% mutate(status = factor(ifelse(step_outcome == "disease_death", "disease", step_outcome),
                levels = c("healthy", "disease", "death")
            ))

            km <- plyr::ddply(young_pop, plyr::.(sex), function(data) {
                if (nrow(data %>% filter(status != "healthy")) == 0) {
                    return(data.frame(time = steps[i], est = 0, var = 0, group = 1, outcome = "disease"))
                }
                fit <- cmprsk::cuminc(data$obsT, data$status, cencode = "healthy")
                ret <- purrr::map2_df(fit, names(fit), ~ as.data.frame(.x[1:3]) %>% mutate(name = .y)) %>%
                    tidyr::separate(name, into = c("group", "outcome"), sep = " (?=[^ ]*$)") %>%
                    filter(time <= steps[i]) %>%
                    arrange(desc(time), desc(est)) %>% # this is the estimate at latest measured time before outcome years.
                    distinct(outcome, .keep_all = TRUE)
                return(ret)
            })
            return(km %>% mutate(mode = state))
        })
        km[[i]] <- k
    }
    names(km) <- names(population)
    return(list(prob = km, required_conditions = required_conditions))
}



#' calculate expected number of disease patients
#' @param index - index of entry in the empirical disase prob of the current age to start propogation from
#' @param population_count - data.frame containing age, sex, and the number of patients available
#' @param edp - empirical disease prob, output of mldp_disease_empirical_prob_for_disease
mldp_disease_expected <- function(index, population_count, edp) {
    plyr::adply(population_count, 1, function(pc) {
        pc %>% mutate(disease_n = mldp_disease_expected_sex(
            pc$sex,
            pc$n,
            edp[1:index],
            mode = 0
        ))
    })
}

mldp_disease_expected_sex <- function(sex, n, disease_prob, mode) {
    if (n == 0 | length(disease_prob) == 0) {
        return(0)
    }
    disease_count <- round(n * disease_prob[[length(disease_prob)]] %>%
        filter(sex == !!sex, mode == !!mode, outcome == "disease") %>%
        pull(est))
    if (identical(disease_count, numeric(0))) {
        disease_count <- 0
    }

    # message(paste(age, " : ", size, " : ", round(disease_count)))
    next_n <- round(n * (1 - (disease_prob[[length(disease_prob)]] %>% filter(sex == !!sex, mode == !!mode) %>% pull(est) %>% sum())))
    return(disease_count + mldp_disease_expected_sex(sex, next_n, head(disease_prob, -1), min(mode + 1, 2)))
}


mldp_disease_assign_expected <- function(index, ordered_patients, empirical_disease_prob) {
    patients_filtered <- ordered_patients %>%
        filter(eval(rlang::parse_expr(empirical_disease_prob$required_conditions))) %>%
        mutate(target_class = 0)
    expected <- mldp_disease_expected(index, patients_filtered %>% count(sex), empirical_disease_prob$prob)
    patients_filtered_class <- plyr::adply(expected, 1, function(a) {
        top_res <- a$disease_n
        ret <- patients_filtered %>% inner_join(a %>% select(sex), by = "sex")
        ret$target_class[1:top_res] <- 1
        return(ret)
    }) %>% select(-n, -disease_n)
    return(patients_filtered_class %>%
        bind_rows(ordered_patients %>% anti_join(patients_filtered %>% select(id), by = "id")))
}
