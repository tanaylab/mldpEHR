
#' Predict scores / probabilities of multi-age models for new data
#'
#' @description Given a multi-age model, predict scores for new data. If a Markov model
#' is provided, life-long probabilities are calculated as well. NA would be returned for patients with age outside the range of the model.
#'
#' @param data A data frame with "id", "age", "sex" and the same features as the training data, or a MldpEHR object. The "age" field should be in years and the "sex" field should contain
#' 1 for male and 2 for female.
#' @param predictors A list of predictors, one for each age group. Output of \code{mldp_mortality_multi_age_predictors} or \code{mldp_disease_multi_age_predictors}.
#' Note that if predictor model is not a single XGBoost model but a list of models (e.g. from cross validation) then average model score will be computed.
#' @param markov_models A Markov model, output of \code{mldp_mortality_markov} or \code{mldp_disease_markov}. Can be NULL, and then only scores are calculated.
#' @param outcome A character vector indicating the outcome to calculate the life-long probabilities for. For mortality models use \code{"death"} and for disease models please use \code{c("disease", "disease_death")}. If NULL, the function would try to infer which outcome to use ("death" or "disease + disease_death") type: and if it fails, it chooses the first outcome in the \code{markov$prob} data frame.
#' @param years_before_model The number of years to allow the age of the patient to be before the age the model was trained on. e.g., if the model was trained on patients 30-80 years old, and \code{years_before_model} is 5, then patients with age 25-80 would be allowed. Defaults to 5.
#'
#' @return a data frame with the following columns:
#' \itemize{
#' \item{id: }{unique patient id}
#' \item{age: }{patient age.}
#' \item{sex: }{patient sex}
#' \item{model_age: }{age of the model used to predict}
#' \item{score: }{the predicted score for the patient}
#' \item{quantile: }{the quantile of the score of the patient}
#' \item{prob: }{the predicted life-long probability for the patient (in case \code{markov} is provided)}
#' }
#'
#'
#'
#' @examples
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
#' predictors <- mldp_model_export(predictors)
#'
#' new_data <- load_mortality_example_data(N = 1e3, num_age_groups = 3)
#' scores <- mldp_predict_multi_age(new_data, predictors, markov, "death")
#'
#' # For diseases
#'
#' disease_data <- load_disease_example_data(N = 100, num_age_groups = 3)
#'
#' # Build predictors
#' predictors_disease <- mldp_disease_multi_age_predictors(
#'     disease_data@patients,
#'     disease_data@features,
#'     target_years = 5,
#'     nfolds = 2
#' )
#' markov_disease <- mldp_disease_markov(predictors_disease, 5, qbins = seq(0, 1, by = 0.1))
#'
#' new_data <- load_disease_example_data(N = 1e3, num_age_groups = 3)
#'
#' scores_disease <- mldp_predict_multi_age(
#'     new_data,
#'     predictors_disease,
#'     markov_disease,
#'     c("disease", "disease_death")
#' )
#'
#' head(scores)
#'
#' @export
mldp_predict_multi_age <- function(data, predictors, markov_models = NULL, outcome = NULL, years_before_model = 5) {
    if (is(data, "MldpEHR")) {
        data <- purrr::map2_dfr(data@patients, data@features, function(p, f) {
            p %>%
                select(id, sex, age) %>%
                left_join(f %>% select(-any_of(c("age", "sex"))), by = "id")
        })
    }

    # make sure that we have "id", "sex" and "age" columns
    purrr::walk(c("id", "sex", "age"), ~ {
        if (!(.x %in% colnames(data))) {
            cli::cli_abort("Column {.field {.x}} not found in data")
        }
    })

    if (any(data$age < 0)) {
        cli::cli_abort("Age must be positive")
    }

    if (!all(data$sex %in% c(1, 2))) {
        cli::cli_abort("The {.field sex} column should contain 1 for male and 2 for female")
    }

    model_ages <- purrr::map_dbl(predictors, "age")
    rev_model_ages <- rev(model_ages)

    # for each patient, find closest age to model_ages after the patient's age
    get_patient_model_age <- function(age) {
        ma <- rev_model_ages[rev_model_ages >= age][1]
        if (is.na(ma)) {
            return(NA)
        }
        if (age + years_before_model < ma) {
            return(NA)
        }
        return(ma)
    }


    data <- data %>%
        mutate(model_age = purrr::map_dbl(age, get_patient_model_age))

    if (any(is.na(data$model_age))) {
        cli::cli_alert_warning("Some patients were not assigned to a model age. This is probably because their age is either too low or too high. Consider increasing the {.field years_before_model} parameter.")
        cli::cli_alert_warning("The range of ages in the data is: {.val {min(data$age)}}-{.val {max(data$age)}}, while in the models it is: {.val {min(model_ages)}}-{.val {max(model_ages)}}. {.field years_before_model} is {.val {years_before_model}}.")
    }


    scores <- plyr::ddply(data, "model_age", function(x) {
        if (is.na(x$model_age[1])) {
            return(tibble(id = x$id, age = x$age, sex = x$sex, model_age = NA, score = NA, quantile = NA))
        }
        x <- as_tibble(x)
        model_age <- as.character(x$model_age[1])
        model <- predictors[[model_age]]
        # make sure that all features are present in the data. If not, add them with NA
        feature_names <- setdiff(model$feature_names, "id")
        for (f in feature_names) {
            if (!f %in% colnames(x)) {
                cli::cli_alert("Feature {.field {f}} is missing in the data at model age {.val {x$model_age[1]}}. Adding it with NA.")
                x[[f]] <- NA
            }
        }

        # predict scores
        features_mat <- as.matrix(x[, feature_names])
        if (class(model$model) == "xgb.Booster") {
            scores <- predict(model$model, features_mat)
        } else { #assuming list of predictors, will compute score for each and then average them out
            scores <- rowMeans(do.call(cbind, purrr::map(model$model, ~ predict(.x, features_mat))))
        }
        quantiles <- rep(NA, length(scores))
        quantiles[x$sex == 1] <- model$score2quantile[[1]](scores[x$sex == 1])
        quantiles[x$sex == 2] <- model$score2quantile[[2]](scores[x$sex == 2])

        return(tibble(id = x$id, age = x$age, sex = x$sex, model_age = model_age, score = scores, quantile = quantiles))
    })

    # calculate life-long probabilities
    if (!is.null(markov_models)) {
        scores <- mldp_score_to_prob(scores, markov_models, outcome)
    }

    return(scores)
}

#' Calculate the probability of a disease or death given model scores
#'
#' @param scores A data frame with model scores. Output of \code{mldp_predict_multi_age}.
#' @return A data frame with the same columns as \code{scores} and an additional column \code{prob} with the probability of the outcome.
#'
#' @examples
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
#' new_data <- load_mortality_example_data(N = 1e3, num_age_groups = 3)
#' scores <- mldp_predict_multi_age(new_data, predictors)
#'
#' markov <- mldp_mortality_markov(predictors, 5, qbins = seq(0, 1, by = 0.1))
#' probs <- mldp_score_to_prob(scores, markov_models = markov, outcome = "death")
#' head(probs)
#'
#' @export
#' @rdname mldp_predict_multi_age
mldp_score_to_prob <- function(scores, markov_models = NULL, outcome = NULL, years_before_model = 5) {
    # make sure all the fields are present
    fields <- c("id", "age", "sex", "model_age", "score", "quantile")
    purrr::walk(fields, ~ {
        if (!(.x %in% colnames(scores))) {
            cli::cli_abort("Column {.field {.x}} not found in scores. Use {.code mldp_predict_multi_age} to get scores.")
        }
    })

    if (is.null(outcome)) {
        if (all(c("disease", "disease_death") %in% colnames(markov_models$prob))) {
            outcome <- c("disease", "disease_death")
        } else if ("death" %in% colnames(markov_models$prob)) {
            outcome <- "death"
        } else {
            outcome <- setdiff(colnames(markov_models$prob), c("age", "sbin", "sex"))[1]
        }
    }

    cli::cli_alert_info("Using {.field {.val {paste(outcome, collapse = ' + ')}}} as outcome")

    probs <- plyr::ddply(scores, "model_age", function(x) {
        if (is.na(x$model_age[1])) {
            return(scores %>%
                filter(is.na(model_age)) %>%
                mutate(prob = NA))
        }
        model_age <- as.character(x$model_age[1])

        sex_probs <- purrr::map_dfr(c(1, 2), function(sex) {
            model_probs <- markov_models$prob %>% filter(as.character(age) == model_age, sex == !!sex)
            sbins <- model_probs$sbin[!(model_probs$sbin %in% colnames(model_probs)) & model_probs$sbin != "no_score"]
            model_probs <- model_probs %>%
                filter(sbin %in% sbins) %>%
                mutate(qmax = 1 / 20 * as.numeric(sbin), qmin = qmax - 1 / 20, q = pmean(qmax, qmin))
            model_probs$p <- rowSums(model_probs[, outcome, drop = FALSE])

            f <- approxfun(c(0, model_probs$q, 1), c(min(model_probs$p), model_probs$p, max(model_probs$p)))

            cur_scores <- scores %>%
                filter(model_age == !!model_age, sex == !!sex) %>%
                mutate(prob = f(quantile))
            return(cur_scores)
        })

        return(sex_probs)
    })

    scores <- probs %>%
        select(id, age, sex, model_age, score, quantile, prob)

    return(scores)
}
