# this file contains methods for building and analyzing prediction model

#' Build a set of mldpEHR prediction models for mortality
#'
#' @description Given a set of patients and features, builds a set of models for predicting mortality up to \code{survival_years} after the oldest age (e.g. 85) by stitching together models for each age group. Each model is build using a cross validation \code{xgboost} k-fold classification model, where the number of folds is determined by the \code{nfolds} parameter. Stitching is done by first modeling the oldest group, and then taking the top \code{q_thresh} quantile of the predicted score for each patient in the next younger group and building a model to classify them. See more in the vignette.
#'
#' @param survival_years the number of years to predict mortality for (after the oldest age). e.g. 5 for predicting mortality up to age 85
#' @param nfolds number of folds used for k-fold cross validation
#' @param required_conditions a string with an expression for any filter to apply to the patients to filter out from training or testing. Can be used to filter out patients with missing data or other conditions. See \code{\link{dplyr::filter}} for more details.
#' @param q_thresh quantile of the predicted score to use for the next age group. Default is 0.05, which means the top 5% of the predicted score will be used for the next age group.
#' @param xgboost_params parameters used for xgboost model training
#' @param nrounds number of training rounds
#' @return a list of predictors, for each age. Each predictor is a list
#' with the following members:
#' \itemize{
#' \item{model: }{a list of xgboost models, for each fold}
#' \item{train: }{data.frame containing the patients id, fold, target class and predicted value in training}
#' (each id was used in nfolds-1 for training)
#' \item{test: }{data.frame containing the patients id, fold, target class and predicted value in testing
#' (each id was tested once in the fold it was not used for training)}
#' \item{xgboost_params: }{the set of parameters used in xgboost}
#' \item{nrounds: }{number of training iterations conducted}
#' }
#'
#'
#' @examples
#'
#' # Load the example data
#' mortality <- load_mortality_example_data()
#'
#' #' Build predictors
#' predictors <- mldp_mortality_multi_age_predictors(mortality@patients, mortality@features, survival_years = 5, nfolds = 3, q_thresh = 0.2)
#'
#' # Plot the ecdf of the score of each model for colored by outcome
#' mldp_plot_multi_age_predictors_ecdf(predictors)
#'
#' @inheritParams MldpEHR
#' @export
mldp_mortality_multi_age_predictors <- function(patients,
                                                features,
                                                survival_years,
                                                nfolds,
                                                required_conditions = "id==id",
                                                q_thresh = 0.05,
                                                xgboost_params = list(
                                                    booster = "gbtree",
                                                    objective = "binary:logistic",
                                                    subsample = 0.7,
                                                    max_depth = 3,
                                                    colsample_bytree = 1,
                                                    eta = 0.05,
                                                    min_child_weight = 1,
                                                    gamma = 0,
                                                    eval_metric = "auc"
                                                ),
                                                nrounds = 1000) {
    mldp <- MldpEHR(patients, features)
    steps <- c(survival_years, mldp_get_age_steps(mldp))

    predictors <- list()
    pop <- purrr::map(seq_along(mldp@patients), ~ mldp_compute_target_mortality(mldp@patients[[.x]], steps[.x], .x == 1)) %>%
        purrr::set_names(names(mldp@patients))

    predictors[[1]] <- mldp_cv_train_outcome(pop[[1]], mldp@features[[1]], nfolds, required_conditions)
    i <- 2
    while (i <= length(pop)) {
        # the predictor score will be used to set the target class for the next predictor (of younger age)
        target <- predictors[[i - 1]]$test %>%
            mutate(target_class = ifelse(qpredict < q_thresh, 0, 1)) %>%
            select(id, fold, target_class)

        # defining the target_class for the current age features
        # there are two options:
        # 1) patients at younger age that die within the step (known target_class),
        step_target <- pop[[i]] %>% # looking at known target_class (will not appear in the older predictor score)
            filter(!is.na(target_class)) %>%
            group_by(sex) %>%
            mutate(fold = sample(1:nfolds, n(), replace = TRUE)) %>%
            ungroup()

        # 2) patients that their target class is defined by the predictor score for the advanced age
        predictor_target <- pop[[i]] %>%
            filter(is.na(target_class)) %>%
            select(-any_of(c("target_class", "fold"))) %>%
            left_join(target, by = "id") # this will leave patients with target_class=NA

        # combining the two:
        source_target <- step_target %>% bind_rows(predictor_target)

        # training the predictor
        predictors[[i]] <- mldp_cv_train_outcome(
            source_target,
            mldp@features[[i]],
            nfolds,
            required_conditions = required_conditions,
            xgboost_params = predictors[[i - 1]]$xgboost_params,
            nrounds = predictors[[i - 1]]$nrounds
        )
        i <- i + 1
    }
    names(predictors) <- names(mldp@patients)
    return(predictors)
}




#' build an xgboost cross validation classification model with k-fold cross-validation for each featureset provided, assumed that the classification is defined by the previous model
#' @param patients - list of data.frames of all the patients in the system going back in time. For example the first data.frame represents age 80, next is 75 and so forth.
#' Each patient data.frame contains the following columns:
#' - patient id
#' - sex
#' - age
#' - death - age at death, NA if unknown
#' - disease - age at disease, NA if unknown
#' - followup - available followup time (in years) for this patient - time until end of database or until patient exists the system (not due to death)
#' and any additional columns required for patient filtering in the future
#' @param features - list of data.frames of features
#' @param step - time between prediction models
#' @param nfolds - number of folds used for k-fold cross validation
#' @param required_conditions - any filter to apply to the features to filter out training/testing samples (e.g. missing data)
#' @param xgboost_params - parameters used for xgboost model training
#' @param nrounds - number of training rounds
#' @return the full list of predictors, according to provided patients, Each predictor is a list with the following members:
#' - model - list of xgboost models, for each fold
#' - train - data.frame containing the patients id, fold, target class and predicted value in training (each id was used in nfolds-1 for training)
#' - test - data.frame containing the patients id, fold, target class and predicted value in testing (each id was tested once in the fold it was not used for training)
#' - xgboost_params - the set of parameters used in xgboost
#' - nrounds - number of training iterations conducted

#' @examples
#'
#' library(dplyr)
#' library(ggplot2)
#' # build base predictor
#' N <- 1000
#' patients <- purrr::map(0:5, ~ data.frame(
#'     id = 1:N,
#'     sex = rep(1, N),
#'     age = 80 - .x * 5,
#'     death = c(rep(NA, 0.4 * N), rep(82, 0.6 * N)),
#'     disease = rep(rep(c(NA, 81), each = N / 4), 2),
#'     followup = .x * 5 + 5
#' )) %>%
#'     setNames(seq(80, by = -5, length.out = 6))
#' features <- purrr::map(0:5, ~ data.frame(
#'     id = 1:N,
#'     a = c(rnorm(0.4 * N), rnorm(0.6 * N, mean = 2, sd = 1)),
#'     b = rep(c(rnorm(N / 4), rnorm(N / 4, mean = 3)), 2)
#' )) %>% setNames(seq(80, by = -5, length.out = 6))
#' predictors <- mldp_disease_multi_age_predictors(patients, features, 5, 3)
#'
#'
#' test <- purrr::map2_df(predictors, names(predictors), ~ .x$test %>%
#'     mutate(n = .y) %>%
#'     arrange(id) %>%
#'     mutate(
#'         outcome =
#'             c(
#'                 rep("healthy", 0.25 * N),
#'                 rep("disease", 0.15 * N),
#'                 rep("disease_death", 0.1 * N),
#'                 rep("death", 0.25 * N),
#'                 rep("disease_death", 0.25 * N)
#'             )
#'     ))
#' ggplot(test, aes(x = predict, colour = factor(outcome))) +
#'     facet_wrap(~n, nrow = 1) +
#'     stat_ecdf() +
#'     theme_bw()
#'
#' @export
mldp_disease_multi_age_predictors <- function(patients,
                                              features,
                                              step,
                                              nfolds,
                                              required_conditions = "id==id",
                                              xgboost_params = list(
                                                  booster = "gbtree",
                                                  objective = "binary:logistic",
                                                  subsample = 0.7,
                                                  max_depth = 3,
                                                  colsample_bytree = 1,
                                                  eta = 0.05,
                                                  min_child_weight = 1,
                                                  gamma = 0,
                                                  eval_metric = "auc"
                                              ),
                                              nrounds = 1000) {
    predictors <- list()

    pop <- purrr::map(1:length(patients), ~ mldp_compute_target_disease(patients[[.x]], step, .x == 1)) %>%
        purrr::set_names(names(patients))
    empirical_disease_prob <- .mldp_disease_empirical_prob_for_disease(pop, step, required_conditions)
    predictors[[1]] <- mldp_cv_train_outcome(pop[[1]], features[[1]], nfolds, required_conditions)

    i <- 2
    while (i <= length(pop)) {
        # the predictor score will be used to set the target class for the next predictor (of younger age)
        target <- predictors[[i - 1]]$test %>%
            filter(!is.na(qpredict)) %>%
            arrange(desc(qpredict)) %>%
            select(id, fold, qpredict)

        # defining the target_class for the current age features
        # there are two options:
        # 1) patients at younger age that get sick
        step_target <- pop[[i]] %>% # looking at known target_class (will not appear in the older predictor score)
            filter(!is.na(target_class)) %>%
            group_by(sex) %>%
            mutate(fold = sample(1:nfolds, n(), replace = TRUE)) %>%
            ungroup() %>%
            mutate(qpredict = 1)

        # 2) patients that their target class is defined by the predictor score for the advanced age
        predictor_target <- pop[[i]] %>%
            filter(is.na(target_class)) %>%
            select(-any_of(c("target_class", "fold"))) %>%
            inner_join(target, by = "id") %>%
            arrange(desc(qpredict)) # %>%
        # select(-qpredict)

        # combining the two:
        source_target <- step_target %>%
            bind_rows(predictor_target)

        # find out how many will have the disease
        source_disease <- .mldp_disease_assign_expected(i, source_target, empirical_disease_prob)

        # need to add the patients with missing target
        source_disease <- source_disease %>% bind_rows(pop[[i]] %>% anti_join(source_disease %>% select(id), by = "id"))

        # training the predictor
        predictors[[i]] <- mldp_cv_train_outcome(
            source_disease,
            features[[i]],
            nfolds,
            required_conditions = required_conditions,
            xgboost_params = predictors[[i - 1]]$xgboost_params,
            nrounds = predictors[[i - 1]]$nrounds
        )
        i <- i + 1
    }
    names(predictors) <- names(patients)
    return(predictors)
}





#' Train an xgboost cross validation classification model with k-fold cross-validation
#'
#' @param target a data.frame containing the patient id, sex, target_class (0/1) and fold
#' (number used to assigne to cross validation folds)
#' @param features a data.frame containing patient id along with all other features to be used in
#' classification model
#' @param folds number of cross-validation folds
#' @param xgboost_params parameters used for xgboost model training
#' @param nrounds number of training rounds
#' @return a predictor, a list with the following elements
#' \itemize{
#' \item{model: }{a list of xgboost models, for each fold}
#' \item{train: }{data.frame containing the patients id, fold, target class and predicted value in training}
#' (each id was used in nfolds-1 for training)
#' \item{test: }{data.frame containing the patients id, fold, target class and predicted value in testing
#' (each id was tested once in the fold it was not used for training)}
#' \item{xgboost_params: }{the set of parameters used in xgboost}
#' \item{nrounds: }{number of training iterations conducted}
#' }
#'
#' @inheritParams mldp_disease_multi_age_predictors
#'
#' @noRd
mldp_cv_train_outcome <- function(target,
                                  features,
                                  folds,
                                  required_conditions = "id==id",
                                  xgboost_params = list(
                                      booster = "gbtree",
                                      objective = "binary:logistic",
                                      subsample = 0.7,
                                      max_depth = 3,
                                      colsample_bytree = 1,
                                      eta = 0.05,
                                      min_child_weight = 1,
                                      gamma = 0,
                                      eval_metric = "auc"
                                  ),
                                  nrounds = 1000) {
    # assign folds to patiens, controlling for sex and target randomly
    if (!"fold" %in% colnames(target)) {
        target$fold <- NA
    }
    target_fold <- target %>%
        filter(is.na(fold)) %>%
        group_by(sex, target_class) %>%
        mutate(fold = sample(1:folds, n(), replace = TRUE)) %>%
        ungroup() %>%
        bind_rows(target %>% filter(!is.na(fold)))


    target_features <- target_fold %>%
        filter(eval(rlang::parse_expr(required_conditions))) %>%
        select(id, target_class, fold) %>%
        inner_join(features, by = "id")

    pb <- progress::progress_bar$new(
        format = "  Training [:bar] :current/:total (:percent) in :elapsed",
        total = folds, clear = FALSE, width = 60, show_after = 0
    )
    invisible(pb$tick(0))
    model_folds <- purrr::map(1:max(target_fold$fold), function(cur_fold) {
        pb$tick()

        dtrain <- xgboost::xgb.DMatrix(
            data = as.matrix(target_features %>% filter(fold != cur_fold, !is.na(target_class)) %>% select(-id, -fold, -target_class)),
            label = target_features %>% filter(fold != cur_fold, !is.na(target_class)) %>% pull(target_class),
            missing = NA
        )
        dtest <- as.matrix(target_features %>% filter(fold == cur_fold) %>% select(-id, -fold, -target_class))

        xgb <- xgboost::xgb.train(data = dtrain, nthread = 24, params = xgboost_params, nrounds = nrounds, verbose = 1)
        train <- target_features %>%
            filter(fold != cur_fold, !is.na(target_class)) %>%
            select(id, fold, target_class) %>%
            mutate(predict = predict(xgb, dtrain))
        test <- target_features %>%
            filter(fold == cur_fold) %>%
            select(id, fold, target_class) %>%
            mutate(predict = predict(xgb, dtest))
        return(list(model = xgb, test = test, train = train))
    })
    model <- purrr::map(model_folds, ~ .x$model)
    train <- purrr::map_df(model_folds, ~ .x$train)
    test <- purrr::map_df(model_folds, ~ .x$test) %>%
        full_join(target_fold, by = c("id", "fold", "target_class")) %>%
        group_by(sex) %>%
        mutate(qpredict = ecdf(predict)(predict)) %>%
        ungroup()
    return(list(model = model, test = test, train = train, target = target_fold, features = features, xgboost_params = xgboost_params, nrounds = nrounds))
}



#' Extract statistics on the features used in a k-fold cross-validation model using SHAP (SHapley Additive exPlanations)
#'
#'
#'
#' @param predictor a single predictor, a single element of the list returned by \code{mldp_disease_multi_age_predictors} or \code{mldp_mortality_multi_age_predictors}
#'
#' @return a list with 3 elements:
#' \itemize{
#' \item{summary: }{a data frame containing for each feature the mean absolute shaply value for the feature across all training data.}
#' \item{shap_by_patient: }{a data frame data frame containing for each patient and feature the feature value and mean shap value across all training folds.}
#' \item{shap_by_fold: }{similar to shap_by_patient, but for each fold seperately.}
#' }
#'
#' @examples
#' # Load the example data
#' mortality <- load_mortality_example_data(100)
#'
#' #' Build predictors
#' predictors <- mldp_mortality_multi_age_predictors(mortality@patients, mortality@features, survival_years = 5, nfolds = 3, q_thresh = 0.05)
#'
#' predictor_features <- mldp_prediction_model_features(predictors[[1]])
#'
#' names(predictor_features)
#' head(predictor_features$summary)
#' head(predictor_features$shap_by_patient)
#' head(predictor_features$shap_by_fold)
#'
#' @export
mldp_prediction_model_features <- function(predictor) {
    if (!"model" %in% names(predictor) | !"features" %in% names(predictor)) {
        cli::cli_abort("predictor must be a list with elements {.field model} and {.field features}")
    }
    # going over all folds
    shap_fold <- plyr::adply(1:length(predictor$model), 1, function(fold) {
        model_fold <- predictor$model[[fold]]
        train_features_fold <- predictor$features %>%
            left_join(predictor$target %>% select(id, fold), by = "id") %>%
            filter(fold != !!fold)
        train_features_fold %>%
            select(id) %>%
            bind_cols(predict(model_fold,
                as.matrix(train_features_fold %>% select(one_of(model_fold$feature_names))),
                predcontrib = TRUE,
                approxcontrib = TRUE
            ) %>%
                as_tibble())
    }) %>%
        rename(fold = X1) %>%
        arrange(id)

    # compute average contribution per patient (avereging folds)
    number_of_contributions_per_patient <- length(predictor$model) - 1

    features <- colnames(predictor$features %>% select(-id))
    shap_id <- shap_fold %>%
        distinct(id) %>%
        bind_cols(
            as.data.frame(do.call(
                cbind,
                purrr::map(features, ~
                    .colMeans(shap_fold[, .x], number_of_contributions_per_patient, nrow(shap_fold) / number_of_contributions_per_patient))
            )) %>% purrr::set_names(colnames(predictor$features %>% select(-id)))
        )
    shap_id_val <- shap_id %>%
        pivot_longer(!id, names_to = "feature", values_to = "shap") %>%
        left_join(
            predictor$features %>% pivot_longer(!id, names_to = "feature", values_to = "value"),
            by = c("id", "feature")
        )
    shap_summary <- data.frame(
        feature = features,
        mean_abs_shap = colMeans(abs(shap_id %>% select(one_of(features))))
    )
    return(list(summary = shap_summary, shap_by_patient = shap_id_val, shap_by_fold = shap_fold))
}


mldp_compute_target_mortality <- function(pop, step, final_outcome) {
    pop %>%
        mutate(
            step_outcome = ifelse(!is.na(death) & death <= age + step & death <= age + followup, "death", "alive"),
            target_class = ifelse(followup < step,
                NA,
                ifelse(step_outcome == "death", 1, ifelse(final_outcome, 0, NA))
            ),
            obsT = pmax(0, pmin(step, death - age, followup, na.rm = TRUE))
        )
}


mldp_compute_target_disease <- function(pop, step, final_outcome) {
    pop %>%
        mutate(
            step_outcome =
                ifelse(!is.na(disease) & disease <= age + step,
                    ifelse(!is.na(death) & death <= age + step & death <= age + followup, "disease_death", "disease"),
                    ifelse(!is.na(death) & death <= age + step & death <= age + followup, "death", "healthy")
                ),
            target_class = ifelse(followup < step,
                NA,
                ifelse(step_outcome == "disease" | step_outcome == "disease_death",
                    1,
                    ifelse(step_outcome == "death", NA, ifelse(final_outcome, 0, NA))
                )
            ),
            obsT = ifelse(!is.na(disease) & disease < age, # already sick
                pmax(0, pmin(step, death - age, followup, na.rm = TRUE)),
                pmax(0, pmin(step, disease - age, death - age, followup, na.rm = TRUE))
            )
        )
}
