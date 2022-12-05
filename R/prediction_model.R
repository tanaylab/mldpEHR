# this file contains methods for building and analyzing prediction model

#' build an xgboost cross validation classification model with k-fold cross-validation for each featureset provided, assumed that the classification is defined by the previous model
#' @param patients - list of data.frames of all the patients in the system going back in time. For example the first data.frame represents age 80, next is 75 and so forth.
#' Each patient data.frame contains the following columns:
#' - patient id
#' - sex
#' - age
#' - death - age at death, NA if unknown
#' - followup - available followup time (in years) for this patient - time until end of database or until patient exists the system (not due to death)
#' and any additional columns required for patient filtering in the future
#' @param features - list of data.frames of features. Muste contain patient id column.
#' @param step - time between prediction models
#' @param nfolds - number of folds used for k-fold cross validation
#' @param required_conditions - any filter to apply to the patients to filter out training/testing
#' samples (e.g. missing data)
#' @param q_thresh - score quantile threshold for target classification of 1
#' @param xgboost_params - parameters used for xgboost model training
#' @param nrounds - number of training rounds
#' @return the full list of predictors, according to provided patients. Each predictor is a list
#' with the following members:
#' - model - list of xgboost models, for each fold
#' - train - data.frame containing the patients id, fold, target class and predicted value in training
#' (each id was used in nfolds-1 for training)
#' - test - data.frame containing the patients id, fold, target class and predicted value in testing
# (each id was tested once in the fold it was not used for training)
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
#'     sex = rep(c(1, 2), N / 2),
#'     age = 80 - .x * 5,
#'     death = c(rep(NA, 0.2 * N), rep(82, 0.8 * N)),
#'     followup = .x * 5 + 5
#' )) %>%
#'     setNames(seq(80, by = -5, length.out = 6))
#' features <- purrr::map(0:5, ~ data.frame(
#'     id = 1:N,
#'     a = c(rnorm(0.2 * N), rnorm(0.8 * N, mean = 2, sd = 0.5))
#' )) %>% setNames(seq(80, by = -5, length.out = 6))
#' predictors <- mldpEHR.mortality_multi_age_predictors(patients, features, 5, 3, q_thresh = 0.2)
#' test <- purrr::map2_df(predictors, names(predictors), ~ .x$test %>%
#'     mutate(n = .y) %>%
#'     arrange(id) %>%
#'     mutate(
#'         outcome =
#'             c(
#'                 rep("alive", 0.2 * N),
#'                 rep("death", 0.8 * N)
#'             )
#'     ))
#' ggplot(test, aes(x = predict, colour = factor(outcome))) +
#'     facet_wrap(~n, nrow = 1) +
#'     stat_ecdf() +
#'     theme_bw()

#' @export


mldpEHR.mortality_multi_age_predictors <- function(patients,
                                                   features,
                                                   step,
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
    predictors <- list()
    pop <- purrr::map(1:length(patients), ~ .mldpEHR.compute_target_mortality(patients[[.x]], step, .x == 1)) %>%
        purrr::set_names(names(patients))

    predictors[[1]] <- .mldpEHR.cv_train_outcome(pop[[1]], features[[1]], nfolds, required_conditions)
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
        predictors[[i]] <- .mldpEHR.cv_train_outcome(
            source_target,
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
#' predictors <- mldpEHR.disease_multi_age_predictors(patients, features, 5, 3)
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


mldpEHR.disease_multi_age_predictors <- function(patients,
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

    pop <- purrr::map(1:length(patients), ~ .mldpEHR.compute_target_disease(patients[[.x]], step, .x == 1)) %>%
        purrr::set_names(names(patients))
    empirical_disease_prob <- .mldpEHR.disease_empirical_prob_for_disease(pop, step, required_conditions)
    predictors[[1]] <- .mldpEHR.cv_train_outcome(pop[[1]], features[[1]], nfolds, required_conditions)

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
        source_disease <- .mldpEHR.disease_assign_expected(i, source_target, empirical_disease_prob)

        # need to add the patients with missing target
        source_disease <- source_disease %>% bind_rows(pop[[i]] %>% anti_join(source_disease %>% select(id), by = "id"))

        # training the predictor
        predictors[[i]] <- .mldpEHR.cv_train_outcome(
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





#' train an xgboost cross validation classification model with k-fold cross-validation
#' @param target - data.frame containing the patient id, sex, target_class (0/1) and fold
#' (number used to assigne to cross validation folds)
#' @param features - data.frame containing patient id along with all other features to be used in
#' classification model
#' @param folds - number of cross-validation folds
#' @param required_conditions - any filter to apply to the features to filter out training/testing samples (e.g. missing data)
#' @param xgboost_params - parameters used for xgboost model training
#' @param nrounds - number of training rounds
#' @return a predictor, a list with the following elements
#' - model - list of xgboost models, for each fold
#' - train - data.frame containing the patients id, fold, target class and predicted value in training
#' (each id was used in nfolds-1 for training)
#' - test - data.frame containing the patients id, fold, target class and predicted value in testing
#' (each id was tested once in the fold it was not used for training)
#' - xgboost_params - the set of parameters used in xgboost
#' - nrounds - number of training iterations conducted
.mldpEHR.cv_train_outcome <- function(target,
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



#' analyze feature significance for a k-fold cross validation model using shaply values
#' @param predictor - classification model trained by mldpEHR.cv_train_outcome
#' @return a list with the 3 objects
#' - summary - data frame containing for each feature the mean absolute shaply value for the feature across all training data
#' - shap_by_patient - data frame containing for each patient and feature the feature value and mean shap value across all training folds
#' - shap_by_fold - similar to shap_by_patient, but for each fold seperately
#' @examples
#' N <- 100
#' patients <- list(data.frame(
#'     id = 1:N,
#'     sex = rep(c(1, 2), N / 2),
#'     age = 80,
#'     death = c(rep(NA, 0.2 * N), rep(82, 0.8 * N)),
#'     followup = 5
#' ))
#' features <- list(data.frame(
#'     id = 1:N,
#'     a = c(rnorm(0.2 * N), rnorm(0.8 * N, mean = 2, sd = 1)),
#'     b = rep(c(rnorm(N / 4), rnorm(N / 4, mean = 3)), 2)
#' ))
#' predictor <- mldpEHR.mortality_multi_age_predictors(patients, features, 5, 3, q_thresh = 0.05)
#' predictor_features <- mldpEHR.prediction_model_features(predictor[[1]])
#' @export

mldpEHR.prediction_model_features <- function(predictor) {
    if (!"model" %in% names(predictor) | !"features" %in% names(predictor)) {
        stop("predictor must contain model and features information")
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


.mldpEHR.compute_target_mortality <- function(pop, step, final_outcome) {
    return(pop %>% mutate(
        step_outcome = ifelse(!is.na(death) & death <= age + step & death <= age + followup, "death", "alive"),
        target_class = ifelse(followup < step,
            NA,
            ifelse(step_outcome == "death", 1, ifelse(final_outcome, 0, NA))
        ),
        obsT = pmax(0, pmin(step, death - age, followup, na.rm = TRUE))
    ))
}


.mldpEHR.compute_target_disease <- function(pop, step, final_outcome) {
    return(pop %>% mutate(
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
    ))
}
