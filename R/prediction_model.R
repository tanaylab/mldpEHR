# this file contains methods for building and analyzing prediction model

#' train an xgboost cross validation classification model with k-fold cross-validation
#' @param target - data.frame containing the patient id, sex, target_class (0/1) and fold (number used to assigne to cross validation folds)
#' @param features - data.frame containing patient id along with all other features to be used in classification model
#' @param folds - number of cross-validation folds
#' @param xgboost_params - parameters used for xgboost model training
#' @param nrounds - number of training rounds
#' @return a predictor, a list with the following elements
#' - model - list of xgboost models, for each fold
#' - train - data.frame containing the patients id, fold, target class and predicted value in training (each id was used in nfolds-1 for training)
#' - test - data.frame containing the patients id, fold, target class and predicted value in testing (each id was tested once in the fold it was not used for training)
#' - xgboost_params - the set of parameters used in xgboost
#' - nrounds - number of training iterations conducted
#' @examples
#' library(dplyr)
#' library(ggplot2)
#' target <- data.frame(id = 1:1000, sex = rep(c("male", "female"), 500), target_class = rep(c(0, 1), each = 500))
#' features <- data.frame(id = 1:500, a = rnorm(500), b = rnorm(500)) %>%
#'     bind_rows(
#'         data.frame(id = 501:1000, a = rnorm(500, mean = 0.5, sd = 2), b = rnorm(500, mean = -0.5, sd = 2))
#'     )
#' predictor <- mldpEHR.cv_train_outcome(target, features, folds = 3)
#' ggplot(predictor$test, aes(x = predict, colour = factor(target_class))) +
#'     geom_density() +
#'     theme_bw()
#' @export

mldpEHR.cv_train_outcome <- function(target,
                                     features,
                                     folds,
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
        select(id, target_class, fold) %>%
        left_join(features, by = "id")
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
        left_join(target %>% select(id, sex), by = "id") %>%
        group_by(sex) %>%
        mutate(qpredict = ecdf(predict)(predict)) %>%
        ungroup()
    return(list(model = model, test = test, train = train, target = target_fold, features = features, xgboost_params = xgboost_params, nrounds = nrounds))
}


#' build an xgboost cross validation classification model with k-fold cross-validation for each featureset provided, assumed that the classification is defined by the previoud model
#' @param base_name - name of predictor for latest time point.
#' @param base_predictor - classification model used as the emprical classifier at latest time point. This contains the target outcome model.
#' @param target_list - list of patients with target classification, going back in time, representing the previous (time) snapshot.
#' For patients with the target outcome known before the next model time, the target_classification should be set to 1.
#' Otherwise, target classificaiton should be set to NA
#' @param feature_list - list of features, each represents the previous (time) snapshot.
#' @param q_thresh - score quantile threshold for target classification of 1
#' @return the full list of predictors, according to provided target_list
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
#' predictors <- mldpEHR.cv_train_stitch_outcome(0, predictor, target_list, feature_list, q_thresh = 0.5)
#' test <- purrr::map2_df(predictors, names(predictors), ~ .x$test %>% mutate(n = .y))
#' ggplot(test, aes(x = predict, colour = factor(target_class))) +
#'     facet_wrap(~n, nrow = 1) +
#'     geom_density() +
#'     theme_bw()
#'
#' @export
mldpEHR.cv_train_stitch_outcome <- function(base_name,
                                            base_predictor,
                                            target_list,
                                            feature_list,
                                            q_thresh = 0.05) {
    predictors <- list()
    predictors[[as.character(base_name)]] <- base_predictor
    predictor <- base_predictor
    nfolds <- length(predictor$model)

    for (age in names(feature_list)) {
        message(age)
        # the predictor score will be used to set the target class for the next predictor (of younger age)
        target <- predictor$test %>%
            mutate(target_class = ifelse(qpredict < q_thresh, 0, 1)) %>%
            select(id, fold, target_class)

        # defining the target_class for the current age features
        # there are two options:
        # 1) patients at younger age that die within the step (known target_class),
        step_target <- target_list[[age]] %>% # looking at known target_class (will not appear in the older predictor score)
            filter(!is.na(target_class)) %>%
            group_by(sex) %>%
            mutate(fold = sample(1:nfolds, n(), replace = TRUE)) %>%
            ungroup() %>%
            select(id, sex, fold, target_class)
        # and
        # 2) patients that their target class is defined by the predictor score for the advanced age
        predictor_target <- target_list[[age]] %>%
            filter(is.na(target_class)) %>%
            select(id, sex) %>%
            left_join(target, by = "id") # this will leave patients with target_class=NA

        # combining the two:
        source_target <- step_target %>% bind_rows(predictor_target)

        # training the predictor
        predictor <- mldpEHR.cv_train_outcome(
            source_target,
            feature_list[[age]],
            nfolds,
            xgboost_params = predictor$xgboost_params,
            nrounds = predictor$nrounds
        )
        predictors[[age]] <- predictor
    }
    return(predictors)
}

#' analyze feature significance for a k-fold cross validation model using shaply values
#' @param predictor - classification model trained by mldpEHR.cv_train_outcome
#' @param nrounds - number of training rounds
#' @return a data.frame containing all features and their statistics:
#' - model - list of xgboost models, for each fold
#' - train - data.frame containing the patients id, fold, target class and predicted value in training (each id was used in nfolds-1 for training)
#' - test - data.frame containing the patients id, fold, target class and predicted value in testing (each id was tested once in the fold it was not used for training)
#' - xgboost_params - the set of parameters used in xgboost
#' - nrounds - number of training iterations conducted
#' @examples
#' library(dplyr)
#' library(ggplot2)
#' target <- data.frame(id = 1:1000, sex = rep(c("male", "female"), 500), target_class = rep(c(0, 1), each = 500))
#' features <- data.frame(id = 1:500, a = rnorm(500), b = rnorm(500)) %>%
#'     bind_rows(
#'         data.frame(id = 501:1000, a = rnorm(500, mean = 0.5, sd = 2), b = rnorm(500, mean = -0.5, sd = 2))
#'     )
#' predictor <- mldpEHR.cv_train_outcome(target, features, folds = 3)
#' @export

mldpEHR.cv_model_features <- function(predictor) {
    if (!"model" %in% names(predictor) | !"features" %in% names(predictor)) {
        stop("predictor must contain model and features information")
    }
    # going over all folds
    shap_fold <- plyr::adply(1:length(predictor$model), 1, function(fold) {
        model_fold <- predictor$model[[fold]]
        train_features_fold <- predictor$features %>%
            left_join(predictor$target %>% select(id, fold)) %>%
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
            )) %>% set_names(colnames(predictor$features %>% select(-id)))
        )
    shap_summary <- data.frame(
        feature = features,
        mean_abs_shap = colMeans(abs(shap_id %>% select(one_of(features))))
    )
    return(list(summary = shap_summary, shap_by_patient = shap_id, shap_by_fold = shap_fold))
}
