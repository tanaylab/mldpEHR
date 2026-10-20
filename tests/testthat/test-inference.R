set.seed(17)
mortality <- load_mortality_example_data(N = 200, num_age_groups = 2)
predictors <- suppressMessages(mldp_mortality_multi_age_predictors(
    mortality@patients,
    mortality@features,
    survival_years = 5,
    nfolds = 2,
    q_thresh = 0.2,
    nthread = 2
))
new_data <- load_mortality_example_data(N = 100, num_age_groups = 2)
new_data <- purrr::map2_dfr(new_data@patients, new_data@features, ~ .x %>%
    dplyr::select(id, sex, age) %>%
    dplyr::left_join(.y, by = "id"))

# The average score of boosters(predictor) for each patient, at the age of their model (the patients are at model ages)
expected_scores <- function(predictors, boosters) {
    boosters <- purrr::as_mapper(boosters)
    purrr::map_dfr(names(predictors), function(age) {
        x <- new_data %>% dplyr::filter(age == !!as.numeric(age))
        scores <- purrr::map(boosters(predictors[[age]]), ~ predict(.x, as.matrix(x[, .x$feature_names])))
        return(x %>% dplyr::transmute(id, model_age = !!age, expected = rowMeans(do.call(cbind, scores))))
    })
}

compare_scores <- function(scores, expected) {
    scores <- scores %>% dplyr::inner_join(expected, by = c("id", "model_age"))
    expect_equal(nrow(scores), nrow(expected))
    expect_equal(scores$score, scores$expected, tolerance = 1e-6)
}

test_that("predictors score by their model trained on all the data", {
    scores <- mldp_predict_multi_age(new_data, predictors)
    compare_scores(scores, expected_scores(predictors, ~ list(.x$model)))
})

test_that("use_cv scores by the average of the cross validation models", {
    scores <- mldp_predict_multi_age(new_data, predictors, use_cv = TRUE)
    compare_scores(scores, expected_scores(predictors, ~ .x$model_cv))
    expect_false(isTRUE(all.equal(scores$score, mldp_predict_multi_age(new_data, predictors)$score)))
})

test_that("a predictor whose model is a list of models scores by their average, without age and features fields", {
    cv_predictors <- purrr::map(predictors, ~ list(model = .x$model_cv, score2quantile = .x$score2quantile))
    scores <- mldp_predict_multi_age(new_data, cv_predictors)
    expect_equal(scores, mldp_predict_multi_age(new_data, predictors, use_cv = TRUE))
})

test_that("each model gets its features in its own order", {
    x <- as.matrix(new_data %>% dplyr::filter(age == 80) %>% dplyr::select(feature1, feature2))
    y <- as.numeric(x[, "feature1"] + rnorm(nrow(x)) > 2)
    train <- function(features) {
        xgboost::xgb.train(data = xgboost::xgb.DMatrix(x[, features], label = y), params = list(objective = "binary:logistic"), nrounds = 5, nthread = 1)
    }
    boosters <- list(train(c("feature1", "feature2")), train(c("feature2", "feature1")))
    reordered <- list("80" = list(model = boosters, score2quantile = predictors[["80"]]$score2quantile))
    scores <- mldp_predict_multi_age(new_data %>% dplyr::filter(age == 80), reordered)
    expected <- (predict(boosters[[1]], x[, c("feature1", "feature2")]) + predict(boosters[[2]], x[, c("feature2", "feature1")])) / 2
    expect_equal(scores$score, expected, tolerance = 1e-6)
})

test_that("exported predictors keep their cross validation models", {
    exported <- mldp_model_export(predictors)
    expect_equal(mldp_predict_multi_age(new_data, exported, use_cv = TRUE), mldp_predict_multi_age(new_data, predictors, use_cv = TRUE))
    expect_equal(mldp_predict_multi_age(new_data, exported), mldp_predict_multi_age(new_data, predictors))
})

test_that("use_cv fails for predictors without cross validation models", {
    single <- purrr::map(predictors, ~ .x[c("model", "score2quantile", "age")])
    expect_error(mldp_predict_multi_age(new_data, single, use_cv = TRUE), "no cross validation models")
})
