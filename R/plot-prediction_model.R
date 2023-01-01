#' Plot an ecdf of the score of each model for colored by outcome
#'
#' @description For the oldest age the outcome is death, for younger ages the outcome is the score quantile of the next age group.
#'
#' @param predictors a list of predictors, output of \code{mldp_mortality_multi_age_predictors}
#'
#' @return a ggplot object
#'
#' @examples
#' # Load the example data
#' data <- load_mortality_example_data()
#' predictors <- mldp_mortality_multi_age_predictors(mortality@patients, mortality@features, survival_years = 5, nfolds = 3, q_thresh = 0.2)
#' mldp_plot_multi_age_predictors_ecdf(predictors)
#'
#' @export
mldp_plot_multi_age_predictors_ecdf <- function(predictors) {
    purrr::map2_dfr(predictors, names(predictors), ~ .x$test %>% mutate(model = .y)) %>%
        mutate(outcome = ifelse(!is.na(death), "outcome", "no-outcome")) %>%
        mutate(outcome = factor(outcome, levels = c("outcome", "no-outcome"))) %>%
        ggplot2::ggplot(ggplot2::aes(x = predict, color = outcome)) +
        ggplot2::scale_color_manual(values = c("darkred", "grey")) +
        ggplot2::stat_ecdf() +
        ggplot2::facet_wrap(~model) +
        ggplot2::theme_bw() +
        ggplot2::xlab("Score") +
        ggplot2::ylab("Cumulative probability")
}

#' Plot SHAP (SHapley Additive exPlanations) values for a model
#'
#' @param predictor_shap output of \code{mldp_prediction_model_features}
#' @param n_features the number of features to plot
#' @param n_patients the number of patients to plot
#'
#' @return a ggplot object
#'
#' @examples
#' # Load the example data
#' data <- load_mortality_example_data()
#' predictors <- mldp_mortality_multi_age_predictors(mortality@patients, mortality@features, survival_years = 5, nfolds = 3, q_thresh = 0.2)
#' predictor_shap <- mldp_prediction_model_features(predictors$`80`)
#' mldp_plot_prediction_model_features(predictor_shap)
#'
#' @export
mldp_plot_multi_age_predictors_shap <- function(predictor_shap) {
    browser()
}
