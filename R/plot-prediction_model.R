#' Plot an ecdf of the score of each model for colored by outcome
#'
#' @description For the oldest age the outcome is death, for younger ages the outcome is the score quantile of the next age group.
#'
#' @param predictors a list of predictors, output of \code{mldp_mortality_multi_age_predictors}
#'
#' @return a ggplot object
#'
#' @examples
#' # Load a small example data
#' data <- load_mortality_example_data(N = 100, num_age_groups = 3)
#' predictors <- mldp_mortality_multi_age_predictors(
#'     data@patients,
#'     data@features,
#'     survival_years = 5,
#'     nfolds = 2,
#'     q_thresh = 0.2,
#'     nthread = 2 # CRAN allows only 2 cores
#' )
#' mldp_plot_multi_age_predictors_ecdf(predictors)
#'
#' @export
mldp_plot_multi_age_predictors_ecdf <- function(predictors) {
    x <- purrr::map2_dfr(predictors, names(predictors), ~ .x$test %>% mutate(model = .y))
    if ("disease" %in% colnames(x)) {
        x <- x %>%
            rename(outcome = step_outcome)
        colors <- c("death" = "black", "disease" = "orange", "disease_death" = "red", "healthy" = "gray")
    } else {
        x <- x %>%
            mutate(outcome = ifelse(!is.na(death), "outcome", "no-outcome")) %>%
            mutate(outcome = factor(outcome, levels = c("outcome", "no-outcome")))
        colors <- c("outcome" = "darkred", "no-outcome" = "gray")
    }
    x %>%
        ggplot2::ggplot(ggplot2::aes(x = predict, color = outcome)) +
        ggplot2::scale_color_manual(values = colors) +
        ggplot2::stat_ecdf() +
        ggplot2::facet_wrap(~model) +
        ggplot2::theme_bw() +
        ggplot2::xlab("Score") +
        ggplot2::ylab("Cumulative probability")
}

#' Plot SHAP (SHapley Additive exPlanations) values for a model
#'
#' @param predictor_shap output of \code{mldp_prediction_model_features}
#' @param n_patients the number of patients to plot. Default is to sample 1000 patients. If \code{n_patients} is greater than the number of patients in the data, all patients are plotted.
#' @param n_features the number of features to plot. Features will be sorted by the average absolute SHAP. If \code{NULL} all features are plotted.
#' @param point_size the size of the points in the plot. Default is 0.01.
#' @param point_alpha the transparency of the points in the plot. Default is 0.3.
#'
#' @return a ggplot object
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
#' predictor_shap <- mldp_model_features(predictors$`80`)
#' mldp_plot_shap(predictor_shap)
#'
#' @inheritDotParams ggplot2::geom_point
#' @export
mldp_plot_shap <- function(predictor_shap, n_patients = 1000, n_features = NULL, point_size = 0.01, point_alpha = 0.3, ...) {
    if (is.null(n_features)) {
        n_features <- nrow(predictor_shap$summary)
    }

    top_features <- predictor_shap$summary %>%
        arrange(desc(mean_abs_shap)) %>%
        slice(1:n_features) %>%
        pull(feature)

    features <- predictor_shap$shap_by_patient %>%
        filter(feature %in% top_features) %>%
        mutate(feature = factor(feature, levels = top_features))

    min_num_patients <- predictor_shap$shap_by_patient %>%
        count(feature) %>%
        pull(n) %>%
        min()

    if (n_patients > min_num_patients) {
        cli::cli_alert_info("{.field n_patients} ({.val {n_patients}}) is greater than the number of patients in the data ({.val {min_num_patients}}). All patients will be plotted.")
        n_patients <- min_num_patients
    }


    features <- features %>%
        group_by(feature) %>%
        sample_n(n_patients)

    p <- features %>%
        ggplot2::ggplot(ggplot2::aes(x = value, y = shap)) +
        ggplot2::geom_point(size = 0.01, alpha = 0.3, ...) +
        ggplot2::facet_wrap(~feature, nrow = 1, scales = "free_y") +
        ggplot2::theme_bw()

    return(p)
}
