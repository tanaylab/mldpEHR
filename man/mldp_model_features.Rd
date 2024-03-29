% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prediction_model.R
\name{mldp_model_features}
\alias{mldp_model_features}
\title{Extract statistics on the features used in a k-fold cross-validation model using SHAP (SHapley Additive exPlanations)}
\usage{
mldp_model_features(predictor)
}
\arguments{
\item{predictor}{a single predictor, a single element of the list returned by \code{mldp_disease_multi_age_predictors} or \code{mldp_mortality_multi_age_predictors}}
}
\value{
a list with 3 elements:
\itemize{
\item{summary: }{a data frame containing for each feature the mean absolute Shapley value for the feature across all training data.}
\item{shap_by_patient: }{a data frame data frame containing for each patient and feature the feature value and mean shap value across all training folds.}
\item{shap_by_fold: }{similar to shap_by_patient, but for each fold separately.}
}
}
\description{
Extract statistics on the features used in a k-fold cross-validation model using SHAP (SHapley Additive exPlanations)
}
\examples{
# Load a small example data
mortality <- load_mortality_example_data(N = 100, num_age_groups = 3)

#' Build predictors
predictors <- mldp_mortality_multi_age_predictors(
    mortality@patients,
    mortality@features,
    survival_years = 5,
    nfolds = 2,
    q_thresh = 0.05,
    nthread = 2 # CRAN allows only 2 cores
)

predictor_features <- mldp_model_features(predictors[[1]])

names(predictor_features)
head(predictor_features$summary)
head(predictor_features$shap_by_patient)
head(predictor_features$shap_by_fold)

}
