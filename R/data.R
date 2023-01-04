#' Download a single dataset
#'
#' @inheritParams mldp_download_example_data
#' @noRd
download_dataset <- function(output_dir, dataset, url, overwrite, timeout = 60 * 60 * 2) {
    # set the 'timeout' option to 2 hours
    withr::with_options(list(timeout = timeout), {
        for (type in c("patients", "features")) {
            file <- paste0(dataset, ".", type, ".rds")
            file_path <- file.path(output_dir, file)
            if (!file.exists(file_path) || overwrite) {
                cli::cli_alert_info("Downloading {.file {file}}")
                utils::download.file(paste0(url, "/", file), file_path)
            } else {
                cli::cli_alert_info("Skipping {.file {file}} (set {.code overwrite = TRUE} to overwrite)")
            }
        }
    })
}

#' Download example data
#'
#' @description Download example data. See description of the datasets below.
#'
#' @section longevity.patients:
#' A dataset for patient survival.
#'
#' This dataset is a list of 9 data frames, representing the population of patients for
#' each age starting at 80 and going down to 40 at 5 year interavls. For each patient,
#' the age at death is provided.
#' \cr
#'
#' Each data frame has the following 5 columns:
#' \itemize{
#'   \item{id: }{unique patient id}
#'   \item{age: }{age of patient}
#'   \item{sex: }{1 for male, 2 for female}
#'   \item{death: }{age at death, NA if patient doesn't die by end of followup}
#'   \item{has_cbc: }{binary indicating if patient has required feature data}
#'   \item{followup: }{available followup time in years}
#' }
#' \cr
#'
#' @section longevity.features:
#' A dataset for patient survival features.
#'
#' This dataset is a list of 9 data frames, representing the features of patients available in longevity.patients.
#' Similar to longevity.patients, each data frame is for a given age population,
#' starting at 80 and going down to 40 at 5 year interavls
#' \cr
#' Each data frame has 86 columns. The features include sex and common labs, quantile normalized values.
#'
#' @section diabetes.patients:
#' A dataset for patient disease outcome.
#'
#' This dataset is a list of 9 data frames, representing the population of patients for
#' each age starting at 80 and going down to 40 at 5 year interavls. For each patient,
#' the age at disease is provided.
#' \cr
#' Each dataframe with the following 5 columns:
#' \itemize{
#'   \item{id: }{unique patient id}
#'   \item{age: }{age of patient}
#'   \item{sex: }{1 for male, 2 for female}
#'   \item{death: }{age at death, NA if patient doesn't die by end of follow-up}
#'   \item{disease: }{age at diabetes disease, NA if patient doesn't get sick by end of follow-up}
#'   \item{has_cbc: }{binary indicating if patient has available data}
#'   \item{followup: }{available follow-up time in years}
#' }
#' \cr
#'
#' @section diabetes.features:
#' A dataset for patient diabetes features.
#'
#' This dataset is a list of 9 data frames, representing the features of patients available in diabetes.patients.
#' Similar to diabetes.patients, each data frame is for a given age population,
#' starting at 80 and going down to 40 at 5 year interavls.
#' \cr
#' Each data frame has 86 columns. The features include sex and common labs, quantile normalized values.
#'
#'
#' @param output_dir Directory to save the data to. Defaults to "examples" in the current working directory.
#' @param datasets A character vector of datasets to download. Defaults to "longevity" and "diabetes".
#' @param overwrite Overwrite existing files. Defaults to FALSE.
#' @param timeout Timeout in seconds for downloading files. Defaults to 2 hours. If you are having trouble downloading the files, try increasing this value.
#'
#' @return Nothing
#'
#' @examples
#' \donttest{
#' mldp_download_example_data()
#' }
#'
#' \dontshow{
#' unlink("examples", recursive = TRUE)
#' }
#'
#' @export
mldp_download_example_data <- function(output_dir = file.path(getwd(), "examples"), datasets = c("longevity", "diabetes"), overwrite = FALSE, timeout = 60 * 60 * 2) {
    datasets_opt <- list(
        "longevity" = "https://mldp-ehr.s3.eu-west-1.amazonaws.com/example",
        "diabetes" = "https://mldp-ehr.s3.eu-west-1.amazonaws.com/example"
    )

    # create output_dir if it doesn't exist
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    }


    for (dataset in names(datasets)) {
        if (!(dataset %in% names(datasets_opt))) {
            cli::cli_abort("Dataset must be one of {.field {names(datasets_opt)}}, not {.field {dataset}}")
        }
        download_dataset(output_dir, dataset, datasets_opt[[dataset]], overwrite)
    }

    cli::cli_alert_success("Successfully downloaded data to {.file {output_dir}}. Use {.code mldpEHR::load_data()} to load the data.")
}

#' Load example data
#'
#' @description Load example data that was downloaded with \code{mldpEHR::download_example_data()}. See description of the datasets below.
#'
#' @param dataset A character vector of the dataset to load, e.g. "longevity" or "diabetes".
#' @param data_dir Directory to load the data from. Defaults to "examples" in the current working directory.
#' You can download the data using \code{mldpEHR::download_example_data()}.
#'
#' @return an MldpEHR object
#'
#' @examples
#' \donttest{
#' mldp_download_example_data(file.path(getwd(), "examples"))
#' longevity_data <- mldp_load_data("longevity", file.path(getwd(), "examples"))
#' diabetes_data <- mldp_load_data("diabetes", file.path(getwd(), "examples"))
#' }
#'
#' \dontshow{
#' unlink("examples", recursive = TRUE)
#' }
#'
#' @inheritSection mldp_download_example_data longevity.patients
#' @inheritSection mldp_download_example_data longevity.features
#' @inheritSection mldp_download_example_data diabetes.patients
#' @inheritSection mldp_download_example_data diabetes.features
#'
#' @export
mldp_load_data <- function(dataset, data_dir = file.path(getwd(), "examples")) {
    if (!dir.exists(data_dir)) {
        cli::cli_abort("Data directory {.file {data_dir}} does not exist. Use {.code mldp_download_example_data()} to download the data.")
    }

    patients <- readRDS(file.path(data_dir, paste0(dataset, ".patients.rds")))
    features <- readRDS(file.path(data_dir, paste0(dataset, ".features.rds")))

    return(MldpEHR(patients = patients, features = features))
}
