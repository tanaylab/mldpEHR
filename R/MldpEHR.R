#' Data for mldpEHR models
#'
#' @slot patients list of data.frames of all the patients in the system going back in time. For example the first data.frame represents age 80, next is 75 and so forth. See the creation function for more details.
#' @slot features a list of data.frames of all the features for the patients in the system going back in time. For example the first data.frame represents age 80, next is 75 and so forth. See the creation function for more details.
#' @slot age_groups (optional) labels for the age groups
#'
#'
#' @exportClass MldpEHR
setClass(
    "MldpEHR",
    representation(patients = "list", features = "list", age_groups = "vector")
)

#' Create MldpEHR object
#'
#' @param patients list of data.frames of all the patients in the system going back in time. For example the first data.frame represents age 80, next is 75 and so forth. Each patient data.frame contains the following columns:
#' \itemize{
#'   \item{id: }{unique patient id}
#'   \item{age: }{age of patient. All patients in the data.frame must be the same age}
#'   \item{sex: }{1 for male, 2 for female}
#'   \item{death: }{age at death, NA if patient doesn't die by end of followup}
#'   \item{followup: }{available followup time (in years) for this patient - time until end of database or until patient exists the system (not due to death)}
#' }\cr
#'
#' The data frame can contain any additional columns required for patient filtering in the future.\cr
#'
#' Note that every age group except the last one must have at least two patients per sex which exist in the next age group. This is to ensure that the model can be trained.
#'
#' @param features list of data.frames of all the features for the patients in the system going back in time. For example the first data.frame represents age 80, next is 75 and so forth. Each feature data.frame must contain an id column that matches the id column in the patient data.frame. The feature data.frame can contain any additional feature columns.
#' @param age_groups (optional) labels for the age groups. If the \code{patients} list is named, the age_groups will be the names of the lists.
#'
#' @return a MldpEHR object
#'
#' @examples
#' # Create a MldpEHR object
#' patients <- list(
#'     tibble::tibble(
#'         id = 1:10,
#'         age = 80,
#'         sex = sample(1:2, 10, replace = TRUE),
#'         followup = sample(0:5, 10, replace = TRUE),
#'         death = pmin(age + followup, sample(c(NA, 80:85), 10, replace = TRUE))
#'     ),
#'     tibble::tibble(
#'         id = 1:10,
#'         age = 75,
#'         sex = sample(1:2, 10, replace = TRUE),
#'         followup = sample(0:5, 10, replace = TRUE),
#'         death = pmin(age + followup, sample(c(NA, 75:80), 10, replace = TRUE))
#'     )
#' )
#'
#' names(patients) <- c("80", "75")
#'
#' features <- list(
#'     tibble::tibble(
#'         id = 1:10,
#'         feature1 = rnorm(10),
#'         feature2 = rnorm(10),
#'         feature3 = sample(1:3, 10, replace = TRUE)
#'     ),
#'     tibble::tibble(
#'         id = 1:10,
#'         feature1 = rnorm(10),
#'         feature2 = rnorm(10),
#'         feature4 = sample(0:1, 10, replace = TRUE)
#'     )
#' )
#'
#' mldp <- MldpEHR(patients = patients, features = features)
#'
#' @rdname MldpEHR
#' @export
MldpEHR <- function(patients, features, age_groups = NULL) {
    if (!all(sapply(patients, function(x) inherits(x, "data.frame")))) {
        cli::cli_abort("All patients must be data.frames")
    }

    required_columns <- c("id", "sex", "age", "death", "followup")

    # check that all the dataframes have the required columns
    purrr::walk(patients, function(x) {
        if (!all(required_columns %in% colnames(x))) {
            cli::cli_abort("All patients must have the following columns: {.field {required_columns}}")
        }
    })

    # validate that "sex" column has only 1 and 2
    purrr::walk(patients, function(x) {
        if (!all(x$sex %in% c(1, 2))) {
            cli::cli_abort("{.field sex} column must only have 1 for 'male' and 2 for 'female")
        }
    })

    # validate that "death" column is NA if it is more than age+followup
    purrr::walk(patients, function(x) {
        if (!all(is.na(x$death) | x$death <= x$age + x$followup)) {
            cli::cli_abort("{.field death} column must be NA if it is more than {.field age} + {.field followup}")
        }
    })

    # validate that the length of the patients and features are the same
    if (length(patients) != length(features)) {
        cli::cli_abort("The number of age groups of {.field patients} and {.field features} must be the same")
    }

    # make sure that all features data frames have an 'id' column
    purrr::walk(features, function(x) {
        if (!"id" %in% colnames(x)) {
            cli::cli_abort("All features must have an {.field id} column")
        }
    })

    purrr::walk2(patients, features, function(x, y) {
        if (!all(y$id %in% x$id)) {
            cli::cli_abort("All ids in {.field features} must be in {.field patients}")
        }
    })

    # make sure that all features are numeric
    purrr::walk(features, function(x) {
        if (!all(sapply(select(x, -id), is.numeric))) {
            cli::cli_abort("All features must be numeric")
        }
    })

    if (!is.null(names(patients))) {
        age_groups <- names(patients)
        if (!is.null(names(features))) {
            if (!all(names(patients) == names(features))) {
                cli::cli_abort("The names of {.field patients} and {.field features} must be the same")
            }
        }
    }

    if (!is.null(age_groups)) {
        if (length(age_groups) != length(patients)) {
            cli::cli_abort("The number of age groups of {.field age_groups} and {.field patients} must be the same")
        }
    } else {
        age_groups <- seq_along(patients)
    }

    # make sure that each age group has at least 2 patients per gender that exists in the next age group
    for (i in 1:(length(patients) - 1)) {
        x <- patients[[i]] %>% filter(id %in% patients[[i + 1]]$id)

        if (nrow(x) < 2) {
            cli::cli_abort("Each age group must have at least 2 patients that exist in the next age group. Age group {.field {age_groups[i]}} does not meet this requirement")
        }

        x <- x %>%
            count(sex)

        if (any(x$n < 2)) {
            cli::cli_abort("Each age group must have at least 2 patients from each sex that exists in the next age group. Age group {.field {age_groups[i]}} does not meet this requirement")
        }
    }

    # make sure that each age group as only a single age
    purrr::walk(patients, function(x) {
        if (length(unique(x$age)) > 1) {
            cli::cli_abort("Each age group must have only a single age. Age group {.field {age_groups[i]}} does not meet this requirement")
        }
    })

    # make sure that age groups are decreasing monotonically
    if (!all(diff(purrr::map_dbl(patients, ~ .x$age[1])) < 0)) {
        cli::cli_abort("Age groups must be decreasing monotonically")
    }

    new("MldpEHR", patients = patients, features = features, age_groups = age_groups)
}

#' @export
#' @noRd
setMethod(
    "show",
    signature = "MldpEHR",
    definition = function(object) {
        total_patients <- sum(sapply(object@patients, nrow))
        cli::cli_text("MldpEHR object with {.val {length(object@patients)}} age groups and {.val {total_patients}} patients")
        if (!is.null(object@age_groups)) {
            cli::cli_text("Age groups: {.val {object@age_groups}}")
        }
        cli::cli_text("Slots are: {.val {slotNames(object)}}")
    }
)


#' Compute the age difference between age groups in a MldpEHR object
#'
#' @param mldp A MldpEHR object
#'
#' @noRd
mldp_get_age_steps <- function(mldp) {
    abs(diff(purrr::map_dbl(mldp@patients, ~ .x$age[1])))
}
