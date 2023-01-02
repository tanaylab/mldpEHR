
#' Load a mock dataset for mortality prediction
#'
#' @param N The number of patients to generate
#' @param num_age_groups The number of age groups to generate
#'
#' @return an MldpEHR object
#'
#' @examples
#' mortality_data <- load_mortality_example_data(1000)
#'
#' @export
load_mortality_example_data <- function(N = 1000, num_age_groups = 5) {
    patients <- purrr::map(0:(num_age_groups - 1), ~ data.frame(
        id = 1:N,
        sex = rep(c(1, 2), N / 2),
        age = 80 - .x * 5,
        death = c(rep(NA, round(0.2 * N)), rep(82, round(0.8 * N))),
        followup = .x * 5 + 5
    )) %>%
        setNames(seq(80, by = -5, length.out = num_age_groups))

    features <- purrr::map(0:(num_age_groups - 1), ~ data.frame(
        id = 1:N,
        feature1 = c(rnorm(round(0.2 * N)), rnorm(round(0.8 * N), mean = 2, sd = 0.5)),
        feature2 = c(rnorm(round(0.5 * N)), rnorm(round(0.5 * N), mean = 3, sd = 0.1))
    )) %>% setNames(seq(80, by = -5, length.out = num_age_groups))

    return(MldpEHR(patients = patients, features = features))
}


#' Load a mock dataset for disease prediction
#'
#' @param N The number of patients to generate
#' @param num_age_groups The number of age groups to generate
#'
#' @return an MldpEHR object
#'
#' @examples
#' disease_data <- load_disease_example_data(1000)
#'
#' @export
load_disease_example_data <- function(N = 1000, num_age_groups = 5) {
    patients <- purrr::map(0:(num_age_groups - 1), ~ data.frame(
        id = 1:N,
        sex = rep(c(1, 2), N / 2),
        age = 80 - .x * 5,
        death = c(rep(NA, round(0.4 * N)), rep(82, round(0.6 * N))),
        disease = sample(c(NA, 81), N, prob = c(0.4, 0.6), replace = TRUE),
        followup = .x * 5 + 5
    )) %>%
        setNames(seq(80, by = -5, length.out = num_age_groups))

    features <- purrr::map(0:(num_age_groups - 1), ~ data.frame(
        id = 1:N,
        feature1 = c(rnorm(round(0.2 * N)), rnorm(round(0.8 * N), mean = 2, sd = 0.5)),
        feature2 = c(rnorm(round(0.5 * N)), rnorm(round(0.5 * N), mean = 3, sd = 0.1))
    )) %>% setNames(seq(80, by = -5, length.out = num_age_groups))

    return(MldpEHR(patients = patients, features = features))
}
