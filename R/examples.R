
#' Load a mock dataset for mortality prediction
#'
#' @param N The number of patients to generate
#'
#' @return an MldpEHR object
#'
#' @examples
#' mortality_data <- load_mortality_data(1000)
#'
#' @export
load_mortality_example_data <- function(N = 1000) {
    patients <- purrr::map(0:5, ~ data.frame(
        id = 1:N,
        sex = rep(c(1, 2), N / 2),
        age = 80 - .x * 5,
        death = c(rep(NA, 0.2 * N), rep(82, 0.8 * N)),
        followup = .x * 5 + 5
    )) %>%
        setNames(seq(80, by = -5, length.out = 6))
    features <- purrr::map(0:5, ~ data.frame(
        id = 1:N,
        feature1 = c(rnorm(0.2 * N), rnorm(0.8 * N, mean = 2, sd = 0.5)),
        feature2 = rep(c(rnorm(N / 4), rnorm(N / 4, mean = 3)), 2)
    )) %>% setNames(seq(80, by = -5, length.out = 6))

    return(MldpEHR(patients = patients, features = features))
}
