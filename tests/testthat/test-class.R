test_that("MldpEHR class works", {
    # Create a MldpEHR object
    patients <- list(
        tibble::tibble(
            id = 1:10,
            age = 80,
            sex = sample(1:2, 10, replace = TRUE),
            followup = sample(0:5, 10, replace = TRUE),
            death = pmin(age + followup, sample(c(NA, 80:85), 10, replace = TRUE))
        ),
        tibble::tibble(
            id = 1:10,
            age = 75,
            sex = sample(1:2, 10, replace = TRUE),
            followup = sample(0:5, 10, replace = TRUE),
            death = pmin(age + followup, sample(c(NA, 75:80), 10, replace = TRUE))
        )
    )

    names(patients) <- c("80", "75")

    features <- list(
        tibble::tibble(
            id = 1:10,
            feature1 = rnorm(10),
            feature2 = rnorm(10),
            feature3 = sample(1:3, 10, replace = TRUE)
        ),
        tibble::tibble(
            id = 1:10,
            feature1 = rnorm(10),
            feature2 = rnorm(10),
            feature4 = sample(0:1, 10, replace = TRUE)
        )
    )

    mldp <- MldpEHR(patients = patients, features = features)

    expect_equal(mldp@patients, patients)
    expect_equal(mldp@features, features)
})
