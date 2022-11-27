
#' To implement this, all patients at a given age will be binned according to their model score (using quantiles).
#' Each bin is assigned a state, and we are computing the probability for traversing from each state to the next model state
#' Patients with missing score are also included for this model to reflect actual population numbers
#' @param models - list of prediction models (output of mldpEHR.cv_train_stitch_outcome)
#' @param outcome - time from oldest model (first) to target outcome
#' @param step - time between prediction models
#' @param qbins - quantile bin size of prediction score for which the markov model will define a state
#' @return a list of markov models (per age), each is a list with the following members:
#' - model - matrix containing the probability for each quantile(score) bin to reach each of the target_classes provided in the oldest model.
#' - local.model - data.frame containing the probability for each quantile(score) bin to reach each of the quantile(score) bins of the next model by age.
#' - qbins -  bins
#' - target - data frame containing the target bin for this age model (to be used as outcome for the younger age model)
#' @examples
#'
#' library(dplyr)
#' library(ggplot2)
#' # build base predictor
#' outcome <- data.frame(id = 1:1000, sex = rep(0:1, 500), age=rep(80, 1000), death = rep(c(NA, 82), each = 500), disease=rep(rep(c(NA, 81), each=250),2), followup=rep(5, 1000))
#' patients <- c(list(outcome), purrr::map(1:4, ~ data.frame(
#'      id = 1:1000, 
#'      sex = rep(0:1, 500),
#'      age = 80 - .x * 5,
#'      death = rep(c(NA, 82), each = 500),
#'      disease=rep(rep(c(NA, 81), each=250),2),
#'      followup = .x*5+5))) %>% setNames(seq(80, by=-5, length.out=5))
#' features <- purrr::map(1:5, ~ data.frame(id = 1:1000, 
#'      a = c(rnorm(500),rnorm(500, mean = 2, sd = 1)) , 
#'      b = c(rnorm(500), rnorm(500, mean = -2, sd=1)),
#'      c = rep(c(rnorm(250), rnorm(250, mean=3)),2)
#'     )) %>% setNames(seq(80, by=-5, length.out=5))
#' predictors <- mldpEHR.disease_multi_age_predictors(patients, features, 5, 3)
#' markov <- mldpEHR.disease_markov(predictors, 5, 5, qbins = seq(0, 1, by = 0.25))
#' prob <- purrr::map2_df(markov, names(markov), ~
#'      as_tibble(.x$model[[1]], rownames = "sbin") %>%
#'      mutate(sex = 0, model = .y) %>%
#'      bind_rows(as_tibble(.x$model[[2]], rownames = "sbin") %>% mutate(sex = 1, model = .y))) %>%
#'      mutate(sbin = factor(sbin, levels = c("death", 1:4, "disease", "disease_death")))
#' ggplot(prob, aes(x = sbin, y = disease+disease_death, colour = factor(sex), group = factor(sex))) + 
#'      geom_point() + geom_line() + facet_wrap(~model, nrow = 1) + theme_bw()
#' @export


mldpEHR.disease_markov <- function(models, outcome, step, qbins = seq(0, 1, by = 0.05), required_conditions="id==id") {
    markov_models <- list()

    # first model is the oldest model, used to compute the actual risk
    markov_models[[1]] <- .disease_markov_model_for_outcome_model(models[[1]], outcome, qbins, required_conditions)
    for (i in 2:length(models)) {
        markov_models[[i]] <- .disease_markov_model_for_stitch_model(markov_models[[i - 1]], models[[i]], step, qbins, required_conditions)
    }
    names(markov_models) <- names(models)
    return(markov_models)
}


#' build an Markov probability model from multi-age prediction models
#' @param markov - the markov model computed for the next age (older). the states in this age will be mapped to the states in this markov layer
#' @param model - prediction model (output of build_cross_validation_time_stitch_classification_models)
#' @param follow -  data.frames  defining for each patient how much follow-up time they had. Dataframe must include the following columns:
#' - id
#' - time_in_system
#' - target_class - the classification if known within step (i.e. patient died). NA for unknown step outcome.
#' @param step - time between prediction models
#' @param qbins - quantile bin size of prediction score for which the markov model will define a state
#' @return a markov model, a list with the following members:
#' - model - matrix containing the probability for each quantile(score) bin to reach each of the target_classes provided in the oldest model.
#' - local.model - data.frame containing the probability for each quantile(score) bin to reach each of the quantile(score) bins of the next model by age.
#' - qbins -  bins
#' - target - data frame containing the target bin for this age model (to be used as outcome for the younger age model)

.disease_markov_model_for_stitch_model <- function(markov, model, step, qbins, required_conditions, min_obs_for_estimate = 10) {
    m <- model$test %>% # contains only patients with score (some were dropped)
        mutate(sbin = as.numeric(cut(qpredict, qbins, right = FALSE, include.lowest = TRUE))) %>%
        mutate(sbin = ifelse(!is.na(disease) & disease < age, "disease", sbin)) %>% 
        mutate(sbin = replace_na(sbin, "no_score")) %>%
        mutate(sbin = factor(sbin, levels = c(1:(length(qbins) - 1), "no_score", "disease")))
    pm <- m %>% 
        left_join(markov$target %>% select(id, target_sbin), by = "id") %>% 
        filter(!is.na(target_sbin) | step_outcome !="healthy", eval(rlang::parse_expr(required_conditions))) %>% 
        mutate(outcome=ifelse(is.na(target_sbin), step_outcome, as.character(target_sbin)))

    km_model <- plyr::ddply(pm, plyr::.(sbin, sex), function(data) {
        cencode <- 'healthy'
        if (data$sbin[1] == 'disease') {
            #setting censor outcome to disease as possible outcomes are disease or death
            cencode <- 'disease'
        }
        #message(data$sbin[1], " :: " , data$sex[1])
        if (nrow(data) < min_obs_for_estimate) {
            # bin is too small to compute probabilities, will use entire pop prob to fill in
            data <- pm %>% filter(sex == data$sex[1])
            warning("Insufficient stats for age:", data$age[1], " for sbin:", data$sbin[1], " for sex:", c("male", "female")[data$sex[1]], " ,using pop")
        }
        #checking if no events occur, just censoring (this will cause an error)
        if (sum(data$outcome == cencode) == nrow(data)) {
            return(data.frame(time=step, est=0, var=0, group=1, outcome="dead"))
        }

        fit <- cmprsk::cuminc(data$obsT, data$outcome, cencode = cencode)
        ret <- purrr::map2_df(fit, names(fit), ~ as.data.frame(.x[1:3]) %>% mutate(name = .y)) %>%
            tidyr::separate(name, into = c("group", "outcome"), sep = " (?=[^ ]*$)") %>%
            filter(time <= step) %>%
            arrange(desc(time), desc(est)) %>% # this is the estimate at latest measured time before outcome years.
            distinct(outcome, .keep_all = TRUE)
        return(ret)
    }) %>%
        select(sbin, sex, est, outcome) %>%
        arrange(sex, sbin, outcome)
    # adding estimates for censored data (patient remains sick)
    censored_estimate <- km_model %>%
        filter(sbin == "disease") %>% 
        group_by(sbin, sex) %>%
        summarize(est = 1 - sum(est), .groups = "drop") %>%
        mutate(outcome = "disease")
    # combining all outcome states
   
    local_model <- km_model %>%
        bind_rows(censored_estimate) %>% # change target to factor?
        #correcting for numerical errors
        group_by(sbin, sex) %>% mutate(s=sum(est)) %>% ungroup %>% mutate(est=est/s) %>% 
        select(-s) %>% 
        mutate(outcome=factor(outcome, levels=c(levels(markov$target$target_sbin), 'disease_death', 'death'))) %>% 
        arrange(outcome) %>% 
        pivot_wider(id_cols = c(sbin, sex), names_from = outcome, values_from = est) %>%
        arrange(sex, sbin)

    local_model_by_sex <- plyr::dlply(local_model, plyr::.(sex), function(s) {
        s %>%
            select(-sex) %>%
            bind_rows(data.frame(sbin = c("disease_death", "death"), disease_death=c(1,0), death = c(0, 1))) %>%
            mutate(sbin = factor(sbin, levels = c(1:(length(qbins) - 1), "no_score", "disease", "disease_death", "death"))) %>%
            replace(is.na(.), 0) %>%
            arrange(sbin) %>%
            column_to_rownames("sbin")
    })
    model_by_sex <- purrr::map2(local_model_by_sex, markov$model, ~ as.matrix(.x) %*% .y)

    return(list(
        model = model_by_sex,
        local_model = local_model_by_sex,
        qbins = qbins,
        target = m %>% mutate(target_sbin = sbin)
    ))
}

.disease_markov_model_for_outcome_model <- function(model, step, qbins, required_conditions, min_obs_for_estimate = 10) {
    m <- model$test %>% # contains only patients with score (some were dropped)
        mutate(sbin = as.numeric(cut(qpredict, qbins, right = FALSE, include.lowest = TRUE))) %>%
        mutate(sbin = ifelse(!is.na(disease) & disease < age, "disease", sbin)) %>% 
        mutate(sbin = replace_na(sbin, "no_score")) %>%
        mutate(sbin = factor(sbin, levels = c(1:(length(qbins) - 1), "no_score", "disease")))
    pm <- m %>% filter(!is.na(step_outcome), eval(rlang::parse_expr(required_conditions))) %>% 
        rename(outcome=step_outcome) %>% 
        filter(outcome != "healthy" | obsT >= step)

    # compute competing risk models for each sex indeptendently according to source bin(sbin)
    km_model <- plyr::ddply(pm, plyr::.(sbin, sex), function(data) {
        cencode <- 'healthy'
        #message(data$sbin[1], " :: " , data$sex[1])
        if (nrow(data) < min_obs_for_estimate) {
            # bin is too small to compute probabilities, will use entire pop prob to fill in
            data <- pm %>% filter(sex == data$sex[1])
            warning("Insufficient stats for age:", data$age[1], " for sbin:", data$sbin[1], " for sex:", c("male", "female")[data$sex[1]], " ,using pop")
        }
        if (data$sbin[1] == 'disease') {
            #setting censor outcome to disease as possible outcomes are disease or death
            cencode <- 'disease'
        }
        #checking if no events occur, just censoring (this will cause an error)
        if (sum(data$outcome == cencode) == nrow(data)) {
            return(data.frame(time=step, est=0, var=0, group=1, outcome="dead"))
        }

        fit <- cmprsk::cuminc(data$obsT, data$outcome, cencode = cencode)
        ret <- purrr::map2_df(fit, names(fit), ~ as.data.frame(.x[1:3]) %>% mutate(name = .y)) %>%
            tidyr::separate(name, into = c("group", "outcome"), sep = " (?=[^ ]*$)") %>%
            filter(time <= step) %>%
            arrange(desc(time), desc(est)) %>% # this is the estimate at latest measured time before outcome years.
            distinct(outcome, .keep_all = TRUE)
        return(ret)
    }) %>%
        select(sbin, sex, est, outcome) %>%
        arrange(sex, sbin, outcome)
    # adding estimates for censored data (patient is healthy)
    censored_estimate <- km_model %>%
        group_by(sbin, sex) %>%
        summarize(est = 1 - sum(est), .groups = "drop") %>%
        mutate(outcome = ifelse(sbin=="disease", "disease", "healthy"))
    # combining all outcome states
   
    local_model <- km_model %>%
        bind_rows(censored_estimate) %>% # change target to factor?
        #correcting for numerical errors
        group_by(sbin, sex) %>% mutate(s=sum(est)) %>% ungroup %>% mutate(est=est/s) %>% 
        select(-s) %>% 
        mutate(outcome=factor(outcome, levels=c('disease', 'disease_death', 'death', 'healthy'))) %>% 
        arrange(outcome) %>% 
        pivot_wider(id_cols = c(sbin, sex), names_from = outcome, values_from = est) %>%
        arrange(sex, sbin)

    local_model_by_sex <- plyr::dlply(local_model, plyr::.(sex), function(s) {
        s %>%
            select(-sex) %>%
            bind_rows(data.frame(sbin = c("disease_death", "death"), disease_death=c(1,0), death = c(0, 1))) %>%
            mutate(sbin = factor(sbin, levels = c(1:(length(qbins) - 1), "no_score", "disease", "disease_death", "death"))) %>%
            replace(is.na(.), 0) %>%
            arrange(sbin) %>%
            column_to_rownames("sbin")
    })
    model_by_sex <- purrr::map(local_model_by_sex, ~ as.matrix(.x))

    return(list(
        model = model_by_sex,
        local_model = local_model_by_sex,
        qbins = qbins,
        target = m %>% mutate(target_sbin = sbin)
    ))
}



#' this method computes the probability for future disease according to the availablitiy of data
#' by basically computing a stitched markov model
#' as score availability is not unbiased (patients with less measurements tend to be healthier), 
#' the probaility will be computed seperately for patients with / without score
#' we will define the following modes:
#' - 0 has required conditions at young age and known score / outcome at older age
#' - 1 has required conditions at young age and no restrictions at older age
#' - 2 no restrictions at young/older age
#' @param population - list of data.frames of all the patients in the system going back in time. For example the first data.frame represents age 80, next is 75 and so forth.
#' Each patient data.frame contains the following columns:
#' - patient id
#' - sex
#' - target_class (does this patient have a known outcome)
#' - disease age of disease
#' - death age of death
#' - followup
#' - any other columns that can be used in required conditions
#' @param step - time between populations
#' @param required_conditions - conditions that will be applied on patients that will be have a prediction score
.mldpEHR.disease_empirical_prob_for_disease <- function(population, step, required_conditions="id==id") {
	km <- list()
	for (i in 1:length(population)) {
		k <- plyr::adply(0:2, 1, function(state) {
			if (state == 2) {
				young_pop <- population[[i]] %>% filter(is.na(disease) | disease > age)#?do we need ot filter out sick?
			} else {
				young_pop <- population[[i]] %>% filter(is.na(disease) | disease > age, eval(rlang::parse_expr(required_conditions))) 
				if (state == 0) {
					if (i == 1) {
						old_pop <- population[[i]] %>% filter(!is.na(target_class)) %>% pull(id)
					} else {
						old_pop <- c(population[[i]] %>% filter(!is.na(target_class)) %>% pull(id),
							population[[i-1]] %>% filter(eval(rlang::parse_expr(required_conditions))) %>% pull(id))
					}
					young_pop <- young_pop %>% filter(id %in% old_pop)
				}
			}
            young_pop<- young_pop %>% mutate(status=factor(ifelse(step_outcome=="disease_death", "disease", step_outcome), 
                levels=c("healthy", "disease", "death")))
    		 
			km <- plyr::ddply(young_pop, plyr::.(sex), function(data) { 
                if (nrow(data %>% filter(status != "healthy")) == 0) {
                    return(data.frame(time=step, est=0, var=0, group=1, outcome="disease"))
                }
				fit <- cmprsk::cuminc(data$obsT, data$status, cencode = "healthy")
        		ret <- purrr::map2_df(fit, names(fit), ~ as.data.frame(.x[1:3]) %>% mutate(name = .y)) %>%
            		tidyr::separate(name, into = c("group", "outcome"), sep = " (?=[^ ]*$)") %>%
            		filter(time <= step) %>%
            		arrange(desc(time), desc(est)) %>% # this is the estimate at latest measured time before outcome years.
            		distinct(outcome, .keep_all = TRUE)
        		return(ret)
			}) 
			return(km %>% mutate(mode=state))
			
		})
		km[[i]] <- k
	}
	names(km) <- names(population)
	return(list(prob=km, required_conditions=required_conditions))
}






#' calculate expected number of disease patients
#' @param population_count - data.frame containing age, sex, and the number of patients available
#' @param empirical_disease_prob - output of .mldpEHR.disease_empirical_prob_for_disease
#' @param step - time between populations
.mldpEHR.disease_expected <- function(index, population_count, edp) {
	plyr::adply(population_count, 1, function(pc) {
		pc %>% mutate(disease_n=.mldpEHR.disease_expected_sex(
			pc$sex, 
			pc$n, 
			edp[1:index],
			mode=0))
		})
}

.mldpEHR.disease_expected_sex <- function(sex, n, disease_prob, mode) {
	if (n == 0 | length(disease_prob) == 0) {
		return(0)
	}
	disease_count= round(n * disease_prob[[length(disease_prob)]] %>% filter(sex == !!sex, mode == !!mode, outcome=="disease") %>% pull(est))
	if (identical(disease_count, numeric(0))) {
			disease_count <- 0
	}

	#message(paste(age, " : ", size, " : ", round(disease_count)))
	next_n <- round(n * (1 - (disease_prob[[length(disease_prob)]] %>% filter(sex == !!sex, mode==!!mode) %>% pull(est) %>% sum)) )
	return(disease_count + .mldpEHR.disease_expected_sex(sex, next_n, head(disease_prob, -1), min(mode+1, 2)))
}


.mldpEHR.disease_assign_expected <- function(index, ordered_patients, empirical_disease_prob) 
{
	patients_filtered <- ordered_patients %>% filter(eval(rlang::parse_expr(empirical_disease_prob$required_conditions))) %>% 
		mutate(target_class = 0)
	expected <- .mldpEHR.disease_expected(index, patients_filtered %>% count(sex), empirical_disease_prob$prob)
	patients_filtered_class <- plyr::adply(expected,1, function(a) {
		top_res <- a$disease_n
		ret <- patients_filtered %>% inner_join(a %>% select(-n, -disease_n))
		ret$target_class[1:top_res] <- 1
		return(ret)
	}) %>% select(-n, -disease_n)
	return(patients_filtered_class %>% bind_rows(ordered_patients %>% anti_join(patients_filtered %>% select(id))))
}