globals <- unique(c(
    "age", "sbin", "status", "sex", "name", "time", "est", "outcome", ".", "step_outcome", "obsT", "s", "target_sbin", "qpredict", "disease", "followup", "death", "target_class", "fold", "disease_n", "value", "shap", "X1", "mean_abs_shap", "feature", "model_age", "prob", "qmax", "qmin", "score"
))
utils::globalVariables(globals)
