# Declare non-standard-evaluation names so R CMD check stays quiet.
utils::globalVariables(c("HR", "HR95L", "HR95H", "lower", "upper", "Pvalue",
                         "col_group", "group", "time", "status", "S", "SE",
                         "Symbol", "FPR", "TPR", "cohort", "table", "row",
                         "label", "marker", "left_out", "order"))
utils::globalVariables(c("time", "cif", "lower", "upper", "cause_lab"))
utils::globalVariables(c("sig", "label", "k"))
utils::globalVariables(c("dataset", "cif", "cause_lab"))
# B5 batch diagnostics: the ggplot2 aesthetics of plot_batch_diagnostics()
utils::globalVariables(c("value", "cut_median", "cut_top_pct", "gene"))
