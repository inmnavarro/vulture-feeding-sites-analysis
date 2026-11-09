# Extract final posterior SD estimates without refitting any Bayesian model.
# Run in the analysis folder containing the corrected .rds files.
suppressPackageStartupMessages(library(brms))
project_dir <- Sys.getenv("SFS_DATA_DIR", unset=".")
model_dir <- Sys.getenv("SFS_MODEL_DIR", unset=file.path(project_dir, "models"))
fits <- c(
  'Facultative-scavenger detection' = 'modelo_riesgo_bayes_revisado.rds',
  'Site safety (fox based)' = 'modelo_seguridad_bayes.rds',
  'First facultative arrival' = 'modelo_primero_bayes.rds',
  'Time to first vulture' = 'modelo_tiempo_buitres_bayes_revisado.rds',
  'Vulture attendance' = 'modelo_abundancia_bayes_revisado.rds'
)
missing <- fits[!file.exists(file.path(model_dir, fits))]
if (length(missing)) stop('Missing these final models: ', paste(missing, collapse=', '))
summary_param <- function(x, parameter, model) {
  if (!parameter %in% colnames(x)) stop('Missing parameter ',parameter,' in model ',model)
  draws <- x[,parameter]
  data.frame(Model=model, Parameter=parameter,
             Mean=mean(draws),
             Lower95=unname(quantile(draws, .025)),
             Upper95=unname(quantile(draws, .975)))
}
rows <- list()
for (i in seq_along(fits)) {
  modelname <- names(fits)[i]
  fit <- readRDS(file.path(model_dir, fits[i]))
  x <- as.matrix(fit)
  rows[[length(rows)+1]] <- summary_param(x,'sd_ID_COTO__Intercept',modelname)
  rows[[length(rows)+1]] <- summary_param(x,'sds_sfecha_sc_1',modelname)
  rm(fit,x);gc()
}
out <- do.call(rbind,rows)
out[c('Mean','Lower95','Upper95')] <- lapply(out[c('Mean','Lower95','Upper95')],round,3)
file <- file.path(project_dir,'Supplementary_Table_S1_FINAL_values.csv')
write.csv(out,file,row.names=FALSE)
print(out,row.names=FALSE)
message('Created: ', file)
