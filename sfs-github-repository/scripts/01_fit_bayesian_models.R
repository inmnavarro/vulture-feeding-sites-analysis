# Bayesian models for supplementary feeding site study (central Spain)
#
# Analysis reconstruction from the original project script.
# This version documents the three convergence-control adjustments made for
# the final models (risk, vulture arrival time, vulture attendance).
# Model 2 "site safety" is FOX-BASED in the original fitted analysis.
# Outputs are local brms objects: NEVER commit RDS files (they can embed data).
# No raw data are distributed in this code repository.
#
# Usage:
#   SFS_DATA_DIR=/path/to/private/data Rscript scripts/01_fit_bayesian_models.R
#
# Required source CSV name: BASE_DATOS_BUITRES.csv (semicolon-separated, Latin1).
# Reproducing results requires the non-public curated data and package versions.
# The original analysis was not a randomized experiment.

needed <- c("readr", "dplyr", "lubridate", "brms", "ggplot2", "bayesplot")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing)) stop("Install R packages: ", paste(missing, collapse=", "))
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(lubridate)
  library(brms)
  library(ggplot2)
  library(bayesplot)
})

project_dir <- Sys.getenv("SFS_DATA_DIR", unset=".")
input_csv <- file.path(project_dir, "BASE_DATOS_BUITRES.csv")
if (!file.exists(input_csv)) stop("Private data not found: ", input_csv)
out_dir <- file.path(project_dir, "models")
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
options(mc.cores = min(parallel::detectCores(), 4L))

# IMPORTAR DATOS
###############################################################

datos <- read_csv2(
  input_csv,
  locale = locale(encoding = "Latin1")
)

datos <- datos %>%
  filter(!is.na(ID_COTO))

###############################################################
# FECHAS
###############################################################

datos$Fecha_aporte <- dmy(datos$Fecha_aporte)
datos$Fecha_evento <- dmy(datos$Fecha_evento)

###############################################################
# FACTORES
###############################################################

factores <- c(
  "Tipo_recinto",
  "Coto_tipo",
  "PRES_GANADO",
  "Valla_acondicionada",
  "Clasificación",
  "Especie_detectada",
  "behav",
  "ID_COTO",
  "ID_aporte"
)

datos[factores] <- lapply(
  datos[factores],
  as.factor
)

###############################################################
# VARIABLES TEMPORALES
###############################################################

datos <- datos %>%
  mutate(
    
    fecha_num =
      as.numeric(
        Fecha_aporte -
          min(Fecha_aporte, na.rm = TRUE)
      ),
    
    mes =
      month(Fecha_aporte),
    
    juliano =
      yday(Fecha_aporte)
  )

###############################################################
# ESCALADO
###############################################################

datos <- datos %>%
  mutate(
    
    altura_sc =
      scale(Altura_cercado_cm)[,1],
    
    peso_sc =
      scale(Peso_estimado_kg)[,1],
    
    temp_sc =
      scale(Temperatura_evento)[,1],
    
    fecha_sc =
      scale(fecha_num)[,1]
  )

###############################################################
# VARIABLES RESPUESTA
###############################################################

datos <- datos %>%
  mutate(
    
    vector_riesgo =
      ifelse(
        Clasificación == "Oportunista",
        1,
        0
      ),
    
    consumo =
      ifelse(
        behav == "Consumo",
        1,
        0
      ),
    
    zorro =
      ifelse(
        Especie_detectada == "Zorro",
        1,
        0
      ),
    
    buitre =
      ifelse(
        Especie_detectada == "Buitres",
        1,
        0
      )
  )

###############################################################
# SOLO EVENTOS POST APORTE
###############################################################

datos_aporte <- datos %>%
  filter(
    POST_aporte == "SI"
  )

###############################################################
# PRIMERA ESPECIE
###############################################################

primeros <- datos_aporte %>%
  filter(
    !is.na(Tiempo_desde_aporte_min)
  ) %>%
  group_by(ID_aporte) %>%
  slice_min(
    Tiempo_desde_aporte_min,
    with_ties = FALSE
  ) %>%
  ungroup()

primeros <- primeros %>%
  mutate(
    primero_vector =
      ifelse(
        Clasificación == "Oportunista",
        1,
        0
      )
  )

###############################################################
# TIEMPOS DE LLEGADA
###############################################################

tiempo_buitres <- datos_aporte %>%
  filter(
    Especie_detectada == "Buitres",
    !is.na(Tiempo_desde_aporte_min)
  ) %>%
  group_by(ID_aporte) %>%
  summarise(
    tiempo_buitres =
      min(Tiempo_desde_aporte_min)
  )

tiempo_zorro <- datos_aporte %>%
  filter(
    Especie_detectada == "Zorro",
    !is.na(Tiempo_desde_aporte_min)
  ) %>%
  group_by(ID_aporte) %>%
  summarise(
    tiempo_zorro =
      min(Tiempo_desde_aporte_min)
  )

###############################################################
# TABLA DE APORTES
###############################################################

aportes <- datos_aporte %>%
  distinct(
    ID_aporte,
    .keep_all = TRUE
  ) %>%
  left_join(
    tiempo_buitres,
    by = "ID_aporte"
  ) %>%
  left_join(
    tiempo_zorro,
    by = "ID_aporte"
  )

aportes <- aportes %>%
  mutate(
    
    zorro_antes =
      case_when(
        
        is.na(tiempo_zorro) ~ 0,
        
        is.na(tiempo_buitres) &
          !is.na(tiempo_zorro) ~ 1,
        
        tiempo_zorro <
          tiempo_buitres ~ 1,
        
        TRUE ~ 0
      ),
    
    seguro =
      ifelse(
        zorro_antes == 1,
        0,
        1
      )
  )

###############################################################
# PRIORS
###############################################################

priors <- c(
  
  prior(
    normal(0,2),
    class = "b"
  ),
  
  prior(
    student_t(3,0,2.5),
    class = "Intercept"
  ),
  
  prior(
    exponential(1),
    class = "sd"
  )
)

###############################################################
# MODELO 1
# RIESGO SANITARIO
###############################################################

modelo_riesgo <- brm(
  
  vector_riesgo ~
    
    Tipo_recinto +
    Valla_acondicionada +
    PRES_GANADO +
    
    altura_sc +
    peso_sc +
    temp_sc +
    
    s(fecha_sc) +
    
    (1|ID_COTO),
  
  family = bernoulli(),
  
  data = datos_aporte,
  
  prior = priors,
  
  chains = 4,
  iter = 4000,
  warmup = 2000,
  
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  seed = 123
)

summary(modelo_riesgo)

bayes_R2(modelo_riesgo)

###############################################################
# MODELO 2
# SEGURIDAD
###############################################################

modelo_seguridad <- brm(
  
  seguro ~
    
    Tipo_recinto +
    Valla_acondicionada +
    PRES_GANADO +
    
    s(fecha_sc) +
    
    (1|ID_COTO),
  
  family = bernoulli(),
  
  data = aportes,
  
  prior = priors,
  
  chains = 4,
  iter = 4000,
  warmup = 2000,
  
  control = list(
    adapt_delta = 0.99
  ),
  
  seed = 123
)

summary(modelo_seguridad)

bayes_R2(modelo_seguridad)

###############################################################
# MODELO 3
# PRIMERA ESPECIE
###############################################################

modelo_primero <- brm(
  
  primero_vector ~
    
    Tipo_recinto +
    Valla_acondicionada +
    PRES_GANADO +
    
    s(fecha_sc)+
    
    (1|ID_COTO)
  ,
  
  family = bernoulli(),
  
  data = primeros,
  
  prior = priors,
  
  chains = 4,
  iter = 4000,
  warmup = 2000,
  
  seed = 123
)

summary(modelo_primero)

bayes_R2(modelo_primero)

###############################################################
# MODELO 4
# TIEMPO HASTA LLEGADA DE BUITRES
###############################################################

aportes_buitres <- aportes %>%
  filter(
    is.finite(tiempo_buitres),
    tiempo_buitres > 0
  )

modelo_tiempo_buitres <- brm(
  
  tiempo_buitres ~
    
    Tipo_recinto +
    Valla_acondicionada +
    PRES_GANADO +
    
    peso_sc +
    temp_sc +
    
    s(fecha_sc)+
    
    (1|ID_COTO)
  ,
  
  family = Gamma(link="log"),
  
  data = aportes_buitres,
  
  prior = priors,
  
  chains = 4,
  iter = 4000,
  warmup = 2000,
  
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  seed = 123
)

summary(modelo_tiempo_buitres)

bayes_R2(modelo_tiempo_buitres)

###############################################################
# MODELO 5
# ABUNDANCIA DE BUITRES
###############################################################

buitres <- datos_aporte %>%
  filter(
    Especie_detectada == "Buitres",
    N > 0
  )

modelo_abundancia <- brm(
  
  N ~
    
    Tipo_recinto +
    Valla_acondicionada +
    PRES_GANADO +
    
    peso_sc +
    temp_sc +
    
    s(fecha_sc) +
    
    (1|ID_COTO),
  
  family = negbinomial(),
  
  data = buitres,
  
  prior = priors,
  
  chains = 4,
  iter = 4000,
  warmup = 2000,
  
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  seed = 123
)

summary(modelo_abundancia)

bayes_R2(modelo_abundancia)

###############################################################
# CHEQUEOS
###############################################################

pp_check(modelo_riesgo)

pp_check(modelo_seguridad)

pp_check(modelo_primero)

pp_check(modelo_tiempo_buitres)

pp_check(modelo_abundancia)

###############################################################


# Final model objects: local storage ONLY; do not upload.
saveRDS(modelo_riesgo, file.path(out_dir, "modelo_riesgo_bayes_revisado.rds"))
saveRDS(modelo_seguridad, file.path(out_dir, "modelo_seguridad_bayes.rds"))
saveRDS(modelo_primero, file.path(out_dir, "modelo_primero_bayes.rds"))
saveRDS(modelo_tiempo_buitres, file.path(out_dir, "modelo_tiempo_buitres_bayes_revisado.rds"))
saveRDS(modelo_abundancia, file.path(out_dir, "modelo_abundancia_bayes_revisado.rds"))

message("Fitted models saved only locally in: ", out_dir)
