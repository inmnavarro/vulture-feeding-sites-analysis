# NOTE: the field camera-trap dataset is private and is NOT included in this repo.
# All CSVs produced by this script remain in the local SFS_DATA_DIR output folder.
# The descriptive consumption rule excludes observations explicitly marked as
# eating vegetation (e.g. dogs grazing), even when coded behav == 'Consumo'.
################################################################################
# SUPPLEMENTARY FEEDING SITES | MANUSCRIPT FIGURES | VERSION 3
# Clean, journal-style graphics | English terminology | subtle vector silhouettes
# V3.1: English month abbreviations regardless of local R/Windows locale
# Reads existing brms RDS fits. Does NOT re-fit models or compile Stan.
# Outputs: 5 main/optional + 3 supplementary figures (PDF, PNG and TIFF at 600 dpi),
#          posterior prediction tables, contrasts, coefficients and model R2.
################################################################################

# 0. PACKAGES AND PATHS ---------------------------------------------------------
needed <- c('brms', 'ggplot2', 'dplyr', 'readr', 'lubridate', 'patchwork')
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing)) {
  stop('Install packages first: install.packages(c(',
       paste(sprintf('"%s"', missing), collapse=', '), '))')
}
suppressPackageStartupMessages({
  library(brms)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(lubridate)
  library(patchwork)
})

project_dir <- Sys.getenv("SFS_DATA_DIR", unset=".")
if (!file.exists(file.path(project_dir, "BASE_DATOS_BUITRES.csv"))) {
  stop("Private data not found: ", file.path(project_dir, "BASE_DATOS_BUITRES.csv"))
}
output_dir <- file.path(project_dir, 'FIGURAS_NATURE_V3')
dir.create(output_dir, showWarnings=FALSE, recursive=TRUE)

fit_dir <- Sys.getenv("SFS_MODEL_DIR", unset=file.path(project_dir, "models"))
fit_files <- c(
  risk       = 'modelo_riesgo_bayes_revisado.rds',
  safety     = 'modelo_seguridad_bayes.rds',
  arrival    = 'modelo_primero_bayes.rds',
  vulture_t  = 'modelo_tiempo_buitres_bayes_revisado.rds',
  vulture_n  = 'modelo_abundancia_bayes_revisado.rds'
)
missing_fits <- fit_files[!file.exists(file.path(fit_dir, fit_files))]
if (length(missing_fits)) {
  stop('Missing corrected RDS files in project folder: ',
       paste(missing_fits, collapse=', '),
       '\nDo not silently use the previous divergent fits.')
}
fits <- lapply(fit_files, function(f) readRDS(file.path(fit_dir, f)))

# 1. DATA AND AGREED ENGLISH TERMINOLOGY ---------------------------------------
raw <- readr::read_csv2(
  file.path(project_dir, 'BASE_DATOS_BUITRES.csv'),
  locale=readr::locale(encoding='Latin1'),
  show_col_types=FALSE
)
raw <- raw %>%
  filter(!is.na(ID_COTO), trimws(as.character(ID_COTO))!='') %>%
  mutate(
    Fecha_aporte=lubridate::dmy(Fecha_aporte),
    Fecha_evento=lubridate::dmy(Fecha_evento)
  )
post <- raw %>% filter(POST_aporte=='SI', !is.na(ID_aporte),
                       trimws(as.character(ID_aporte))!='')

site_names <- c(MULADAR='Big feeding station', PUNTO='Carcass deposition point')
taxon_names <- c(
  'Buitres'='Vultures', 'Zorro'='Red fox', 'Perro'='Dog',
  'Ciervo'='Red deer', 'Liebre'='Hare',
  'Mustélido'='Mustelid', 'Mustelidos'='Mustelid', 'Mustelido'='Mustelid'
)
# Site-safety fitted response is fox-based; first-arrival is guild-based.
model_names <- c(
  risk='Facultative-scavenger detections',
  safety='Site safety',
  arrival='Facultative scavenger first',
  vulture_t='Time to vulture arrival',
  vulture_n='Vultures per detection'
)
predictor_names <- c(
  'Tipo_recintoPUNTO'='Carcass deposition point',
  'Valla_acondicionadaSI'='Adapted fence',
  'PRES_GANADOSI'='Livestock presence',
  'altura_sc'='Fence height (1 SD)',
  'peso_sc'='Carcass biomass (1 SD)',
  'temp_sc'='Temperature (1 SD)'
)
relabel <- function(x, lookup) {
  x <- as.character(x)
  v <- unname(lookup[x])
  ifelse(is.na(v), x, v)
}

# Cross-tabulation of combinations ACTUALLY in the observed dataset.
observed <- post %>%
  distinct(ID_aporte, Tipo_recinto, Valla_acondicionada, PRES_GANADO) %>%
  count(Tipo_recinto, Valla_acondicionada, PRES_GANADO, name='n_provisionings')
write.csv(observed, file.path(output_dir,'Observed_design_combinations.csv'),row.names=FALSE)
print(observed)

# 2. REPRODUCIBLE JOURNAL THEME -------------------------------------------------
INK <- '#253640'
TEAL <- '#315A65'
OCHRE <- '#B87950'
PALE <- '#D4DBDE'
PALETTE <- c('Big feeding station'=TEAL, 'Carcass deposition point'=OCHRE)

theme_nature <- function(base_size=11) {
  theme_classic(base_size=base_size, base_family='sans') +
    theme(
      plot.title=element_text(face='bold',colour=INK,size=base_size+2,
                              margin=margin(b=6)),
      plot.subtitle=element_text(colour='#52636B',size=base_size-1,
                                 margin=margin(b=9)),
      plot.caption=element_text(colour='#56666C',size=base_size-2,
                                hjust=0,margin=margin(t=9)),
      axis.title=element_text(colour=INK),
      axis.text=element_text(colour=INK),
      axis.line=element_line(linewidth=.35,colour=INK),
      axis.ticks=element_line(linewidth=.3,colour=INK),
      strip.background=element_blank(),
      strip.text=element_text(colour=INK,face='bold',size=base_size),
      panel.grid=element_blank(),
      legend.title=element_text(face='bold'),
      legend.position='none',
      plot.margin=margin(10,13,11,12)
    )
}

export_figure <- function(g, name, width=12, height=6) {
  path <- file.path(output_dir,name)
  ggplot2::ggsave(paste0(path,'.pdf'),plot=g,width=width,height=height,
                  units='in',device='pdf',bg='white')
  ggplot2::ggsave(paste0(path,'.png'),plot=g,width=width,height=height,
                  units='in',dpi=600,device='png',bg='white')
  ggplot2::ggsave(paste0(path,'.tiff'),plot=g,width=width,height=height,
                  units='in',dpi=600,device='tiff',compression='lzw',bg='white')
  message('Saved: ', name)
}

# 3. POSTERIOR PREDICTIONS -------------------------------------------------------
# Every contrast is between two configurations represented by observed data.
# We do NOT compare livestock present in big feeding stations, since absent.
comparisons <- data.frame(
  comparison=c('Site type','Site type',
               'Fence adaptation','Fence adaptation',
               'Livestock presence','Livestock presence'),
  label=c('Big feeding station','Carcass deposition point',
          'No','Yes','Absent','Present'),
  Tipo_recinto=c('MULADAR','PUNTO','PUNTO','PUNTO','PUNTO','PUNTO'),
  Valla_acondicionada=c('NO','NO','NO','SI','SI','SI'),
  PRES_GANADO=c('NO','NO','NO','NO','NO','SI'),
  stringsAsFactors=FALSE
)
comparisons$comparison <- factor(
  comparisons$comparison,
  levels=c('Site type','Fence adaptation','Livestock presence')
)

make_newdata <- function(model, des, dates_z=NULL) {
  md <- model$data
  if (is.null(md) || nrow(md)<1L) stop('Saved model does not contain data.')
  nd <- md[rep(1L,nrow(des)), , drop=FALSE]
  for (v in c('Tipo_recinto','Valla_acondicionada','PRES_GANADO')) {
    if (v %in% names(md)) {
      nd[[v]] <- if (is.factor(md[[v]])) {
        factor(des[[v]],levels=levels(md[[v]]))
      } else as.character(des[[v]])
    }
  }
  for (v in c('altura_sc','peso_sc','temp_sc','fecha_sc')) {
    if (v %in% names(nd)) nd[[v]] <- 0
  }
  if (!is.null(dates_z) && 'fecha_sc' %in% names(nd)) nd$fecha_sc <- dates_z
  if ('ID_COTO' %in% names(md)) {
    coto <- as.character(md$ID_COTO[which(!is.na(md$ID_COTO))[1]])
    nd$ID_COTO <- if (is.factor(md$ID_COTO)) {
      factor(rep(coto,nrow(nd)),levels=levels(md$ID_COTO))
    } else rep(coto,nrow(nd))
  }
  nd
}

summarise_draws <- function(x) {
  if (!is.matrix(x)) stop('Expected a posterior matrix [draws x observations].')
  data.frame(
    Estimate=colMeans(x),
    Q2.5=apply(x,2,quantile,probs=.025,names=FALSE),
    Q97.5=apply(x,2,quantile,probs=.975,names=FALSE)
  )
}

prediction_data <- dplyr::bind_rows(lapply(names(fits),function(nm) {
  nd <- make_newdata(fits[[nm]], comparisons)
  dr <- brms::posterior_epred(fits[[nm]],newdata=nd,re_formula=NA)
  z <- dplyr::bind_cols(comparisons,summarise_draws(dr))
  z$model <- nm
  z
}))
write.csv(prediction_data,file.path(output_dir,'Posterior_predictions.csv'),row.names=FALSE)

prediction_plot <- function(df, title=NULL, y_label='Predicted probability',probability=TRUE) {
  # each comparison has TWO observations in source and these are not paired samples
  df$label <- factor(df$label,
    levels=c('Big feeding station','Carcass deposition point','No','Yes','Absent','Present'))
  p <- ggplot(df,aes(x=label,y=Estimate)) +
    geom_linerange(aes(ymin=Q2.5,ymax=Q97.5),colour=TEAL,linewidth=.75) +
    geom_point(size=2.8,shape=21,stroke=.7,fill='white',colour=INK) +
    facet_wrap(~comparison,scales='free_x',nrow=1) +
    labs(title=title,x=NULL,y=y_label) +
    theme_nature(10.5) +
    theme(axis.text.x=element_text(angle=18,hjust=1))
  if (probability) p <- p+coord_cartesian(ylim=c(0,1))
  p
}

# 4. FIGURE 1 — SCHEMATIC + SAMPLING COUNTS (VECTOR SILHOUETTES) -----------------
# Stylised animal silhouettes made from embedded vector polygon coordinates.
# These are illustrative design elements, NOT photographs or identification aids.
fox_xy <- rbind(
  c(.02,.51),c(.07,.65),c(.18,.70),c(.32,.64),c(.43,.55),
  c(.60,.56),c(.67,.65),c(.71,.85),c(.76,.72),c(.81,.85),
  c(.84,.68),c(.98,.56),c(.88,.50),c(.81,.51),
  c(.79,.22),c(.72,.22),c(.70,.49),c(.48,.49),
  c(.43,.19),c(.36,.19),c(.38,.51),c(.21,.54),c(.10,.47)
)
vulture_xy <- rbind(
  c(.02,.61),c(.13,.69),c(.24,.77),c(.35,.72),c(.42,.64),
  c(.47,.65),c(.49,.78),c(.52,.85),c(.56,.76),c(.57,.63),
  c(.65,.74),c(.76,.78),c(.87,.69),c(.98,.59),c(.84,.55),
  c(.73,.60),c(.63,.58),c(.57,.49),c(.54,.34),c(.51,.23),
  c(.48,.35),c(.46,.49),c(.39,.58),c(.27,.60),c(.16,.56)
)
icon_df <- function(xy,x0,y0,w,h,kind,id) {
  data.frame(x=x0+xy[,1]*w,y=y0+xy[,2]*h,kind=kind,id=id)
}
icons <- dplyr::bind_rows(
  icon_df(vulture_xy,2.4,2.35,1.6,1.1,'vulture','v1'),
  icon_df(vulture_xy,7.45,2.35,1.35,1.1,'vulture','v2'),
  icon_df(fox_xy,.15,1.1,1.35,1.1,'fox','f1'),
  icon_df(fox_xy,5.95,1.1,1.20,1.1,'fox','f2')
)

f1a <- ggplot() +
  annotate('rect',xmin=1.55,xmax=5.45,ymin=1.05,ymax=3.75,
           fill='#ECF1F2',colour=TEAL,linewidth=.65) +
  annotate('rect',xmin=7.15,xmax=9.75,ymin=1.05,ymax=3.75,
           fill='#F6F0EA',colour=OCHRE,linewidth=.65) +
  geom_segment(data=data.frame(x=c(1.55,2.85,4.15,5.45,7.15,8.45,9.75)),
               aes(x=x,xend=x,y=1.05,yend=3.85),linewidth=.6,colour='#6F848B') +
  geom_polygon(data=icons,aes(x=x,y=y,group=id,fill=kind),colour=NA) +
  scale_fill_manual(values=c(vulture=INK,fox='#6A7378')) +
  annotate('point',x=c(3.0,8.1),y=c(1.65,1.65),size=5.0,
           shape=23,fill='#958273',colour='#958273') +
  annotate('text',x=3.5,y=4.45,label='Big feeding station',size=4.15,
           fontface='bold',colour=INK) +
  annotate('text',x=8.4,y=4.45,label='Carcass deposition point',size=3.65,
           fontface='bold',colour=INK) +
  annotate('text',x=5.25,y=.44,label='Illustrative design; not to scale',
           colour='#69777C',size=2.8) +
  coord_fixed(xlim=c(0,10.5),ylim=c(0,5),expand=FALSE) +
  theme_void(base_family='sans')+
  theme(legend.position='none',plot.margin=margin(7,5,7,8))

estate_n <- raw %>% distinct(ID_COTO,Tipo_recinto) %>%
  count(Tipo_recinto,name='Count') %>% mutate(Metric='Hunting estates')
provision_n <- post %>% distinct(ID_aporte,Tipo_recinto) %>%
  count(Tipo_recinto,name='Count') %>% mutate(Metric='Carcass provisionings')
counts <- bind_rows(estate_n,provision_n) %>%
  mutate(Site=factor(relabel(Tipo_recinto,site_names),levels=names(PALETTE)),
         Metric=factor(Metric,levels=c('Hunting estates','Carcass provisionings')))
f1b <- ggplot(counts,aes(x=Site,y=Count,fill=Site)) +
  geom_col(width=.6) +
  geom_text(aes(label=Count),vjust=-.4,size=3.8,colour=INK) +
  facet_wrap(~Metric,scales='free_y',ncol=2) +
  scale_fill_manual(values=PALETTE) +
  scale_y_continuous(expand=expansion(mult=c(0,.18))) +
  labs(x=NULL,y='Number')+
  theme_nature(9.5) +
  theme(axis.text.x=element_text(angle=25,hjust=1))

f1 <- (f1a / f1b) + patchwork::plot_layout(heights=c(1.05,1)) +
  plot_annotation(title='Supplementary feeding sites: study design',
                  subtitle='Schematic comparison and number of independent sampling units',
                  tag_levels='A',
                  theme=theme(plot.title=element_text(face='bold',size=14,colour=INK),
                              plot.subtitle=element_text(size=10,colour='#52636B')))
export_figure(f1,'Figure_1_study_design',width=13,height=9.1)

# 5. FIGURE 2 — DETECTIONS OF FACULTATIVE SCAVENGERS ----------------------------
f2 <- prediction_plot(filter(prediction_data,model=='risk'),
  title='Facultative scavenger detections',
  y_label='Pr(recorded detection is facultative)') +
  labs(caption=paste('Conditional prediction for observed site configurations;',
                     'continuous predictors fixed at their standardized means.',
                     '95% credible intervals.'))
export_figure(f2,'Figure_2_facultative_detection',width=13,height=5.3)

# 6. FIGURE 3 — SITE SAFETY + FIRST ARRIVAL -------------------------------------
f3data <- prediction_data %>%
  filter(model %in% c('safety','arrival')) %>%
  mutate(Response=factor(ifelse(model=='safety','Site safety',
                                  'Facultative scavenger first'),
                         levels=c('Site safety','Facultative scavenger first')),
         label=factor(label,levels=c('Big feeding station','Carcass deposition point',
                                     'No','Yes','Absent','Present')))
f3 <- ggplot(f3data,aes(x=label,y=Estimate)) +
  geom_linerange(aes(ymin=Q2.5,ymax=Q97.5),linewidth=.8,colour=TEAL) +
  geom_point(shape=21,fill='white',colour=INK,stroke=.7,size=2.5) +
  facet_grid(Response~comparison,scales='free_x',space='free_x') +
  coord_cartesian(ylim=c(0,1)) +
  labs(title='Site safety and recorded arrival order',x=NULL,
       y='Predicted probability',
       caption=paste('Site safety is the study\'s record-based index:',
                     'no red fox detected before vultures.',
                     'It does not imply measured pathogen absence.')) +
  theme_nature(10.5) +
  theme(axis.text.x=element_text(angle=17,hjust=1),
        strip.text.y=element_text(angle=0))
export_figure(f3,'Figure_3_site_safety_and_arrival',width=14.2,height=7.5)

# 7. FIGURE 4 — VULTURE RESPONSE ------------------------------------------------
ft <- prediction_data %>% filter(model=='vulture_t') %>%
  mutate(Estimate=Estimate/60,Q2.5=Q2.5/60,Q97.5=Q97.5/60)
fn <- prediction_data %>% filter(model=='vulture_n')
f4a <- prediction_plot(ft,title='Time to first recorded vulture',
                       y_label='Predicted time (hours)',probability=FALSE)
f4b <- prediction_plot(fn,title='Vultures per recorded detection',
                       y_label='Predicted number of vultures',probability=FALSE)
f4 <- (f4a / f4b) + patchwork::plot_annotation(
  title='Vulture responses to supplementary feeding sites',tag_levels='A',
  theme=theme(plot.title=element_text(face='bold',colour=INK,size=14)))
export_figure(f4,'Figure_4_vulture_responses',width=13.5,height=10.2)

# 8. FIGURE 5 — DOCUMENTED DETECTION AND CARCASS CONSUMPTION ---------------------
# Observed proportions per provisioning ID (descriptive, not detection-adjusted).
base_provision <- post %>%
  mutate(facultative=ifelse(is.na(`Clasificación`),FALSE,
                            `Clasificación`=='Oportunista'),
         # The field log marked some dog grazing (grass eating) as 'Consumo'.
         # These are NOT documented consumption of hunting carrion.
         grass_only=grepl('hierba|grass', coalesce(as.character(Notas), ''),
                          ignore.case=TRUE),
         consumed=facultative & !is.na(behav) & behav=='Consumo' & !grass_only) %>%
  group_by(ID_aporte,Tipo_recinto,Valla_acondicionada) %>%
  summarise(facultative_recorded=any(facultative),
            facultative_consumed=any(consumed),.groups='drop')
if (anyDuplicated(base_provision$ID_aporte)) {
  stop('One ID_aporte has inconsistent site/fence coding. Check your database.')
}

calculate_fractions <- function(column,facet_label) {
  z <- base_provision %>%
    group_by(Category=.data[[column]]) %>%
    summarise(n=n(),recorded=sum(facultative_recorded),
              consumed=sum(facultative_consumed),.groups='drop')
  if (identical(column,'Tipo_recinto')) z$Category <- relabel(z$Category,site_names)
  if (identical(column,'Valla_acondicionada')) {
    z$Category <- relabel(z$Category,c(NO='No',SI='Yes'))
  }
  dplyr::bind_rows(
    z %>% transmute(Group=facet_label,Category,n,
                    Outcome='Facultative scavenger recorded',positive=recorded),
    z %>% transmute(Group=facet_label,Category,n,
                    Outcome='Carcass consumption documented',positive=consumed)
  )
}
observed_rates <- bind_rows(
  calculate_fractions('Tipo_recinto','Site type'),
  calculate_fractions('Valla_acondicionada','Fence adaptation')
) %>%
  mutate(Rate=positive/n,
         Annotation=paste0(positive,'/',n),
         Outcome=factor(Outcome,levels=c('Facultative scavenger recorded',
                                       'Carcass consumption documented')),
         Group=factor(Group,levels=c('Site type','Fence adaptation')))
write.csv(observed_rates,file.path(output_dir,'Observed_detection_and_consumption.csv'),
          row.names=FALSE)

f5 <- ggplot(observed_rates,aes(x=Category,y=Rate)) +
  geom_col(width=.60,fill=TEAL) +
  geom_text(aes(label=Annotation),vjust=-.45,colour=INK,size=3.7) +
  facet_grid(Outcome~Group,scales='free_x',space='free_x')+
  scale_y_continuous(labels=function(z)paste0(round(100*z),'%'),
                     limits=c(0,1.12),breaks=seq(0,1,.25),
                     expand=expansion(mult=c(0,0)))+
  labs(title='Observed facultative scavenger activity',x=NULL,
       y='Proportion of carcass provisionings',
       caption=paste('Counts shown as positive/identified provisionings.',
                     'Descriptive proportions; repeated provisionings are nested within estates.',
                     'Non-detection does not imply absence.'))+
  theme_nature(10.5)+
  theme(axis.text.x=element_text(angle=15,hjust=1),
        strip.text.y=element_text(angle=0))
export_figure(f5,'Figure_5_documented_consumption',width=12,height=7.3)

# 8. SUPPLEMENTARY S1 — ALL TAXA TRANSLATED -----------------------------------
sp <- post %>% filter(!is.na(Especie_detectada)) %>%
  count(Especie_detectada,name='n') %>%
  mutate(Taxon=relabel(Especie_detectada,taxon_names),
         Guild=case_when(
           Especie_detectada=='Buitres' ~ 'Obligate scavengers',
           Especie_detectada %in% c('Zorro','Perro','Mustélido','Mustelidos','Mustelido') ~
             'Facultative scavengers',
           TRUE ~ 'Other wildlife')) %>%
  arrange(n)
sp$Taxon <- factor(sp$Taxon,levels=unique(sp$Taxon))
write.csv(sp,file.path(output_dir,'Recorded_taxa_English.csv'),row.names=FALSE)
fs1 <- ggplot(sp,aes(x=Taxon,y=n,fill=Guild)) +
  geom_col(width=.67) +
  geom_text(aes(label=n),hjust=-.20,size=3.7,colour=INK) +
  coord_flip() +
  scale_fill_manual(values=c('Obligate scavengers'=TEAL,
                              'Facultative scavengers'=OCHRE,
                              'Other wildlife'='#9DA9AE')) +
  scale_y_continuous(expand=expansion(mult=c(0,.15))) +
  labs(title='Wildlife detected following carcass provisioning',
       x=NULL,y='Number of recorded detections',fill=NULL,
       caption='Detection records do not represent independent individuals.') +
  theme_nature(11) + theme(legend.position='bottom')
export_figure(fs1,'Supplementary_Figure_S1_taxa',width=8.4,height=5.5)

# 9. SUPPLEMENTARY S2 — MODEL COEFFICIENTS ------------------------------------
coef_df <- bind_rows(lapply(names(fits),function(nm) {
  z <- as.data.frame(brms::fixef(fits[[nm]],probs=c(.025,.975)))
  z$term <- rownames(z)
  z$model <- nm
  z
})) %>%
  filter(!grepl('^Intercept$|^sfecha_sc',term)) %>%
  mutate(Predictor=relabel(term,predictor_names),
         Model=factor(relabel(model,model_names),levels=unname(model_names)))
write.csv(coef_df,file.path(output_dir,'Bayesian_fixed_effects.csv'),row.names=FALSE)
fs2 <- ggplot(coef_df,aes(x=Estimate,y=Predictor)) +
  geom_vline(xintercept=0,linetype='dashed',linewidth=.4,colour='#8E9DA4')+
  geom_segment(aes(x=Q2.5,xend=Q97.5,yend=Predictor),
               colour=TEAL,linewidth=.8)+
  geom_point(shape=21,fill='white',colour=INK,stroke=.65,size=2.3)+
  facet_wrap(~Model,scales='free_x',ncol=2)+
  labs(title='Posterior fixed effects across models',
       x='Coefficient on model-specific link scale (95% CrI)',y=NULL,
       caption='Different likelihood families use different link scales; compare within models only.')+
  theme_nature(10)
export_figure(fs2,'Supplementary_Figure_S2_coefficients',width=11.8,height=10.3)

# 10. SUPPLEMENTARY S3 — TEMPORAL PREDICTIONS ----------------------------------
# Scale date separately from each model's saved training data, when available.
all_dates <- sort(unique(post$Fecha_aporte[!is.na(post$Fecha_aporte)]))
if (length(all_dates)<2) stop('Not enough provisioning dates for temporal predictions.')
date_grid <- seq(min(all_dates),max(all_dates),length.out=75)
# Robust date mapping, directly from known points in the fit's model data.
# (The coefficient is the same scale as when the fit was constructed.)
model_date_z <- function(model,dates) {
  md <- model$data
  if ('Fecha_aporte' %in% names(md) && 'fecha_sc' %in% names(md)) {
    d0 <- as.numeric(as.Date(md$Fecha_aporte))
    z0 <- as.numeric(md$fecha_sc)
    ok <- is.finite(d0)&is.finite(z0)
    if (length(unique(d0[ok]))>=2) {
      b <- stats::coef(stats::lm(z0[ok] ~ d0[ok]))
      return(unname(b[1]+b[2]*as.numeric(dates)))
    }
  }
  d0 <- as.numeric(raw$Fecha_aporte)
  warning('Date scaling approximated from full dataset for a saved model.')
  (as.numeric(dates)-mean(d0,na.rm=TRUE))/sd(d0,na.rm=TRUE)
}
ref_design <- data.frame(Tipo_recinto='PUNTO',Valla_acondicionada='SI',PRES_GANADO='NO')
temporal <- bind_rows(lapply(names(fits),function(nm) {
  des <- ref_design[rep(1,length(date_grid)),,drop=FALSE]
  nd <- make_newdata(fits[[nm]],des,dates_z=model_date_z(fits[[nm]],date_grid))
  draws <- brms::posterior_epred(fits[[nm]],newdata=nd,re_formula=NA)
  data.frame(Date=date_grid,model=nm,summarise_draws(draws))
})) %>% mutate(Model=factor(relabel(model,model_names),levels=unname(model_names)))
write.csv(temporal,file.path(output_dir,'Temporal_predictions.csv'),row.names=FALSE)
fs3 <- ggplot(temporal,aes(x=Date,y=Estimate)) +
  geom_ribbon(aes(ymin=Q2.5,ymax=Q97.5),fill=TEAL,alpha=.16)+
  geom_line(colour=TEAL,linewidth=.9)+
  facet_wrap(~Model,scales='free_y',ncol=2)+
  # Month labels from numeric month numbers: independent of Catalan/Spanish Windows locale.
  scale_x_date(date_breaks='1 month',
    labels=function(x) paste0(month.abb[as.integer(format(x,'%m'))], ' ', format(x,'%Y')))+
  labs(title='Temporal profiles of recorded scavenger responses',
       x='Carcass provisioning date',y='Posterior expected value',
       caption=paste('Carcass deposition point, adapted fence, no livestock;',
                     'other predictors at standardized mean. Units vary by panel.'))+
  theme_nature(10.5)
export_figure(fs3,'Supplementary_Figure_S3_temporal',width=11.5,height=9.2)

# 11. ADDITIONAL TABLES --------------------------------------------------------
r2 <- bind_rows(lapply(names(fits),function(nm) {
  z <- as.data.frame(brms::bayes_R2(fits[[nm]]))
  z$Model <- model_names[[nm]]
  z
}))
write.csv(r2,file.path(output_dir,'Bayesian_R2.csv'),row.names=FALSE)

# Within-model posterior contrasts on the response scale.
contrast_defs <- data.frame(
  Contrast=c('Deposition point minus big feeding station',
             'Adapted minus unadapted fence',
             'Livestock present minus absent'),
  i0=c(1L,3L,5L),i1=c(2L,4L,6L)
)
contrasts <- bind_rows(lapply(names(fits),function(nm) {
  nd <- make_newdata(fits[[nm]],comparisons)
  dr <- brms::posterior_epred(fits[[nm]],newdata=nd,re_formula=NA)
  bind_rows(lapply(seq_len(nrow(contrast_defs)),function(j) {
    v <- dr[,contrast_defs$i1[j]]-dr[,contrast_defs$i0[j]]
    data.frame(Model=model_names[[nm]],Contrast=contrast_defs$Contrast[j],
               Median=median(v),CrI_low=unname(quantile(v,.025)),
               CrI_high=unname(quantile(v,.975)),Pr_positive=mean(v>0))
  }))
}))
write.csv(contrasts,file.path(output_dir,'Posterior_pairwise_contrasts.csv'),row.names=FALSE)

capture.output(sessionInfo(),file=file.path(output_dir,'Session_info.txt'))
fig_names <- c('Figure_1_study_design','Figure_2_facultative_detection',
               'Figure_3_site_safety_and_arrival','Figure_4_vulture_responses',
               'Figure_5_documented_consumption',
               'Supplementary_Figure_S1_taxa',
               'Supplementary_Figure_S2_coefficients',
               'Supplementary_Figure_S3_temporal')
expected <- unlist(lapply(fig_names,function(nm)
  file.path(output_dir,paste0(nm,'.',c('pdf','png','tiff')))))
missing_output <- expected[!file.exists(expected)]
if (length(missing_output)) stop('These outputs were NOT generated: ',
                                paste(basename(missing_output),collapse=', '))
message('\nSUCCESS: 8 English-language figures exported to: ',normalizePath(output_dir))
message('All RDS models are unchanged. Manuscript labels: Big feeding station,',
        ' Carcass deposition point, Red deer, Livestock presence, Site safety.')
