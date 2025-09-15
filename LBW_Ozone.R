# --- Load Packages ---
library(dplyr)
library(tidyr)
library(purrr)
library(lubridate)
library(flextable)
library(officer)

# --- Load Data ---
O3_data <- readRDS("C:/Users/smoss/Downloads/O3_8hTWA_maxDiaria.RDS")
load("C:/Users/smoss/Downloads/base_nacimientos_O3.RData")

hospital_data <- data.frame(base.muni.O3)
rm(base.muni.O3)

# --- Prepare O3 Data ---
O3_data <- O3_data %>%
  mutate(
    PROVINCIA = as.integer(PROVINCIA),
    MUNICIPIO = as.integer(MUNICIPIO),
    MUNI_INE = PROVINCIA * 1000 + MUNICIPIO
  )

# Add missing_7day_stretch and clean

O3_data_clean <- O3_data %>%
  arrange(MUNI_INE, FECHA) %>%
  group_by(MUNI_INE) %>%
  mutate(
    is_na = is.na(value),
    run_id = with(rle(is_na), rep(seq_along(lengths), lengths)),
    run_length = ave(is_na, run_id, FUN = length),
    missing_7day_stretch = is_na & run_length >= 7
  ) %>%
  ungroup() %>%
  filter(FECHA <= as.Date("2017-12-31")) %>%    # FIX: was DATE, should be FECHA
  filter(!missing_7day_stretch)

# Keep only municipalities that overlap
overlapping_muni_ine <- intersect(
  unique(O3_data_clean$MUNI_INE),
  unique(hospital_data$MUNI_INE)
)

O3_data_clean <- O3_data_clean %>%
  filter(MUNI_INE %in% overlapping_muni_ine)

# --- Prepare Hospital Data ---
hospital_data <- hospital_data %>%
  mutate(
    MUNI_INE = as.integer(MUNI_INE),
    date = make_date(year = ANOPAR, month = MESPAR, day = DIAPAR),
    low_birth_weight = ifelse(PESON < 2500, 1, 0)
  )

hospital_filtered <- hospital_data %>%
  filter(MUNI_INE %in% overlapping_muni_ine)

# --- Weekly Ozone Matrix ---
create_weekly_matrix <- function(df, value_col) {
  df %>%
    arrange(FECHA) %>%
    mutate(week = floor_date(FECHA, "week")) %>%
    group_by(week) %>%
    summarise(mean_val = mean(.data[[value_col]], na.rm = TRUE), .groups = "drop")
}

weekly_df <- O3_data_clean %>%
  group_by(MUNI_INE) %>%
  group_modify(~ create_weekly_matrix(.x, "value")) %>%
  ungroup()

# --- Exposure Extractor Function ---
get_exposure_fast <- function(muni, dob, weeks) {
  df <- weekly_df[weekly_df$MUNI_INE == muni, ]
  if (nrow(df) == 0) return(rep(NA_real_, weeks))
  df <- df[df$week <= dob, ]
  if (nrow(df) == 0) return(rep(NA_real_, weeks))
  df <- df[order(df$week, decreasing = TRUE), ]
  vals <- head(df$mean_val, weeks)
  if (length(vals) < weeks) vals <- c(vals, rep(NA_real_, weeks - length(vals)))
  rev(vals)
}

# --- Apply in Chunks ---
chunk_size <- 5000
n <- nrow(hospital_filtered)
max_weeks <- max(hospital_filtered$SEMGES, na.rm = TRUE)

exposure_list <- vector("list", n)

for (i in seq(1, n, by = chunk_size)) {
  cat("Processing rows", i, "to", min(i + chunk_size - 1, n), "...\n")
  idx <- i:min(i + chunk_size - 1, n)
  exposures <- map2(
    hospital_filtered$MUNI_INE[idx],
    hospital_filtered$date[idx],
    ~ get_exposure_fast(.x, .y, hospital_filtered$SEMGES[idx[1]])
  )
  exposure_list[idx] <- exposures
  gc()
}

# --- Convert to DataFrame ---
exposure_df <- do.call(rbind, lapply(exposure_list, function(x) {
  x <- as.numeric(x)
  c(x, rep(NA_real_, max_weeks - length(x)))
}))

exposure_df <- as.data.frame(exposure_df)
colnames(exposure_df) <- paste0("V", 1:max_weeks)

hospital_exposure <- bind_cols(hospital_filtered, exposure_df)

# --- Compute Trimester Exposures ---
compute_Ts <- function(df_batch) {
  df_batch %>%
    rowwise() %>%
    mutate(
      T1 = mean(c_across(V1:V12), na.rm = TRUE),
      T2 = mean(c_across(V13:V27), na.rm = TRUE),
      T3 = if (SEMGES >= 28) {
        vals <- unlist(across(starts_with("V")))
        mean(vals[28:SEMGES], na.rm = TRUE)
      } else {
        NA_real_
      }
    ) %>%
    ungroup()
}


batch_size <- 5000
n <- nrow(hospital_exposure)
batches <- split(hospital_exposure, (seq_len(n) - 1) %/% batch_size)

results_list <- lapply(batches, compute_Ts)
hospital_exposure <- bind_rows(results_list)

# --- Final cleaned dataset ---
cleaned_data <- hospital_exposure %>%
  filter(!is.na(T1) & !is.na(T2) & !is.na(T3))


# analysis
cleaned_data$mean<-(cleaned_data$T1 + cleaned_data$T2 + cleaned_data$T3)/3

model_trimester_cont <- glm(
  low_birth_weight ~ mean+ T1 + T2 + T3,
  family = binomial,
  data = cleaned_data
)



summary(model_trimester)


cuartiles_t1<- quantile(cleaned_data$T1,c(0,0.25,0.5,0.75,1))

cleaned_data$T1_4c<-1

cleaned_data[cleaned_data$T1 > cuartiles_t1[2], "T1_4c"]<-2
cleaned_data[cleaned_data$T1 > cuartiles_t1[3], "T1_4c"]<-3
cleaned_data[cleaned_data$T1 > cuartiles_t1[4], "T1_4c"]<-4
cleaned_data$T1_4c<-as.factor(cleaned_data$T1_4c)


cuartiles_t2<-quantile(cleaned_data$T2,c(0,0.25,0.5,0.75,1))

cleaned_data$T2_4c<-1

cleaned_data[cleaned_data$T2> cuartiles_t2[2], "T2_4c"]<-2
cleaned_data[cleaned_data$T2> cuartiles_t2[3], "T2_4c"]<-3
cleaned_data[cleaned_data$T2> cuartiles_t2[4], "T2_4c"]<-4
cleaned_data$T2_4c<-as.factor(cleaned_data$T2_4c)

cuartiles_t3<-quantile(cleaned_data$T3,c(0,0.25,0.5,0.75,1))

cleaned_data$T3_4c<-1

cleaned_data[cleaned_data$T3> cuartiles_t3[2], "T3_4c"]<-2
cleaned_data[cleaned_data$T3> cuartiles_t3[3], "T3_4c"]<-3
cleaned_data[cleaned_data$T3> cuartiles_t3[4], "T3_4c"]<-4
cleaned_data$T3_4c<-as.factor(cleaned_data$T3_4c)

cuartiles_mean<- quantile(cleaned_data$mean,c(0,0.25,0.5,0.75,1))

cleaned_data$mean_4c<-1
cleaned_data[cleaned_data$mean > cuartiles_mean[2], "mean_4c"]<-2
cleaned_data[cleaned_data$mean > cuartiles_mean[3], "mean_4c"]<-3
cleaned_data[cleaned_data$mean > cuartiles_mean[4], "mean_4c"]<-4
cleaned_data$mean_4c<-as.factor(cleaned_data$mean_4c)

cleaned_data <- cleaned_data %>%
  mutate(cat_date_birth = case_when(
    ANOPAR < 2005 ~ 1,
    ANOPAR >= 2005 & ANOPAR <= 2008 ~ 2,
    ANOPAR >=2009 & ANOPAR <= 2012 ~3,
    ANOPAR >=2013 & ANOPAR <= 2017 ~ 4,
    TRUE ~ NA_real_
  ))
cleaned_data$cat_date_birth<-as.factor(cleaned_data$cat_date_birth)
summary(cleaned_data$cat_date_birth)

cleaned_data <- cleaned_data %>%
  mutate(cat_edad_mother = case_when(
    EDADM < 20 ~1,
    EDADM >= 20 & EDADM <= 34 ~ 2,
    EDADM >=35 & EDADM <= 40 ~3,
    EDADM >40 ~ 4,
    TRUE ~ NA_real_
  ))
cleaned_data$cat_edad_mother<-as.factor(cleaned_data$cat_edad_mother)
summary(cleaned_data$cat_edad_mother)


cleaned_data <- cleaned_data %>%
  mutate(cat_edad_mother = case_when(
    EDADM < 20 ~1,
    EDADM >= 20 & EDADM <= 34 ~ 2,
    EDADM >=35 & EDADM <= 40 ~3,
    EDADM >40 ~ 4,
    TRUE ~ NA_real_
  ))
cleaned_data$cat_edad_mother<-as.factor(cleaned_data$cat_edad_mother)
summary(cleaned_data$cat_edad_mother)

cleaned_data <- cleaned_data %>%
  mutate(cat_socioec = case_when(
    cond_socio_media>= 0.31 & cond_socio_media <= 0.57 ~1,
    cond_socio_media >= 0.58 & cond_socio_media <= 0.83 ~ 2,
    cond_socio_media >=0.84 & cond_socio_media <= 1.09 ~3,
    cond_socio_media >= 1.10 & cond_socio_media<= 1.35~4,
    EDADM >1.36 ~ 5,
    TRUE ~ NA_real_
  ))

cleaned_data$cat_socioec<-as.factor(cleaned_data$cat_socioec)
summary(cleaned_data$cat_socioec)

cleaned_data <- cleaned_data %>%
  mutate(cat_tasa_act = case_when(
    tasa_actividad< 41 ~1,
    tasa_actividad < 61 ~ 2,
    tasa_actividad <81 ~3,
    tasa_actividad >=81 ~ 4,
    TRUE ~ NA_real_
  ))

cleaned_data$cat_tasa_act<-as.factor(cleaned_data$cat_tasa_act)
summary(cleaned_data$cat_tasa_act)
# single association

model_trim_T1_4c<-glm(low_birth_weight ~ T1, family = binomial, data=cleaned_data)
model_trim_T2_4c<-glm(low_birth_weight ~ T2, family = binomial, data=cleaned_data)
model_trim_T3_4c<-glm(low_birth_weight ~ T3, family = binomial, data=cleaned_data)
model_mean_4c<-glm(low_birth_weight ~  mean, family = binomial, data=cleaned_data)


# Get summary of the model
summary_model <- summary(model_trim_T2_4c)

# Extract coefficients, SE, z, p
coef_table <- summary_model$coefficients
# exp() for ORs
OR_CI <- exp(cbind(OR = coef(model_trim_T2_4c),
                   confint(model_trim_T2_4c)))

# Combine into one data frame
result_table <- cbind(
  OR_CI,
  p_value = coef_table[, "Pr(>|z|)"]
)

# Round nicely
result_table <- round(result_table, 3)

result_table



# association and adjusted GLM CAT
model_trim_t1_4c_adj<-glm(low_birth_weight ~ factor(T1_4c) +SEXO + factor(cat_edad_mother) + EDADP + factor(cat_date_birth) + factor(cat_socioec) + cat_tasa_act + SEMGES + factor(CCAA) + factor(TMUNR),family = binomial, data=cleaned_data)
model_trim_t2_4c_adj<-glm(low_birth_weight ~ factor(T2_4c) +SEXO + factor(cat_edad_mother) + EDADP + factor(cat_date_birth) + factor(cat_socioec) + cat_tasa_act + SEMGES + factor(CCAA) + factor(TMUNR),family = binomial, data=cleaned_data)
model_trim_t3_4c_adj<-glm(low_birth_weight ~ factor(T3_4c) +SEXO + factor(cat_edad_mother) + EDADP + factor(cat_date_birth) + factor(cat_socioec) + cat_tasa_act + SEMGES + factor(CCAA) + factor(TMUNR),family = binomial, data=cleaned_data)
model_mean_4c_adj<-glm(low_birth_weight ~ factor(mean_4c) +SEXO + factor(cat_edad_mother) + EDADP + factor(cat_date_birth) + factor(cat_socioec) + cat_tasa_act + SEMGES + factor(CCAA) + factor(TMUNR),family = binomial, data=cleaned_data)

# association and adjusted GLM cont
model_trim_t1_4c_adj_glm<-glm(low_birth_weight ~ T1 +SEXO + factor(cat_edad_mother) + EDADP + factor(cat_date_birth) + factor(cat_socioec) + cat_tasa_act + SEMGES + factor(CCAA) + factor(TMUNR),family = binomial, data=cleaned_data)
model_trim_t2_4c_adj_glm<-glm(low_birth_weight ~ T2 +SEXO + factor(cat_edad_mother) + EDADP + factor(cat_date_birth) + factor(cat_socioec) + cat_tasa_act + SEMGES + factor(CCAA) + factor(TMUNR),family = binomial, data=cleaned_data)
model_trim_t3_4c_adj_glm<-glm(low_birth_weight ~ T3 +SEXO + factor(cat_edad_mother) + EDADP + factor(cat_date_birth) + factor(cat_socioec) + cat_tasa_act + SEMGES + factor(CCAA) + factor(TMUNR),family = binomial, data=cleaned_data)
model_mean_4c_adj_glm<-glm(low_birth_weight ~ mean +SEXO + factor(cat_edad_mother) + EDADP + factor(cat_date_birth) + factor(cat_socioec) + cat_tasa_act + SEMGES + factor(CCAA) + factor(TMUNR),family = binomial, data=cleaned_data)


# lm continous adjusted

model_trim_t1_adj_lm<-lm(PESON ~ T1 +SEXO + factor(cat_edad_mother) + EDADP + factor(cat_date_birth) + factor(cat_socioec) + cat_tasa_act + SEMGES + factor(CCAA) + factor(TMUNR),family = binomial, data=cleaned_data)
model_trim_t2_adj_lm<-lm(PESON ~ T2 +SEXO + factor(cat_edad_mother) + EDADP + factor(cat_date_birth) + factor(cat_socioec) + cat_tasa_act + SEMGES + factor(CCAA) + factor(TMUNR),family = binomial, data=cleaned_data)
model_trim_t3_adj_lm<-lm(PESON ~ T3 +SEXO + factor(cat_edad_mother) + EDADP + factor(cat_date_birth) + factor(cat_socioec) + cat_tasa_act + SEMGES + factor(CCAA) + factor(TMUNR),family = binomial, data=cleaned_data)
model_mean_adj_lm<-lm(PESON ~ mean +SEXO + factor(cat_edad_mother) + EDADP + factor(cat_date_birth) + factor(cat_socioec) + cat_tasa_act + SEMGES + factor(CCAA) + factor(TMUNR),family = binomial, data=cleaned_data)

# lm cat 
model_trim_t1_4c_lm<-lm(PESON ~ factor(T1_4c), data=cleaned_data)
model_trim_t2_4c_lm<-lm(PESON ~ factor(T2_4c), data=cleaned_data)
model_trim_t3_4c_lm<-lm(PESON ~ factor(T3_4c), data=cleaned_data)
model_trim_mean_4c_lm<-lm(PESON ~ factor(mean_4c), data=cleaned_data)


model_trim_t1_4c_adj_lm<-lm(PESON ~ factor(T1_4c) +SEXO + factor(cat_edad_mother) + EDADP + factor(cat_date_birth) + factor(cat_socioec) + cat_tasa_act + SEMGES + factor(CCAA) + factor(TMUNR),family = binomial, data=cleaned_data)
model_trim_t2_4c_adj_lm<-lm(PESON ~ factor(T2_4c) +SEXO + factor(cat_edad_mother) + EDADP + factor(cat_date_birth) + factor(cat_socioec) + cat_tasa_act + SEMGES + factor(CCAA) + factor(TMUNR),family = binomial, data=cleaned_data)
model_trim_t3_4c_adj_lm<-lm(PESON ~ factor(T3_4c) +SEXO + factor(cat_edad_mother) + EDADP + factor(cat_date_birth) + factor(cat_socioec) + cat_tasa_act + SEMGES + factor(CCAA) + factor(TMUNR),family = binomial, data=cleaned_data)
model_mean_4c_adj_lm<-lm(PESON ~ factor(mean_4c) +SEXO + factor(cat_edad_mother) + EDADP + factor(cat_date_birth) + factor(cat_socioec) + cat_tasa_act + SEMGES + factor(CCAA) + factor(TMUNR),family = binomial, data=cleaned_data)





# Example with your model
model <- model_mean_4c_adj

# Get coefficients
coef_summary <- summary(model)$coefficients 

# Odds Ratios and CIs
OR <- exp(coef(model))
CI <- exp(confint(model))   # profile likelihood CI
p_values <- coef_summary[, "Pr(>|z|)"]

# Build table
OR_table <- data.frame(
  Variable = names(OR),
  OR = OR,
  CI_lower = CI[,1],
  CI_upper = CI[,2],
  p_value = p_values
)


OR_table <- OR_table %>%
  mutate(
    OR = round(OR, 2),
    CI_lower = round(CI_lower, 2),
    CI_upper = round(CI_upper, 2),
    CI = paste0(CI_lower, " – ", CI_upper),
    p_value = formatC(p_value, format = "f", digits = 3)  # 3 decimal places
  ) %>%
  select(Variable, OR, CI, p_value)


library(flextable)
library(officer)

ft <- flextable(OR_table) %>%
  set_caption("Table 7. Trimester 3 Odds Ratios (OR) with 95% CI and p-values for Low Birth Weight") %>%
  autofit() %>%
  theme_vanilla()


read_docx() %>%
  body_add_flextable(ft) %>%
  print(target = "OR_Ajd_trim2.docx")




library(gtsummary)
library(flextable)

# --- Create regression tables for each trimester ---

install.packages("gtsummary")
library(gtsummary)


t1 <- tbl_regression(
  model_trim_T1_4c,
  exponentiate = TRUE
) %>%
  bold_p(t = 0.05) %>%
  italicize_labels()

t2 <- tbl_regression(
  model_trim_T2_4c,
  exponentiate = TRUE
) %>%
  bold_p(t = 0.05) %>%
  italicize_labels()

t3 <- tbl_regression(
  model_trim_T3_4c,
  exponentiate = TRUE
) %>%
  bold_p(t = 0.05) %>%
  italicize_labels()

tmean<-tbl_regression(
  model_mean_4c,
  exponentiate = TRUE
) %>%
  bold_p(t = 0.05) %>%
  italicize_labels()

# --- Merge into one table ---
combined_tbl <- tbl_merge(
  tbls = list(t1, t2, t3,tmean),
  tab_spanner = c("**Trimester 1**", "**Trimester 2**", "**Trimester 3**","**Mean**")
)

# --- Save as Word docx ---
combined_tbl %>%
  as_flex_table() %>%
  save_as_docx(path = "regression_results_all.docx")


# plot 
mod.gam.bp <- gam(low_birth_weight~ s(mean, bs="cr"), family=binomial, data=cleaned_data )
plot(mod.gam.bp, col="blue", xlab = "10 µg/m³  Ozone", ylab = "s(Mean Trimetser)", main = "LBW risk")
abline(h=0, col="grey")
mod_gam.M <- gam(PESON ~ s(mean, bs="cr"), data=cleaned_data) ## este es el bueno

plot(mod_gam.M, col="blue", xlab = "10 µg/m³ Ozones",ylab = "s(AADT)", main = "Birth Weight")#, shade=TRUE)#,seWithMean=TRUE,pch=19,1,cex=.55 ) #### este es el bueno
abline(h=0, col="grey")
