
observed_cases_bahia <- as.matrix(bra_bh_22[,4])
observed_cases_bahia <- matrix(observed_cases_bahia, nrow = 18, ncol = 52, byrow = FALSE)

observed_cases_ce <- as.matrix(bra_ce_22[,4])
observed_cases_ce <- matrix(observed_cases_ce, nrow = 18, ncol = 52, byrow = FALSE)

observed_cases_mg <- as.matrix(bra_mg_22[,4])
observed_cases_mg <- matrix(observed_cases_mg, nrow = 18, ncol = 52, byrow = FALSE)

observed_cases_pn <- as.matrix(bra_pn_22[,4])
observed_cases_pn <- matrix(observed_cases_pn, nrow = 18, ncol = 52, byrow = FALSE)

observed_cases_pa <- as.matrix(bra_pa_22[,4])
observed_cases_pa <- matrix(observed_cases_pa, nrow = 18, ncol = 52, byrow = FALSE)

observed_cases_rg <- as.matrix(bra_rg_22[,4])
observed_cases_rg <- matrix(observed_cases_rg, nrow = 18, ncol = 52, byrow = FALSE)

observed_cases_pi <- as.matrix(bra_pi_22[,4])
observed_cases_pi <- matrix(observed_cases_pi, nrow = 18, ncol = 52, byrow = FALSE)

observed_cases_ag <- as.matrix(bra_ag_22[,4])
observed_cases_ag <- matrix(observed_cases_ag, nrow = 18, ncol = 52, byrow = FALSE)

observed_cases_tc <- as.matrix(bra_tc_22[,4])
observed_cases_tc <- matrix(observed_cases_tc, nrow = 18, ncol = 52, byrow = FALSE)

stan_model_age <- stan_model("01_Script/1_2_SIR_models/age_struc_vacc_bra_v4.1.stan")

# state specific foi
bra_foi <- allfoi %>% filter(country == "Brazil")
bra_foi_sf <- st_as_sf(bra_foi, coords = c("x", "y"), crs = 4326)
br_states <- st_read("00_Data/0_1_Raw/country_shape/gadm41_BRA_shp/gadm41_BRA_1.shp")
br_states <- st_transform(br_states, crs = st_crs(bra_foi_sf))
bra_foi_states <- st_join(bra_foi_sf, br_states, join = st_intersects)

bra_foi_state_summ <- bra_foi_states %>%
  group_by(NAME_1) %>%
  summarise(avg_foi = mean(foi_mid, na.rm = TRUE))%>%
  filter(!is.na(NAME_1))

# age groups
age_groups <- c(mean(0:4),
                mean(5:9),
                mean(10:14),
                mean(15:19),
                mean(20:24),
                mean(25:29),
                mean(30:34),
                mean(35:39),
                mean(40:44),
                mean(45:49),
                mean(50:54),
                mean(55:59),
                mean(60:64),
                mean(65:69),
                mean(70:74),
                mean(75:79),
                mean(80:84),
                mean(85:89)
)

## prevacc data 
stan_data_prevacc_bh <- list(
  T = 52,
  A = 18,
  N = as.vector(N_bahia), 
  observed_cases_by_age = round(observed_cases_bahia),
  r = rep(0, 18),
  delay = 53,
  VE_block = 0,
  vaccination_rate = 0,
  vaccine_target_age = rep(0,18),
  indexP = rep(0,52),
  eta = 0,
  prior_I0 = round(observed_cases_bahia[,1]),
  prior_sd_I0 = 100,
  sero = 1 - exp(- bra_foi_state_summ$avg_foi[bra_foi_state_summ$NAME_1 == "Bahia"] * age_groups)
)

stan_data_prevacc_ce <- list(
  T = 52,
  A = 18,
  N = as.vector(N_ceara), 
  observed_cases_by_age = round(observed_cases_ce),
  r = rep(0, 18),
  delay = 53,
  VE_block = 0,
  vaccination_rate = 0,
  vaccine_target_age = rep(0,18),
  indexP = rep(0,52),
  eta = 0,
  prior_I0 = round(observed_cases_ce[,1]),
  prior_sd_I0 = 100,
  sero = 1 - exp(- bra_foi_state_summ$avg_foi[bra_foi_state_summ$NAME_1 == "Ceará"] * age_groups)
)

stan_data_prevacc_mg <- list(
  T = 30,
  A = 18,
  N = as.vector(N_mg), 
  observed_cases_by_age = round(observed_cases_mg),
  r = rep(0, 18),
  delay = 53,
  VE_block = 0,
  vaccination_rate = 0,
  vaccine_target_age = rep(0,18),
  indexP = rep(0,30),
  eta = 0,
  prior_I0 = round(observed_cases_mg[,1]),
  prior_sd_I0 = 100,
  sero = 1 - exp(- bra_foi_state_summ$avg_foi[bra_foi_state_summ$NAME_1 == "Minas Gerais"] * age_groups)
)

stan_data_prevacc_pn <- list(
  T = 52,
  A = 18,
  N = as.vector(N_pemam), 
  observed_cases_by_age = round(observed_cases_pn),
  r = rep(0, 18),
  delay = 53,
  VE_block = 0,
  vaccination_rate = 0,
  vaccine_target_age = rep(0,18),
  indexP = rep(0,52),
  eta = 0,
  prior_I0 = round(observed_cases_pn[,1]),
  prior_sd_I0 = 100,
  sero = 1 - exp(- bra_foi_state_summ$avg_foi[bra_foi_state_summ$NAME_1 == "Pernambuco"] * age_groups)
)

stan_data_prevacc_pa <- list(
  T = 52,
  A = 18,
  N = as.vector(N_pa), 
  observed_cases_by_age = round(observed_cases_pa),
  r = rep(0, 18),
  delay = 53,
  VE_block = 0,
  vaccination_rate = 0,
  vaccine_target_age = rep(0,18),
  indexP = rep(0,52),
  eta = 0,
  prior_I0 = round(observed_cases_pa[,1]),
  prior_sd_I0 = 100,
  sero = 1 - exp(- bra_foi_state_summ$avg_foi[bra_foi_state_summ$NAME_1 == "Paraíba"] * age_groups)
)

stan_data_prevacc_rg <- list(
  T = 52,
  A = 18,
  N = as.vector(N_rg), 
  observed_cases_by_age = round(observed_cases_rg),
  r = rep(0, 18),
  delay = 53,
  VE_block = 0,
  vaccination_rate = 0,
  vaccine_target_age = rep(0,18),
  indexP = rep(0,52),
  eta = 0,
  prior_I0 = round(observed_cases_rg[,1]),
  prior_sd_I0 = 100,
  sero = 1 - exp(- bra_foi_state_summ$avg_foi[bra_foi_state_summ$NAME_1 == "Rio Grande do Norte"] * age_groups)
)

stan_data_prevacc_pi <- list(
  T = 52,
  A = 18,
  N = as.vector(N_pi), 
  observed_cases_by_age = round(observed_cases_pi),
  r = rep(0, 18),
  delay = 53,
  VE_block = 0,
  vaccination_rate = 0,
  vaccine_target_age = rep(0,18),
  indexP = rep(0,52),
  eta = 0,
  prior_I0 = round(observed_cases_pi[,1]),
  prior_sd_I0 = 100,
  sero = 1 - exp(- bra_foi_state_summ$avg_foi[bra_foi_state_summ$NAME_1 == "Piauí"] * age_groups)
)

stan_data_prevacc_ag <- list(
  T = 52,
  A = 18,
  N = as.vector(N_ag), 
  observed_cases_by_age = round(observed_cases_ag),
  r = rep(0, 18),
  delay = 53,
  VE_block = 0,
  vaccination_rate = 0,
  vaccine_target_age = rep(0,18),
  indexP = rep(0,52),
  eta = 0,
  prior_I0 = round(observed_cases_ag[,1]),
  prior_sd_I0 = 100,
  sero = 1 - exp(- bra_foi_state_summ$avg_foi[bra_foi_state_summ$NAME_1 == "Alagoas"] * age_groups)
)

stan_data_prevacc_tc <- list(
  T = 52,
  A = 18,
  N = as.vector(N_tc), 
  observed_cases_by_age = round(observed_cases_tc),
  r = rep(0, 18),
  delay = 53,
  VE_block = 0,
  vaccination_rate = 0,
  vaccine_target_age = rep(0,18),
  indexP = rep(0,52),
  eta = 0,
  prior_I0 = round(observed_cases_tc[,1]),
  prior_sd_I0 = 100,
  sero = 1 - exp(- bra_foi_state_summ$avg_foi[bra_foi_state_summ$NAME_1 == "Tocantins"] * age_groups)
)

### fitting 
fit_prevacc_bh <- sampling(
  object = stan_model_age,
  data = stan_data_prevacc_bh,
  iter = 1000,            # Reduced from 4000
  chains = 1,             # Reduced from 4
  warmup = 100,           # Specify warmup period
  thin = 2,               # Add thinning
  seed = 123,
  control = list(
    adapt_delta = 0.8,    # Start with lower adaptation
    max_treedepth = 15    # Reduced tree depth
  )
)

fit_prevacc_ce <- sampling(
  object = stan_model_age,
  data = stan_data_prevacc_ce,
  iter = 1000,            # Reduced from 4000
  chains = 1,             # Reduced from 4
  warmup = 100,           # Specify warmup period
  thin = 2,               # Add thinning
  seed = 123,
  control = list(
    adapt_delta = 0.8,    # Start with lower adaptation
    max_treedepth = 15    # Reduced tree depth
  )
)

fit_prevacc_mg <- sampling(
  object = stan_model_age,
  data = stan_data_prevacc_mg,
  iter = 1000,            # Reduced from 4000
  chains = 1,             # Reduced from 4
  warmup = 100,           # Specify warmup period
  thin = 2,               # Add thinning
  seed = 123,
  control = list(
    adapt_delta = 0.8,    # Start with lower adaptation
    max_treedepth = 15    # Reduced tree depth
  )
)

fit_prevacc_pn <- sampling(
  object = stan_model_age,
  data = stan_data_prevacc_pn,
  iter = 1000,            # Reduced from 4000
  chains = 1,             # Reduced from 4
  warmup = 100,           # Specify warmup period
  thin = 2,               # Add thinning
  seed = 123,
  control = list(
    adapt_delta = 0.8,    # Start with lower adaptation
    max_treedepth = 15    # Reduced tree depth
  )
)

fit_prevacc_pa <- sampling(
  object = stan_model_age,
  data = stan_data_prevacc_pa,
  iter = 1000,            # Reduced from 4000
  chains = 1,             # Reduced from 4
  warmup = 100,           # Specify warmup period
  thin = 2,               # Add thinning
  seed = 123,
  control = list(
    adapt_delta = 0.8,    # Start with lower adaptation
    max_treedepth = 15    # Reduced tree depth
  )
)

fit_prevacc_rg <- sampling(
  object = stan_model_age,
  data = stan_data_prevacc_rg,
  iter = 1000,            # Reduced from 4000
  chains = 1,             # Reduced from 4
  warmup = 100,           # Specify warmup period
  thin = 2,               # Add thinning
  seed = 123,
  control = list(
    adapt_delta = 0.8,    # Start with lower adaptation
    max_treedepth = 15    # Reduced tree depth
  )
)

fit_prevacc_pi <- sampling(
  object = stan_model_age,
  data = stan_data_prevacc_pi,
  iter = 1000,            # Reduced from 4000
  chains = 1,             # Reduced from 4
  warmup = 100,           # Specify warmup period
  thin = 2,               # Add thinning
  seed = 123,
  control = list(
    adapt_delta = 0.8,    # Start with lower adaptation
    max_treedepth = 15    # Reduced tree depth
  )
)

fit_prevacc_ag <- sampling(
  object = stan_model_age,
  data = stan_data_prevacc_ag,
  iter = 1000,            # Reduced from 4000
  chains = 1,             # Reduced from 4
  warmup = 100,           # Specify warmup period
  thin = 2,               # Add thinning
  seed = 123,
  control = list(
    adapt_delta = 0.8,    # Start with lower adaptation
    max_treedepth = 15    # Reduced tree depth
  )
)

fit_prevacc_tc <- sampling(
  object = stan_model_age,
  data = stan_data_prevacc_tc,
  iter = 1000,            # Reduced from 4000
  chains = 1,             # Reduced from 4
  warmup = 100,           # Specify warmup period
  thin = 2,               # Add thinning
  seed = 123,
  control = list(
    adapt_delta = 0.8,    # Start with lower adaptation
    max_treedepth = 15    # Reduced tree depth
  )
)

save(fit_prevacc_bh, file = "00_Data/0_2_Processed/fit_prevacc_bh.RData")
save(fit_prevacc_ce, file = "00_Data/0_2_Processed/fit_prevacc_ce.RData")
save(fit_prevacc_mg, file = "00_Data/0_2_Processed/fit_prevacc_mg.RData")
save(fit_prevacc_pn, file = "00_Data/0_2_Processed/fit_prevacc_pn.RData")
save(fit_prevacc_pa, file = "00_Data/0_2_Processed/fit_prevacc_pa.RData")
save(fit_prevacc_rg, file = "00_Data/0_2_Processed/fit_prevacc_rg.RData")
save(fit_prevacc_pi, file = "00_Data/0_2_Processed/fit_prevacc_pi.RData")
save(fit_prevacc_ag, file = "00_Data/0_2_Processed/fit_prevacc_ag.RData")
save(fit_prevacc_tc, file = "00_Data/0_2_Processed/fit_prevacc_tc.RData")

## post processing 
list_bh <- create_summary_df(fit_prevacc_bh,
                             bra_sum_bh,
                             region = "Bahia")

list_ce <- create_summary_df(fit_prevacc_ce,
                             bra_sum_ce,
                             region = "Ceará")

list_mg <- create_summary_df(fit_prevacc_mg,
                             bra_sum_mg,
                             region = "Minas Gerais")

list_pn <- create_summary_df(fit_prevacc_pn,
                             bra_sum_pn,
                             region = "Pernambuco")

list_pa <- create_summary_df(fit_prevacc_pa,
                             bra_sum_pa,
                             region = "Paraíba")

list_rg <- create_summary_df(fit_prevacc_rg,
                             bra_sum_rg,
                             region = "Rio Grande do Norte")

list_pi <- create_summary_df(fit_prevacc_pi,
                             bra_sum_pi,
                             region = "Piauí")

list_ag <- create_summary_df(fit_prevacc_ag,
                             bra_sum_ag,
                             region = "Alagoas")

list_tc <- create_summary_df(fit_prevacc_tc,
                             bra_sum_tc,
                             region = "Tocantins")

## bahia
df_bh <- list_bh$df_out

df_bh_summ <- list_bh$df_summ

observed_bh <- list_bh$observed

age_strat_gg(df_bh)

i <- overall_fit_gg(observed_bh, df_bh_summ)

## ceara
df_ce <- list_ce$df_out

df_ce_summ <- list_ce$df_summ

observed_ce <- list_ce$observed

age_strat_gg(df_ce)

a <- overall_fit_gg(observed_ce, df_ce_summ)

combined_bh <- age_struc_fit(fit_prevacc_bh,
                             model_age_bins,
                             bra_bh_22,
                             bra_sum_bh)
agestrat_ui_gg(combined_bh)

## minas gerais
df_mg <- list_mg$df_out

df_mg_summ <- list_mg$df_summ

observed_mg <- list_mg$observed

age_strat_gg(df_mg)

b <- overall_fit_gg(observed_mg, df_mg_summ)

combined_mg <- age_struc_fit(fit_prevacc_mg,
                             model_age_bins,
                             bra_mg_22,
                             bra_sum_mg)
agestrat_ui_gg(combined_mg)

## pernambuco
df_pn <- list_pn$df_out

df_pn_summ <- list_pn$df_summ

observed_pn <- list_pn$observed

age_strat_gg(df_pn)

c <- overall_fit_gg(observed_pn, df_pn_summ)

combined_pn <- age_struc_fit(fit_prevacc_pn,
                             model_age_bins,
                             bra_pn_22,
                             bra_sum_pn)
agestrat_ui_gg(combined_pn)


## paraiba
df_pa <- list_pa$df_out

df_pa_summ <- list_pa$df_summ

observed_pa <- list_pa$observed

age_strat_gg(df_pa)

d <- overall_fit_gg(observed_pa, df_pa_summ)

combined_pa <- age_struc_fit(fit_prevacc_pa,
                             model_age_bins,
                             bra_pa_22,
                             bra_sum_pa)
agestrat_ui_gg(combined_pa)


## rio grande norte
df_rg <- list_rg$df_out

df_rg_summ <- list_rg$df_summ

observed_rg <- list_rg$observed

age_strat_gg(df_rg)

e <- overall_fit_gg(observed_rg, df_rg_summ)

combined_rg <- age_struc_fit(fit_prevacc_rg,
                             model_age_bins,
                             bra_rg_22,
                             bra_sum_rg)
agestrat_ui_gg(combined_rg)

## Piuai
df_pi <- list_pi$df_out

df_pi_summ <- list_pi$df_summ

observed_pi <- list_pi$observed

age_strat_gg(df_pi)

f <- overall_fit_gg(observed_pi, df_pi_summ)

combined_pi <- age_struc_fit(fit_prevacc_pi,
                             model_age_bins,
                             bra_pi_22,
                             bra_sum_pi)
agestrat_ui_gg(combined_pi)

## Alagoas
df_ag <- list_ag$df_out

df_ag_summ <- list_ag$df_summ

observed_ag <- list_ag$observed

age_strat_gg(df_ag)

g <- overall_fit_gg(observed_ag, df_ag_summ)

combined_ag <- age_struc_fit(fit_prevacc_ag,
                             model_age_bins,
                             bra_ag_22,
                             bra_sum_ag)
agestrat_ui_gg(combined_ag)

## tocantins
df_tc <- list_tc$df_out

df_tc_summ <- list_tc$df_summ

observed_tc <- list_tc$observed

age_strat_gg(df_tc)

h <- overall_fit_gg(observed_tc, df_tc_summ)

combined_tc <- age_struc_fit(fit_prevacc_tc,
                             model_age_bins,
                             bra_tc_22,
                             bra_sum_tc)
agestrat_ui_gg(combined_tc)

## overall graph 
observed_ce$region <- "Ceará" 
observed_bh$region <- "Bahia" 
observed_ag$region <- "Alagoas" 
observed_mg$region <- "Minas Gerais" 
observed_tc$region <- "Tocantins"
observed_pa$region <- "Paraíba" 
observed_pi$region <- "Piauí"
observed_pn$region <- "Pernambuco"
observed_rg$region <- "Rio Grande do Norte"

observed_all <- bind_rows(observed_ce, observed_bh, observed_ag,
                          observed_mg, observed_tc, observed_pa,
                          observed_pi, observed_pn, observed_rg)

pred_all <- bind_rows(
  df_ag_summ, df_bh_summ, df_ce_summ,
  df_mg_summ, df_pa_summ, df_pi_summ,
  df_rg_summ, df_pn_summ, df_tc_summ
)

p <- 

ggplot()+
  geom_point(data = observed_all, aes(x = Week, y = Observed), size = 0.8)+
  facet_wrap(~region)+
  # Predicted line
  geom_line(data = pred_all, aes(x = Week, y = Median, color = Type), size = 1) +
  
  # Prediction interval (ribbon)
  geom_ribbon(data = pred_all, aes(x = Week, ymin = Lower, ymax = Upper, fill = Type), 
              alpha = 0.2)+
  theme_light()+
  ylab("Predicted and observed reported symptomatic cases")+
  scale_y_continuous(labels = comma)

ggsave(filename = "02_Outputs/2_1_Figures/fig_2022_fit.jpg", p, width = 10, height = 8, dpi = 1200)


save("observed_ce", "observed_bh", "observed_ag", 
     "observed_mg", "observed_tc", "observed_pa",
     "observed_pi", "observed_pn", "observed_rg",
     file = "00_Data/0_2_Processed/observed_2022.RData")


