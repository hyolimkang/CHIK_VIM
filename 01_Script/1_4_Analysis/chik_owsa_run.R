
run_owsa_for_region <- function(region_name) {
  # Lookup abbreviation for observed/N object names
  region_abbr <- c(
    "Ceará" = "ce",
    "Pernambuco" = "pn",
    "Minas Gerais" = "mg",
    "Bahia" = "bh",
    "Paraíba" = "pa",
    "Rio Grande do Norte" = "rg",
    "Piauí" = "pi",
    "Alagoas" = "ag",
    "Tocantins" = "tc"
  )
  
  abbr <- region_abbr[[region_name]]
  observed <- get(paste0("observed_", abbr))
  N_region <- get(paste0("N_", abbr))
  
  # Step 1: Set regional parameters (use full region name)
  params_list <- make_params(
    region = region_name,
    bra_foi_state_summ = bra_foi_state_summ,
    posterior = get(paste0("posterior_", abbr)),
    lhs_sample_young = lhs_sample_young,
    hosp = hosp,
    fatal = fatal,
    nh_fatal = nh_fatal
  )
  
  baseline_params    <- params_list$baseline_params
  sensitivity_ranges <- params_list$sensitivity_ranges
  
  # Step 2: Run sensitivity analysis simulations
  sensitivity_results <- list()
  for (param in names(sensitivity_ranges)) {
    param_values <- sensitivity_ranges[[param]]
    impacts <- sapply(param_values, function(val) {
      run_simulation_with_params(
        target_age_list = target_age_list,
        param_overrides = setNames(list(val), param),
        baseline_params = baseline_params,
        posterior = get(paste0("posterior_", abbr)),
        bra_foi_state_summ = bra_foi_state_summ,
        age_groups = age_groups,
        N = N_region,
        region = region_name,  # full name used here
        observed = observed,
        age_gr_levels = age_gr_levels,
        lhs_sample_young = lhs_sample_young,
        lhs_old = lhs_old,
        le_sample = le_sample
      )
    })
    
    sensitivity_results[[param]] <- data.frame(
      Parameter = rep(param, length(param_values)),
      Value = param_values,
      Outcome = impacts
    )
  }
  
  # Step 3: Extract outcomes
  extract_owsa <- lapply(sensitivity_results, function(x) {
    list(
      foi = x$foi,
      lo = list(
        scenario_week_sums = x$Outcome.1$scenario_week_sums,
        pre_tot_cases_owsa = x$Outcome.1$pre_tot_cases_owsa
      ),
      pre_tot_cases = x$Outcome.1$pre_tot_cases,
      scenario_week_sums = x$Outcome.1$scenario_week_sums,
      mid_week_sums = x$Outcome.1$mid_week_sums,
      pre_cases_sum = x$Outcome.1$pre_cases_sum,
      hi = list(
        pre_tot_cases = x$Outcome.2$pre_tot_cases,
        scenario_week_sums = x$Outcome.2$scenario_week_sums,
        mid_week_sums = x$Outcome.2$mid_week_sums,
        pre_cases_sum = x$Outcome.2$pre_cases_sum,
        pre_tot_cases_owsa = x$Outcome.2$pre_tot_cases_owsa
      ),
      VE_block = x$VE_block,
      weekly_delivery_speed = x$weekly_delivery_speed
    )
  })
  
  # Step 4: Process into tornado_data
  tornado_data <- bind_rows(lapply(names(extract_owsa), function(param_name) {
    lo_vec <- extract_owsa[[param_name]]$lo$scenario_week_sums
    hi_vec <- extract_owsa[[param_name]]$hi$scenario_week_sums
    mid_vec <- extract_owsa[[param_name]]$mid_week_sums
    pre_cases_sum <- extract_owsa[[param_name]]$pre_cases_sum
    pre_tot_cases_owsa_lo <- extract_owsa[[param_name]]$lo$pre_tot_cases_owsa
    pre_tot_cases_owsa_hi <- extract_owsa[[param_name]]$hi$pre_tot_cases_owsa
    
    tibble(
      Region = region_name,
      Parameter = param_name,
      Scenario = paste0("Scenario_", seq_along(mid_vec)),
      Low = lo_vec,
      High = hi_vec,
      Mid = mid_vec,
      Pre_Cases_Sum = rep(pre_cases_sum, length(mid_vec)),
      Pre_Tot_Cases_OWSA_Lo = rep(pre_tot_cases_owsa_lo, length(mid_vec)),
      Pre_Tot_Cases_OWSA_Hi = rep(pre_tot_cases_owsa_hi, length(mid_vec))
    )
  })) %>%
    mutate(
      Mid_Impact = Pre_Cases_Sum - Mid,
      Lo_Impact = Pre_Tot_Cases_OWSA_Lo - Low,
      Hi_Impact = Pre_Tot_Cases_OWSA_Hi - High,
      Mid_Impact_pct = (Mid_Impact / Pre_Cases_Sum) * 100,
      Lo_Impact_pct = (Lo_Impact / Pre_Tot_Cases_OWSA_Lo) * 100,
      Hi_Impact_pct = (Hi_Impact / Pre_Tot_Cases_OWSA_Hi) * 100
    )
  
  return(tornado_data)
}

regions <- c("Ceará", "Pernambuco",  "Minas Gerais","Bahia","Paraíba",            
             "Rio Grande do Norte", "Piauí", "Alagoas", "Tocantins")

regions <- names(region_abbr)

tornado_results_all <- bind_rows(lapply(regions, run_owsa_for_region))

save("tornado_results_all", file = "00_Data/0_2_Processed/tornado_results_all.RData")

tornado_bar_pct_all <- tornado_results_all %>% group_by(Region) %>%
  select(Parameter, Scenario, Mid_Impact_pct, Lo_Impact_pct, Hi_Impact_pct) %>%
  pivot_longer(cols = c(Lo_Impact_pct, Hi_Impact_pct), names_to = "Type", values_to = "Value") %>%
  mutate(
    Direction = ifelse(Type == "Lo_Impact_pct", "Lower", "Upper"),
    Bar_start = Mid_Impact_pct,
    Bar_end = Value,
    Scenario = as.factor(Scenario)
  )%>%
  mutate(Parameter = dplyr::recode(Parameter,
                                   VE_block = "Vaccine efficacy",            # <-- Capital "VE"
                                   weekly_delivery_speed = "Weekly delivery speed",
                                   foi = "Force of infection",
                                   .default = Parameter)) %>%
  group_by(Region) %>%
  mutate(TotalImpact = sum(abs(Bar_end - Bar_start), na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(Region = fct_reorder(Region, TotalImpact, .desc = TRUE))


owsa_subnat <- 
ggplot(tornado_bar_pct_all, aes(y = fct_rev(Parameter))) +
  geom_segment(aes(x = Bar_start, xend = Bar_end, yend = fct_rev(Parameter),
                   color = Direction), linewidth = 5) +
  geom_point(aes(x = Bar_start), shape = 21, fill = "white",
             color = "black", size = 2.5) +
  facet_grid(Region ~ Scenario, scales = "free_y", switch = "y",
             labeller = labeller(
               Scenario = as_labeller(c(
                 "Scenario_1" = "<20 years only", 
                 "Scenario_2" = "20–59 years only", 
                 "Scenario_3" = ">60 years only"
               ))
             )) +
  scale_color_manual(values = c("Lower" = "#F8766D", "Upper" = "#00BFC4")) +
  labs(
    x = "% Cases averted",
    y = "Parameter",
    color = "Bounds"  
  ) +
  theme_light(base_size = 13) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.major.y = element_blank(),
    legend.position = "bottom",
    strip.text.y = element_text(size = 6)
  )

ggsave(filename = "02_Outputs/2_1_Figures/owsa_subnat.jpg", owsa_subnat, width = 10, height = 9, dpi = 1200)

# national level---------------------------------------------------------------

tornado_nat <- tornado_results_all %>%
  group_by(Scenario, Parameter) %>%
  summarise(
    Low = sum(Low, na.rm = TRUE),
    High = sum(High, na.rm = TRUE),
    Mid = sum(Mid, na.rm = TRUE),
    Pre_Cases_Sum = sum(Pre_Cases_Sum, na.rm = TRUE),
    Pre_Tot_Cases_OWSA_Lo = sum(Pre_Tot_Cases_OWSA_Lo, na.rm = TRUE),
    Pre_Tot_Cases_OWSA_Hi = sum(Pre_Tot_Cases_OWSA_Hi, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Mid_Impact = Pre_Cases_Sum - Mid,
    Lo_Impact = Pre_Tot_Cases_OWSA_Lo - Low,
    Hi_Impact = Pre_Tot_Cases_OWSA_Hi - High,
    
    Mid_Impact_pct = (Mid_Impact / Pre_Cases_Sum) * 100,
    Lo_Impact_pct  = (Lo_Impact / Pre_Tot_Cases_OWSA_Lo) * 100,
    Hi_Impact_pct  = (Hi_Impact / Pre_Tot_Cases_OWSA_Hi) * 100
  )

tornado_bar_pct_nat <- tornado_nat %>% 
  select(Parameter, Scenario, Mid_Impact_pct, Lo_Impact_pct, Hi_Impact_pct) %>%
  pivot_longer(cols = c(Lo_Impact_pct, Hi_Impact_pct), names_to = "Type", values_to = "Value") %>%
  mutate(
    Direction = ifelse(Type == "Lo_Impact_pct", "Lower", "Upper"),
    Bar_start = Mid_Impact_pct,
    Bar_end = Value,
    Scenario = as.factor(Scenario)
  ) %>%
  mutate(Parameter = dplyr::recode(Parameter,
                            VE_block = "Vaccine efficacy",            # <-- Capital "VE"
                            weekly_delivery_speed = "Weekly delivery speed",
                            foi = "Force of infection",
                            .default = Parameter)) 

owsa_all <- 
ggplot(tornado_bar_pct_nat, aes(y = fct_rev(Parameter))) +
  geom_segment(aes(x = Bar_start, xend = Bar_end, yend = fct_rev(Parameter),
                   color = Direction), linewidth = 5) +
  geom_point(aes(x = Bar_start), shape = 21, fill = "white",
             color = "black", size = 2.5) +
  facet_grid(. ~ Scenario, scales = "free_y", switch = "y",
             labeller = as_labeller(c(
               "Scenario_1" = "<20 years only", 
               "Scenario_2" = "20-59 years only", 
               "Scenario_3" = ">60 years only"
             )))  +
  scale_color_manual(values = c("Lower" = "#F8766D", "Upper" = "#00BFC4")) +
  labs(
    x = "% Cases averted",
    y = "Parameter",
    color = "Bounds"  
  ) +
  theme_light(base_size = 13) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.major.y = element_blank(),
    legend.position = "bottom"
  )

