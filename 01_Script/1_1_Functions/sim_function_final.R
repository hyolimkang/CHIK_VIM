sirv_sim_coverageSwitch <- function(
    # Model dimensions
  T, A,          # T = number of weeks, A = number of age groups
  N,             # population by age
  r,             # aging rates
  
  # Epidemic parameters
  base_beta,     # length T
  I0_draw,       # initial infected by age
  R0,            # initial underlying immunity fraction by age
  rho,           # detection probability
  gamma,         # recovery rate
  
  # Vaccination parameters
  delay,                 # the week index to START physical vaccination
  coverage_threshold,    # e.g., 0.7 => 70%
  target_age,            # 0/1 vector of length A
  total_coverage,        # fraction of total population eventually vaccinated
  weekly_delivery_speed, # fraction of total supply used per week
  VE_block,              # vaccine efficacy
  
  # NEW: toggle coverage switch
  use_coverage_switch = FALSE  # set FALSE to disable Phase 2 reallocation
){
  # ------------------------------------------------------------------------
  # 1) Initialize compartments and data structures
  # ------------------------------------------------------------------------
  
  # S, I, R, V each dimension [A x T]
  S  <- matrix(0, nrow = A, ncol = T)
  I  <- matrix(0, nrow = A, ncol = T)
  R_ <- matrix(0, nrow = A, ncol = T)
  V  <- matrix(0, nrow = A, ncol = T)
  
  # Delayed vaccination storage: doses administered at time t that become effective 2 weeks later.
  vacc_delayed <- matrix(0, nrow = A, ncol = T)
  
  # For reporting
  age_stratified_cases <- matrix(0, nrow = A, ncol = T)  # weekly detected cases by age
  raw_allocation_age   <- matrix(0, nrow = A, ncol = T)  # how many doses allocated each week
  coverage_frac        <- matrix(0, nrow = A, ncol = T)  # fraction of each age group in V
  coverage_target      <- numeric(T)                     # coverage fraction in target group
  phi_vec              <- numeric(T)                     # force of infection each week
  
  # Identify target group
  target_indices <- which(target_age == 1)
  target_pop     <- sum(N[target_indices])
  no_target      <- (length(target_indices) == 0 || target_pop == 0)
  
  # Total supply (when we define coverage as a fraction of total pop, or fraction of targeted pop the definition changes)
  total_supply <- target_pop  * total_coverage # target pop or sum(N)
  
  # ------------------------------------------------------------------------
  # 1A) Phase 1 supply: target only
  # ------------------------------------------------------------------------
  
  weekly_dose_total <- total_supply * weekly_delivery_speed
  
  # In Phase 1, only target ages get a share of weekly_dose_total
  weekly_supply_age <- numeric(A)
  if(!no_target) {
    for(a in seq_len(A)) {
      if(target_age[a] == 1) {
        weekly_supply_age[a] <- weekly_dose_total * (N[a] / target_pop)
      } else {
        weekly_supply_age[a] <- 0
      }
    }
  } else {
    weekly_supply_age <- rep(0, A)
  }
  
  # total_available_supply_age: maximum doses that can be allocated to each age a
  total_available_supply_age <- numeric(A)
  for(a in seq_len(A)) {
    if(target_age[a] == 1 && !no_target) {
      total_available_supply_age[a] <- total_supply * (N[a] / target_pop)
    } else {
      total_available_supply_age[a] <- 0
    }
  }
  
  # Track usage vs leftover supply
  total_supply_age  <- numeric(A) # how many doses have actually been used for age [a]
  unused_supply_age <- total_available_supply_age
  
  # coverage_switch becomes TRUE once coverage_threshold is reached in the target group
  coverage_switch <- FALSE
  switch_day      <- NA
  
  # ------------------------------------------------------------------------
  # 1B) Vaccination duration tracking
  # ------------------------------------------------------------------------
  vacc_start_week <- rep(NA, A)
  vacc_end_week   <- rep(NA, A)
  
  # ------------------------------------------------------------------------
  # 1C) Set initial conditions at t = 1
  # ------------------------------------------------------------------------
  for(a in seq_len(A)) {
    R_[a, 1] <- R0[a] * N[a]
    I[a, 1]  <- I0_draw[a]
    S[a, 1]  <- N[a] - I[a,1] - R_[a,1]
    V[a, 1]  <- 0
    
    age_stratified_cases[a, 1] <- rho * I[a,1]
    coverage_frac[a, 1]       <- V[a, 1] / N[a]
  }
  
  coverage_target[1] <- 0
  
  # ------------------------------------------------------------------------
  # 2) Main loop, t = 2..T
  # ------------------------------------------------------------------------
  for(t_ in 2:T) {
    
    # ============================================================
    # New Step 1: Apply 2-week vaccination effect BEFORE computing FOI
    # ============================================================
    # Vaccines administered at time (t_ - 2) become effective now.
    if(t_ - 2 >= 1) {
      for(a in seq_len(A)) {
        newly_effective <- VE_block * vacc_delayed[a, t_ - 2]
        # Move at most what is available in susceptibles from last time step.
        actual_move <- min(S[a, t_ - 1], newly_effective)
        S[a, t_ - 1] <- S[a, t_ - 1] - actual_move
        V[a, t_ - 1] <- V[a, t_ - 1] + actual_move
      }
    }
    
    # ============================================================
    # 2A) Update coverage in target group (cumulative doses administered)
    # ============================================================
    if(!no_target) {
      target_vacc_cumulative <- sum(raw_allocation_age[target_indices, 1:(t_ - 1)])
      coverage_now <- target_vacc_cumulative / target_pop
      
      if(use_coverage_switch) {
        if(!coverage_switch && coverage_now >= coverage_threshold) {
          coverage_switch <- TRUE
          switch_day <- t_
          
          # Compute the leftover supply from phase 1:
          leftover_total <- sum(unused_supply_age)
          
          # Reallocate leftover supply: now, for the non-target groups only.
          non_target_indices <- which(target_age == 0)
          total_non_target_pop <- sum(N[non_target_indices])
          
          if(total_non_target_pop > 0 && leftover_total > 0) {
            # For each non-target age group, set its maximum available supply
            for(a in non_target_indices) {
              total_available_supply_age[a] <- leftover_total * (N[a] / total_non_target_pop)
              unused_supply_age[a] <- total_available_supply_age[a]  # reset usage for these groups
            }
            
            # Now, recalculate the weekly dose total for Phase 2 based on leftover supply:
            new_weekly_dose_total <- leftover_total * weekly_delivery_speed
            
            # Set the weekly supply for non-target groups only:
            for(a in seq_len(A)) {
              if(a %in% non_target_indices) {
                weekly_supply_age[a] <- new_weekly_dose_total * (N[a] / total_non_target_pop)
              } else {
                weekly_supply_age[a] <- 0  # no further allocation to already-targeted groups
              }
            }
          }
        }
      }
      
    } else {
      coverage_now <- 0
    }
    
    # ============================================================
    # 2B) Compute Force of Infection (using updated state from t_-1)
    # ============================================================
    infected_prev <- sum(I[, t_ - 1])
    region_pop <- sum(N)
    phi <- base_beta[t_ - 1] * (infected_prev / region_pop)
    phi_vec[t_] <- phi
    
    # ============================================================
    # 2C) Vaccination allocation for week t_
    # ============================================================
    for(a in seq_len(A)) {
      # Vaccination starts exactly at t = delay.
      if(t_ < delay) {
        doses_planned <- 0
      } else {
        doses_planned <- weekly_supply_age[a]
      }
      
      susceptibles_now <- S[a, t_ - 1]
      remaining_supply <- unused_supply_age[a]
      
      feasible <- min(doses_planned, remaining_supply, susceptibles_now)
      feasible <- max(feasible, 0)
      
      raw_allocation_age[a, t_] <- feasible
      
      if(feasible > 0) {
        if(is.na(vacc_start_week[a])) {
          vacc_start_week[a] <- t_
        }
        vacc_end_week[a] <- t_
        
        total_supply_age[a] <- total_supply_age[a] + feasible
        unused_supply_age[a] <- total_available_supply_age[a] - total_supply_age[a]
        
        # Record the doses so that they become effective 2 weeks later.
        vacc_delayed[a, t_] <- feasible
      }
    }
    
    # ============================================================
    # 2D) S-I-R transitions for week t_ (using state from t_-1)
    # ============================================================
    for(a in seq_len(A)) {
      if(a == 1) {
        S[a, t_] <- S[a, t_ - 1] - phi * S[a, t_ - 1] - r[a] * S[a, t_ - 1]
        I[a, t_] <- I[a, t_ - 1] + phi * S[a, t_ - 1] - gamma * I[a, t_ - 1] - r[a] * I[a, t_ - 1]
        R_[a, t_] <- R_[a, t_ - 1] + gamma * I[a, t_ - 1] - r[a] * R_[a, t_ - 1]
        V[a, t_] <- V[a, t_ - 1] - r[a] * V[a, t_ - 1]
      } else {
        S[a, t_] <- S[a, t_ - 1] - phi * S[a, t_ - 1] - r[a] * S[a, t_ - 1] + r[a - 1] * S[a - 1, t_ - 1]
        I[a, t_] <- I[a, t_ - 1] + phi * S[a, t_ - 1] - gamma * I[a, t_ - 1] - r[a] * I[a, t_ - 1] + r[a - 1] * I[a - 1, t_ - 1]
        R_[a, t_] <- R_[a, t_ - 1] + gamma * I[a, t_ - 1] - r[a] * R_[a, t_ - 1] + r[a - 1] * R_[a - 1, t_ - 1]
        V[a, t_] <- V[a, t_ - 1] - r[a] * V[a, t_ - 1] + r[a - 1] * V[a - 1, t_ - 1]
      }
      
      # Ensure no compartment goes negative.
      S[a, t_] <- pmax(0, S[a, t_])
      I[a, t_] <- pmax(0, I[a, t_])
      R_[a, t_] <- pmax(0, R_[a, t_])
      V[a, t_] <- pmax(0, V[a, t_])
    }
    
    # ============================================================
    # 2F) Post-update: update coverage and recorded cases for week t_
    # ============================================================
    for(a in seq_len(A)) {
      coverage_frac[a, t_] <- V[a, t_] / N[a]
      age_stratified_cases[a, t_] <- rho * I[a, t_]
    }
    
    if(!no_target) {
      target_vacc_cumulative <- sum(raw_allocation_age[target_indices, 1:t_])
      coverage_target[t_] <- target_vacc_cumulative / target_pop
    } else {
      coverage_target[t_] <- 0
    }
    
  }  # end of main loop
  
  # ------------------------------------------------------------------------
  # 3) Compute total vaccination duration by age
  # ------------------------------------------------------------------------
  real_vaccine_duration_age <- numeric(A)
  for(a in seq_len(A)) {
    if(!is.na(vacc_start_week[a]) && !is.na(vacc_end_week[a])) {
      real_vaccine_duration_age[a] <- vacc_end_week[a] - vacc_start_week[a] + 1
    } else {
      real_vaccine_duration_age[a] <- 0
    }
  }
  
  # ------------------------------------------------------------------------
  # 4) Return results
  # ------------------------------------------------------------------------
  return(list(
    # compartments
    S = S,
    I = I,
    R = R_,
    V = V,
    
    # delayed vaccination (to become effective 2 weeks later)
    vacc_delayed = vacc_delayed,
    
    # supply usage
    total_supply               = total_supply,
    total_available_supply_age = total_available_supply_age,
    total_supply_age           = total_supply_age,
    unused_supply_age          = unused_supply_age,
    
    # coverage switch info
    coverage_switch = coverage_switch,
    switch_day      = switch_day,
    
    # epidemiological reporting
    age_stratified_cases = age_stratified_cases,
    coverage_frac        = coverage_frac,
    coverage_target      = coverage_target,
    phi                  = phi_vec,
    
    # vaccination window info
    vacc_start_week           = vacc_start_week,
    vacc_end_week             = vacc_end_week,
    real_vaccine_duration_age = real_vaccine_duration_age,
    
    # raw allocation by age each week
    raw_allocation_age = raw_allocation_age
  ))
}
