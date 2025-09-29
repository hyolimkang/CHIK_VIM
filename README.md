# CHIK_VIM
Chikungunya outbreak response vaccine impact modelling

# 🦠 CHIK_VIM: Chikungunya Vaccine Impact Modelling

A mathematical modelling project to estimate the impact of chikungunya vaccines in Brazil using age-structured transmission models and outbreak scenarios.

![License: MIT](https://img.shields.io/badge/license-MIT-green)
![R](https://img.shields.io/badge/R-4.3.3-blue)

---

## 📁 Project Structure
```markdown

sim_functions_final.R
1. Ixchiq
2. Vimkunya

FOI = base_beta * Sum(I) / N assuming homogenous mixing across ages
target_ages = indexses all eligibile age groups (20 age groups <1 to >80 years)
total_supply = target_pop * total_coverage
how fast = weekly_dose_total = total_supply * weekly_delivery_speed

protection mechanism ---
1. infection blocking: successful vaccinations moved people out of S and adds to V
2. disease-only blocking: infections still occur but progression to symptomatic diseases is reduced by the coverage (VE * VC)
