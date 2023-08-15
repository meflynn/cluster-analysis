library(tidyverse)
library(brms)
library(here)


data <- readxl::read_xlsx(here("data/anchor regions analysis.xlsx")) |>
  mutate(LatentGroup = as.numeric(NA),
         across(c("totpop_19", "real_gdp_21", "net_mig"),
                ~arm::rescale(.x)))

bf1 <- bf(totpop_19 ~ 0 + mi(LatentGroup),
   decomp = "QR")

bf2 <- bf(real_gdp_21 ~ 0 + mi(LatentGroup),
          decomp = "QR")

bf3 <- bf(net_mig ~ 0 + mi(LatentGroup),
          decomp = "QR")

bf4 <- bf(LatentGroup | mi() ~ 0)

get_prior(formula = bf1 + bf2 + bf3 + bf4,
          family = gaussian(),
          data = data)

PRIORS <- c(prior(constant(1), coef = miLatentGroup, resp = totpop19),
            prior(normal(1, 1), coef = miLatentGroup, resp = realgdp21),
            prior(normal(1, 1), coef = miLatentGroup, resp = netmig))

test_mod <- brm(
  bf1 + bf2 + bf3 + bf4 + set_rescor(FALSE),
  data = data,
  prior = PRIORS,
  family = gaussian(),
  iter = 4000,
  warmup = 2000,
  thin = 1,
  cores = 4,
  chains = 4,
  threads = threading(2),
  backend = "cmdstanr")


out <-  tidybayes::gather_draws(test_mod,
                                            "LatentGroup")
