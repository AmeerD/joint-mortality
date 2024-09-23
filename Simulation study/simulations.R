library(tidyverse)
library(tidybayes)
library(rstan)
library(magrittr)
library(usmap) #use fips_info to decode fips codes
library(bayesplot)

options(mc.cores=4)
rstan_options(auto_write=TRUE)

#Input data for the model should contain the following columns:
# - state
# - county and/or county code (code is preferable to prevent duplicates across states)
# - year
# - age category
# - deaths
# - population
load("pcs.rda")
su <- pcs_tot$u
sing <- pcs_tot$d
pcs <- pcs_tot$v

#Specify the number of principal components
n_pcs <- 2

set.seed(1234)
dat.sim <- expand_grid(state="SIM", 
                       code=1:25, 
                       year=2001:2010, 
                       race=letters[1:5],
                       ages=dat %>% select(ages) %>% distinct %>% pull) %>% 
  mutate(prop = case_when(race == "a" ~ 0.5, race == "b" ~ 0.2, TRUE ~ 0.1)) %>% 
  mutate(pop = 100000*code*prop*(1+0.01*(year-2001))) %>% ########100000
  select(-prop) %>%
  left_join(
    dat %>% filter(code == "06037") %>% group_by(state, code, year, ages) %>% 
      summarise(deaths=sum(deaths), pop=sum(pop)) %>% mutate(prop = pop/sum(pop)) %>%
      ungroup %>% select(ages, year, prop)
  ) %>%
  mutate(prop = pmax(rnorm(nrow(.), mean=prop, sd=0.005), 0.01)) %>%
  mutate(pop2 = pop*prop) %>%
  select(state, code, year, race, ages, pop=pop2)

base_mort <- mean((su %*% diag(sing))[,1]) * pcs[, 1]

corr1 <- array(0, dim=c(5,5,10))
corr1[,,1] <- corr1[,,2] <- corr1[,,3] <- diag(rep(1,5))
tempmat <- matrix(0.5, nrow=5, ncol=5) 
diag(tempmat) <- rep(1,5)
corr1[,,4] <- corr1[,,5] <- corr1[,,6] <- tempmat
tempmat <- matrix(c(1,0.8,0.6,0.9,0.5,
                    0.8,1,0.9,0.9,0.5,
                    0.6,0.9,1,0.8,0.5,
                    0.9,0.9,0.8,1,0.5,
                    0.5,0.5,0.5,0.5,1), nrow=5, ncol=5) 
corr1[,,7] <- corr1[,,8] <- corr1[,,9] <- corr1[,,10] <- tempmat

var1 <- rep(0.1, 5)
cov1 <- array(0, dim=c(5,5,10))
for (i in 1:10) {
  cov1[,,i] <- diag(var1) %*% corr1[,,i] %*% diag(var1)
}



corr2 <- array(0, dim=c(5,5,10))
corr2[,,1] <- corr2[,,2] <- corr2[,,3] <- diag(rep(1,5))
tempmat <- matrix(0.5, nrow=5, ncol=5) 
diag(tempmat) <- rep(1,5)
corr2[,,4] <- corr2[,,5] <- corr2[,,6] <- tempmat
tempmat <- matrix(c(1,0.9,0,0,0,
                    0.9,1,0,0,0,
                    0,0,1,0.5,0.5,
                    0,0,0.5,1,0.5,
                    0,0,0.5,0.5,1), nrow=5, ncol=5) 
corr2[,,7] <- corr2[,,8] <- corr2[,,9] <- corr2[,,10] <- tempmat

var2 <- rep(0.5, 5)
cov2 <- array(0, dim=c(5,5,10))
for (i in 1:10) {
  cov2[,,i] <- diag(var2) %*% corr2[,,i] %*% diag(var2)
}

basemat <- expand_grid(state="SIM",
                       code=pull(distinct(select(dat.sim, code))),
                       year=pull(distinct(select(dat.sim, year))))
stretch <- matrix(0, nrow=nrow(basemat), ncol=5)
bump <- matrix(0, nrow=nrow(basemat), ncol=5)
for (i in 1:nrow(basemat)) {
  stretch[i,] <- MASS::mvrnorm(1, mu=rep(1,5), Sigma=cov1[,,basemat$year[i]-2000])
  bump[i,] <- MASS::mvrnorm(1, mu=rep(0,5), Sigma=cov2[,,basemat$year[i]-2000])
}

logmx <- bind_cols(
  basemat,
  stretch %>% set_colnames(pull(distinct(select(dat.sim, race)))) %>% as.tibble
) %>%
  pivot_longer(cols=c(-state, -code, -year), names_to="race", values_to="stretch") %>%
  left_join(
    bind_cols(
      basemat,
      bump %>% set_colnames(pull(distinct(select(dat.sim, race)))) %>% as.tibble
    ) %>%
      pivot_longer(cols=c(-state, -code, -year), names_to="race", values_to="bump")
  ) %>%
  expand_grid(ages = sort(pull(distinct(select(dat.sim, ages))))) %>%
  left_join(tibble(ages = sort(pull(distinct(select(dat.sim, ages)))),
                   base = base_mort,
                   acc = c(0,0,0,0,0.25,0.5,1,1,0.5,0.25,0,0,0,0,0,0,0,0,0))) %>%
  mutate(lmx = stretch * base + bump*acc)

dat.sim2 <- dat.sim %>%
  inner_join(logmx) %>%
  mutate(deaths = rpois(nrow(.), pop*exp(lmx))) %>%
  select(-lmx)

dat.sim2 %>% 
  filter(year == 2001) %>% 
  ggplot(aes(x=ages, y=log(deaths/pop), colour=race)) + 
  geom_point() + 
  facet_wrap(~code)

#Set the baseyear
baseyear <- dat.sim2 %>% select(year) %>% pull %>% min

pcs_adj <- rep(1, n_pcs)
pcs_sim <- matrix(c(base_mort, acc), ncol=2)

#Build the model input list
stan_data <- dat.sim2 %>%
  mutate(year = year - baseyear + 1,
         code = as.factor(code)) %>% 
  compose_data(state = x_at_y(state, code)) %>%
  within({
    n_year <- length(unique(year))
    pcs <- pcs_sim 
    #pcs <- pcs[,1:n_pcs] %*% diag(pcs_adj)
    n_pcs <- n_pcs
  }) 

model_file <- "joint_model_multi.stan"

model <- stan_model(model_file, model_name = "base mortality")

print(proc.time())

fit.sim2 = sampling(model,
               data = stan_data,
               #init = intercepts,
               pars = c('u', 'eps', 'beta'),
               #pars = c('logmx', 'beta'),
               include = FALSE,
               iter = 3000,
               save_warmup = FALSE,
               warmup = 1500,
               chains = 2,
               #thin = 2,
               control = list(adapt_delta = 0.95,
                              max_treedepth = 12)
               )

print(proc.time())

check_hmc_diagnostics(fit.sim2)
get_elapsed_time(fit.sim2)

simcov <- summary(fit.sim2, pars="beta_cor", probs=c(0.025, 0.05, 0.1, 0.9, 0.95, 0.975))$summary %>% 
  as.data.frame %>% 
  rownames_to_column %>% 
  as_tibble %>% 
  mutate(rowname=gsub("\\]", "", gsub("beta_cor\\[", "", rowname))) %>% 
  separate(rowname, into=c("PC", "year", "race1", "race2"), sep="\\,") %>% 
  mutate(across(PC:race2, ~as.numeric(.))) %>%
  select(PC, year, race1, race2, `2.5%`:`97.5%`) %>% filter(race1 < race2) %>%
  mutate(year = year + 2000, race1 = letters[race1], race2 = letters[race2]) %>%
  left_join(bind_rows(
    bind_cols(expand_grid(PC = 1, year = 2001:2010, race1 = letters[1:5]), 
              matrix(corr1, ncol=5, byrow=T) %>% set_colnames(letters[1:5])) %>% 
      pivot_longer(cols=letters[1:5], names_to="race2", values_to="truth"),
    bind_cols(expand_grid(PC = 2, year = 2001:2010, race1 = letters[1:5]), 
              matrix(corr2, ncol=5, byrow=T) %>% set_colnames(letters[1:5])) %>% 
      pivot_longer(cols=letters[1:5], names_to="race2", values_to="truth")
  )) %>%
  mutate(c80 = ifelse(truth >= `10%` & truth <= `90%`, 1, 0),
         c90 = ifelse(truth >= `5%` & truth <= `95%`, 1, 0),
         c95 = ifelse(truth >= `2.5%` & truth <= `97.5%`, 1, 0)) %>%
  summarise(cov80 = sum(c80)/n(),
            cov90 = sum(c90)/n(),
            cov95 = sum(c95)/n())

summary(fit.sim2, pars="beta_cor")$summary %>% 
  as.data.frame %>% 
  rownames_to_column %>% 
  as_tibble %>% 
  mutate(rowname=gsub("\\]", "", gsub("beta_cor\\[", "", rowname))) %>% 
  separate(rowname, into=c("PC", "year", "race1", "race2"), sep="\\,") %>% 
  mutate(across(PC:race2, ~as.numeric(.))) %>% 
  select(PC, year, race1, race2, med = `50%`) %>%
  mutate(year = year + 2000, race1 = letters[race1], race2 = letters[race2]) %>%
  pivot_wider(names_from=race2, values_from=med)

summary(fit.sim2, pars="beta_cor")$summary %>% 
  as.data.frame %>% 
  rownames_to_column %>% 
  as_tibble %>% 
  mutate(rowname=gsub("\\]", "", gsub("beta_cor\\[", "", rowname))) %>% 
  separate(rowname, into=c("PC", "year", "race1", "race2"), sep="\\,") %>% 
  mutate(across(PC:race2, ~as.numeric(.)),
         source = "model") %>% 
  select(source, PC, year, race1, race2, med = `50%`) %>%
  mutate(year = year + 2000, 
         race1 = letters[race1], 
         race2 = letters[race2]) %>%
  bind_rows(
    bind_rows(
      bind_cols(expand_grid(PC = 1, year = 2001:2010, race1 = letters[1:5]), 
                matrix(corr1, ncol=5, byrow=T) %>% set_colnames(letters[1:5])) %>% 
        pivot_longer(cols=letters[1:5], names_to="race2", values_to="med") %>% 
        mutate(source = "truth"),
      bind_cols(expand_grid(PC = 2, year = 2001:2010, race1 = letters[1:5]), 
                matrix(corr2, ncol=5, byrow=T) %>% set_colnames(letters[1:5])) %>% 
        pivot_longer(cols=letters[1:5], names_to="race2", values_to="med") %>% 
        mutate(source = "truth")
    )
  ) %>%
  mutate(race1 = factor(race1, levels=letters[1:10]),
         race2 = factor(race2, levels=rev(letters[1:10]))) %>%
  filter(year == 2001) %>%
  ggplot(aes(x=race1, y=race2, fill=med)) + 
  geom_tile() +
  scale_fill_gradient2(low="blue", high="red", mid="white", midpoint=0, limits=c(-1,1)) + 
  facet_grid(PC~source)

simcordat <- summary(fit.sim2, pars="beta_cor")$summary %>%
  as.data.frame %>%
  tibble::rownames_to_column() %>%
  tibble::as_tibble() %>%
  mutate(rowname=gsub("\\]", "", gsub("beta_cor\\[", "", rowname))) %>%
  separate(rowname, into=c("PC", "year", "race1", "race2"), sep="\\,") %>%
  mutate(across(PC:race2, ~as.numeric(.)),
         source = "model") %>%
  select(source, PC, year, race1, race2, med = `50%`) %>%
  mutate(year = year + 2000,
         race1 = letters[race1],
         race2 = letters[race2]) %>%
  bind_rows(
    bind_rows(
      bind_cols(expand_grid(PC = 1, year = 2001:2010, race1 = letters[1:5]), 
                matrix(corr1, ncol=5, byrow=T) %>% set_colnames(letters[1:5])) %>% 
        pivot_longer(cols=letters[1:5], names_to="race2", values_to="med") %>% 
        mutate(source = "truth"),
      bind_cols(expand_grid(PC = 2, year = 2001:2010, race1 = letters[1:5]), 
                matrix(corr2, ncol=5, byrow=T) %>% set_colnames(letters[1:5])) %>% 
        pivot_longer(cols=letters[1:5], names_to="race2", values_to="med") %>% 
        mutate(source = "truth")
    )
  ) %>%
  mutate(race1 = factor(race1, levels=letters[1:10]),
         race2 = factor(race2, levels=rev(letters[1:10]))) %>%
  mutate(source = factor(source, levels=c("truth", "model")))

acc = c(0,0,0,0,0.25,0.5,1,1,0.5,0.25,0,0,0,0,0,0,0,0,0)

lmxcov <- summary(fit.sim2, pars="logmx", probs=c(0.025, 0.05, 0.1, 0.9, 0.95, 0.975))$summary %>% 
  as.data.frame %>% 
  rownames_to_column %>% 
  as_tibble %>% 
  mutate(rowname=gsub("\\]", "", gsub("logmx\\[", "", rowname))) %>% 
  separate(rowname, into=c("code", "year", "ages", "race"), sep="\\,") %>% 
  mutate(across(code:race, ~as.numeric(.))) %>%
  select(code, year, ages, race, `2.5%`:`97.5%`) %>% 
  mutate(year = year + 2000, ages=sort(pull(dat %>% select(ages) %>% distinct))[ages], 
         race = letters[race]) %>%
  left_join(logmx %>% select(code, year, ages, race, truth=lmx)) %>%
  mutate(c80 = ifelse(truth >= `10%` & truth <= `90%`, 1, 0),
         c90 = ifelse(truth >= `5%` & truth <= `95%`, 1, 0),
         c95 = ifelse(truth >= `2.5%` & truth <= `97.5%`, 1, 0)) %>%
  summarise(cov80 = sum(c80)/n(),
            cov90 = sum(c90)/n(),
            cov95 = sum(c95)/n())

save(corr1, corr2, dat.sim2, base_mort, acc, simcordat, simcov, lmxcov, file="simm6.rda")

# saveRDS(dat, "datCAjt.rds")
# saveRDS(fit, "fitCAjt.rds")

summary(fit.sim2, pars="beta_cor", probs=c(0.025, 0.05, 0.1, 0.9, 0.95, 0.975))$summary %>% 
  as.data.frame %>% 
  rownames_to_column %>% 
  as_tibble %>% 
  mutate(rowname=gsub("\\]", "", gsub("beta_cor\\[", "", rowname))) %>% 
  separate(rowname, into=c("PC", "year", "race1", "race2"), sep="\\,") %>% 
  mutate(across(PC:race2, ~as.numeric(.))) %>%
  select(PC, year, race1, race2, `2.5%`:`97.5%`) %>% filter(race1 < race2) %>%
  mutate(year = year + 2000, race1 = letters[race1], race2 = letters[race2]) %>%
  left_join(bind_rows(
    bind_cols(expand_grid(PC = 1, year = 2001:2010, race1 = letters[1:5]), 
              matrix(corr1, ncol=5, byrow=T) %>% set_colnames(letters[1:5])) %>% 
      pivot_longer(cols=letters[1:5], names_to="race2", values_to="truth"),
    bind_cols(expand_grid(PC = 2, year = 2001:2010, race1 = letters[1:5]), 
              matrix(corr2, ncol=5, byrow=T) %>% set_colnames(letters[1:5])) %>% 
      pivot_longer(cols=letters[1:5], names_to="race2", values_to="truth")
  )) %>%
  mutate(c80 = ifelse(truth >= `10%` & truth <= `90%`, 1, 0),
         c90 = ifelse(truth >= `5%` & truth <= `95%`, 1, 0),
         c95 = ifelse(truth >= `2.5%` & truth <= `97.5%`, 1, 0)) %>%
  summarise(cov80 = sum(c80)/n(),
            cov90 = sum(c90)/n(),
            cov95 = sum(c95)/n())

