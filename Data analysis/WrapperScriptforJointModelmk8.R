# WrapperScriptforJointModelmk6.R
#
# An R-script which estimates US county mortality rates using Bayesian Inference 
# using a joint sex model fitting principal age-specific mortality components, 
# pooling data at the county level, smoothing over time.
#
# The model code 'joint_model_mk6.stan'
#
# Script authors: Monica Alexander, Ameer Dharamshi, Celeste Winant
# Model authors: Monica Alexander, Ameer Dharamshi
# 2021
#
# Dependencies:
# 1. input file of county-level deaths for a given state (e.g. "CAdeaths.csv") (format CSV)
# 2. input file of county-level populations for for a given state ("CApop.csv") (format CSV)
# 3. input file of principal components from reference population ("pcs.rda") (format R binary)
# 4. Stan model file ('joint_model_mk6.stan')
#
# It is recommended that you do *not* run this script interactively (i.e. with rstudio)
# Instead, run from a dedicated 'R' session (launched through a dedicated x-terminal)
# Or from a batch command
#####################
rm(list=ls()) # Clear memory as this script will consume large amounts of memory


library(tidyverse)
library(tidybayes) # not available on quigley
library(rstan)
library(magrittr)
library(usmap) # For US counties: use fips_info to decode fips codes

#sink(paste0("log",Sys.Date(),".txt"),append=FALSE,split=FALSE)

options(mc.cores=5)
rstan_options(auto_write=TRUE)

## The user must modify these next few lines

# This is the pathway for the working directory from which input files are retrived and the the output of the 
# Bayesian estimation are written
# The script expects separate folders for each state (or equivalent geographic region (such as a country) 
# which is further divided into sub-regions (e.g. counties, departments, cantons)
# For this example, we expect all dependencies as well as outputs to be in the same working directory
# that you run the code from
# But you can edit this path to reflect your working directory structure
## wpath <- "/naphsis/cwinant/BayesianEstimation/WORKING" 
wpath <- "~/Data/BayesianEstimation/FRANCE/USCountyBayesianCodePackage/"

#This is the pathway for the program dependencies: the model file and the principal-component file
# For this example, we expect all dependencies as well as outputs to be in the same working directory
# that you run the code from
# But you can edit this path to reflect your working directory structure
## codepath <- "/naphsis/cwinant/BayesianEstimation/codebase"
codepath <- "."

states <- c("AL","AK","AZ","AR","CA","CO","CT","DE","FL","GA","HI","ID","IL","IN","IA","KS","KY","LA",
            "ME","MD","MA","MI","MN","MS","MO","MT","NE","NV","NH","NJ","NM","NY","NC","ND",
            "OH",'OK',"OR","PA","RI","SC","SD","TN","TX","UT","VT","VA","WA","WV","WI","WY")

# The script re-formats the age coding in the inputs using this age-encoder tribble
# Change age encoder depending on input and desired output
age_encoder = 
  tribble(
    ~age_group, ~age_floor, ~ages,  
    '< 1 year'   , 0,   '< 1',
    '1-4 years'  , 1,   '1-4',
    '5-9 years'  , 5,   '5-9',
    '10-14 years', 10,  '10-14',
    '15-19 years', 15,  '15-19',
    '20-24 years', 20,  '20-24',
    '25-29 years', 25,  '25-29',
    '30-34 years', 30,  '30-34',
    '35-39 years', 35,  '35-39',
    '40-44 years', 40,  '40-44',
    '45-49 years', 45,  '45-49',
    '50-54 years', 50,  '50-54',
    '55-59 years', 55,  '55-59',
    '60-64 years', 60,  '60-64',
    '65-69 years', 65,  '65-69',
    '70-74 years', 70,  '70-74',
    '75-79 years', 75,  '75-79',
    '80-84 years', 80,  '80-84',
    '85+ years'  , 85,  '85+'
  ) %>%
  mutate(ages = factor(ages, levels = unique(ages)))

# The principal components are encoded in a file containing three list-objects (e.g. pcs_mf) 
# that are each the output of of the function svd(),
# We use the right singular vectors ( list element 'v') from 
# the principal components constructed from the stacked male and female (so "mf") US State life tables
pcsfile <- "HMD_pcs.rda"
load(file.path(codepath,pcsfile))
pcs <- pcs_mf$v

# Specify the number of principal components
n_pcs <- 4
# Adjust pcs if desired
pcs_adj <- rep(1, n_pcs)

# Main loop 
### If running all states, use this first line to paramaterize the loop
# for (s in states) { # use this line to run all 50 states (serially)

### If running one state (or country), use this line to paramaterize the loop and hard-code the chosen state
for (s in "FRANCE1968_2022") { #Just run with California in this example

  print(s)

## import death data - change file name to target file
death_file <- file.path(wpath,s,"InputDB",paste0(s,"deaths.csv"))
deaths <- read.csv(death_file,header=TRUE) %>%
  rename(state = State,
         code = PopCode,
         age_floor = Age,
         sex = Sex,
         year = Year,
         deaths = Deaths) %>%
  left_join(age_encoder, by='age_floor') %>%
  select(state, code, sex, year, ages, age_floor, deaths) %>%
  #mutate(code = as.numeric(code)) %>%
  mutate(deaths = as.integer(round(deaths))) # We round the deaths since
                                              # the Bayesian estimator only accept deaths in integer format
                                              # Since deaths are modeled as the outcome of a Poisson process

## import population data - change file name to target file
pop_file <- file.path(wpath,s,"InputDB",paste0(s,"pop.csv"))
pop <- read.csv(pop_file,header=TRUE) %>%
  rename(code = PopCode,
         age_floor = Age,
         sex = Sex,
         year = Year) %>% 
  left_join(age_encoder, by='age_floor') %>%
  group_by(code, sex, year, age_floor, ages) %>%
 mutate(code = ifelse(str_length(code)==4,paste0("0",code),code)) %>% # for US county codes
  #mutate(code = as.numeric(code)) %>%
  summarise(pop = sum(Population)) %>%
  ungroup

## Note that with the right_join(), the code will impute deaths = 0 for each region/age/sex/year 
## combination represented
## in the population series where data are missing in the death series.  This accommodates the 
## fact that detailed mortality files will not explicity enumerate cases where deaths = 0.
## The script also accepts explicit entries of deaths = 0
dat <- right_join(deaths, pop, 
                  by = c('code', 'sex', 'year', 'ages', 'age_floor')) %>%
  # mutate( code = as.numeric(code)) %>%
  mutate(deaths = ifelse(is.na(deaths),0,deaths)) %>% # assume deaths = NA means deaths = 0
  mutate(state = s) %>%
  mutate(code = as.character(ifelse(str_length(code)==1,paste0("0",code),code))) %>% # for US county codes
  arrange(code, sex, year, age_floor) %>%
  select(!age_floor)

## Write joint deaths/population tibble to file
save(dat, file=file.path(wpath,s,"InputDB",paste0(s,"input_data.rda")))

#dat <- dat %>%
#  mutate(codei = as.numeric(dense_rank(code)))

# Set the baseyear in case it differs from the earliest year in your input data series
baseyear <- dat %>% select(year) %>% pull %>% min

# Build the model input list
stan_data <- dat %>%
  mutate(year = year - baseyear + 1) %>% 
  compose_data(state = x_at_y(state, code)) %>%
  within({
    n_year <- length(unique(year))
    pcs <- pcs[,1:n_pcs] %*% diag(pcs_adj)
    n_pcs <- n_pcs
  })

model_file <- file.path(codepath,"joint_model_mk8.stan")

if (file.exists(file.path(codepath,"joint_model_mk8.rds"))) {
  unlink(file.path(codepath, "joint_model_mk8.rds"))
}
model <- rstan::stan_model(model_file, model_name = "base mortality joint model",verbose=TRUE)

print(proc.time())

iters <- 3000
warmups <- iters - 2500

fit = sampling(model,
                 data = stan_data,
                 #pars = c('u', 'eps'),
                 #include = FALSE,
                 iter = iters,
                 warmup = warmups,
                 chains = 4,
                 control = list(adapt_delta = 0.9,
                                max_treedepth = 12)
)

sampleTag <- ifelse(iters==5000,"BigWarmup","")
pcsTag <- ifelse(pcsfile=="HMD_pcs.rda","HMDpcs",
                 ifelse(pcsfile=="FRATNP_pcs.rda","FRApcs",""))
modelTag <- ifelse(model_file==file.path(codepath,"joint_model_mk8.stan"),"mk8","")

print(proc.time())

# Write output (fit) to file (state+'fitJointModel.rds') under the subdirectory 'fit'
# Note that the file is saved as a single-object R-binary
saveRDS(fit, file=file.path(wpath,s,"fit",paste0("fitJointModelNewClone",pcsTag,sampleTag,s,modelTag,".rds")))
# rm(fit) # if you are using the loop to launch serial Bayesian estimations for different states
          # it's good practice to remove the output object 'fit' before the next loop iteration
}

#sink()


