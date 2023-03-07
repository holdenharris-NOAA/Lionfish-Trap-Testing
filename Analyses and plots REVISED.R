########################################################################
## ---------------------------------------------------------------------
##
## R CODE FOR 'Three Trap Designs Evaluated for 
##             a Deepwater Lionfish Fishery'
##
## Analyses:  Holden Earl Harris
## Contact:   holden.harris@noaa.gov
##            holdenharris@ufl.edu
## Last edit: March 2023
##            Revised based on reviewer suggestions
##
## ---------------------------------------------------------------------
########################################################################

## Push test

## WORKING LIBRARIES----------------------------------------------------
rm(list=ls()); gc(); windows()
library(dplyr)
library(car)
library(MASS)
library(ggplot2) 
library(lme4)
library(cowplot)
library(pscl)
library(glmmTMB)
library(car)
library(MASS)
library(lme4)


########################################################################
## FUNCTIONS TO WRITE OUT MODEL

## FUNCTION TO NOTE SIGNIFICANCE----------------------------------------
sig_stars <- function(P){
  stars = P
  stars[P>0.1] = "" 
  stars[P<=0.1] = "." 
  stars[P<=0.05] = "*" 
  stars[P<=0.01] = "**" 
  stars[P<=0.001] = "***" 
  return(stars)
}

## FUNCTION TO OUTPUT PARAMETERS AND SIGNIFICANCE FROM LOG-TRANSFORMED GLM
coeftab_glm <- function(mod, response = "", model = "", dec_places = 3){
  #mod = lf_nb
  coefs = as.data.frame(coef(summary(mod)))
  out = data.frame(Effects = rownames(coefs))
  out$Est = round(exp(coefs$Estimate), dec_places)
  out$Est[1] = round(coefs$Estimate[1], dec_places)
  upCI  = round(exp(coefs$Estimate + 1.96 * coefs$`Std. Error`), dec_places)
  lowCI = round(exp(coefs$Estimate - 1.96 * coefs$`Std. Error`), dec_places)
  upCI[1]  = round(coefs$Estimate[1] + 1.96 * coefs$`Std. Error`[1], dec_places)
  lowCI[1] = round(coefs$Estimate[1] - 1.96 * coefs$`Std. Error`[1], dec_places)
  out$'95% CI' = paste(sprintf("%.3f", lowCI), "-", sprintf("%.3f", upCI), sep = "")
  out$test_stat = round(coefs[3], dec_places)
  out$P = coefs[4]
  names(out$test_stat)=names(coefs[3])
  P = coefs[,4]
  P = sprintf("%.3f", P)
  P[P<0.001] = "<0.001" 
  out$P = P
  out$Sig = sig_stars(P)
  out$Estimate = NULL; out$`Std. Error`=NULL
  if (response != "") out$Response = response
  if (model    != "") out$Model    = model
  return(out)
}


###########################################################################
## ------------------------------------------------------------------------
##
## PART I - CATCHES AND CPUE
##
## ------------------------------------------------------------------------
###########################################################################

fish       <- read.csv("./data/trap_catches.csv")
retrievals <- read.csv("./data/trap_retrievals.csv")

## ---------------------------------------------------------------------
## COMPILE CATCHES

## PULL OUT SITE COVARIATES
site_covars <- data.frame(
  Replicate = retrievals$Replicate,
  Site_ID = retrievals$Site_ID,
  Trip = retrievals$Trip,
  Type = retrievals$Type,
  Array_num = factor(retrievals$Array_num),
  Num_traps = retrievals$Drop_Succ,
  Region = retrievals$Region,
  Soakdays = retrievals$Soakdays,
  Depth = retrievals$Dep_m,
  Relief_m = retrievals$Relief_m,
  Relief_cat = retrievals$Relief_cat,
  LF_dens = retrievals$Mean_dens,
  Rep_id = retrievals$Rep_id
)

site_covars <- site_covars[order(site_covars$Rep_id), ]; site_covars$Rep_id = NULL

## SUMMARIZE LIONFISH CATCHES AND ADD COLUMN
summ_catch_lionfish <-   
  fish %>% 
  filter(!is.na(Fishery_cat)) %>% 
  filter(!is.na(WT)) %>%
  filter(Fishery_cat == "lionfish") %>% 
  group_by(Replicate) %>% 
  summarise(biom_lf = sum(WT)) %>%  
  as.data.frame(); summ_catch_lionfish


## SUMMARIZE FISHERY CATCHES AND ADD COLUMN
summ_catch_fishery <-   
  fish %>% 
  filter(!is.na(Fishery_cat)) %>% 
  filter(!is.na(WT)) %>%
  filter(Fishery_cat == "fishery") %>% 
  group_by(Replicate) %>% 
  summarise(biom_fish = sum(WT)) %>%  
  as.data.frame(); summ_catch_fishery


## SUMMARIZE NONFISHERY CATCHES AND ADD COLUMN
summ_catch_nonfishery <-   
  fish %>% 
  filter(!is.na(Fishery_cat)) %>% 
  filter(!is.na(WT)) %>%
  filter(Fishery_cat == "non-fishery") %>% 
  group_by(Replicate) %>% 
  summarise(biom_sdrf = sum(WT)) %>%  
  as.data.frame(); summ_catch_nonfishery

## LIONFISH
catch_mod_tab <- merge(site_covars, summ_catch_lionfish, by = "Replicate", all = TRUE)
catch_mod_tab$biom_lf[is.na(catch_mod_tab$biom_lf)] <- 0
catch_mod_tab$cpue_lf = catch_mod_tab$biom_lf / catch_mod_tab$Num_traps
catch_mod_tab$bin_lf = ifelse(catch_mod_tab$biom_lf > 0, 1, 0)

## FISHERY
catch_mod_tab <- merge(catch_mod_tab, summ_catch_fishery, by = "Replicate", all = TRUE)
catch_mod_tab$biom_fish[is.na(catch_mod_tab$biom_fish)] <- 0
catch_mod_tab$cpue_fish = catch_mod_tab$biom_fish / catch_mod_tab$Num_traps
catch_mod_tab$bin_fish = ifelse(catch_mod_tab$biom_fish > 0, 1, 0)

## NONFISHERY
catch_mod_tab <- merge(catch_mod_tab, summ_catch_nonfishery, by = "Replicate", all = TRUE)
catch_mod_tab$biom_sdrf[is.na(catch_mod_tab$biom_sdrf)] <- 0
catch_mod_tab$cpue_sdrf = catch_mod_tab$biom_sdrf / catch_mod_tab$Num_traps
catch_mod_tab$bin_sdrf = ifelse(catch_mod_tab$biom_sdrf > 0, 1, 0)

## QUICK LOOK AT BOXPLOTS
catch_mod_tab$Type <- as.factor(catch_mod_tab$Type)
par(mfrow=c(2,2))
plot(catch_mod_tab$cpue_lf   ~ catch_mod_tab$Type, col = 1:3, xlab = 'Trap type', ylab = 'Lionfish CPUE')
plot(catch_mod_tab$cpue_fish ~ catch_mod_tab$Type, col = 1:3, xlab = 'Trap type', ylab = 'Fisheries spp. CPUE')
plot(catch_mod_tab$cpue_sdrf ~ catch_mod_tab$Type, col = 1:3, xlab = 'Trap type', ylab = 'Non-fisheries spp. CPUE')
par(mfrow=c(1,1))


###########################################################################
## CONVERT BIOMASSES FROM KILOGRAMS TO GRAMS
catch_mod_tab$biom_fish = round(catch_mod_tab$biom_fish * 1000)
catch_mod_tab$cpue_fish = round(catch_mod_tab$cpue_fish * 1000)

catch_mod_tab$biom_sdrf = round(catch_mod_tab$biom_sdrf * 1000)
catch_mod_tab$cpue_sdrf = round(catch_mod_tab$cpue_sdrf * 1000)

catch_mod_tab$biom_lf = round(catch_mod_tab$biom_lf * 1000)
catch_mod_tab$cpue_lf = round(catch_mod_tab$cpue_lf * 1000)

## REORDER: LOBSTER, SEA BASS, GITTINGS
catch_mod_tab$Type = factor(catch_mod_tab$Type, levels = c("L", "S", "G"))


###########################################################################
## ------------------------------------------------------------------------
## PLOT CPUE 
## ------------------------------------------------------------------------

## SUMMARIZE (SUMM) MEANS AND SE FOR BOXPLOT---------------------------------------
summ_lf <- catch_mod_tab %>% 
  group_by(Type) %>% 
  summarise(
    mu_cpue = mean(cpue_lf),
    se_cpue = sd(cpue_lf)/sqrt(n()),
  ) %>% 
  as.data.frame(); summ_lf$Category = "Lionfish"

summ_fish <- catch_mod_tab %>% 
  group_by(Type) %>% 
  summarise(
    mu_cpue = mean(cpue_fish),
    se_cpue = sd(cpue_fish)/sqrt(n()),
  ) %>% 
  as.data.frame(); summ_fish$Category = "Fishery species"

summ_sdrf <- catch_mod_tab %>% 
  group_by(Type) %>% 
  summarise(
    mu_cpue = mean(cpue_sdrf),
    se_cpue = sd(cpue_sdrf)/sqrt(n()),
  ) %>% 
  as.data.frame(); summ_sdrf$Category = "Non-fishery species"

summ_bytype <- rbind(summ_lf, summ_fish, summ_sdrf); summ_bytype

## ORDER: LOBSTER, SEABASS, GITTINGS
summ_bytype$Type = factor(summ_bytype$Type, levels = c("L", "S", "G"))

## PLOT--------------------------------------------------------------------
ax.txt = 11
ax.ttl = 14
leg.ttl = 13
leg.txt = 13
strp.txt = 12
ax.tit = 12

cols = c("lightcoral", "palegreen3", "steelblue1") 
s = 0.5
w = 0.2

plot_cpue <- 
  ggplot(summ_bytype) +
  geom_bar(aes(fill = factor(Type, levels = c("L", "S", "G")),
               x = Type,
               y = mu_cpue), 
           stat = "identity",
           color = 'black',
           position= position_dodge(width=0.90),
           width = 0.8) +
  facet_wrap(~factor(Category, levels = c("Lionfish", "Fishery species", "Non-fishery species")), 
             scales = "free_y", strip.position = "top") +
  geom_errorbar(aes(x = Type, 
                    ymin = pmax(mu_cpue - (1.96 * se_cpue), 0.00),   
                    ymax = mu_cpue + (1.96 * se_cpue)),
                size = s, 
                position=position_dodge(width=0.90),
                width = w) +
  xlab ("") +
  ylab ("CPUE (g per trap)") +
  scale_fill_manual(values = cols, guide = guide_legend(title = "Trap type"), 
                    labels = c("Lobster", "Sea bass","Gittings")) +
  scale_y_continuous(expand = c(0,0,.05,0)) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    panel.border = element_rect(fill = NA, colour = "white"),
    axis.title = element_text(size = ax.ttl, color = "black"),
    axis.text = element_text(size = ax.txt, color = "black"),
    axis.title.y = element_text(size = ax.tit),
    legend.title = element_text(size = leg.ttl),
    legend.text  = element_text(size = leg.txt),
    axis.line = element_line(),
    panel.background = element_rect(fill = "grey100"),
    legend.position = "right",
    legend.background = element_rect(size = 0.5, color = "black", linetype="solid")
  ); plot_cpue

## WRITE OUT PLOT
# (Remove comments to make plot)
tiff(filename = "./CPUE_plots.tif", 
     units = "in", width = 7, height = 3.5, res = 400)
plot(plot_cpue) 
dev.off()


###########################################################################
## ------------------------------------------------------------------------
## STATISTICAL MODELS FOR CPUE
## ------------------------------------------------------------------------


## First, look at the linear model of lionfish catch against site density
par(mfrow=c(1,1))
dens_lm <- lm(biom_lf ~ LF_dens, data = catch_mod_tab); summary(dens_lm)
plot(cpue_lf ~ LF_dens, data = catch_mod_tab,
     ylab = "Lionfish CPUE (kg)", xlab = "Lionfish density (fish per 100 m2)")
abline(dens_lm) 
## --> Catch not correlated to LF density
## --> Looks zero-inflated
## --> A GLM appears needed


## Examine ERROR STRUCTURE------------------------------------------------
## FOR THE THREE CATCH TYPES
pos_lf   =  subset(catch_mod_tab, catch_mod_tab$bin_lf == 1)
pos_fish = subset(catch_mod_tab, catch_mod_tab$bin_fish == 1)
pos_sdrf = subset(catch_mod_tab, catch_mod_tab$bin_sdrf == 1)

## LIONFISH
## All catches
par(mfrow=c(2,2))
qqp(catch_mod_tab$cpue_lf)
qqp(catch_mod_tab$cpue_lf, distribution = "lnorm")
poisson <- fitdistr(as.integer(catch_mod_tab$cpue_lf), "Poisson"); qqp(catch_mod_tab$cpue_lf, "pois", lambda = poisson$estimate) 
nbinom <- fitdistr(as.integer(catch_mod_tab$cpue_lf), "Negative Binomial"); qqp(catch_mod_tab$cpue_lf, "nbinom", size = nbinom$estimate[[1]], mu = nbinom$estimate[[2]]) ## BEST FIT DUE TO OVERDISPERSION
## --> Zero inflated

## Exclude zeros
par(mfrow=c(2,2))
qqp(pos_lf$cpue_lf)
qqp(pos_lf$cpue_lf, distribution = "lnorm")
poisson <- fitdistr(as.integer(pos_lf$cpue_lf), "Poisson"); qqp(pos_lf$cpue_lf, "pois", lambda = poisson$estimate) 
nbinom <- fitdistr(as.integer(pos_lf$cpue_lf), "Negative Binomial"); qqp(pos_lf$cpue_lf, "nbinom", size = nbinom$estimate[[1]], mu = nbinom$estimate[[2]]) ## BEST FIT DUE TO OVERDISPERSION
## --> Negative binomial best fit

## FISHERY
## All catches
par(mfrow=c(2,2))
qqp(catch_mod_tab$cpue_fish)
qqp(catch_mod_tab$cpue_fish + 1, distribution = "lnorm") ## BEST FIT
poisson <- fitdistr(as.integer(catch_mod_tab$cpue_fish), "Poisson"); qqp(catch_mod_tab$cpue_fish, "pois", lambda = poisson$estimate) 
nbinom <- fitdistr(as.integer(catch_mod_tab$cpue_fish), "Negative Binomial"); qqp(catch_mod_tab$cpue_fish, "nbinom", size = nbinom$estimate[[1]], mu = nbinom$estimate[[2]]) ## BEST FIT DUE TO OVERDISPERSION
## --> Less than lionfish, but still very zero inflated

## Exclude zeros
par(mfrow=c(2,2))
qqp(pos_fish$cpue_fish)
qqp(pos_fish$cpue_fish + 1, distribution = "lnorm") ## BEST FIT
poisson <- fitdistr(as.integer(pos_fish$cpue_fish), "Poisson"); qqp(pos_fish$cpue_fish, "pois", lambda = poisson$estimate) 
nbinom <- fitdistr(as.integer(pos_fish$cpue_fish), "Negative Binomial"); qqp(pos_fish$cpue_fish, "nbinom", size = nbinom$estimate[[1]], mu = nbinom$estimate[[2]]) ## BEST FIT DUE TO OVERDISPERSION
## --> Log normal and negative binomial best fits


## NON-FISHERY
## All catches
par(mfrow=c(2,2))
qqp(catch_mod_tab$cpue_sdrf, distribution = "norm")
qqp(catch_mod_tab$cpue_sdrf + 1, distribution = "lnorm")
poisson <- fitdistr(as.integer(catch_mod_tab$cpue_sdrf), "Poisson"); qqp(catch_mod_tab$cpue_sdrf, "pois", lambda = poisson$estimate) 
nbinom <- fitdistr(as.integer(catch_mod_tab$cpue_sdrf), "Negative Binomial"); qqp(catch_mod_tab$cpue_sdrf, "nbinom", size = nbinom$estimate[[1]], mu = nbinom$estimate[[2]]) ## BEST FIT DUE TO OVERDISPERSION
par(mfrow=c(1,1))
## --> Again, zero inflated

## Exclude zeros
par(mfrow=c(2,2))
qqp(pos_sdrf$cpue_fish)
qqp(pos_sdrf$cpue_fish + 1, distribution = "lnorm") ## BEST FIT
poisson <- fitdistr(as.integer(pos_sdrf$cpue_fish), "Poisson"); qqp(pos_sdrf$cpue_fish, "pois", lambda = poisson$estimate) 
nbinom <- fitdistr(as.integer(pos_sdrf$cpue_fish), "Negative Binomial"); qqp(pos_sdrf$cpue_fish, "nbinom", size = nbinom$estimate[[1]], mu = nbinom$estimate[[2]]) ## BEST FIT DUE TO OVERDISPERSION
## --> Again, log normal and negative binomial best fits



###########################################################################
## CPUE ANALYSES AND MODEL SIMPLIFICATION
library(pscl)
library(car)

## -----------------------------------------------------------------------------
## Lionfish 
catch_mod_tab$response     = catch_mod_tab$cpue_lf
catch_mod_tab$response_bin = catch_mod_tab$bin_lf

## ZINB (response = CPUE)
lf_zi <- zeroinfl(response ~ Type + Num_traps + Soakdays + LF_dens + Depth + Relief_m + Region ## Positive component
                           | Type + Num_traps, ## Zero component
                  data = catch_mod_tab, 
                  dist = "negbin")
summary(lf_zi)
exp(lf_zi$coefficients$count)
exp(lf_zi$coefficients$zero)
Anova(lf_zi) ##--> Deviance table only gives deviance table for positive component, so look at the individual models
## -->
## Logistic regression: Gittings trap type and configuration (number of traps) significant
##   Prob. of no catch is 83% lower with Gittings traps and 53% lower per additional trap. 
## Neg. binom. pos. count regression: Config. is the only significant factor, with about 45% lower per additional trap. 

## Examine constituent positive and binomial models
## Negative binomial (response = CPUE)
lf_nb <- glm.nb(response ~ Type + Num_traps + Soakdays + LF_dens + Depth + Relief_m + Region,
                data = subset(catch_mod_tab, catch_mod_tab$bin_lf == 1))
summary(lf_nb) ## Estimates same as ZINB positive component
exp(lf_nb$coefficients)
Anova(lf_nb)
coeftab_glm(lf_nb)

## Binomial (response = binary catch)
lf_bi <- glm(response_bin ~ Type + Num_traps,
             data = catch_mod_tab, family = "binomial")
summary(lf_bi) ## Similar as ZINB zero component except estimates multiplied by -1
exp(lf_bi$coefficients) ## These equal the inverse of the ZINB zero component, i.e., 1 / exp(lf_zi$coefficients$zero)
Anova(lf_bi)
coeftab_glm(lf_bi)

catch_mod_tab$response = catch_mod_tab$response_bin = NULL ## Clear response

## ------------------------------------------------------------------------
## FISHERY CATCHES
catch_mod_tab$response     = catch_mod_tab$cpue_fish
catch_mod_tab$response_bin = catch_mod_tab$bin_fish

fish_zi <- zeroinfl(response ~ Type + Num_traps + Soakdays + LF_dens + Depth + Relief_m + Region
                              | Type + Num_traps,
                    data = subset(catch_mod_tab, catch_mod_tab$Type != "G"), 
                                  dist = "negbin")
summary(fish_zi)
exp(fish_zi$coefficients$count)
exp(fish_zi$coefficients$zero)
Anova(fish_zi)
## --> Note that Gittings traps were removed due to zero catches, which was causing the model to have trouble (return: "system is computationally singular"). 
##     Colineary between depth and zero Gittings catches appeared to cause singularities (maybe division by zero?) and not allow the model to converge.
## -->
## Logistic regression: Trap type (sea bass) and soak time significant. 
##   Probability of no catch 83% lower with sea bass traps and increases with increasing soak time. 
## Neg. binom. pos. count regression: Sea bass trap type highly significant with about 105% higher catch. 
##   Depth significant with about 3% lower catch amount per increasing m of depth. Theta (i.e., dispersion) very high. 

## Examine constituent positive and binomial models
## Negative binomial (response = CPUE)
fish_nb <- glm.nb(response ~ Type + Num_traps + Soakdays + LF_dens + Depth + Relief_m + Region,
                data = subset(catch_mod_tab, catch_mod_tab$bin_fish == 1))
summary(fish_nb) ## Estimates same as ZINB positive component
exp(fish_nb$coefficients)
Anova(fish_nb)
coeftab_glm(fish_nb)

## Binomial (response = binary catch)
fish_bi <- glm(response_bin ~ Type + Num_traps,
             data = catch_mod_tab, family = "binomial")
summary(fish_bi) ## Same as ZINB zero component except estimates multiplied by -1
exp(fish_bi$coefficients) ## These equal the inverse of the ZINB zero component, i.e., 1 / exp(fish_zi$coefficients$zero)
Anova(fish_bi)
coeftab_glm(fish_bi)

catch_mod_tab$response = catch_mod_tab$response_bin = NULL ## Clear response

## ------------------------------------------------------------------------
## NON-FISHERY CATCHES
catch_mod_tab$response     = catch_mod_tab$cpue_sdrf
catch_mod_tab$response_bin = catch_mod_tab$bin_sdrf

sdrf_zi <- zeroinfl(response ~ Type + Num_traps + Soakdays + LF_dens + Depth + Relief_m + Region
                              | Type + Num_traps,
                    data = catch_mod_tab, dist = "negbin")
summary(sdrf_zi)
exp(sdrf_zi$coefficients$count)
exp(sdrf_zi$coefficients$zero)
Anova(sdrf_zi)
## Logistic regression: Seabass trap type and configuration (number of traps) significant.
##   Probability of no catch 84% lower with sea bass traps and 58% lower per additional trap. 
## Neg. binom. pos. count regression: Sea bass trap type highly significant with about 505% higher catch. 
##   Number of traps also significant, with a 21% decrease in CPUE per additional trap. 

## Examine constituent positive and binomial models
## Negative binomial (response = CPUE)
sdrf_nb <- glm.nb(response ~ Type + Num_traps + Soakdays + LF_dens + Depth + Relief_m + Region,
                data = subset(catch_mod_tab, catch_mod_tab$bin_sdrf == 1))
summary(sdrf_nb) ## Estimates same as ZINB positive component
exp(sdrf_nb$coefficients)
Anova(sdrf_nb)
coeftab_glm(sdrf_nb)

## Binomial (response = binary catch)
sdrf_bi <- glm(response_bin ~ Type + Num_traps,
             data = catch_mod_tab, family = "binomial")
summary(sdrf_bi) ## Same as ZINB zero component except estimates multiplied by -1
exp(sdrf_bi$coefficients) ## These equal the inverse of the ZINB zero component, i.e., 1 / exp(sdrf_zi$coefficients$zero)
Anova(sdrf_bi)
coeftab_glm(sdrf_bi)

catch_mod_tab$response = catch_mod_tab$response_bin = NULL ## Clear response

## Examine histogram of non-fishery catches by seabass traps
sb = subset(catch_mod_tab, catch_mod_tab$Type=="S")
hist(sb$cpue_sdrf)
hist(catch_mod_tab$cpue_sdrf)
## --> dominated by two outliers




###########################################################################
## ------------------------------------------------------------------------
##
## PART II - IN SITU RECRUITMENT
##           from REMOTE TIMELAPSE CAMERA VIDEO DATA
##
## ------------------------------------------------------------------------
###########################################################################

vids <- read.csv('./data/timelapse_fish_counts.csv')

## Merge lionfish densities
vids <- merge(vids, catch_mod_tab[ , c("Site_ID", "LF_dens")], by = "Site_ID", all.x = TRUE)

## Aggregate individual species counts into "lionfish", "fishery", or "non-fishery"
byvid <- vids %>% 
  filter(!is.na(Min_count)) %>% 
  group_by(Vid_stamp, Replicate_video, Trip, Type, TOD, Site_ID, Trap_footprint, soak_d) %>% 
  summarise(
    lf_in = sum(Min_count[Fishery_cat=="lionfish" & In.Out=="in"], na.rm = TRUE),
    lf_out = sum(Min_count[Fishery_cat=="lionfish" & In.Out=="out"], na.rm = TRUE),
    fish_in = sum(Min_count[Fishery_cat=="fishery" & In.Out=="in"], na.rm = TRUE),
    fish_out = sum(Min_count[Fishery_cat=="fishery" & In.Out=="out"], na.rm = TRUE),
    nonfish_in = sum(Min_count[Fishery_cat=="non-fishery" & In.Out=="in"], na.rm = TRUE),
    nonfish_out = sum(Min_count[Fishery_cat=="non-fishery" & In.Out=="out"], na.rm = TRUE),
    nvid = n()); head(byvid)

## Merge with site data to get lionfish density per site
byvid <-merge(byvid, catch_mod_tab[,c(2,12)], by="Site_ID")

## For repeat reads at same site and time of day, aggregate into Dawn, Midday, or Dusk
bytod <- byvid %>% 
  group_by(Replicate_video, Trip, Type, TOD, Site_ID, Trap_footprint) %>% 
  summarise(
    soak_d      = mean(soak_d),
    lf_in       = mean(lf_in),
    lf_out      = mean(lf_out),
    fish_in     = mean(fish_in),
    fish_out    = mean(fish_out),
    nonfish_in  = mean(nonfish_in),
    nonfish_out = mean(nonfish_out),
    nvid = sum(nvid),
    ntod = n()) %>% 
  mutate(across(where(is.numeric), round, 0))

## Merge with site data to get lionfish density per site
bytod <-merge(bytod, catch_mod_tab[,c(2,12)], by="Site_ID")

## Order factor levels
bytod$Type = factor(bytod$Type, levels = c("L", "S", "G"))
bytod$TOD = factor(bytod$TOD, levels = c("A", "M", "U"))


# SUMMARIZE DATA----------------------------------------------------------------
nrow(bytod)
length(subset(bytod$lf_in, bytod$Type == "G" & bytod$lf_in == 0)) 
length(subset(bytod$lf_in, bytod$Type == "L" & bytod$lf_in == 0)) 
length(subset(bytod$lf_in, bytod$Type == "S" & bytod$lf_in == 0)) 

summ_bytype <- bytod %>% 
  filter(!is.na(TOD)) %>% 
  group_by(Type) %>% 
  summarise(
    soak_d = mean(soak_d),
    mu_lf_in = mean(lf_in),
    se_lf_in = sd(lf_in)/sqrt(n()),
    cnt_zero = sum(.==0),
    mu_lf_out = mean(lf_out),
    mu_fish_in = mean(fish_in),
    se_fish_in = sd(fish_in)/sqrt(n()),
    mu_fish_out = mean(fish_out),
    mu_nonfish_in = mean(nonfish_in),
    se_nonfish_in = sd(nonfish_in)/sqrt(n()),
    mu_nonfish_out = mean(nonfish_out)); summ_bytype

summ_bytod <- bytod %>% 
  filter(!is.na(TOD)) %>% 
  group_by(Type, TOD) %>% 
  summarise(
    soak_d = mean(soak_d),
    mu_lf_in = mean(lf_in),
    se_lf_in = sd(lf_in)/sqrt(n()),
    cnt_zero = sum(.==0),
    mu_lf_out = mean(lf_out),
    mu_fish_in = mean(fish_in),
    se_fish_in = sd(fish_in)/sqrt(n()),
    mu_fish_out = mean(fish_out),
    mu_nonfish_in = mean(nonfish_in),
    se_nonfish_in = sd(nonfish_in)/sqrt(n()),
    mu_nonfish_out = mean(nonfish_out)) %>% 
  as.data.frame(); summ_bytod
summ_bytod$Type = factor(summ_bytod$Type, levels = c("L", "S", "G"))

summ_soak <- bytod %>% 
  filter(!is.na(TOD)) %>% 
  group_by(Soak_int = round(soak_d, 0), Type) %>% 
  summarise(
    mu_lf_in = mean(lf_in),
    se_lf_in = sd(lf_in)/sqrt(n()),
  ); summ_soak

summ_bydens <- vids %>% 
  filter(!is.na(Min_count)) %>% 
  filter(!is.na(TOD)) %>% 
  filter(Fishery_cat == "lionfish") %>% 
  group_by(Site_ID, Type, In.Out) %>% 
  summarise(
    mu = mean(Min_count, na.rm = TRUE),
    dens = mean(LF_dens),
    nvid = n()
  ) ; summ_bydens
summ_bydens$Type = factor(summ_bydens$Type, levels = c("L", "S", "G"))

###########################################################################
## ------------------------------------------------------------------------
## RECRUITMENT PLOTS
## ------------------------------------------------------------------------

## PLOT SETTINGS-----------------------------------------------------------
dielcols = c("sienna2", "darkslategray3", "violet")
trapcols = c("lightcoral", "palegreen3", "steelblue1") 
s = 0.9
w = 0.2
ax.title = 11
ax.txt = 9
leg.txt = 10
leg.ttl = 10
leg.pos.x = 0.2
leg.pos.y = 0.85
endday = 9
ylab_lf = "Lionfish count"
ylab_fish = "Fishery species count"
ylab_nonfish = "Non-fishery species count"

## PLOT LIONFISH RECRUITMENT---------------------------------------------- 
## by SOAK TIME
endday = 10
plot_lf_by_soak <- 
  ggplot(bytod, aes(x = soak_d,
                    color = factor(bytod$Type, levels = c("L", "S", "G")),
                    y = lf_in)) +
  geom_smooth(stat= "smooth", method = "loess", size = 1.5, alpha = 0.25, span = 0.55) +
  scale_color_manual(values = trapcols, guide = guide_legend(title = "Trap type"), 
                     labels = c("Lobster", "Seabass", "Gittings")) +
  scale_y_continuous(expand = c(0,0,0,0)) +
  scale_x_continuous(limits = c(0,endday), breaks=seq(0,endday,1)) +
  xlab ("Days deployed") +
  ylab (ylab_lf) +
  coord_cartesian(xlim = c(0.5, endday), ylim = c(0, 0.8), expand = TRUE) + 
  theme(
    legend.position = c(0.85, 0.8),
    axis.title = element_text(size = ax.title, color = "black"),
    axis.text = element_text(size = ax.txt, color = "black"),
    legend.title = element_text(size = leg.ttl),
    legend.text  = element_text(size = leg.txt),
    axis.line = element_line(),
    panel.background = element_rect(fill = "grey100"),
    legend.background = element_rect(size = 0.5, color = "black", linetype="solid")
  ); plot_lf_by_soak

## PLOT LIONFISH RECRUITMENT---------------------------------------------- 
## by TRAP TYPE and by TIME OF DAY
plot_lf_tod <- 
  ggplot(summ_bytod) +
  geom_bar(aes(x = Type,
               y = mu_lf_in, 
               fill = TOD),
           stat = "identity",
           col = "black",
           position= position_dodge(width=0.90)) +
  geom_errorbar(aes(x = Type, fill = TOD, 
                    ymin = pmax(mu_lf_in - (1.96 * se_lf_in), 0.00),   
                    ymax = mu_lf_in + (1.96 * se_lf_in)),
                size = s, 
                position=position_dodge(width=0.9),
                width = w) +
  scale_x_discrete(labels = c("Lobster", "Seabass", "Gittings")) +
  scale_fill_manual(values = dielcols, guide = guide_legend(title = "Diel period"), labels = c("Dawn", "Midday", "Dusk")) +
  scale_y_continuous(expand = c(0,0,0,0.02)) +
  xlab ("Trap type") +
  ylab (ylab_lf) +
  theme(
    axis.title = element_text(size = ax.title, color = "black"),
    axis.text = element_text(size = ax.txt, color = "black"),
    legend.title = element_text(size = leg.ttl),
    legend.text  = element_text(size = leg.txt),
    axis.line = element_line(),
    panel.background = element_rect(fill = "grey100"),
    legend.position='right',
    legend.background = element_rect(size = 0.5, color = "black", linetype="solid")
  ); plot_lf_tod

## PLOT LIONFISH RECRUITMENT---------------------------------------------- 
## by SITE DENSITY OF LIONFISH
plot_vid_bydens <- 
  summ_bydens %>% 
  filter (In.Out == "in") %>% 
  ggplot(aes(x = dens, y = mu, shape = Type, color = Type)) +
  geom_point(size = 2) +
  geom_smooth(method = "glm", span = 0.95, alpha = 0.25, size = 1.5) +
  scale_color_manual(values = trapcols, 
                     name = "Trap type", 
                     labels = c("Lobster", "Seabass", "Gittings")) +
  scale_shape_manual(values = c(15, 19, 17),                     
                     name = "Trap type", 
                     labels = c("Lobster", "Seabass", "Gittings")) +
  scale_y_continuous(expand = c(0,0,0,0)) +
  scale_x_continuous(limits = c(0,1.5), breaks=seq(0,1.5,0.5)) +
  coord_cartesian(ylim = c(0, 1.9), xlim = c(0.07, 1.5)) + 
  xlab ("Lionfish density") +
  ylab ("Lionfish count") +
  theme(
    legend.title = element_text(size = leg.ttl),
    legend.text  = element_text(size = leg.txt),
    axis.line = element_line(),
    panel.background = element_rect(fill = "grey100"),
    legend.position = c(0.85, 0.8),
    legend.background = element_rect(size = 0.5, color = "black", linetype="solid")
  ); plot_vid_bydens

## COMPILE AND WRITE OUT PLOTS-------------------------------------------
library(gridExtra)
vidplots <- 
  grid.arrange(plot_lf_by_soak, plot_vid_bydens, plot_lf_tod,
             widths = c(1, 2, 2, 1),
             layout_matrix = rbind(c(1, 1, 2, 2), c(NA, 3, 3, NA)))

png(filename = "./Video_plots.png", 
     units = "in", width = 8, height = 9, res = 2000)
plot(vidplots) 
dev.off()

###########################################################################
## ------------------------------------------------------------------------
## STATISTICAL MODELS FOR FOR LIONFISH RECRUITMENT 
## ------------------------------------------------------------------------

bytod$TOD = factor(bytod$TOD, levels = c("M", "A", "U"))

## DETERMINE ERROR STRUCTURE-----------------------------------------------
par(mfrow=c(2,2))
qqp(bytod$lf_in)
qqp(bytod$lf_in, distribution = "lnorm")
poisson <- fitdistr(as.integer(bytod$lf_in), "Poisson"); qqp(bytod$lf_in, "pois", lambda = poisson$estimate) 
nbinom <- fitdistr(as.integer(bytod$lf_in), "Negative Binomial"); qqp(bytod$lf_in, "nbinom", size = nbinom$estimate[[1]], mu = nbinom$estimate[[2]]) ## BEST FIT DUE TO OVERDISPERSION
par(mfrow=c(1,1))



## ------------------------------------------------------------------------
## ZINB 
## --> Better approach given zero-inflation, althouh findings similiar to 
##     lognormal GLMM.



## DETERMINE ERROR STRUCTURE-----------------------------------------------
par(mfrow=c(2,2))
qqp(byvid$lf_in)
qqp(byvid$lf_in, distribution = "lnorm")
poisson <- fitdistr(as.integer(byvid$lf_in), "Poisson"); qqp(byvid$lf_in, "pois", lambda = poisson$estimate) 
nbinom <- fitdistr(as.integer(byvid$lf_in), "Negative Binomial"); qqp(byvid$lf_in, "nbinom", size = nbinom$estimate[[1]], mu = nbinom$estimate[[2]]) ## BEST FIT DUE TO OVERDISPERSION
par(mfrow=c(1,1))

lf_in_zi <- glmmTMB(lf_in ~ (Type * TOD) + (Type * LF_dens) + (1|Site_ID),  ## Positive component
                      zi =~ Type + LF_dens + (1|Site_ID), ## Zero component
                    data = bytod, 
                    family = nbinom1(link = "log"),
                    control = glmmTMBControl(optimizer = optim, optArgs = list(method="BFGS")))
summary(lf_in_zi)
effects <- fixef(lf_in_zi)
exp(effects$cond)
exp(effects$zi)
Anova(lf_in_zi) ##--> Analysis of deviance table only given for positive component

## Binomial (response = lionfish counts)
bytod$lf_binary = ifelse(bytod$lf_in > 0, 1 , 0)

lf_bi <- glmer(lf_binary ~ Type + LF_dens + (1|Site_ID),
             data = bytod, family = "binomial")
summary(lf_bi) ## Similar as ZINB zero component except estimates multiplied by -1
exp(lf_bi$coefficients) 
Anova(lf_bi)

##_END_OF_R_SCRIPT_________________________________________________________
