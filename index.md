[Home](https://paho-ghe.github.io/PAHO/)

## 1. Background 
Excess death = expected death - observed death 
{% include example1.html %}

## 2. Method 
Expected death forecasting  : Using GHE data from 2000 - 2019 (WHO) 
### Call dataset
```r
setwd("./")
GHE_PAHO<-read.csv("./GHE_PAHO.csv")
GHE_PAHO_NDeaths_AgeGroup     <- GHE_PAHO[GHE_PAHO$ghecause==0,-c(5,6,8,9)]
GHE_PAHO_NDeaths              <- data.frame(xtabs(dths~iso3+year, GHE_PAHO_NDeaths_AgeGroup))
names(GHE_PAHO_NDeaths)[3]    <- "dths"
#GHE_PAHO_NDeaths_wide         <- tidyr::spread(GHE_PAHO_NDeaths, year, dths)
write.csv(GHE_PAHO_NDeaths,"./GHE20002019.csv")

#GHE_PAHO_population_AgeGroup  <- GHE_PAHO[GHE_PAHO$ghecause==0,-c(5,7:9)]
#GHE_PAHO_population           <- data.frame(xtabs(pop~iso3+year, GHE_PAHO_population_AgeGroup))
#names(GHE_PAHO_population)[3] <- "pop"
#GHE_PAHO_population_wide      <- tidyr::spread(GHE_PAHO_population, year, pop)

GHE_PAHO_NDeaths<-GHE_PAHO_NDeaths%>%arrange(iso3, year)
```
### Random walk model to predict expected number of deaths in 2020 
```r

year <- c(2020)
iso3 <- c( "ARG" ,"ATG", "BHS", "BLZ", "BOL", "BRA" ,"BRB", "CAN", "CHL", "COL", "CRI", "CUB", "DOM" ,"ECU" ,"GRD" ,"GTM" ,"GUY" ,"HND" ,"HTI" ,"JAM" ,"LCA" ,"MEX", "NIC", "PAN" ,"PER" ,"PRY", "SLV", "SUR", "TTO", "URY", "USA" ,"VCT", "VEN")
dths<-c(NA)
df_fill = expand.grid(iso3 = iso3,year = year,dths=dths)
total <- rbind(df_fill, GHE_PAHO_NDeaths)%>%arrange(iso3, year)

# Code referenced from https://msemburi.github.io/#data-input
expected.ghe <- function(){
  set.seed(12345)
  iso.list   <- sort(unique(total$iso3))  
  ex.list    <- list()                        

for (i in 1:33){
#  print(paste0(round(100*i/33,1), "%"))
  group3.subset <- total %>% filter(iso3 == iso.list[i]) %>% 
                    arrange(year) %>% select(iso3, year, dths) %>% mutate(lall = log(dths))
 fit <- inla(lall ~ f(year, model = "rw2"), data = group3.subset, control.predictor= list(compute=TRUE))
 fit.out <- data.table(expected = exp(fit$summary.fitted.values$mean),
                              fit.l    = exp(fit$summary.fitted.values$"0.025quant"),
                              fit.u    = exp(fit$summary.fitted.values$"0.975quant"))
    ex.list[[i]]<- cbind(group3.subset %>% mutate(lall = NULL), fit.out)
}
rbindlist(ex.list)
}

expectedout <- expected.ghe()
write.csv(expectedout, file = "./expectedout.csv", row.names = F)

```
![download (1)](https://user-images.githubusercontent.com/81782228/126556506-7c46e283-fb0e-4408-939c-95b3b14578a3.png)



Retrieve observed death in year 2020. If missing, replace with value driven from the model using extra dataset with explanatory variables as below 
  - SDI index data : [source](http://ghdx.healthdata.org/sites/default/files/record-attached-files/IHME_GBD_2019_SDI_1990_2019_Y2020M10D15.XLSX)
  - WHO COVID19 data: [source](https://covid19.who.int/WHO-COVID-19-global-data.csv)
  - Oxford goverment policy data: [source](https://www.nature.com/articles/s41562-021-01079-8)
  - Karolinsky World mortality data: [source](https://raw.github.com/akarlinsky/world_mortality/)
  - Our world data: [source](https://covid.ourworldindata.org/data/owid-covid-data.csv)
  - World bank income group data: [source](http://databank.worldbank.org/data/download/site-content/CLASS.xls)



```r
# Aggregating from quarters to the sum for the entire year

df.lin      <- GHE_PAHO_NDeaths %>% 
  group_by(Country, iso3, WHO_region, Nx, ncdr, sdi) %>%
  summarize(covid = sum(cov_deaths), 
            expected = sum(expected), 
            observed = sum(observed), 
            .groups = "drop") %>%
  ungroup() %>% 
  mutate(excess    = observed - expected, 
         expectedc = expected + covid, 
         ratio     = observed/expectedc, 
         covidr    = covid/Nx, 
         covsdi    = covidr * sdi) %>%
  arrange(iso3)

# Indices to identify rows
df.lin      <- df.lin %>% mutate(row = 1:nrow(df.lin))
ind.obsl    <- df.lin %>% filter(!is.na(excess)) %>% pull(row)


reg.mod     <- lm(log(ratio) ~ . - 1, data = df.lin %>% 
                    dplyr::select(ratio, sdi, covidr, covsdi, ncdr))
pe          <- reg.mod$coefficients  # point estimates
vc          <- vcov(reg.mod)         # variance-covariance

set.seed(1234)
draws       <- 4000

# Distribution of coefficients
simbetas    <- MASS::mvrnorm(draws, pe, vc)

# Input data to scale the coefficients
df.lin.res  <- df.lin %>% dplyr::select(dimnames(simbetas)[[2]]) %>% as.matrix()

# Distribution of fitted values (these are scaling factors from the ratio of observed to covid plus expected)
replic      <- exp(df.lin.res %*% t(simbetas))

# Use the scaling factors to obtain distribution of excess deaths
fitdist1e   <- replic*replicate(draws, df.lin$expectedc) - replicate(draws, df.lin$expected)
obsvals     <- df.lin$excess # vector of observed values

# For the non-missing, replace predictions with the actual observed
for (i in 1:draws){
  fitdist1e[ind.obsl,] <- obsvals[ind.obsl]
}
# Basic functions to summarise 80% bootstrap intervals or sample from gaussian

getsum <- function(x){
  tibble(mean = quantile(x, 0.5, na.rm = T), 
         lwr  = quantile(x, 0.025, na.rm = T), 
         uppr = quantile(x, 0.975, na.rm = T)) %>% round()
}

getsum2 <- function(x){
  c(quantile(x, 0.5, na.rm = T), 
    quantile(x, 0.025, na.rm = T), 
    quantile(x, 0.975, na.rm = T))
}

get.sample <- function(x){
  rnorm(draws, x[1], x[2])
}

df.lin.out   <- data.table::data.table(t(apply(fitdist1e, 1, getsum2))) %>%
  rename(y = "50%", low = "2.5%", high = "97.5%") %>%
  cbind(df.lin)  %>% 
  mutate(y    = ifelse(!is.na(excess), excess, y), #mean value
         low  = ifelse(!is.na(excess), NA, low),
         high = ifelse(!is.na(excess), NA, high), 
         source = ifelse(is.na(excess), "Predicted", "Observed")) %>% arrange(-y) %>%
  data.frame()
```

```r
# Identify the rows for each region to aggregate the matrices accordingly
analysis_mod_paho<-analysis_mod %>% filter(WHO_region == "AMRO")
analysis_mod_paho<-analysis_mod_paho%>% mutate(row = 1:nrow(analysis_mod_paho))

# Get a summary of each aggregate, by draw and by across aggregate

dlply(analysis_mod_paho$Country, function(a) a.ind   <- analysis_mod_paho %>% filter(Country == a) %>% pull(row))
dlply(analysis_mod_paho$Country, function(a) a.sum.exc   <- getsum(apply(as.data.frame(t(as.matrix(fit.dist3e[a.ind,]))), 2, sum)) %>%
                    mutate(Country = a,
                           covid = sum(analysis_mod_paho[a.ind,]$cov_deaths), 
                            Nx = .5*sum(analysis_mod_paho[a.ind,]$Nx))  #For missing and non-missing countries with observed values
                            
            
)

sum.exc.exc <- rbind(Argentina.sum.exc  ,AntiguaandBarbudaa.sum.exc   ,Bahamas.sum.exc  ,Belize.sum.exc  ,Barbados.sum.exc ,Dominica.sum.exc  ,DominicanRepublic.sum.exc  ,Grenada.sum.exc   ,Guyana.sum.exc  ,Honduras.sum.exc  ,Haiti.sum.exc  ,SaintKittsandNevis.sum.exc  ,SaintLucia.sum.exc  ,Nicaragua.sum.exc  ,ElSalvador.sum.exc   ,Suriname.sum.exc  ,TrinidadandTobago.sum.exc   ,Uruguay.sum.exc  ,SaintVincentandtheGrenadines.sum.exc  ,VenezuelaBolivarianRepublicof.sum.exc  ,BoliviaPlurinationalStateof.sum.exc   ,Brazil.sum.exc   ,Canada.sum.exc  ,Chile.sum.exc  ,Colombia.sum.exc   ,CostaRica.sum.exc  ,Cuba.sum.exc  ,Ecuador.sum.exc   ,Guatemala.sum.exc   ,Jamaica.sum.exc   ,Mexico.sum.exc  ,Panama.sum.exc  ,Peru.sum.exc   ,Paraguay.sum.exc  ,UnitedStatesofAmerica.sum.exc  
) %>% 
  rename(low = lwr, high = uppr) %>%
  mutate(excess = round(mean - covid), 
         mean = round(mean), low = round(low), high = round(high),
         region = factor(Country, levels = c( "Argentina" ,"Antigua and Barbuda","Bahamas" ,"Belize","Bolivia (Plurinational State of)","Brazil","Barbados","Canada","Chile","Colombia","Costa Rica","Cuba","Dominica","Dominican Republic","Ecuador","Grenada","Guatemala","Guyana","Honduras","Haiti","Jamaica","Saint Kitts and Nevis","Saint Lucia","Mexico","Nicaragua","Panama","Peru","Paraguay","El Salvador","Suriname","Trinidad and Tobago","Uruguay","United States of America","Saint Vincent and the Grenadines","Venezuela (Bolivarian Republic of)")), 
         model = "Mixed Effects") %>%
  dplyr::select(model, region, Nx, covid, excess, low, mean, high)
```
## 3. Result 



This page's data and code are forked from https://msemburi.github.io/#background and edited. 
