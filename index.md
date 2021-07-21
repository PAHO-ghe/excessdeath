[Home](https://paho-ghe.github.io/PAHO/)

## 1. Background 
Excess death = expected death - observed death 
{% include example1.html %}

![download](https://user-images.githubusercontent.com/81782228/126429763-6c879790-408d-4e57-8f7d-0f50c687c73a.png)


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
  - SDI index data  
  - GBD data
  - Oxford national policy data
  - Karolinsky baseline data 
  - Our world data
  - World bank income group data

## 3. Result 

