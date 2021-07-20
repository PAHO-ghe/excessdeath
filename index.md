[Home](https://paho-ghe.github.io/PAHO/)

## 1. Background 
Excess death = expected death - observed death 

## 2. Method 
Expected death forecasting (GAN or ARIMA) : Using GHE data from 2000 - 2019 (WHO) 
Retrieve observed death in year 2020. If missing, replace with value driven from the model using extra dataset with explanatory variables as below 
  - SDI index data  
  - GBD data
  - Oxford national policy data
  - Karolinsky baseline data 
  - Our world data
  - World bank income group data

## 3. Result 
