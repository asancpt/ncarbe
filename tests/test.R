library(tidyverse)

concData <- read.csv('inst/example/beConc.csv', stringsAsFactors = FALSE) %>% 
  rename('Period' = 'PRD')

head(concData)

#PRD = 'Period'

ncarbe::beNCA(concData, PRD = 'Period')
