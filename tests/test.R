library(tidyverse)

concData <- read.csv('inst/example/beConc.csv', stringsAsFactors = FALSE) %>%
  rename('Period' = 'PRD') %>% filter(SUBJ %in% c(5:10))



head(concData)

#PRD = 'Period'

ncarbe::beNCA(concData, PRD = 'Period')

# sas



dataBioeq <- read_csv('beNCAdataset.csv') %>%
  print()

model.Cmax <- nlme::lme(fixed = log(Cmax) ~ TRT + PRD + GRP,
                        random = ~1|SUBJ,
                        data = dataBioeq)

ci.Cmax <- nlme::intervals(model.Cmax, level=0.9, which="fixed")

model.AUClast <- nlme::lme(fixed = log(AUClast) ~ TRT + PRD + GRP,
                           random = ~1|SUBJ,
                           data = dataBioeq)

ci.AUClast <- nlme::intervals(model.AUClast, level=0.9, which="fixed")

confidenceInterval <- bind_rows(exp(ci.Cmax$fixed["TRTT", ]),
                                exp(ci.AUClast$fixed["TRTT", ])) %>%
  mutate(parameter = c('Cmax', 'AUClast')) %>%
  select(parameter,
         `Lower limit of 90% CI` = lower,
         `T/R ratio` = est.,
         `Upper limit of 90% CI` = upper)

ncarbestat <- list(Data = as.data.frame(dataBioeq),
                      Cmax = model.Cmax,
                      Summary.Cmax = summary(model.Cmax),
                      AUClast = model.AUClast,
                      Summary.AUClast = summary(model.AUClast),
                      Confidence.Interval = confidenceInterval)
print(ncarbestat)

#capture.output(yyNca$tadStat, file = 'docs/reports/nlme-report-Tad.txt')
