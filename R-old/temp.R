
print(betest(bedata, "lnAUClast", logtransformed=T))





bedata = bedata[order(bedata$GRP, bedata$PRD, bedata$SUBJ),];
if(!assert(bedata)) {
  cat("\n Bad Data Format !\n");
  return(NULL);
}



# misc ----


ncares[c(ncaKey, 'lnAUClast', 'lnCmax', 'Tmax')]


bedata$lnAUClast = log(bedata$AUClast);
bedata$lnCmax = log(bedata$Cmax);


kbe()

left_join(as_tibble(ncares) %>% gather(PPTESTCD, PPORRES, ),
          tibble(PPTESTCD = attr(ncares, 'dimnames')[[2]], UNIT = attr(ncares, 'units')),
          by = 'PPTESTCD') %>% 
  arrange(Subject, PPTESTCD) %>% 
  head(20)


bedata$lnAUClast = log(bedata$AUClast);
bedata$lnCmax = log(bedata$Cmax);

cat("\n\n[AUClast]\n\n");
print(betest(bedata, "lnAUClast", logtransformed=T));




library(tidyverse)
Theoph %>% 
  as_tibble() %>% 
  mutate_at(vars(Subject), as.numeric) %>% 
  select(Subject, Time, conc) %>% 
  arrange(Subject)

NonCompart::tblNCA()
function(concData = Theoph, PRD, TRT, GRP, key = , ...){
  
}
  
read.csv()
  
tblNCA(concData, key = "Subject", colTime = "Time", colConc = "conc", dose = 0, 
         adm = "Extravascular", dur = 0, doseUnit = "mg", timeUnit = "h", 
         concUnit = "ug/L", down = "Linear", MW = 0)