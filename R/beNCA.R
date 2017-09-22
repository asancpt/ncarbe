#' Performs a statistical analysis of 2x2 bioequivalence study
#' 
#' \code{beNCA} returns text output of a statistical analysis of 2x2 bioequivalence study. Analysis of variance, least square means, confidence interval, and sample size will be calculated for AUClast and Cmax. Wilcoxon Signed-Rank Test and Hodges-Lehmann Estimate will be performed for Tmax. 
#' 
#' * Reference: Shein-Chung Chow, Jen-pei Liu. Design and Analysis of Bioavailability and Bioequivalence Studies, 3rd ed. 2008. (ISBN:9781584886686)
#' 
#' @param SUBJ Subject ID, any data type
#' @param GRP column name in which information of "RT" or "TR" exists.
#' @param PRD column name in which information of 1 or 2 exists.
#' @param TRT column name in which information of "R" or "T" exists.
#' @return Text results of statistical analysis of 2x2 bioequivalence study including 
#' @import NonCompart
#' @import dplyr
#' @export
#' @examples
#' file <- system.file('example', 'beConc.csv', package = 'ncarbe')
#' concData <- read.csv(file, as.is = TRUE)
#' beNCA(concData)
#'  

beNCA <- function(concData, SUBJ = 'SUBJ', GRP = 'GRP', PRD = 'PRD', TRT = 'TRT', ...){
  # SUBJ = 'SUBJ'; GRP = 'GRP';  PRD = 'PRD'; TRT = 'TRT'
  
  betestKey <- c('SUBJ', 'GRP', 'PRD', 'TRT')
  ncaKey <- setNames(c(SUBJ, GRP, PRD, TRT), betestKey)
  colnames(concData) <- sub(ncaKey['SUBJ'], 'SUBJ', colnames(concData))
  colnames(concData) <- sub(ncaKey['GRP'], 'GRP', colnames(concData))
  colnames(concData) <- sub(ncaKey['PRD'], 'PRD', colnames(concData))
  colnames(concData) <- sub(ncaKey['TRT'], 'TRT', colnames(concData))
  colnames(concData)
  
  bedata <- NonCompart::tblNCA(as.data.frame(concData), 
                     key= betestKey, 
                     colTime="TIME", colConc="CONC", dose=100000, ...) %>% 
    as_tibble() %>% 
    mutate_at(vars(AUCLST, CMAX, TMAX), as.numeric) %>% 
    mutate(lnAUClast = log(AUCLST),
           lnCmax = log(CMAX),
           Tmax = TMAX) %>% 
    arrange(GRP, PRD, SUBJ) %>% 
    as.data.frame()
  
  if(!assert(bedata)) {
    cat("\n Bad Data Format !\n");
    return(NULL);
  }
  cat("\n\n[AUClast]\n\n");
  print(betest(bedata, "lnAUClast", logtransformed=T));
  cat("\n\n[Cmax]\n\n");
  print(betest(bedata, "lnCmax", logtransformed=T));
  cat("\n\n[Tmax]\n\n");
  print(hodges(bedata, "Tmax"))
}
