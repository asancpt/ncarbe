#' Performs a statistical analysis of 2x2 bioequivalence study
#'
#' \code{beNCA} returns text output of a statistical analysis of 2x2 bioequivalence study. Analysis of variance, least square means, confidence interval, and sample size will be calculated for AUClast and Cmax. Wilcoxon Signed-Rank Test and Hodges-Lehmann Estimate will be performed for Tmax.
#'
#' Reference: Shein-Chung Chow, Jen-pei Liu. Design and Analysis of Bioavailability and Bioequivalence Studies, 3rd ed. 2008. (ISBN:9781584886686)
#'
#' @param SUBJ Subject ID, any data type
#' @param GRP column name in which information of "RT" or "TR" exists.
#' @param PRD column name in which information of 1 or 2 exists.
#' @param TRT column name in which information of "R" or "T" exists.
#' @param method \code{kbe} by authors or \code{nlme} package uploaded on CRAN
#' @return returns text results of statistical analysis of 2x2 bioequivalence study including and \code{beNCAdataset.csv} file in the working directory which can be run in SAS.
#' @import NonCompart
#' @import dplyr
#' @importFrom nlme intervals lme
#' @export
#' @examples
#' file <- system.file('example', 'beConc.csv', package = 'ncarbe')
#' concData <- read.csv(file, as.is = TRUE)
#' beNCA(concData)

beNCA <- function(concData, SUBJ = 'SUBJ', GRP = 'GRP', PRD = 'PRD', TRT = 'TRT', method = 'kbe', ...){
  # SUBJ = 'SUBJ'; GRP = 'GRP';  PRD = 'PRD'; TRT = 'TRT'

  betestKey <- c('SUBJ', 'GRP', 'PRD', 'TRT')
  ncaKey <- setNames(c(SUBJ, GRP, PRD, TRT), betestKey)
  colnames(concData) <- sub(ncaKey['SUBJ'], 'SUBJ', colnames(concData))
  colnames(concData) <- sub(ncaKey['GRP'], 'GRP', colnames(concData))
  colnames(concData) <- sub(ncaKey['PRD'], 'PRD', colnames(concData))
  colnames(concData) <- sub(ncaKey['TRT'], 'TRT', colnames(concData))
  colnames(concData)

  bedataRaw <- NonCompart::tblNCA(as.data.frame(concData),
                                  key= betestKey,
                                  colTime="TIME", colConc="CONC", dose=100000, ..., R2ADJ = 0.5) %>%
    as.data.frame() %>%
    mutate_at(vars(AUCLST, CMAX, TMAX), as.numeric) %>%
    arrange(GRP, PRD, SUBJ) %>%
    select(SUBJ, TRT, GRP, PRD, AUClast = AUCLST, Cmax = CMAX, Tmax = TMAX)

  bedata <- bedataRaw %>%
    mutate(lnAUClast = log(AUClast),
           lnCmax = log(Cmax)) %>%
    as.data.frame()

  write.csv(bedataRaw, 'beNCAdataset.csv', row.names = FALSE)

  if(!assert(bedata)) {
    cat("\n Bad Data Format !\n");
    return(NULL);
  }

  if (method == 'kbe'){
    cat("\n\n[AUClast]\n\n");
    print(betest(bedata, "lnAUClast", logtransformed=TRUE));
    cat("\n\n[Cmax]\n\n");
    print(betest(bedata, "lnCmax", logtransformed=TRUE));
    cat("\n\n[Tmax]\n\n");
    print(hodges(bedata, "Tmax"))
  } else if (method == 'nlme') {

    dataBioeq <- bedata

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

    # ncarbestat <- list(Data = as.data.frame(dataBioeq),
    #                    Cmax = model.Cmax,
    #                    Summary.Cmax = summary(model.Cmax),
    #                    AUClast = model.AUClast,
    #                    Summary.AUClast = summary(model.AUClast),
    #                    Confidence.Interval = confidenceInterval)
    ncarbestat <- list(Confidence.Interval = confidenceInterval)
    print(ncarbestat)
  }
}

