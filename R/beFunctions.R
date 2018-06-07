assert = function(bedata)
{
  Si11 = bedata[bedata$GRP=="RT" & bedata$PRD==1, "SUBJ"]
  Si21 = bedata[bedata$GRP=="RT" & bedata$PRD==2, "SUBJ"]
  Si12 = bedata[bedata$GRP=="TR" & bedata$PRD==1, "SUBJ"]
  Si22 = bedata[bedata$GRP=="TR" & bedata$PRD==2, "SUBJ"]

  return(identical(Si11, Si21) & identical(Si12, Si22))
}

betest = function(bedata, var, logtransformed)
{
# GRP, PRD
  if(!assert(bedata)) {
    cat("\n Bad Data Format !\n");
    return(NULL);
  }

  Yijk = bedata[, var]
  Yi11 = bedata[bedata$GRP=="RT" & bedata$PRD==1, var]
  Yi21 = bedata[bedata$GRP=="RT" & bedata$PRD==2, var]
  Yi12 = bedata[bedata$GRP=="TR" & bedata$PRD==1, var]
  Yi22 = bedata[bedata$GRP=="TR" & bedata$PRD==2, var]

  n1 = length(Yi11)
  n2 = length(Yi12)

  Y... = mean(Yijk)
  SStotal = sum((Yijk - Y...)^2)

  Y.11 = mean(Yi11)
  Y.21 = mean(Yi21)
  Y.12 = mean(Yi12)
  Y.22 = mean(Yi22)

  Yi.1 = (Yi11 + Yi21) / 2
  Yi.2 = (Yi12 + Yi22) / 2

  Y..1 = mean(Yi.1)
  Y..2 = mean(Yi.2)

  mu.r = (Y.11 + Y.22) / 2
  mu.t = (Y.21 + Y.12) / 2

  di1 = (Yi21 - Yi11) / 2
  di2 = (Yi22 - Yi12) / 2

  d.1 = mean(di1)
  d.2 = mean(di2)

#  sig2d = 1/(n1 + n2 - 2)*(sum((di1 - d.1)^2) + sum((di2 - d.2)^2))

#  Chat = Y.12 + Y.22 - Y.11 - Y.21
#  Fhat = mu.t - mu.r
#  Phat = (Y.21 - Y.11 - Y.12 + Y.22)/2

  SScarry   = 2*n1*n2/(n1 + n2)*(Y.12 + Y.22 - Y.11 - Y.21)^2 / 4
  SSinter   = (sum((Yi.1 - Y..1)^2) + sum((Yi.2 - Y..2)^2)) * 2
  SSbetween = SScarry + SSinter

  SSperiod  = 2*n1*n2/(n1+n2)*(Y.21 + Y.22 - Y.11 - Y.12)^2 / 4
  SSdrug    = 2*n1*n2/(n1+n2)*(Y.21 + Y.12 - Y.11 - Y.22)^2 / 4

  SSintra  = 2*(sum((di1 - d.1)^2) + sum((di2 - d.2)^2))
#  SSmodel = SStotal - SSintra

  Source = c("SUBJECT", "GROUP", "SUBJECT(GROUP)", "PERIOD", "DRUG", "ERROR", "TOTAL");
  SS     = c(SSbetween, SScarry, SSinter, SSperiod, SSdrug, SSintra, SStotal);
  DF     = c(n1 + n2 - 1, 1, n1 + n2 - 2, 1, 1, n1 + n2 - 2, 2*n1 + 2*n2 - 1);
  MS     = SS / DF
  mse    = SSintra/(n1 + n2 - 2)
  F      = MS / c(mse, MS[3], mse, mse, mse, mse, mse);
  p1 = 1 - pf(F[1], n1 + n2 - 1, n1 + n2 - 2)
  p2 = 1 - pf(F[2], 1, n1 + n2 - 2);
  p3 = 1 - pf(F[3], n1 + n2 - 2, n1 + n2 - 2);
  p4 = 1 - pf(F[4], 1, n1 + n2 - 2);
  p5 = 1 - pf(F[5], 1, n1 + n2 - 2);
  p  = c(p1, p2, p3, p4, p5, NA, NA)
  F[6] = F[7] = MS[7] = NA

  ANOVA = cbind(SS, DF, MS, F, p)
  dimnames(ANOVA) = list(Source,c("SS", "DF", "MS", "F", "p"))

  pe = mu.t - mu.r
  se = sqrt(mse / 2 * (1/n1 + 1/n2))   # See pp 62-63
  t0 = qt(0.95, n1 + n2 - 2);
  ci0 = cbind(pe - t0 * se, pe, pe + t0 * se)

  sig2b = (MS[3] - MS[6])/2
  sig2w = MS[6]

  if (logtransformed == TRUE) {
    cvs = cbind(sqrt(exp(sig2b) - 1), sqrt(exp(sig2w) - 1)) * 100

    lsm = cbind(exp(mu.r), exp(mu.t))
    dimnames(lsm) = list("Geometric Means", cbind("Reference Drug", "Test Drug"))

    ci = exp(ci0);
    dimnames(ci) = list("90% CI for Ratio", c("Lower Limit", "Point Estimate", "Upper Limit"));

    sampsize1 = bss.r.mse(mse);
    sampsize2 = bss.r.mse(mse, true.r=exp(pe));
    ss = cbind(sampsize1, sampsize2)
    dimnames(ss) = list("80% Power Sample Size", c("True Ratio=1", "True Ratio=Point Estimate"));
  } else {
    cvs = cbind(sqrt(sig2b)/mu.r, sqrt(sig2w)/mu.r) * 100

    lsm = cbind(mu.r, mu.t)
    dimnames(lsm) = list("Arithmetic Means", cbind("Reference Drug", "Test Drug"))

    ci1 = (1 + ci0 / mu.r) * 100
    ci = rbind(ci0, ci1)
    dimnames(ci) = list(c("90% CI for Difference", "90% CI for Difference(%)"), c("Lower Limit", "Point Estimate", "Upper Limit"));

    sampsize1 = bss.d.mse(mu.r, mse);
    sampsize2 = bss.d.mse(mu.r, mse, true.d=pe);
    ss = cbind(sampsize1, sampsize2)
    dimnames(ss) = list("80% Power Sample Size", c("True Difference=0", "True Difference=Point Estimate"));
  }
  cvs = rbind(cbind(sig2b, sig2w), cvs)
  dimnames(cvs) = list(cbind("Variance Estimate", "Coefficient of Variation, CV(%)"), cbind("Between Subject", "Within Subject"))

  result = list(ANOVA, cvs, lsm, ci, ss);
  names(result) = c("Analysis of Variance", "Between and Within Subject Variability", "Least Square Means", "90% Confidence Interval", "Sample Size")

  return(result);
}

hodges = function(bedata, var)
{
  if(!assert(bedata)) {
    cat("\n Bad Data Format !\n");
    return(NULL);
  }

  Yi11 = bedata[bedata$GRP=="RT" & bedata$PRD==1, var]
  Yi21 = bedata[bedata$GRP=="RT" & bedata$PRD==2, var]
  Yi12 = bedata[bedata$GRP=="TR" & bedata$PRD==1, var]
  Yi22 = bedata[bedata$GRP=="TR" & bedata$PRD==2, var]

  n1 = length(Yi11)
  n2 = length(Yi12)

  if(n1 * n2 < 12) {
    cat("\n Too Small Sample Size for 90% Confidence Interval !\n");
    return(NULL);
  }

  mu.r = (mean(Yi11) + mean(Yi22)) / 2;

  G1D = (Yi21 - Yi11) / 2
  G2D = (Yi22 - Yi12) / 2
  D = sort(outer(G1D, G2D, "-"));

  pval = pwilcox(min(length(D[D>0]), length(D[D<0])), n1, n2)
  w05 = qwilcox(0.05, n1, n2)
  w95 = qwilcox(0.95, n1, n2)

  names(pval) = list(c("p-value"));

  est1 = cbind(D[w05 - 1], median(D), D[w95])
  est2 = (1 + est1 / mu.r ) * 100
  est.a = rbind(est1, est2)
  dimnames(est.a) = list(c("90% Confidence Interval", "90% Confidence Interval(%)"), c("Lower Limit", "Point Estimate", "Upper Limit"));

#  est3 = cbind(D[w05], median(D), D[w95+1])
#  est4 = (1 + est3 / mu.r ) * 100
#  est.b = rbind(est3, est4)
#  dimnames(est.b) = list(c("90% Confidence Interval", "90% Confidence Interval(%)"), c("Lower Limit", "Point Estimate", "Upper Limit"));

#  result = list(pval, est.a, est.b);
#  names(result) = c("Wilcoxon Signed-Rank Test", "Hodges-Lehmann Estimate", "Hodges-Lehmann Estimate Old")

  result = list(pval, est.a);
  names(result) = c("Wilcoxon Signed-Rank Test", "Hodges-Lehmann Estimate")
  return(result);
}


bss.r.mse = function(mse, true.r=1, alpha=0.1, beta=0.2, thetaL=0.8, thetaU=1.25, nmax=999999)
{
  if (true.r <= thetaL | true.r >= thetaU) return(Inf)
  for (i in 2:nmax) {
    power = bpow.r.mse(i, i, mse, true.r, alpha, thetaL, thetaU)
    if (power > 1 - beta) return(i)
  }
  return(paste(">", nmax))
}

bpow.r.mse = function(n1, n2, mse, true.r=1, alpha=0.1, thetaL=0.8, thetaU=1.25)
{
# n1, n2: sample size for each group
  if (true.r <= thetaL | true.r >= thetaU) return(0)
  t0 = qt(1 - alpha/2, n1 + n2 - 2)
  p1 = pt(-1*t0, n1 + n2 - 2, ncp = log(thetaL/true.r) / sqrt(mse/2*(1/n1 + 1/n2)) ) # 1st beta error: LL > 0.8 = 1 - p1
  p2 = pt(t0, n1 + n2 - 2, ncp = log(thetaU/true.r) / sqrt(mse/2*(1/n1 + 1/n2)) )    # 2nd beta error: UL < 1.25 = p2
  power = p1 - p2  # 1 - (1 - p1 + p2) = 1 - both beta errors
  if (power < 0) power = 0
  return(power)
}
