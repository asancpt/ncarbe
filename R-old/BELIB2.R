# R Libraries for Equivalence Test

# Diletti E, Hauscheke D, Steinjans VW. Sample size determination for bioequivalence assessment by means of confidence intervals. Int J Clin Pharmacology, Therapy and Toxicology. 1992;29:S51-58
# Chow SC, Liu JP. Design and Analysis of Bioavailability and Bioequivalence Study. 2nd ed. p158. 2000, Marcel Dekker Inc

##################################
# Sample Size for BE : Difference

bpow.d.cv = function(n1, n2, cv, delta=0, alpha=0.1, thetaL=-20, thetaU=20)
{
# n1, n2: Sample size of each group
# cv, delta in % scale (0 - 100%)
  if (delta <= thetaL | delta >= thetaU) return(0)
  t0 = qt(1 - alpha/2, n1 + n2 - 2)
  p1 = pt(-1*t0, n1 + n2 - 2, ncp=(thetaL - delta)/(cv*sqrt((1/n1 + 1/n2)/2)) )
  p2 = pt(t0, n1 + n2 - 2, ncp=(thetaU - delta)/(cv*sqrt((1/n1 + 1/n2)/2)) )
  power = p1 - p2
  if (power < 0) power = 0
  return(power)
}


bpow.d.mse = function(n1, n2, mu.r, mse, true.d=0, alpha=0.1, thetaL=-20, thetaU=20)
{
  cv = 100*sqrt(mse) / mu.r
  delta = 100*true.d / mu.r
  return(bpow.d.cv(n1, n2, cv, delta, alpha, thetaL, thetaU))
}


bss.d.cv = function(cv, delta=0, alpha=0.1, beta=0.2, thetaL=-20, thetaU=20, nmax=999999)
{
  if (delta <= thetaL | delta >= thetaU) return(Inf)
  for(i in 2:nmax) {
    power = bpow.d.cv(i, i, cv, delta, alpha, thetaL, thetaU)
    if (power > 1 - beta) return(i)
  }
  return(paste(">", nmax))
}


bss.d.mse = function(mu.r, mse, true.d=0, alpha=0.1, beta=0.2, thetaL=-20, thetaU=20, nmax=999999)
{
  cv = 100*sqrt(mse) / mu.r
  delta = 100*true.d / mu.r
  return(bss.d.cv(cv, delta, alpha, beta, thetaL, thetaU, nmax))
}

# LL, UL like 85% - 115%
# N: total N, 2 * n per group
bss.d.ci = function(n1, n2, LL, UL, alpha=0.1)
{
  pe = (LL + UL)/2
  t0 = qt(1 - alpha/2, n1 + n2 - 2)
  cv = (UL - LL)/(2*t0)/sqrt((1/n1 + 1/n2)/2)
  s1 = bss.d.cv(cv)
  s2 = bss.d.cv(cv, delta=(pe - 100))
  sampsize = cbind(s1, s2)
  p1 = round(100 * bpow.d.cv(n1, n2, cv))
  p2 = round(100 * bpow.d.cv(n1, n2, cv, delta=(pe - 100)))
  power = cbind(p1, p2)
  result = rbind(sampsize, power)
  dimnames(result) = list(c("80% Power Sample Size", paste("Power at n1 =", n1, ", n2 =", n2)), c("True Percent=100", sprintf("True Percent=%.2f", pe)))
  return(result)
}

############################
# Sample Size for BE : Ratio

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

# For the explanation of last number, see Chow & Liu 3e. p290
bpow.r.mse0 = function(n1, n2, mse, true.r=1, alpha=0.1, thetaL=0.8, thetaU=1.25) # 2x2 design
{
  if (true.r <= thetaL | true.r >= thetaU) return(0)
  nu = n1 + n2 - 2
  bx = sqrt(mse/2*(1/n1 + 1/n2))
  t0 = qt(1 - alpha/2, nu)
  p1 = pt(-1*t0, nu, ncp = log(thetaL/true.r)/bx) # 1st beta error: LL > 0.8 = 1 - p1
  p2 = pt(t0, nu, ncp = log(thetaU/true.r)/bx)    # 2nd beta error: UL < 1.25 = p2
  power = p1 - p2  # 1 - (1 - p1 + p2) = 1 - both beta errors
  if (power < 0) power = 0
  return(power)
}

bpow.r.mse1 = function(n1, n2, n3, n4, mse, true.r=1, alpha=0.1, thetaL=0.8, thetaU=1.25) # 4x2 design
{
  if (true.r <= thetaL | true.r >= thetaU) return(0)
  nu = n1 + n2 + n3 + n4 - 3
  bx = sqrt(mse/2*(1/n1 + 1/n2 + 1/n3 + 1/n4))
  t0 = qt(1 - alpha/2, nu)
  p1 = pt(-1*t0, nu, ncp = log(thetaL/true.r)/bx) # 1st beta error: LL > 0.8 = 1 - p1
  p2 = pt(t0, nu, ncp = log(thetaU/true.r)/bx)    # 2nd beta error: UL < 1.25 = p2
  power = p1 - p2  # 1 - (1 - p1 + p2) = 1 - both beta errors
  if (power < 0) power = 0
  return(power)
}

bpow.r.mse2 = function(n1, n2, mse, true.r=1, alpha=0.1, thetaL=0.8, thetaU=1.25) # 2x3 design
{
  if (true.r <= thetaL | true.r >= thetaU) return(0)
  nu = 2*(n1 + n2 - 2)
  bx = sqrt(mse/2*3/4*(1/n1 + 1/n2))
  t0 = qt(1 - alpha/2, nu)
  p1 = pt(-1*t0, nu, ncp = log(thetaL/true.r)/bx) # 1st beta error: LL > 0.8 = 1 - p1
  p2 = pt(t0, nu, ncp = log(thetaU/true.r)/bx)    # 2nd beta error: UL < 1.25 = p2
  power = p1 - p2  # 1 - (1 - p1 + p2) = 1 - both beta errors
  if (power < 0) power = 0
  return(power)
}

bpow.r.mse3 = function(n1, n2, mse, true.r=1, alpha=0.1, thetaL=0.8, thetaU=1.25) # 2x4 design
{
  if (true.r <= thetaL | true.r >= thetaU) return(0)
  nu = 3*(n1 + n2) - 5
  bx = sqrt(mse/2*11/20*(1/n1 + 1/n2))
  t0 = qt(1 - alpha/2, nu)
  p1 = pt(-1*t0, nu, ncp = log(thetaL/true.r)/bx) # 1st beta error: LL > 0.8 = 1 - p1
  p2 = pt(t0, nu, ncp = log(thetaU/true.r)/bx)    # 2nd beta error: UL < 1.25 = p2
  power = p1 - p2  # 1 - (1 - p1 + p2) = 1 - both beta errors
  if (power < 0) power = 0
  return(power)
}

bpow.r.mse4 = function(n1, n2, n3, n4, mse, true.r=1, alpha=0.1, thetaL=0.8, thetaU=1.25)  # 2x4x4 design (two treatment)
{
  if (true.r <= thetaL | true.r >= thetaU) return(0)
  nu = 3*(n1 + n2 + n3 + n4) - 5
  bx = sqrt(mse/2*1/8*(1/n1 + 1/n2 + 1/n3 + 1/n4))
  t0 = qt(1 - alpha/2, nu)
  p1 = pt(-1*t0, nu, ncp = log(thetaL/true.r)/bx) # 1st beta error: LL > 0.8 = 1 - p1
  p2 = pt(t0, nu, ncp = log(thetaU/true.r)/bx)    # 2nd beta error: UL < 1.25 = p2
  power = p1 - p2  # 1 - (1 - p1 + p2) = 1 - both beta errors
  if (power < 0) power = 0
  return(power)
}

bpow.r.mse5a = function(n1, n2, n3, n4, n5, n6, mse, true.r=1, alpha=0.1, thetaL=0.8, thetaU=1.25)  # 6x3 design, including carry-over effect
{
  if (true.r <= thetaL | true.r >= thetaU) return(0)
  nu = 2*(n1 + n2 + n3 + n4 + n5 + n6 - 3)
  bx = sqrt(mse/2*5/36*(1/n1 + 1/n2 + 1/n3 + 1/n4 + 1/n5 + 1/n6))
  t0 = qt(1 - alpha/2, nu)
  p1 = pt(-1*t0, nu, ncp = log(thetaL/true.r)/bx) # 1st beta error: LL > 0.8 = 1 - p1
  p2 = pt(t0, nu, ncp = log(thetaU/true.r)/bx)    # 2nd beta error: UL < 1.25 = p2
  power = p1 - p2  # 1 - (1 - p1 + p2) = 1 - both beta errors
  if (power < 0) power = 0
  return(power)
}

bpow.r.mse5b = function(n1, n2, n3, n4, n5, n6, mse, true.r=1, alpha=0.1, thetaL=0.8, thetaU=1.25)  # 6x3 design, excluding carry-over effect
{
  if (true.r <= thetaL | true.r >= thetaU) return(0)
  nu = 2*(n1 + n2 + n3 + n4 + n5 + n6 - 2)
  bx = sqrt(mse/2*1/9*(1/n1 + 1/n2 + 1/n3 + 1/n4 + 1/n5 + 1/n6))
  t0 = qt(1 - alpha/2, nu)
  p1 = pt(-1*t0, nu, ncp = log(thetaL/true.r)/bx) # 1st beta error: LL > 0.8 = 1 - p1
  p2 = pt(t0, nu, ncp = log(thetaU/true.r)/bx)    # 2nd beta error: UL < 1.25 = p2
  power = p1 - p2  # 1 - (1 - p1 + p2) = 1 - both beta errors
  if (power < 0) power = 0
  return(power)
}

bpow.r.mse6a = function(n1, n2, n3, n4, mse, true.r=1, alpha=0.1, thetaL=0.8, thetaU=1.25)  # 4x4x4 design (four treatment), including carry-over effect
{
  if (true.r <= thetaL | true.r >= thetaU) return(0)
  nu = 3*(n1 + n2 + n3 + n4 - 3)
  bx = sqrt(mse/2*11/40*(1/n1 + 1/n2 + 1/n3 + 1/n4))
  t0 = qt(1 - alpha/2, nu)
  p1 = pt(-1*t0, nu, ncp = log(thetaL/true.r)/bx) # 1st beta error: LL > 0.8 = 1 - p1
  p2 = pt(t0, nu, ncp = log(thetaU/true.r)/bx)    # 2nd beta error: UL < 1.25 = p2
  power = p1 - p2  # 1 - (1 - p1 + p2) = 1 - both beta errors
  if (power < 0) power = 0
  return(power)
}

bpow.r.mse6b = function(n1, n2, n3, n4, mse, true.r=1, alpha=0.1, thetaL=0.8, thetaU=1.25)  # 4x4x4 design (four treatment), excluding carry-over effect
{
  if (true.r <= thetaL | true.r >= thetaU) return(0)
  nu = 3*(n1 + n2 + n3 + n4 - 2)
  bx = sqrt(mse/2*11/40*(1/n1 + 1/n2 + 1/n3 + 1/n4))
  t0 = qt(1 - alpha/2, nu)
  p1 = pt(-1*t0, nu, ncp = log(thetaL/true.r)/bx) # 1st beta error: LL > 0.8 = 1 - p1
  p2 = pt(t0, nu, ncp = log(thetaU/true.r)/bx)    # 2nd beta error: UL < 1.25 = p2
  power = p1 - p2  # 1 - (1 - p1 + p2) = 1 - both beta errors
  if (power < 0) power = 0
  return(power)
}

bpow.r.cv = function(n1, n2, cv, true.r=1, alpha=0.1, thetaL=0.8, thetaU=1.25)
{
  mse = log(1 + (cv/100)^2)
  return(bpow.r.mse(n1, n2, mse, true.r, alpha, thetaL, thetaU))
}

bpow.r.cv0 = function(n1, n2, cv, true.r=1, alpha=0.1, thetaL=0.8, thetaU=1.25)
{
  mse = log(1 + (cv/100)^2)
  return(bpow.r.mse0(n1, n2, mse, true.r, alpha, thetaL, thetaU))
}

bpow.r.cv1 = function(n1, n2, n3, n4, cv, true.r=1, alpha=0.1, thetaL=0.8, thetaU=1.25)
{
  mse = log(1 + (cv/100)^2)
  return(bpow.r.mse1(n1, n2, n3, n4, mse, true.r, alpha, thetaL, thetaU))
}

bpow.r.cv2 = function(n1, n2, cv, true.r=1, alpha=0.1, thetaL=0.8, thetaU=1.25)
{
  mse = log(1 + (cv/100)^2)
  return(bpow.r.mse2(n1, n2, mse, true.r, alpha, thetaL, thetaU))
}

bpow.r.cv3 = function(n1, n2, cv, true.r=1, alpha=0.1, thetaL=0.8, thetaU=1.25)
{
  mse = log(1 + (cv/100)^2)
  return(bpow.r.mse3(n1, n2, mse, true.r, alpha, thetaL, thetaU))
}

bpow.r.cv4 = function(n1, n2, n3, n4, cv, true.r=1, alpha=0.1, thetaL=0.8, thetaU=1.25)
{
  mse = log(1 + (cv/100)^2)
  return(bpow.r.mse4(n1, n2, n3, n4, mse, true.r, alpha, thetaL, thetaU))
}

bpow.r.cv5a = function(n1, n2, n3, n4, n5, n6, cv, true.r=1, alpha=0.1, thetaL=0.8, thetaU=1.25)
{
  mse = log(1 + (cv/100)^2)
  return(bpow.r.mse5a(n1, n2, mse, true.r, alpha, thetaL, thetaU))
}

bpow.r.cv5b = function(n1, n2, n3, n4, n5, n6, cv, true.r=1, alpha=0.1, thetaL=0.8, thetaU=1.25)
{
  mse = log(1 + (cv/100)^2)
  return(bpow.r.mse5b(n1, n2, n3, n4, n5, n6, mse, true.r, alpha, thetaL, thetaU))
}

bpow.r.cv6a = function(n1, n2, n3, n4, cv, true.r=1, alpha=0.1, thetaL=0.8, thetaU=1.25)
{
  mse = log(1 + (cv/100)^2)
  return(bpow.r.mse6a(n1, n2, n3, n4, mse, true.r, alpha, thetaL, thetaU))
}

bpow.r.cv6b = function(n1, n2, n3, n4, cv, true.r=1, alpha=0.1, thetaL=0.8, thetaU=1.25)
{
  mse = log(1 + (cv/100)^2)
  return(bpow.r.mse6b(n1, n2, n3, n4, mse, true.r, alpha, thetaL, thetaU))
}

bpow.r.ci = function(n, n1, n2, LL, UL)
{
# n : sample size per group in new study
# n1, n2 : sample size of each group in old study
  pe = exp((log(UL) + log(LL))/2)
  t0 = qt(0.95, n1 + n2 - 2)
  sd = (log(UL) - log(LL))/(2*t0)
  mse = 2*sd^2/(1/n1 + 1/n2)
  pow1 = bpow.r.mse(n, n, mse)
  pow2 = bpow.r.mse(n, n, mse, true.r=pe)
  row1 = cbind(1, pow1)
  row2 = cbind(pe, pow2)
  result = rbind(row1, row2)
  colnames(result) = c("True Ratio", "Power")
  return(result)
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

bss.r.mse1 = function(mse, true.r=1, alpha=0.1, beta=0.2, thetaL=0.8, thetaU=1.25, nmax=999999)
{
  if (true.r <= thetaL | true.r >= thetaU) return(Inf)
  for (i in 2:nmax) {
    power = bpow.r.mse1(i, i, i, i, mse, true.r, alpha, thetaL, thetaU)
    if (power > 1 - beta) return(i)
  }
  return(paste(">", nmax))
}

bss.r.mse2 = function(mse, true.r=1, alpha=0.1, beta=0.2, thetaL=0.8, thetaU=1.25, nmax=999999)
{
  if (true.r <= thetaL | true.r >= thetaU) return(Inf)
  for (i in 2:nmax) {
    power = bpow.r.mse2(i, i, mse, true.r, alpha, thetaL, thetaU)
    if (power > 1 - beta) return(i)
  }
  return(paste(">", nmax))
}

bss.r.mse3 = function(mse, true.r=1, alpha=0.1, beta=0.2, thetaL=0.8, thetaU=1.25, nmax=999999)
{
  if (true.r <= thetaL | true.r >= thetaU) return(Inf)
  for (i in 2:nmax) {
    power = bpow.r.mse3(i, i, mse, true.r, alpha, thetaL, thetaU)
    if (power > 1 - beta) return(i)
  }
  return(paste(">", nmax))
}

bss.r.mse4 = function(mse, true.r=1, alpha=0.1, beta=0.2, thetaL=0.8, thetaU=1.25, nmax=999999)
{
  if (true.r <= thetaL | true.r >= thetaU) return(Inf)
  for (i in 2:nmax) {
    power = bpow.r.mse4(i, i, i, i, mse, true.r, alpha, thetaL, thetaU)
    if (power > 1 - beta) return(i)
  }
  return(paste(">", nmax))
}

bss.r.mse5a = function(mse, true.r=1, alpha=0.1, beta=0.2, thetaL=0.8, thetaU=1.25, nmax=999999)
{
  if (true.r <= thetaL | true.r >= thetaU) return(Inf)
  for (i in 2:nmax) {
    power = bpow.r.mse5a(i, i, i, i, i, i, mse, true.r, alpha, thetaL, thetaU)
    if (power > 1 - beta) return(i)
  }
  return(paste(">", nmax))
}

bss.r.mse5b = function(mse, true.r=1, alpha=0.1, beta=0.2, thetaL=0.8, thetaU=1.25, nmax=999999)
{
  if (true.r <= thetaL | true.r >= thetaU) return(Inf)
  for (i in 2:nmax) {
    power = bpow.r.mse5b(i, i, i, i, i, i, mse, true.r, alpha, thetaL, thetaU)
    if (power > 1 - beta) return(i)
  }
  return(paste(">", nmax))
}

bss.r.mse6a = function(mse, true.r=1, alpha=0.1, beta=0.2, thetaL=0.8, thetaU=1.25, nmax=999999)
{
  if (true.r <= thetaL | true.r >= thetaU) return(Inf)
  for (i in 2:nmax) {
    power = bpow.r.mse6a(i, i, i, i, mse, true.r, alpha, thetaL, thetaU)
    if (power > 1 - beta) return(i)
  }
  return(paste(">", nmax))
}

bss.r.mse6b = function(mse, true.r=1, alpha=0.1, beta=0.2, thetaL=0.8, thetaU=1.25, nmax=999999)
{
  if (true.r <= thetaL | true.r >= thetaU) return(Inf)
  for (i in 2:nmax) {
    power = bpow.r.mse6b(i, i, i, i, mse, true.r, alpha, thetaL, thetaU)
    if (power > 1 - beta) return(i)
  }
  return(paste(">", nmax))
}


bss.r.cv = function(cv, true.r=1, alpha=0.1, beta=0.2, thetaL=0.8, thetaU=1.25)
{
  mse = log(1 + (cv/100)^2)
  return(bss.r.mse(mse, true.r, alpha, beta, thetaL, thetaU))
}

bss.r.cv1 = function(cv, true.r=1, alpha=0.1, beta=0.2, thetaL=0.8, thetaU=1.25)
{
  mse = log(1 + (cv/100)^2)
  return(bss.r.mse1(mse, true.r, alpha, beta, thetaL, thetaU))
}

bss.r.cv2 = function(cv, true.r=1, alpha=0.1, beta=0.2, thetaL=0.8, thetaU=1.25)
{
  mse = log(1 + (cv/100)^2)
  return(bss.r.mse2(mse, true.r, alpha, beta, thetaL, thetaU))
}

bss.r.cv3 = function(cv, true.r=1, alpha=0.1, beta=0.2, thetaL=0.8, thetaU=1.25)
{
  mse = log(1 + (cv/100)^2)
  return(bss.r.mse3(mse, true.r, alpha, beta, thetaL, thetaU))
}

bss.r.cv4 = function(cv, true.r=1, alpha=0.1, beta=0.2, thetaL=0.8, thetaU=1.25)
{
  mse = log(1 + (cv/100)^2)
  return(bss.r.mse4(mse, true.r, alpha, beta, thetaL, thetaU))
}

bss.r.cv5a = function(cv, true.r=1, alpha=0.1, beta=0.2, thetaL=0.8, thetaU=1.25)
{
  mse = log(1 + (cv/100)^2)
  return(bss.r.mse5a(mse, true.r, alpha, beta, thetaL, thetaU))
}

bss.r.cv5b = function(cv, true.r=1, alpha=0.1, beta=0.2, thetaL=0.8, thetaU=1.25)
{
  mse = log(1 + (cv/100)^2)
  return(bss.r.mse5b(mse, true.r, alpha, beta, thetaL, thetaU))
}

bss.r.cv6a = function(cv, true.r=1, alpha=0.1, beta=0.2, thetaL=0.8, thetaU=1.25)
{
  mse = log(1 + (cv/100)^2)
  return(bss.r.mse6a(mse, true.r, alpha, beta, thetaL, thetaU))
}

bss.r.cv6b = function(cv, true.r=1, alpha=0.1, beta=0.2, thetaL=0.8, thetaU=1.25)
{
  mse = log(1 + (cv/100)^2)
  return(bss.r.mse6b(mse, true.r, alpha, beta, thetaL, thetaU))
}

ci2cv = function(n1, n2, LL, UL, alpha=0.1)
{
  pe = exp((log(UL) + log(LL))/2)
  t0 = qt(1 - alpha/2, n1 + n2 - 2)
  sd = (log(UL) - log(LL))/(2*t0)
  mse = 2*sd^2/(1/n1 + 1/n2)
  cv = sqrt(exp(mse) - 1)
  return(100*cv)
}

# LL, UL like 0.85 ~ 1.15
# N: toal N, 2 * n per group
bss.r.ci = function(n1, n2, LL, UL, alpha=0.1)
{
  pe = exp((log(UL) + log(LL))/2)
  t0 = qt(1 - alpha/2, n1 + n2 - 2)
  sd = (log(UL) - log(LL))/(2*t0)
  mse = 2*sd^2/(1/n1 + 1/n2)
  s1 = bss.r.mse(mse)
  s2 = bss.r.mse(mse, true.r=pe)
  sampsize = cbind(s1, s2)
  p1 = round(100 * bpow.r.mse(n1, n2, mse))
  p2 = round(100 * bpow.r.mse(n1, n2, mse, true.r=pe))
  power = cbind(p1, p2)
  result = rbind(sampsize, power)
  dimnames(result) = list(c("80% Power Sample Size", paste("Power at n1 =", n1, ", n2 =", n2)), c("True Ratio=1", sprintf("True Ratio=%.4f", pe)))
  return(result)
}

#########################
# 2x2 BE Test

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

  Y..1 = mean(Yi.1);
  Y..2 = mean(Yi.2);

  mu.r = (Y.11 + Y.22) / 2
  mu.t = (Y.21 + Y.12) / 2

  SScarry   = 2*n1*n2/(n1 + n2)*(Y.12 + Y.22 - Y.11 - Y.21)^2 / 4
#  SSinter   = (sum(Yi.1^2) + sum(Yi.2^2) - Y..1^2 * n1  - Y..2^2 * n2) * 2
  SSinter   = (sum((Yi.1 - Y..1)^2) + sum((Yi.2 - Y..2)^2)) * 2
  SSbetween = SScarry + SSinter

  SSperiod  = 2*n1*n2/(n1+n2)*(Y.21 + Y.22 - Y.11 - Y.12)^2 / 4
  SSdrug    = 2*n1*n2/(n1+n2)*(Y.21 + Y.12 - Y.11 - Y.22)^2 / 4
  SSintra   = SStotal - SScarry - SSinter - SSdrug - SSperiod

  Source = c("SUBJECT", "GROUP", "SUBJECT(GROUP)", "PERIOD", "DRUG", "ERROR", "TOTAL");
  SS     = c(SSbetween, SScarry, SSinter, SSperiod, SSdrug, SSintra, SStotal);
  DF     = c(n1 + n2 - 1, 1, n1 + n2 - 2, 1, 1, n1 + n2 - 2, 2*n1 + 2*n2 - 1);
  MS     = SS / DF
  mse    = SSintra / (n1 + n2 - 2);
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



########################################
# BE Plot


drawind = function(g1l, g1r, g2l, g2r, g1s, g2s)
{
  for (i in 1:length(g1l)) {
    x = jitter(c(1, 2), factor=0.3)
    y = c(g1l[i], g1r[i])
    lines(x, y, type="l", lty=1, col="red")
    text(x[1]-0.05, y[1], paste(g1s[i]), cex=0.6, col="red")
  }

  for (i in 1:length(g2l)) {
    x = jitter(c(1, 2), factor=0.3)
    y = c(g2l[i], g2r[i])
    lines(x, y, type="l", lty=2, col="blue")
    text(x[2]+0.05, y[2], paste(g2s[i]), cex=0.6, col="blue")
  }
}

drawmeansd = function(ma, sa, mb, sb, mc, sc, md, sd, y.max)
{
  sft = 0.03
  delta = mean(ma, mc) - mean(mb, md)
  y.RT = mean(ma, mc) + sign(delta) * y.max * 0.05
  y.TR = mean(mb, md) - sign(delta) * y.max * 0.05

  lines(c(1-sft, 2-sft), c(ma, mc), type="l", lty=1, col="red")
  text(1.5-sft, y.RT, "RT", col="red")
  if (sa > 0) arrows(1-sft, ma-sa, 1-sft, ma+sa, length=0.1, code=3, angle=90, col="red")
  if (sc > 0) arrows(2-sft, mc-sc, 2-sft, mc+sc, length=0.1, code=3, angle=90, col="red")

  lines(c(1+sft, 2+sft), c(mb, md), type="l", lty=2, col="blue")
  text(1.5+sft, y.TR, "TR", col="blue")
  if (sb > 0) arrows(1+sft, mb-sb, 1+sft, mb+sd, length=0.1, code=3, angle=90, col="blue")
  if (sd > 0) arrows(2+sft, md-sd, 2+sft, md+sd, length=0.1, code=3, angle=90, col="blue")
}

beplot = function(bedata, var, id = 'beplot')
{
  if(!assert(bedata)) {
    cat("\n Bad Data Format !\n");
    return(NULL);
  }

  Si11 = bedata[bedata$GRP=="RT" & bedata$PRD==1, "SUBJ"]
  Si21 = bedata[bedata$GRP=="RT" & bedata$PRD==2, "SUBJ"]
  Si12 = bedata[bedata$GRP=="TR" & bedata$PRD==1, "SUBJ"]
  Si22 = bedata[bedata$GRP=="TR" & bedata$PRD==2, "SUBJ"]

  Yi11 = bedata[bedata$GRP=="RT" & bedata$PRD==1, var]
  Yi21 = bedata[bedata$GRP=="RT" & bedata$PRD==2, var]
  Yi12 = bedata[bedata$GRP=="TR" & bedata$PRD==1, var]
  Yi22 = bedata[bedata$GRP=="TR" & bedata$PRD==2, var]

  n1 = length(Yi11)
  n2 = length(Yi12)

  Y.11 = mean(Yi11)
  Y.21 = mean(Yi21)
  Y.12 = mean(Yi12)
  Y.22 = mean(Yi22)

  sY.11 = sd(Yi11)
  sY.21 = sd(Yi21)
  sY.12 = sd(Yi12)
  sY.22 = sd(Yi22)

  y.max = max(Y.11 + sY.11, Y.21 + sY.21, Y.12 + sY.12, Y.22 + sY.22, max(bedata[,var])) * 1.2

  #windows()
  png(paste0('assets/beplot/', id, '-', var, '-equivalence.png'),
      width = 8, height = 6, res = 600, units = 'in')
  par(oma=c(1,1,3,1), mfrow=c(2,2))

  plot(0, 0, type="n", ylim=c(0, y.max), xlim=c(0.5, 2.5), axes=FALSE, xlab="Period",  ylab=var, main="(a) Individual Plot for Period")
  axis(2)
  axis(1, at=c(1,2))
  drawind(Yi11, Yi21, Yi12, Yi22, Si11, Si12)

  plot(0, 0, type="n", ylim=c(0, y.max), xlim=c(0.5, 2.5), axes=FALSE, xlab="Treatment",  ylab=var, main="(b) Individual Plot for Treatment")
  axis(2)
  axis(1, at=c(1,2), labels=c("Test", "Reference"))
  drawind(Yi21, Yi11, Yi12, Yi22, Si11, Si12)

  plot(0, 0, type="n", ylim=c(0, y.max), xlim=c(0.5, 2.5), axes=FALSE, xlab="Period",  ylab=var, main="(c) Mean and SD by Period")
  axis(2)
  axis(1, at=c(1,2))
  drawmeansd(Y.11, sY.11, Y.12, sY.12, Y.21, sY.21, Y.22, sY.22, y.max)

  plot(0, 0, type="n", ylim=c(0, y.max), xlim=c(0.5, 2.5), axes=FALSE, xlab="Treatment",  ylab=var, main="(d) Mean and SD by Treatment")
  axis(2)
  axis(1, at=c(1,2), labels=c("Test", "Reference"))
  drawmeansd(Y.21, sY.21, Y.12, sY.12, Y.11, sY.11, Y.22, sY.22, y.max)

  mtext(outer=T, side=3, paste("Equivalence Plot for", var), cex=1.5)
  dev.off()

  #windows()
  png(paste0('assets/beplot/', id, '-', var, '-box.png'),
      width = 8, height = 6, res = 600, units = 'in')
  par(oma=c(1,1,3,1), mfrow=c(2,2))

  boxplot(Yi11, Yi21, Yi12, Yi22, names=c("PRD=1", "PRD=2", "PRD=1", "PRD=2"), main="(a) By Sequence and Period", sub="SEQ=RT           SEQ=TR")
  boxplot(c(Yi11, Yi21), c(Yi12, Yi22), names=c("Sequence=RT", "Sequence=TR"), main="(b) By Sequence")
  boxplot(c(Yi11, Yi12), c(Yi21, Yi22), names=c("Period=1", "Period=2"), main="(c) By Period")
  boxplot(c(Yi12, Yi21), c(Yi11, Yi22), names=c("Treatment=T", "Treatment=R"), main="(d) By Treatment")
  mtext(outer=T, side=3, paste("Box Plots for", var), cex=1.5)
  dev.off()

}
# bedata = read.csv("d:/csv/propofolbe.csv")
# windows()
# par(mfrow=c(2,2),oma=c(1,1,3,1))
# boxplot(Cmax ~ GRP + PRD, data=bedata)
# boxplot(Cmax ~ GRP, data=bedata)
# boxplot(Cmax ~ PRD, data=bedata)
# boxplot(Cmax ~ TRT, data=bedata)

# options(digits=3)

be = function(filename)
{
  bedata = read.csv(filename);
# File should have the following columns
# SUBJ : Subject ID, any data type
# GRP: "RT" or "TR"
# PRD: 1 or 2
# TRT: "R" or "T"
# AUClast: numeric data type
# AUCinf: numeric data type
# Cmax: numeric data type
# Tmax: numeric data type
# Other columns as you wish

  bedata = bedata[order(bedata$GRP, bedata$PRD, bedata$SUBJ),];
  if(!assert(bedata)) {
    cat("\n Bad Data Format !\n");
    return(NULL);
  }

  beplot(bedata, "AUClast")
  beplot(bedata, "AUCinf")
  beplot(bedata, "Cmax")
  beplot(bedata, "Tmax")

  bedata$lnAUClast = log(bedata$AUClast);
  bedata$lnAUCinf = log(bedata$AUCinf);
  bedata$lnCmax = log(bedata$Cmax);

  cat("\n\n[AUClast]\n\n");
  print(betest(bedata, "lnAUClast", logtransformed=T));

  cat("\n\n[AUCinf]\n\n");
  print(betest(bedata, "lnAUCinf", logtransformed=T));

  cat("\n\n[Cmax]\n\n");
  print(betest(bedata, "lnCmax", logtransformed=T));

  cat("\n\n[Tmax]\n\n");
  print(hodges(bedata, "Tmax"));
}

kbe = function(filename)
{
  bedata = read.csv(filename);
# File should have the following columns
# SUBJ : Subject ID, any data type
# GRP: "RT" or "TR"
# PRD: 1 or 2
# TRT: "R" or "T"
# AUClast: numeric data type
# Cmax: numeric data type
# Tmax: numeric data type
# Other columns as you wish

  bedata = bedata[order(bedata$GRP, bedata$PRD, bedata$SUBJ),];
  if(!assert(bedata)) {
    cat("\n Bad Data Format !\n");
    return(NULL);
  }

  beplot(bedata, "AUClast")
  beplot(bedata, "Cmax")
  beplot(bedata, "Tmax")

  bedata$lnAUClast = log(bedata$AUClast);
  bedata$lnCmax = log(bedata$Cmax);

  cat("\n\n[AUClast]\n\n");
  print(betest(bedata, "lnAUClast", logtransformed=T));

  cat("\n\n[Cmax]\n\n");
  print(betest(bedata, "lnCmax", logtransformed=T));

  cat("\n\n[Tmax]\n\n");
  print(hodges(bedata, "Tmax"));
}

be.bs = function(filename, n=24, N=100, out="")
{
  all.data = read.csv(filename)
  all.subj.id = unique(all.data[,"SUBJ"])

  x.ncol = length(out) + 1
  x = matrix(nrow=length(all.subj.id), ncol=x.ncol)
  x[,x.ncol] = FALSE

  for (i in 1:length(out)) {
    x[,i] = (as.character(all.subj.id) == as.character(out[i]))
    x[,x.ncol] = x[,x.ncol] | x[,i]
  }

  ss.subj = all.subj.id[!x[,x.ncol]]

  result = matrix(ncol=6, nrow=0)
  colnames(result) = c("AUClast LB", "AUClast PE", "AUClast UB", "Cmax LB", "Cmax PE", "Cmax UB")

  for (r in 1:N) {
    sel.id = sample(ss.subj, n, replace=T)
    sel.data = matrix(nrow=0, ncol=ncol(all.data))
    colnames(sel.data) = colnames(all.data)

    for (i in 1:length(sel.id)) {
      y = all.data[all.data[,"SUBJ"] == sel.id[i], ]
      y[,"SUBJ"] = i
      sel.data = rbind(sel.data, y)
    }
    sel.data$lnAUClast = log(sel.data[,"AUClast"])
    sel.data$lnCmax = log(sel.data[,"Cmax"])

    AUClast.R = betest2(sel.data, "lnAUClast", logtransformed=T)
    Cmax.R = betest2(sel.data, "lnCmax", logtransformed=T)
    result = rbind(result, cbind(AUClast.R, Cmax.R))
  }

  return(result)
}

# be.bs("Bedata.csv", n=24, N=100, out=c("13"))

# x = be.bs("Bedata.csv", n=38, N=100, out="13")
# x

# nrow(x[x[,1] > 0.8 & x[,3] < 1.25 & x[,4] > 0.8 & x[,6] <1.25,])
