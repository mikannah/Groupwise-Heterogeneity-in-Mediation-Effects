### R code for following the analysis in section 4
### of Kim & Kim (2025), IMPS 2024 Proceedings
### after gaining access to the HSIS data set.
### The data set is not attached because of confidentiality.


## 1. Load libraries ====

library(psych)
library(foreign)
library(mice)
library(mitools)
library(boot)


## 2. Read and clean HSIS data ====

# Read HSIS data
Dat_copy <- read.spss("D://ICPSR_29462/DS0006/29462.sav",
                      to.data.frame = TRUE)  # HSIS data externally stored

# Re-code variables
Dat_copy$CRgroup <- 
  ifelse(Dat_copy$CHILDRESULTGROUP == "Treatment", 1, 0)
table(Dat_copy$CHILDRESULTGROUP, Dat_copy$CRgroup)  # ok

Dat_copy$bbio <- ifelse(Dat_copy$D_BBIO_WITHCH == "Yes", 1, 0)
table(Dat_copy$D_BBIO_WITHCH, Dat_copy$bbio)  # ok

Dat_copy$cds.num <- as.numeric(Dat_copy$D_CDS) - 1
table(Dat_copy$D_CDS, Dat_copy$cds.num)  # ok

Dat_copy$female <- 
  ifelse(Dat_copy$D_CF_CHSEX01 == "Female", 1, 0)
table(Dat_copy$D_CF_CHSEX01, Dat_copy$female)  # ok

Dat_copy$home.eng <- 
  ifelse(Dat_copy$D_LANG_HM == "English", 1, 0)
table(Dat_copy$D_LANG_HM, Dat_copy$home.eng)  # ok

Dat_copy$mbio <- ifelse(Dat_copy$D_MBIO_RECIMM == "Yes", 1, 0)
table(Dat_copy$D_MBIO_RECIMM, Dat_copy$mbio)  # ok.

Dat_copy$special <- 
  ifelse(Dat_copy$D_PI_SPECNEEDS_F02 == "Not special needs", 
         0, 1)
table(Dat_copy$D_PI_SPECNEEDS_F02, Dat_copy$special)  # ok

Dat_copy$urban <- ifelse(Dat_copy$D_URBAN == "Urban", 1, 0)
table(Dat_copy$D_URBAN, Dat_copy$urban)  # ok

Dat_copy$lowquart <- 
  ifelse(Dat_copy$D_LQUART == "Lowest quartile", 1, 0)
table(Dat_copy$D_LQUART, Dat_copy$lowquart)  # ok

# Data set to be used in multiple imputation
Dat_imp <- Dat_copy[, c("CRgroup", "Black", 
                        "Hispanic", "FOCARR_Y1", "FOCARR_Y2", 
                        "IRT_PPVTML_DS02", "IRT_PPVTML_DS04",
                        "IRT_PPVTML_DS13", "D_1CGVAGE", "bbio",
                        "cds.num", "female", "home.eng", 
                        "mbio", "Mo_LessHS", "Mo_HS",
                        "Mo_NeverMarried", 
                        "Mo_CurrentlyMarried", "special",
                        "D_RISK3_resc", "urban", "lowquart")]
# Note: FOCARR_Y1 - Early HS, FOCARR_Y2 - Regular HS
describe(Dat_imp)  # Table 1 - Observed data set


## 3. Multiple imputation of missing data ====

# Run mice
set.seed(2023)
Dat_mice <- mice(data = Dat_imp, m = 20, method = "pmm", 
                 printFlag = FALSE)  # imputed data sets

# Multiple imputation data descriptives
mice.imp <- NULL
for(i in 1:20) { mice.imp[[i]] <- complete(Dat_mice, 
                                           action = i,
                                           inc = FALSE) }
describe(mice.imp[[1]])  # Table 1 - Analysis data set

Dat_imputed <- mice.imp[[1]]  # Analysis data set


## 4. Groupwise NDE/NIE estimation with bootstrap SEs ====

# 4.1 Function for fitting models to bootstrap samples ----

b.stat2 <- function(data, i) {
  # draw a bootstrap sample (done by the boot() function)
  bdat <- data[i, ]
  
  # mediator model
  mreg.out <- glm(FOCARR_Y2 ~ FOCARR_Y1 * home.eng + 
                    IRT_PPVTML_DS02 + Black + Hispanic + female +
                    lowquart + special + urban + 
                    D_1CGVAGE + bbio + cds.num + mbio + Mo_LessHS + Mo_HS +
                    Mo_NeverMarried + Mo_CurrentlyMarried + D_RISK3_resc, 
                  family = binomial, data = bdat)
  
  # outcome model
  yreg.out <- lm(IRT_PPVTML_DS13 ~ FOCARR_Y1 * FOCARR_Y2 * home.eng +
                   IRT_PPVTML_DS02 + Black + Hispanic + female +
                   lowquart + special + urban +
                   D_1CGVAGE + bbio + cds.num + mbio + Mo_LessHS + Mo_HS +
                   Mo_NeverMarried + Mo_CurrentlyMarried + D_RISK3_resc,
                 data = bdat)
  
  # Extract regression coefficients
  mreg <- mreg.out$coefficients
  yreg <- yreg.out$coefficients
  
  beta0 <- as.numeric(mreg[1])
  beta1 <- as.numeric(mreg[2])
  beta2.vec <- matrix(mreg[3], nrow = 1, ncol = 1, 
                      byrow = TRUE)
  beta3.vec <- matrix(mreg[4:19], nrow = 1, ncol = 16, 
                      byrow = TRUE)
  beta4.vec <- matrix(mreg[20], nrow = 1, ncol = 1, 
                      byrow = TRUE)
  
  theta0 <- as.numeric(yreg[1])
  theta1 <- as.numeric(yreg[2])
  theta2 <- as.numeric(yreg[3])
  theta3.vec <- matrix(yreg[4], nrow = 1, ncol = 1, 
                       byrow = TRUE)
  theta4.vec <- matrix(yreg[5:20], nrow = 1, ncol = 16, 
                       byrow = TRUE)
  theta5 <- as.numeric(yreg[21])
  theta6.vec <- matrix(yreg[22], nrow = 1, ncol = 1, 
                       byrow = TRUE)
  theta7.vec <- matrix(yreg[23], nrow = 1, ncol = 1, 
                       byrow = TRUE)
  theta8.vec <- matrix(yreg[24], nrow = 1, ncol = 1, 
                       byrow = TRUE)
  
  # Conditioning values
  G1 <- matrix(c(1), nrow = 1, ncol = 1, 
               byrow = TRUE)  # English
  G0 <- matrix(c(0), nrow = 1, ncol = 1, byrow = TRUE)  # DLL
  CE <- matrix(sapply(bdat[bdat$home.eng == 1, 
                           c("IRT_PPVTML_DS02", "Black", "Hispanic", 
                             "female", "lowquart", "special", 
                             "urban", "D_1CGVAGE", "bbio", "cds.num",
                             "mbio", "Mo_LessHS", "Mo_HS", "Mo_NeverMarried",
                             "Mo_CurrentlyMarried", "D_RISK3_resc")], mean),
               nrow = 16, ncol = 1, 
               byrow = TRUE)  # covariate means for ENG
  CD <- matrix(sapply(bdat[bdat$home.eng == 0, 
                           c("IRT_PPVTML_DS02", "Black", "Hispanic", 
                             "female", "lowquart", "special", 
                             "urban", "D_1CGVAGE", "bbio", "cds.num",
                             "mbio", "Mo_LessHS", "Mo_HS", "Mo_NeverMarried",
                             "Mo_CurrentlyMarried", "D_RISK3_resc")], mean),
               nrow = 16, ncol = 1, 
               byrow = TRUE)  # covariate means for DLL
  
  # NDE(0): Pure NDE
  a0 <- matrix(c(0), nrow = 1, ncol = 1,
               byrow = TRUE)  # For NDE(A = a0 = 0)
  # mediator prob for NDE(0) calculation
  P.m0.eng <- exp(beta0 + beta1 * a0 + beta2.vec %*% G1 + 
                    beta3.vec %*% CE + 
                    beta4.vec %*% a0 %*% G1) / 
    (1 + exp(beta0 + beta1 * a0 + beta2.vec %*% G1 + 
               beta3.vec %*% CE + beta4.vec %*% a0 %*% G1))
  NDE.0.eng <- (theta1 + theta6.vec %*% G1) + 
    (theta5 + theta8.vec %*% G1) * P.m0.eng
  
  # repeat the steps for the DLL group
  P.m0.DLL <- exp(beta0 + beta1 * a0 + beta2.vec %*% G0 + 
                    beta3.vec %*% CD + 
                    beta4.vec %*% a0 %*% G0) / 
    (1 + exp(beta0 + beta1 * a0 + beta2.vec %*% G0 + 
               beta3.vec %*% CD + beta4.vec %*% a0 %*% G0))
  NDE.0.DLL <- (theta1 + theta6.vec %*% G0) + 
    (theta5 + theta8.vec %*% G0) * P.m0.DLL
  
  # NDE(1): Total NDE
  a1 <- matrix(c(1), nrow = 1, ncol = 1,
               byrow = TRUE)  # For NDE(A = a1 = 1)
  # mediator prob for NDE(1) calculation
  P.m1.eng <- exp(beta0 + beta1 * a1 + beta2.vec %*% G1 + 
                    beta3.vec %*% CE + 
                    beta4.vec %*% a1 %*% G1) / 
    (1 + exp(beta0 + beta1 * a1 + beta2.vec %*% G1 + 
               beta3.vec %*% CE + beta4.vec %*% a1 %*% G1))
  NDE.1.eng <- (theta1 + theta6.vec %*% G1) + 
    (theta5 + theta8.vec %*% G1) * P.m1.eng
  
  # repeat the procedure for the DLL group
  P.m1.DLL <- exp(beta0 + beta1 * a1 + beta2.vec %*% G0 + 
                    beta3.vec %*% CD + 
                    beta4.vec %*% a1 %*% G0) / 
    (1 + exp(beta0 + beta1 * a1 + beta2.vec %*% G0 + 
               beta3.vec %*% CD + beta4.vec %*% a1 %*% G0))
  NDE.1.DLL <- (theta1 + theta6.vec %*% G0) + 
    (theta5 + theta8.vec %*% G0) * P.m1.DLL
  
  # NIE(0): Pure NIE
  P.m_diff.eng <- P.m1.eng - P.m0.eng
  NIE.0.eng <- (theta2 + theta5 %*% a0 + theta7.vec %*% G1 + 
                  theta8.vec %*% a0 %*% G1) * P.m_diff.eng
  P.m_diff.DLL <- P.m1.DLL - P.m0.DLL
  NIE.0.DLL <- (theta2 + theta5 %*% a0 + theta7.vec %*% G0 + 
                  theta8.vec %*% a0 %*% G0) * P.m_diff.DLL
  
  # NIE(1): Total NIE
  NIE.1.eng <- (theta2 + theta5 %*% a1 + theta7.vec %*% G1 + 
                  theta8.vec %*% a1 %*% G1) * P.m_diff.eng
  NIE.1.DLL <- (theta2 + theta5 %*% a1 + theta7.vec %*% G0 + 
                  theta8.vec %*% a1 %*% G0) * P.m_diff.DLL
  
  # TE
  TE.eng <- NDE.1.eng + NIE.0.eng
  TE.DLL <- NIE.1.DLL + NDE.0.DLL
  
  # a vector of estimated values
  est.vec <- c('TE_eng' = TE.eng, 
               'NDE(0)_eng' = NDE.0.eng, 
               'NIE(1)_eng' = NIE.1.eng, 
               'NDE(1)_eng' = NDE.1.eng, 
               'NIE(0)_eng' = NIE.0.eng,
               'TE_DLL' = TE.DLL, 
               'NDE(0)_DLL' = NDE.0.DLL, 
               'NIE(1)_DLL' = NIE.1.DLL, 
               'NDE(1)_DLL' = NDE.1.DLL, 
               'NIE(0)_DLL' = NIE.0.DLL)
  return(est.vec)
}

# 4.2 Get bootstrap outcomes ----
set.seed(2023)  # random seed number for replication

b.out <- boot(Dat_imputed, statistic = b.stat2, R = 10000, 
              stype = 'i')
se <- apply(b.out$t, 2, sd)
ci <- 
  rbind('TE_eng' = boot.ci(b.out, index = 1, 
                           type = 'perc')$percent[4:5],
        'NDE(0)_eng' = boot.ci(b.out, index = 2, 
                               type = 'perc')$percent[4:5],
        'NIE(1)_eng' = boot.ci(b.out, index = 3, 
                               type = 'perc')$percent[4:5],
        'NDE(1)_eng' = boot.ci(b.out, index = 4, 
                               type = 'perc')$percent[4:5],
        'NIE(0)_eng' = boot.ci(b.out, index = 5, 
                               type = 'perc')$percent[4:5],
        'TE_DLL' = boot.ci(b.out, index = 6, 
                           type = 'perc')$percent[4:5],
        'NDE(0)_DLL' = boot.ci(b.out, index = 7, 
                               type = 'perc')$percent[4:5],
        'NIE(1)_DLL' = boot.ci(b.out, index = 8, 
                               type = 'perc')$percent[4:5],
        'NDE(1)_DLL' = boot.ci(b.out, index = 9, 
                               type = 'perc')$percent[4:5],
        'NIE(0)_DLL' = boot.ci(b.out, index = 10, 
                               type = 'perc')$percent[4:5])
colnames(ci) <- c('nparCI(.025)', 'nparCI(.975)')
boot.res <- cbind(b.out$t0, se, ci)
colnames(boot.res) <- c('Estimate', 'SE', 'nparCI(.025)', 
                        'nparCI(.975)')
round(boot.res, 3)  # Table 4


## 5. Group-constant NDE/NIE point estimation ====

# 5.1 Function for fitting models to bootstrap samples ----

b.stat21 <- function(data, i) {
  # draw bootstrap sample (done by the boot() function)
  bdat <- data[i, ]
  
  # mediator model
  mreg.out <- glm(FOCARR_Y2 ~ FOCARR_Y1 + IRT_PPVTML_DS02 + Black + Hispanic + female +
                    lowquart + home.eng + special + urban + 
                    D_1CGVAGE + bbio + cds.num + mbio + Mo_LessHS + Mo_HS +
                    Mo_NeverMarried + Mo_CurrentlyMarried + D_RISK3_resc, 
                  family = binomial, data = bdat)
  # outcome model
  yreg.out <- lm(IRT_PPVTML_DS13 ~ FOCARR_Y1 * FOCARR_Y2 +
                   IRT_PPVTML_DS02 + Black + Hispanic + female +
                   lowquart + home.eng + special + urban +
                   D_1CGVAGE + bbio + cds.num + mbio + Mo_LessHS + Mo_HS +
                   Mo_NeverMarried + Mo_CurrentlyMarried + D_RISK3_resc, 
                 data = bdat)
  
  # Extract regression coefficients
  mreg <- mreg.out$coefficients
  yreg <- yreg.out$coefficients
  
  beta0 <- as.numeric(mreg[1])
  beta1 <- as.numeric(mreg[2])
  beta3.vec <- matrix(mreg[3:19], nrow = 1, ncol = 17, 
                      byrow = TRUE)
  
  theta0 <- as.numeric(yreg[1])
  theta1 <- as.numeric(yreg[2])
  theta2 <- as.numeric(yreg[3])
  theta4.vec <- matrix(yreg[4:20], nrow = 1, ncol = 17, 
                       byrow = TRUE)
  theta5 <- as.numeric(yreg[21])
  
  # Conditioning values
  C <- 
    matrix(sapply(bdat[, c("IRT_PPVTML_DS02", "Black", "Hispanic", 
                           "female", "lowquart", "home.eng", "special", 
                           "urban", "D_1CGVAGE", "bbio", "cds.num",
                           "mbio", "Mo_LessHS", "Mo_HS", "Mo_NeverMarried",
                           "Mo_CurrentlyMarried", "D_RISK3_resc")], mean),
           nrow = 17, ncol = 1, 
           byrow = TRUE)  # covariate means
  
  # NDE(0): Pure NDE
  a0 <- matrix(c(0), nrow = 1, ncol = 1,
               byrow = TRUE)  # For NDE(A = a0 = 0)
  # mediator prob for NDE(0)
  P.m0 <- exp(beta0 + beta1 * a0 + beta3.vec %*% C) / 
    (1 + exp(beta0 + beta1 * a0 + beta3.vec %*% C))
  NDE.0 <- (theta1) + (theta5) * P.m0
  
  # NDE(1): Total NDE
  a1 <- matrix(c(1), nrow = 1, ncol = 1,
               byrow = TRUE)  # For NDE(A = a1 = 1)
  # mediator prob for NDE(1)
  P.m1 <- exp(beta0 + beta1 * a1 + beta3.vec %*% C) / 
    (1 + exp(beta0 + beta1 * a1 + beta3.vec %*% C))
  NDE.1 <- (theta1) + (theta5) * P.m1
  
  # NIE(0): Pure NIE
  P.m_diff <- P.m1 - P.m0
  NIE.0 <- (theta2 + theta5 * a0) * P.m_diff
  
  # NIE(1): Total NIE
  NIE.1 <- (theta2 + theta5 * a1) * P.m_diff
  
  # TE
  TE <- NDE.1 + NIE.0
  
  # a vector of estimated values
  est.vec <- c('TE' = TE, 'NDE(0)' = NDE.0, 'NIE(1)' = NIE.1,
               'NDE(1)' = NDE.1, 'NIE(0)' = NIE.0)
  return(est.vec)
}

# 5.2 Get bootstrap outcomes ----
set.seed(1220)  # random seed for replication

b.out <- boot(Dat_imputed, statistic = b.stat21, R = 10000, 
              stype = 'i')
se <- apply(b.out$t, 2, sd)
ci <- rbind('TE' = boot.ci(b.out, index = 1, 
                           type = 'perc')$percent[4:5],
            'NDE(0)' = boot.ci(b.out, index = 2, 
                               type = 'perc')$percent[4:5],
            'NIE(1)' = boot.ci(b.out, index = 3, 
                               type = 'perc')$percent[4:5],
            'NDE(1)' = boot.ci(b.out, index = 4, 
                               type = 'perc')$percent[4:5],
            'NIE(0)' = boot.ci(b.out, index = 5, 
                               type = 'perc')$percent[4:5])
colnames(ci) <- c('nparCI(.025)', 'nparCI(.975)')
boot.res <- cbind(b.out$t0, se, ci)
colnames(boot.res) <- c('Estimate', 'SE', 'nparCI(.025)', 
                        'nparCI(.975)')
round(boot.res, 3)  # Table 5
