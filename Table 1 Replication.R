data<-read.csv("Data_Estimation.csv")
data

##double hashtag means Diana wrote this in. Singular hashtag is a note from ChatGPT when translating the code from Stata to R. 
##Install correct packages (data.table and fixest)(ex.: install,packages(data.table)), then open them with the library() function.
library(data.table)
library(fixest)

##First step: This step is preparing the corn yield data set to be able to be made into a table in R. The function "dt[]" allows us to view and modify sections of the dataset (Data_Estimation.csv) in order to make it work--both from Stata to R, and in order to make a table. 
setorder(dt,fips,year)
dt[, stateid := as.integer(factor(stateid, levels = sort(unique(stateid))))]
dt[, log_cornyield := log(vAllGrainYield)]
dt[, prec_sq := prec^2]
dt[, isob := as.integer(!is.na(log_cornyield) & !is.na(kdd) & east_meridian == 1)]
dt[, iobsum := sum(isob), by = fips]
dt[, tobsum := sum(isob), by = year]
dt[, sample_corn := as.integer(iobsum > 29 & tobsum > 29 & !is.na(log_cornyield) & !is.na(kdd))]
dt[, corn_weight := mean(vAllGrainHarvest, na.rm = TRUE), by = fips]
dt[!is.finite(corn_weight), corn_weight := NA_real_]


##Second step: this step begins creating the regression. The first part is selecting particular observations: in this case, it contains observations where sample_corn=1).The next part is the actual regression. The function "feols" is for fixed-effects regression. A fixed-effects regression means that we are controlling for certain variables. This one controls for county fixed effects (an example of this is soil quality, as stated by the authors) and common time (year, in this model). In the code, we see the corn yield data being used, the particular specifications they are calculating (GDD, KDD, Precipitation, Precipitation^2, and the controlled variables: county and time). 
s <- dt[sample_corn == 1]
m1 <- feols(log_cornyield ~ gdd + kdd + prec + prec_sq | fips + year,
            data = s, vcov = ~stateid)

##Third step: this step sets the time frame and the counter...as stated below, the translation from Stata to R required modifications. 
# State-specific linear and quadratic time trends (Stata quadtime_1...43).
# Stata's trend counter is year - 1945, even for observations before 1950.
s[, trend := year - 1945]
s[, trend2 := trend^2]

##Fourth step: this step includes the coding of another variable the authors account for which is the U.S. state the sample is from. They use "quadtime" which refers to a quadratic time trend. For Table 1, the second column "resists time effects to be quadratic trends at the State level" (Keane and Neal, 2020).
# Generate explicitly to retain the original Stata design matrix.
for (i in sort(unique(dt$stateid))) {
  v1 <- paste0('quadtime_', i)
  v2 <- paste0('quadtime2_', i)
  s[, (v1) := fifelse(stateid == i & year >= 1946 & year <= 2017, trend, 0)]
  s[, (v2) := get(v1)^2]
}
trend_vars <- grep('^quadtime', names(s), value = TRUE)
form2 <- as.formula(paste('log_cornyield ~ gdd + kdd + prec + prec_sq +',
                          paste(trend_vars, collapse = ' + '), '| fips'))
m2 <- feols(form2, data = s, vcov = ~stateid)

##Fifth step: This step creates the third column of the table, which examines extreme heat (kdd) and its relationship with corn yields (nonlinear).
# Adaptation specification; zeros of kdd become one before log transform.
s[, kdd1 := fifelse(kdd == 0, 1, kdd)]
s[, log_kdd1 := log(kdd1)]
s[, kdd1_integraloflogx := log_kdd1 * kdd1 - kdd1]
m3 <- feols(log_cornyield ~ gdd + kdd1 + kdd1_integraloflogx + prec + prec_sq |
              fips + year, data = s, vcov = ~stateid)

##Sixth step (table 1): This step actually puts together the table after organizing all of the data and variables in the previous steps. 
# Table 1: coefficients and standard errors, with the original six-row layout.
row_labels <- c('gdd','kdd','integraloflogx','prec','prec2','cons')
coeffs <- matrix(0, 6, 3, dimnames = list(row_labels, paste0('Model',1:3)))
stderrs <- coeffs
put <- function(model, column, variables, positions) {
  b <- coef(model); se <- sqrt(diag(vcov(model)))
  for (j in seq_along(variables)) {
    v <- variables[j]; p <- positions[j]
    if (v %in% names(b)) {
      coeffs[p, column] <<- unname(b[v])
      stderrs[p, column] <<- unname(se[v])
    }
  }
}
put(m1, 1, c('gdd','kdd','prec','prec_sq'), c(1,2,4,5))
put(m2, 2, c('gdd','kdd','prec','prec_sq'), c(1,2,4,5))
put(m3, 3, c('gdd','kdd1','kdd1_integraloflogx','prec','prec_sq'), 1:5)

# In fixed-effects estimation, the Stata _cons normalization is not
# interchangeable with fixest's absorbed-effect normalization.
coeffs['cons', ] <- NA_real_
stderrs['cons', ] <- NA_real_

##Seventh step: This step continues to put together the table and displays the results. 
# Stata computes squared correlation of xbu prediction with observed y.
# fixest fitted values include absorbed effects; this is not the same xbu
# normalization, so diagnostics below use within-model R2 as an explicitly
# labeled alternative rather than misrepresenting Stata's diagnostic.
diag <- rbind(within_R2 = c(fitstat(m1,'wr2')[[1]], fitstat(m2,'wr2')[[1]], fitstat(m3,'wr2')[[1]]),
              obs = c(nobs(m1), nobs(m2), nobs(m3)))
colnames(diag) <- paste0('Model',1:3)
cat('\nTable 1 coefficients (FE intercept not comparable):\n'); print(coeffs)
cat('\nTable 1 standard errors:\n'); print(stderrs)
cat('\nAlternative diagnostics (within R2, not Stata corr(xbu,y)^2):\n'); print(diag)

fwrite(as.data.table(coeffs, keep.rownames = 'term'), 'Table1_coefficients_R.csv')
fwrite(as.data.table(stderrs, keep.rownames = 'term'), 'Table1_standard_errors_R.csv')
saveRDS(list(m1=m1, m2=m2, m3=m3), 'Table1_models_R.rds')
cat('\nSaved Table 1 results to the current working directory.\n')