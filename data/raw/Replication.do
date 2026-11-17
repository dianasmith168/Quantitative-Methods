// 
// Replication file for "Climate Change and U.S. Agriculture: Accounting for
// Multidimensional Slope Heterogeneity in Panel Data" by Michael Keane and
// Timothy Neal
// 
// This file replicates the Tables and Figures of the paper. 
// 
// Requirements: 
// Stata 12+ 
// User-written command called "xtmo", included in the replication package.
// 

//
// Import Dataset and Set up Variables
//

qui {
local path = "C:/Keane and Neal (2020) - Replication Package/"

import delimited "`path'Data/Data_Estimation.csv", case(preserve) clear
xtset fips year
set more off
encode stateid, g(stateid2)
drop stateid
rename stateid2 stateid
set matsize 10000
gen log_cornyield = log(vAllGrainYield)
gen prec_sq = prec*prec
gen isob = 0
replace isob = 1 if !missing(log_cornyield) & !missing(kdd) & east_meridian == 1
bysort fips: egen iobsum = sum(isob)
bysort year: egen tobsum = sum(isob) 
gen sample_corn = 0
replace sample_corn = 1 if iobsum > 29 & tobsum > 29 & !missing(log_cornyield) & !missing(kdd)
bysort fips: egen corn_weight = mean(vAllGrainHarvest)

}

//
// Table 1 Replication
//

qui {
matrix Table1_coeffs = J(6,3,0)
matrix Table1_stderr = J(6,3,0)
matrix Table1_diag = J(2,3,0)
matrix rownames Table1_coeffs = gdd kdd integraloflogx prec prec2 cons
matrix rownames Table1_stderr = gdd kdd integraloflogx prec prec2 cons
matrix rownames Table1_diag = r2 obs

xtreg log_cornyield gdd kdd prec prec_sq i.year if sample_corn == 1, fe vce(cluster stateid) // Cty/Year fixed effects
matrix b = e(b)
predict yhat_1, xbu
predict pooled_fixedeffect, u
correlate yhat_1 log_cornyield if sample_corn == 1
matrix corr = r(C)
matrix Table1_diag[1,1] = corr[2,1]*corr[2,1]
matrix Table1_diag[2,1] = e(N)

local varlist = "gdd kdd prec prec_sq _cons"
local v = 1
foreach var in `varlist' {
	if `v' < 3 {
		matrix Table1_coeffs[`v',1] = _b[`var']
		matrix Table1_stderr[`v',1] = _se[`var']
	}
	else {
		matrix Table1_coeffs[`v'+1,1] = _b[`var']
		matrix Table1_stderr[`v'+1,1] = _se[`var']
	}
	local v = `v' + 1
}

qui forvalues i = 1/43 {
	gen quadtime_`i' = 0
	local j = 1
	forvalues t = 1946/2017 {
		replace quadtime_`i' = `j' if year == `t' & stateid == `i'
		local j = `j' + 1
	}
	gen quadtime2_`i' = quadtime_`i' * quadtime_`i'
}

xtreg log_cornyield gdd kdd prec prec_sq quad* if sample_corn == 1, fe vce(cluster stateid) // Cty/State-Yr fixed effects
matrix b = e(b)
predict yhat_2, xbu
correlate yhat_2 log_cornyield if sample_corn == 1
matrix corr = r(C)
matrix Table1_diag[1,2] = corr[2,1]*corr[2,1]
matrix Table1_diag[2,2] = e(N)
drop quadtime*

local varlist = "gdd kdd prec prec_sq _cons"
local v = 1
foreach var in `varlist' {
	if `v' < 3 {
		matrix Table1_coeffs[`v',2] = _b[`var']
		matrix Table1_stderr[`v',2] = _se[`var']
	}
	else {
		matrix Table1_coeffs[`v'+1,2] = _b[`var']
		matrix Table1_stderr[`v'+1,2] = _se[`var']
	}
	local v = `v' + 1
}

gen kdd1 = kdd
replace kdd1 = 1 if kdd == 0
gen log_kdd1 = log(kdd1)
gen kdd1_integraloflogx = log_kdd1*kdd1 - kdd1
xtreg log_cornyield gdd kdd1 kdd1_integraloflogx prec prec_sq i.year if sample_corn == 1, fe vce(cluster stateid) // Cty/Year with adaptation
matrix b = e(b)
predict yhat_3, xbu
predict modols_fixedeffect, u
correlate yhat_3 log_cornyield if sample_corn == 1
matrix corr = r(C)
matrix Table1_diag[1,3] = corr[2,1]*corr[2,1]
matrix Table1_diag[2,3] = e(N)

local varlist = "gdd kdd1 kdd1_integraloflogx prec prec_sq _cons"
local v = 1
foreach var in `varlist' {
	matrix Table1_coeffs[`v',3] = _b[`var']
	matrix Table1_stderr[`v',3] = _se[`var']
	local v = `v' + 1
}
noi display "Table 1 Coefficients:"
noi matlist Table1_coeffs

noi display "Table 1 Standard Errors:"
noi matlist Table1_stderr

noi display "Table 1 Diagnostics:"
noi matlist Table1_diag
}

//
// Table 2 Replication
//

qui {
matrix Table2_values = J(5,6,0)
matrix Table2_stderr = J(5,1,0)
matrix Table2_diag = J(2,1,0)

xtmo log_cornyield gdd kdd prec prec_sq if sample_corn == 1
gen yhatcorn_mo = bit_cons + bit_gdd*gdd + bit_kdd*kdd + bit_prec*prec + bit_prec_sq*prec_sq
correlate yhatcorn_mo log_cornyield if sample_corn == 1
matrix corr = r(C)
matrix Table2_diag[1,1] = corr[2,1]*corr[2,1]
matrix Table2_diag[2,1] = e(N)

local varlist = "gdd kdd prec prec_sq cons"
local v = 1
foreach var in `varlist' {
	su bit_`var', detail
	matrix Table2_values[`v',1] = r(mean)
	matrix Table2_values[`v',3] = r(p50)
	matrix Table2_values[`v',4] = r(sd)
	matrix Table2_values[`v',5] = r(p10)
	matrix Table2_values[`v',6] = r(p90)
	su bit_`var' [aweight=corn_weight], meanonly
	matrix Table2_values[`v',2] = r(mean)
	matrix Table2_stderr[`v',1] = _se[`var']

	local v = `v' + 1
}

matrix rownames Table2_values = GDD KDD Prec Prec2 Cons
matrix rownames Table2_stderr = GDD KDD Prec Prec2 Cons
matrix colnames Table2_values = Mean WMean Median StdDev 10thPct 90thPct
matrix rownames Table2_diag = R2 Obs

noi display "Table 2 values:"
noi matlist Table2_values

noi display "Table 2 Standard Errors:"
noi matlist Table2_stderr

noi display "Table 2 Diagnostics:"
noi matlist Table2_diag

save "`path'/Data/Replication_Temp.dta", replace
}

// 
// Figure 2 Replication
// 

qui{
matrix stats_1 = J(2015-1950+1,10,0)

local i = 1
forvalues t = 1950/2015 {
	qui su bit_kdd if year == `t', detail
	
	matrix stats_1[`i',1] = r(N)
	matrix stats_1[`i',2] = r(mean)
	matrix stats_1[`i',3] = r(min)
	matrix stats_1[`i',4] = r(max)
	matrix stats_1[`i',5] = r(p10)
	matrix stats_1[`i',6] = r(p25)
	matrix stats_1[`i',7] = r(p50)
	matrix stats_1[`i',8] = r(p75)
	matrix stats_1[`i',9] = r(p90)
	matrix stats_1[`i',10] = r(sd)

	local i = `i' + 1
} 
gen pct10_1 = . 
gen pct25_1 = . 
gen pct50_1 = . 
gen pct75_1 = . 
gen pct90_1 = . 
local j = 1
qui forvalues t = 1950/2015 {
	replace pct10_1 = stats_1[`j',5] if year == `t' & fips == 1001
	replace pct25_1 = stats_1[`j',6] if year == `t' & fips == 1001
	replace pct50_1 = stats_1[`j',7] if year == `t' & fips == 1001
	replace pct75_1 = stats_1[`j',8] if year == `t' & fips == 1001
	replace pct90_1 = stats_1[`j',9] if year == `t' & fips == 1001
	local j = `j' + 1
}
}
twoway rarea pct10_1 pct25_1 year if fips == 1001 & year > 1949, fcolor(gs10) lcolor(gs10) lstyle(none) || rarea pct25_1 pct50_1 year if fips == 1001 & year > 1949, fcolor(gs4) lcolor(gs10) lstyle(none) || rarea pct50_1 pct75_1 year if fips == 1001 & year > 1949, fcolor(gs4) lcolor(gs10) lstyle(none)  || rarea pct75_1 pct90_1 year if fips == 1001 & year > 1949, fcolor(gs10) lcolor(gs10) lstyle(none) ///
|| line pct50_1 year if fips == 1001 & year > 1949, graphregion(color(white))  ytitle("Sensitivity to KDD") xtitle("Year") color(1 gs4) color(2 gs0) lwidth(medthick) xlabel(1950 1960 1970 1980 1990 2000 2010) ylabel(-0.02 -0.015 -0.01 -0.005 0) name(Figure2, replace) legend(off)

//
// Table 3 Replication
//

qui {
matrix Table3_upper = J(4,4,0)

gen t = (year - 1949)/10
regress pct50_1 t, vce(hc3)
matrix Table3_upper[1,1] = _b[t]
matrix Table3_upper[2,1] = _b[_cons]
matrix Table3_upper[1,2] = _se[t]
matrix Table3_upper[2,2] = _se[_cons]

matrix rownames Table3_upper = t cons tbreak break
matrix colnames Table3_upper = beta stderr beta stderr

estat sbsingle
local breakdate = r(breakdate)
gen tempbreak = 0 if fips == 1001
replace tempbreak = 1 if fips == 1001 & year >= 1989
gen tempbreak_int = tempbreak * t if fips == 1001
reg pct50_1 t tempbreak tempbreak_int, vce(hc3)

matrix Table3_upper[1,3] = _b[t]
matrix Table3_upper[2,3] = _b[_cons]
matrix Table3_upper[1,4] = _se[t]
matrix Table3_upper[2,4] = _se[_cons]
matrix Table3_upper[3,3] = _b[tempbreak_int]
matrix Table3_upper[4,3] = _b[tempbreak]
matrix Table3_upper[3,4] = _se[tempbreak_int]
matrix Table3_upper[4,4] = _se[tempbreak]

noi display "Table 3 Upper Panel:"
noi matlist Table3_upper

noi display "Table 3 Lower Panel:"
regress pct50_1 t, vce(hc3)
noi estat sbsingle, all

noi display "Estimated Break Date:"
noi display "`breakdate'"

drop temp* t 
}

//
// Figure 3 Replication
// 

qui{
collapse (mean) bit_kdd, by(fips)

merge 1:1 fips using "`path'Data/US_County_LowRes_2013data_Stata11.dta"
drop _merge
}
spmap bit_kdd using "`path'Data/US_County_LowRes_2013coord_Stata11.dta" if continental == 1, id(_ID) fcolor(RdYlBu) legstyle(2) legjunction(" to ") clmethod(custom) clbreaks(-0.08 -0.015 -0.01 -0.005 0 0.005) legend(symy(*1) symx(*2) size(*2)) legorder(lohi) ndlabel("") legtitle("KDD Coefficient:") name(Figure3, replace)

// 
// Figure 4 Replication
//
qui {
use "`path'/Data/Replication_Temp.dta", clear
bysort fips: egen tempavg_kdd = mean(kdd) if sample_corn == 1
bysort fips: egen tempavg_bitkdd = mean(bit_kdd) if sample_corn == 1
gen templog_avgkdd = log(tempavg_kdd)
regress bit_kdd log_kdd, vce(cluster stateid)
predict temppred_bitkdd
predict temppred_bitkdd_stderr, stdp 
gen temppred_bitkdd_high = temppred_bitkdd + 1.96*temppred_bitkdd_stderr
gen temppred_bitkdd_low = temppred_bitkdd - 1.96*temppred_bitkdd_stderr

set seed 123456789
gen rand = runiform()
}
twoway rarea temppred_bitkdd_low temppred_bitkdd_high kdd if year == 1980 & kdd < 300 & temppred_bitkdd > -0.04  & tempavg_bitkdd < 0.002, sort color(gs14) || scatter bit_kdd kdd if rand > 0.97 & kdd < 300 & bit_kdd > -0.04 & bit_kdd < 0.002, graphregion(color(white))  ytitle("Sensitivity to KDD") xtitle("KDD") ylabel(-0.04 -0.03 -0.02 -0.01 0) xlabel(0[100]300) name(Figure4, replace) legend(off) msize(vsmall) mcolor(gs7) || line temppred_bitkdd kdd if year == 1980 & kdd < 300 & bit_kdd > -0.04  & bit_kdd < 0.002, sort lwidth(medthick) color(black)  
qui drop temp* rand

// 
// Figure 5 Replication
//

qui {
gen temp_kdd = 0 in 1/300
replace temp_kdd = temp_kdd[_n-1] + 1 in 2/300
gen marginal_effect_pooled = -.0158332 + .0021673*log(temp_kdd) 
gen marginal_effect_pooled_stdp = sqrt(5.60e-06+2.44e-07*log(temp_kdd)*log(temp_kdd)+2*log(temp_kdd)*(-1.15e-06))
gen marginal_effect_pooled_low = marginal_effect_pooled - 1.96*marginal_effect_pooled_stdp
gen marginal_effect_pooled_high = marginal_effect_pooled + 1.96*marginal_effect_pooled_stdp

gen marginal_effect_mocorn = -.0183599 + .002505*log(temp_kdd)
} 
twoway rarea marginal_effect_pooled_low marginal_effect_pooled_high temp_kdd if marginal_effect_pooled > -0.02, color(gs14) ///
|| line marginal_effect_pooled temp_kdd if marginal_effect_pooled > -0.03, color(black) ///
||  line marginal_effect_mocorn temp_kdd if marginal_effect_pooled > -0.03, lwidth(medthick) color(black) lpattern(dash)  graphregion(color(white))  ytitle("Sensitivity to KDD") xtitle("KDD") ylabel(-0.02 -0.01 0) xlabel(0[100]300) name(Figure5, replace) legend(label(1 "CI 95%") label(2 "FE-OLS") label(3 "MO-OLS") position(4) bmargin(none) ring(0))

qui drop marginal_effect* temp*

// 
// Save estimates for Projections
// 

collapse (mean) pooled_fixedeffect modols_fixedeffect corn_weight, by(fips)
save "`path'/Data/Replication_Temp2.dta", replace


qui {
import delimited "`path'Data/Data_Projections.csv", case(preserve) clear
xtset fips year
}

// 
// Calculate Projections for each climate model
// 

qui{
local rcp26list = "BNUESM_26_r1i1p1 CANESM_26_r1i1p1 CCSM4_26_r1i1p1 CESM1CAM5_26_r1i1p1 CSIROMK3_26_r1i1p1 ECEARTH_26_r1i1p1 FGOALS_G2_26_r1i1p1 FIOESM_26_r1i1p1 GFDL2G_26_r1i1p1 GFDLCM3_26_r1i1p1 GISSE2R_26_r1i1p1 HadGSM_26_r1i1p1 IPSL_CM5_26_r1i1p1 IPSL_MR_26_r1i1p1 MIROC_CH_26_r1i1p1 MIROC_26_r1i1p1 MPI_ESM_LR_26_r1i1p1 NorESM_26_r1i1p1" // 
local rcp45list = "BNUESM_45_r1i1p1 CANESM_45_r1i1p1 CCSM4_45_r1i1p1 CESM1CAM5_45_r1i1p1 CSIROMK3_45_r1i1p1 ECEARTH_45_r1i1p1 FGOALS_G2_45_r1i1p1 FIOESM_45_r1i1p1 GFDL2G_45_r1i1p1 GFDLCM3_45_r1i1p1 GISSE2R_45_r1i1p1 HadGSM_45_r1i1p1 IPSL_CM5_45_r1i1p1 IPSL_MR_45_r1i1p1 MIROC_CH_45_r1i1p1 MIROC_45_r1i1p1 MPI_ESM_LR_45_r1i1p1 NorESM_45_r1i1p1 ACCESS_45_r1i1p1"
local rcp85list = "BNUESM_85_r1i1p1 CANESM_85_r1i1p1 CCSM4_85_r1i1p1 CESM1CAM5_85_r1i1p1 CSIROMK3_85_r1i1p1 ECEARTH_85_r1i1p1 FGOALS_G2_85_r1i1p1 FIOESM_85_r1i1p1 GFDL2G_85_r1i1p1 GFDLCM3_85_r1i1p1 GISSE2R_85_r1i1p1 HadGSM_85_r1i1p1 IPSL_CM5_85_r1i1p1 IPSL_MR_85_r1i1p1 MIROC_CH_85_r1i1p1 MIROC_85_r1i1p1 MPI_ESM_LR_85_r1i1p1 NorESM_85_r1i1p1 ACCESS_85_r1i1p1"
local modellist = "`rcp85list' `rcp45list' `rcp26list'"  

merge 1:1 fips year using "`path'Data/Replication_Temp.dta", keepusing(bit_cons bit_gdd bit_kdd bit_prec bit_prec_sq)
drop _merge

merge m:1 fips using "`path'Data/Replication_Temp2.dta", keepusing(pooled_fixedeffect modols_fixedeffect corn_weight)
drop _merge
gen hist_cornyield = exp(log_cornyield)
sort fips year

// Loop over models to generate forecasts of crop yield
foreach m in `modellist' {
	// Generate index variables for gdd, kdd, and prec
	egen temp1 = mean(gdd_`m') if year > 2005 & year < 2016, by(fips)
	egen temp2 = mean(kdd_`m') if year > 2005 & year < 2016, by(fips)
	egen temp3 = mean(prec_`m') if year > 2005 & year < 2016, by(fips)
	egen temp12 = mean(temp1), by(fips)
	egen temp22 = mean(temp2), by(fips)
	egen temp32 = mean(temp3), by(fips)
	gen index_gdd = gdd_`m' - temp12 if year > 2005
	gen index_kdd = kdd_`m' - temp22 if year > 2005
	gen index_prec = prec_`m' - temp32 if year > 2005
	
	drop temp*
	
	// Apply the index to the historical series
	egen temp1 = mean(gdd) if year > 2005 & year < 2016, by(fips)
	egen temp2 = mean(kdd) if year > 2005 & year < 2016, by(fips)
	egen temp3 = mean(prec) if year > 2005 & year < 2016, by(fips)
	egen temp12 = mean(temp1), by(fips)
	egen temp22 = mean(temp2), by(fips)
	egen temp32 = mean(temp3), by(fips)
	gen extrap_gdd = index_gdd + temp12 if year > 2005
	gen extrap_kdd = index_kdd + temp22 if year > 2005
	gen extrap_prec = index_prec + temp32 if year > 2005
	
	drop temp*
		
	// Now generate new yieldhats
	
	// Standard FE-OLS (no adaptation)
	gen hat1_logcornyield = .0002375*extrap_gdd -.0052091*extrap_kdd + .0007648*extrap_prec -6.19e-07*extrap_prec*extrap_prec + 2.675514 + 1.3486145 + pooled_fixedeffect if year > 2005 // 1.338556 is the 2006-2015 average of the time dummy
	
	// Modified FE-OLS (adaptation)
	gen hat2_logcornyield = 0.0004266*extrap_gdd -.0159265*extrap_kdd + .0021851*(log(extrap_kdd)*extrap_kdd - extrap_kdd) + .0007202*extrap_prec -6.3e-07*extrap_prec*extrap_prec + 2.188991 + 1.3573626 + modols_fixedeffect if year > 2005 & extrap_kdd <=400 // 1.3573626 is the 2006-2015 average of the time dummy
	replace hat2_logcornyield = 0.0004266*extrap_gdd -.0159265*400 + .0021851*(log(400)*400-400) -0.0028*(extrap_kdd-400) + .0007202*extrap_prec -6.3e-07*extrap_prec*extrap_prec + 2.188991 + 1.3573626 + modols_fixedeffect if year > 2005 & extrap_kdd > 400 // Hack to ensure no positive marginal effects of KDD
		
	// MO-OLS (no adaptation)
	egen temp1 = mean(bit_gdd) if year > 2005 & year < 2016, by(fips)
	egen temp2 = mean(bit_kdd) if year > 2005 & year < 2016, by(fips)
	egen temp3 = mean(bit_prec) if year > 2005 & year < 2016, by(fips)
	egen temp4 = mean(bit_prec_sq) if year > 2005 & year < 2016, by(fips)
	egen temp5 = mean(bit_cons) if year > 2005 & year < 2016, by(fips)
	egen temp12 = mean(temp1), by(fips)
	egen temp22 = mean(temp2), by(fips)
	egen temp32 = mean(temp3), by(fips)
	egen temp42 = mean(temp4), by(fips)
	egen temp52 = mean(temp5), by(fips)
	
	gen hat3_logcornyield = temp12*extrap_gdd + temp22*extrap_kdd + temp32*extrap_prec + temp42*extrap_prec*extrap_prec + temp52 if year > 2005 

	// MO-OLS (adaptation)
	gen hat4_logcornyield = temp12*extrap_gdd + (log(extrap_kdd)*0.002505 -.0183599)*extrap_kdd + temp32*extrap_prec + temp42*extrap_prec*extrap_prec + temp52 if year > 2005 & extrap_kdd < 400
	replace hat4_logcornyield = temp12*extrap_gdd + (log(400)*0.002505 -.0183599)*extrap_kdd + temp32*extrap_prec + temp42*extrap_prec*extrap_prec + temp52 if year > 2005 & extrap_kdd >= 400

	drop temp*
	
	// Exponentate the hat variables
	gen hat1_cornyield = exp(hat1_logcornyield)
	gen hat2_cornyield = exp(hat2_logcornyield)
	gen hat3_cornyield = exp(hat3_logcornyield)
	gen hat4_cornyield = exp(hat4_logcornyield)
	
	// Calculate in terms of percentage change from base
	egen temp1 = mean(hist_cornyield) if year > 2005 & year < 2016, by(fips)
	egen temp12 = mean(temp1), by(fips)

	gen predfup_1_`m' = (hat1_cornyield - temp12)/temp12
	gen predfup_2_`m'  = (hat2_cornyield - temp12)/temp12
	gen predfup_3_`m'  = (hat3_cornyield - temp12)/temp12
	gen predfup_4_`m'  = (hat4_cornyield - temp12)/temp12
	
	drop temp*
	
	// Create average forecast series
	gen predm1_`m' = .
	gen predm1p_`m' = .
	gen predm2_`m' = .
	gen predm2p_`m' = .
	gen predm3_`m' = .
	gen predm3p_`m' = .
	gen predm4_`m' = .
	gen predm4p_`m' = .

	forvalues y = 2006/2100 {
		su hat1_cornyield if year == `y' [aweight=corn_weight], meanonly
		replace predm1_`m' = r(mean) if year == `y' & fips == 1001
		su predfup_1_`m' if year == `y' [aweight=corn_weight], meanonly
		replace predm1p_`m' = r(mean) if year == `y' & fips == 1001

		su hat2_cornyield if year == `y' [aweight=corn_weight], meanonly
		replace predm2_`m' = r(mean) if year == `y' & fips == 1001
		su predfup_2_`m' if year == `y' [aweight=corn_weight], meanonly
		replace predm2p_`m' = r(mean) if year == `y' & fips == 1001

		su hat3_cornyield if year == `y' [aweight=corn_weight], meanonly
		replace predm3_`m' = r(mean) if year == `y' & fips == 1001
		su predfup_3_`m' if year == `y' [aweight=corn_weight], meanonly
		replace predm3p_`m' = r(mean) if year == `y' & fips == 1001

		su hat4_cornyield if year == `y' [aweight=corn_weight], meanonly
		replace predm4_`m' = r(mean) if year == `y' & fips == 1001
		su predfup_4_`m' if year == `y' [aweight=corn_weight], meanonly
		replace predm4p_`m' = r(mean) if year == `y' & fips == 1001
	}
			
	drop extrap*
	drop index*
	drop hat*
}

// Average the results by rcp then graph average lines
local rcp26list = "BNUESM_26_r1i1p1 CANESM_26_r1i1p1 CCSM4_26_r1i1p1 CESM1CAM5_26_r1i1p1 CSIROMK3_26_r1i1p1 ECEARTH_26_r1i1p1 FGOALS_G2_26_r1i1p1 FIOESM_26_r1i1p1 GFDL2G_26_r1i1p1 GFDLCM3_26_r1i1p1 GISSE2R_26_r1i1p1 HadGSM_26_r1i1p1 IPSL_CM5_26_r1i1p1 IPSL_MR_26_r1i1p1 MIROC_CH_26_r1i1p1 MIROC_26_r1i1p1 MPI_ESM_LR_26_r1i1p1 NorESM_26_r1i1p1" // 
local rcp45list = "BNUESM_45_r1i1p1 CANESM_45_r1i1p1 CCSM4_45_r1i1p1 CESM1CAM5_45_r1i1p1 CSIROMK3_45_r1i1p1 ECEARTH_45_r1i1p1 FGOALS_G2_45_r1i1p1 FIOESM_45_r1i1p1 GFDL2G_45_r1i1p1 GFDLCM3_45_r1i1p1 GISSE2R_45_r1i1p1 HadGSM_45_r1i1p1 IPSL_CM5_45_r1i1p1 IPSL_MR_45_r1i1p1 MIROC_CH_45_r1i1p1 MIROC_45_r1i1p1 MPI_ESM_LR_45_r1i1p1 NorESM_45_r1i1p1 ACCESS_45_r1i1p1"
local rcp85list = "BNUESM_85_r1i1p1 CANESM_85_r1i1p1 CCSM4_85_r1i1p1 CESM1CAM5_85_r1i1p1 CSIROMK3_85_r1i1p1 ECEARTH_85_r1i1p1 FGOALS_G2_85_r1i1p1 FIOESM_85_r1i1p1 GFDL2G_85_r1i1p1 GFDLCM3_85_r1i1p1 GISSE2R_85_r1i1p1 HadGSM_85_r1i1p1 IPSL_CM5_85_r1i1p1 IPSL_MR_85_r1i1p1 MIROC_CH_85_r1i1p1 MIROC_85_r1i1p1 MPI_ESM_LR_85_r1i1p1 NorESM_85_r1i1p1 ACCESS_85_r1i1p1"
local rcplist = "85 45 26"
qui foreach r in `rcplist' {
	local i = 0
	gen pred1_rcp`r'avg = 0
	gen pred2_rcp`r'avg = 0
	gen pred3_rcp`r'avg = 0
	gen pred4_rcp`r'avg = 0
	gen pred1p_rcp`r'avg = 0
	gen pred2p_rcp`r'avg = 0
	gen pred3p_rcp`r'avg = 0
	gen pred4p_rcp`r'avg = 0
	gen pred1_rcp`r'tsum = 0
	gen pred2_rcp`r'tsum = 0
	gen pred3_rcp`r'tsum = 0
	gen pred4_rcp`r'tsum = 0
	gen pred1p_rcp`r'tsum = 0
	gen pred2p_rcp`r'tsum = 0
	gen pred3p_rcp`r'tsum = 0
	gen pred4p_rcp`r'tsum = 0
	
	foreach m in `rcp`r'list' {
		replace pred1_rcp`r'avg = pred1_rcp`r'avg + predm1_`m'
		replace pred2_rcp`r'avg = pred2_rcp`r'avg + predm2_`m'
		replace pred3_rcp`r'avg = pred3_rcp`r'avg + predm3_`m'
		replace pred4_rcp`r'avg = pred4_rcp`r'avg + predm4_`m'

		replace pred1p_rcp`r'avg = pred1p_rcp`r'avg + predm1p_`m'
		replace pred2p_rcp`r'avg = pred2p_rcp`r'avg + predm2p_`m'
		replace pred3p_rcp`r'avg = pred3p_rcp`r'avg + predm3p_`m'
		replace pred4p_rcp`r'avg = pred4p_rcp`r'avg + predm4p_`m'
		
		local i = `i' + 1
	}
	replace pred1_rcp`r'avg = pred1_rcp`r'avg / `i'
	replace pred2_rcp`r'avg = pred2_rcp`r'avg / `i'
	replace pred3_rcp`r'avg = pred3_rcp`r'avg / `i'
	replace pred4_rcp`r'avg = pred4_rcp`r'avg / `i'
	
	replace pred1p_rcp`r'avg = pred1p_rcp`r'avg / `i'
	replace pred2p_rcp`r'avg = pred2p_rcp`r'avg / `i'
	replace pred3p_rcp`r'avg = pred3p_rcp`r'avg / `i'
	replace pred4p_rcp`r'avg = pred4p_rcp`r'avg / `i'

	// Generate series that say what proportion of counties experience losses
	gen pctloss_1_rcp`r'avg = 0
	gen pctloss_2_rcp`r'avg = 0
	gen pctloss_3_rcp`r'avg = 0
	gen pctloss_4_rcp`r'avg = 0

	local i = 0
	
	foreach m in `rcp`r'list' {
		gen pcttemp1 = 0 
		replace pcttemp1 = 1 if predfup_1_`m' < 0
		gen pcttemp2= 0 
		replace pcttemp2 = 1 if predfup_2_`m' < 0
		gen pcttemp3 = 0 
		replace pcttemp3 = 1 if predfup_3_`m' < 0
		gen pcttemp4 = 0 
		replace pcttemp4 = 1 if predfup_4_`m' < 0
		
		gen pcttemp12 = 0
		gen pcttemp22 = 0
		gen pcttemp32 = 0
		gen pcttemp42 = 0
		
		forvalues y = 2006/2100 {
			su pcttemp1 if year == `y' & corn_weight > 0 & !missing(predfup_1_`m'), meanonly
			replace pcttemp12 = r(mean) if year == `y' & fips == 1001
			su pcttemp2 if year == `y' & corn_weight > 0 & !missing(predfup_2_`m'), meanonly
			replace pcttemp22 = r(mean) if year == `y' & fips == 1001
			su pcttemp3 if year == `y' & corn_weight > 0 & !missing(predfup_3_`m'), meanonly
			replace pcttemp32 = r(mean) if year == `y' & fips == 1001
			su pcttemp4 if year == `y' & corn_weight > 0 & !missing(predfup_4_`m'), meanonly
			replace pcttemp42 = r(mean) if year == `y' & fips == 1001
		}
		
		replace pctloss_1_rcp`r'avg = pctloss_1_rcp`r'avg + pcttemp12
		replace pctloss_2_rcp`r'avg = pctloss_2_rcp`r'avg + pcttemp22
		replace pctloss_3_rcp`r'avg = pctloss_3_rcp`r'avg + pcttemp32
		replace pctloss_4_rcp`r'avg = pctloss_4_rcp`r'avg + pcttemp42
		drop pcttemp*
		local i = `i' + 1
	}
	replace pctloss_1_rcp`r'avg = pctloss_1_rcp`r'avg / `i'
	replace pctloss_2_rcp`r'avg = pctloss_2_rcp`r'avg / `i'
	replace pctloss_3_rcp`r'avg = pctloss_3_rcp`r'avg / `i'
	replace pctloss_4_rcp`r'avg = pctloss_4_rcp`r'avg / `i'

	// Calculate standard errors 
	local i = 0
	foreach m in `rcp`r'list' {
		gen temp1`m' = .
		gen temp2`m' = .
		gen temp3`m' = .
		gen temp4`m' = .	
		replace temp1`m' = ((predm1_`m' + L1.predm1_`m' + L2.predm1_`m' + L3.predm1_`m' + L4.predm1_`m')/5 - pred1_rcp`r'avg)^2 if fips == 1001 
		replace temp2`m' = ((predm2_`m' + L1.predm2_`m' + L2.predm2_`m' + L3.predm2_`m' + L4.predm2_`m' )/5 - pred2_rcp`r'avg)^2 if fips == 1001 
		replace temp3`m' = ((predm3_`m' + L1.predm3_`m' + L2.predm3_`m' + L3.predm3_`m' + L4.predm3_`m')/5 - pred3_rcp`r'avg)^2 if fips == 1001 
		replace temp4`m' = ((predm4_`m' + L1.predm4_`m' + L2.predm4_`m' + L3.predm4_`m' + L4.predm4_`m')/5 - pred4_rcp`r'avg)^2 if fips == 1001
		
		replace pred1_rcp`r'tsum = pred1_rcp`r'tsum + temp1`m' 
		replace pred2_rcp`r'tsum = pred2_rcp`r'tsum + temp2`m' 
		replace pred3_rcp`r'tsum = pred3_rcp`r'tsum + temp3`m' 
		replace pred4_rcp`r'tsum = pred4_rcp`r'tsum + temp4`m' 

		gen temp1p`m' = .
		gen temp2p`m' = .
		gen temp3p`m' = .
		gen temp4p`m' = .	
		replace temp1p`m' = ((predm1p_`m' + L1.predm1p_`m' + L2.predm1p_`m' + L3.predm1p_`m' + L4.predm1p_`m')/5 - pred1p_rcp`r'avg)^2 if fips == 1001 
		replace temp2p`m' = ((predm2p_`m' + L1.predm2p_`m' + L2.predm2p_`m' + L3.predm2p_`m' + L4.predm2p_`m' )/5 - pred2p_rcp`r'avg)^2 if fips == 1001 
		replace temp3p`m' = ((predm3p_`m' + L1.predm3p_`m' + L2.predm3p_`m' + L3.predm3p_`m' + L4.predm3p_`m')/5 - pred3p_rcp`r'avg)^2 if fips == 1001 
		replace temp4p`m' = ((predm4p_`m' + L1.predm4p_`m' + L2.predm4p_`m' + L3.predm4p_`m' + L4.predm4p_`m')/5 - pred4p_rcp`r'avg)^2 if fips == 1001
		
		replace pred1p_rcp`r'tsum = pred1p_rcp`r'tsum + temp1p`m' 
		replace pred2p_rcp`r'tsum = pred2p_rcp`r'tsum + temp2p`m' 
		replace pred3p_rcp`r'tsum = pred3p_rcp`r'tsum + temp3p`m' 
		replace pred4p_rcp`r'tsum = pred4p_rcp`r'tsum + temp4p`m' 
							
		local i = `i' + 1
	}
	drop temp*
	
	gen pred1_rcp`r'std = sqrt(pred1_rcp`r'tsum/(`i'-1))
	gen pred2_rcp`r'std = sqrt(pred2_rcp`r'tsum/(`i'-1))
	gen pred3_rcp`r'std = sqrt(pred3_rcp`r'tsum/(`i'-1))
	gen pred4_rcp`r'std = sqrt(pred4_rcp`r'tsum/(`i'-1))

	gen pred1p_rcp`r'std = sqrt(pred1p_rcp`r'tsum/(`i'-1))
	gen pred2p_rcp`r'std = sqrt(pred2p_rcp`r'tsum/(`i'-1))
	gen pred3p_rcp`r'std = sqrt(pred3p_rcp`r'tsum/(`i'-1))
	gen pred4p_rcp`r'std = sqrt(pred4p_rcp`r'tsum/(`i'-1))

	sort fips year
	tsset fips year
	
	// Create moving averages of the average corn yield and also
	gen pred1_ma_rcp`r'avg = (pred1_rcp`r'avg + L1.pred1_rcp`r'avg + L2.pred1_rcp`r'avg + L3.pred1_rcp`r'avg + L4.pred1_rcp`r'avg)/5 if fips == 1001 & year > 2014
	gen pred2_ma_rcp`r'avg = (pred2_rcp`r'avg + L1.pred2_rcp`r'avg + L2.pred2_rcp`r'avg + L3.pred2_rcp`r'avg + L4.pred2_rcp`r'avg)/5 if fips == 1001 & year > 2014
	gen pred3_ma_rcp`r'avg = (pred3_rcp`r'avg + L1.pred3_rcp`r'avg + L2.pred3_rcp`r'avg + L3.pred3_rcp`r'avg + L4.pred3_rcp`r'avg)/5 if fips == 1001 & year > 2014
	gen pred4_ma_rcp`r'avg = (pred4_rcp`r'avg + L1.pred4_rcp`r'avg + L2.pred4_rcp`r'avg + L3.pred4_rcp`r'avg + L4.pred4_rcp`r'avg)/5 if fips == 1001 & year > 2014
	
	gen pred1p_ma_rcp`r'avg = (pred1p_rcp`r'avg + L1.pred1p_rcp`r'avg + L2.pred1p_rcp`r'avg + L3.pred1p_rcp`r'avg + L4.pred1p_rcp`r'avg)/5 if fips == 1001 & year > 2014
	gen pred2p_ma_rcp`r'avg = (pred2p_rcp`r'avg + L1.pred2p_rcp`r'avg + L2.pred2p_rcp`r'avg + L3.pred2p_rcp`r'avg + L4.pred2p_rcp`r'avg)/5 if fips == 1001 & year > 2014
	gen pred3p_ma_rcp`r'avg = (pred3p_rcp`r'avg + L1.pred3p_rcp`r'avg + L2.pred3p_rcp`r'avg + L3.pred3p_rcp`r'avg + L4.pred3p_rcp`r'avg)/5 if fips == 1001 & year > 2014
	gen pred4p_ma_rcp`r'avg = (pred4p_rcp`r'avg + L1.pred4p_rcp`r'avg + L2.pred4p_rcp`r'avg + L3.pred4p_rcp`r'avg + L4.pred4p_rcp`r'avg)/5 if fips == 1001 & year > 2014

	gen pctloss1_ma_rcp`r'avg = (pctloss_1_rcp`r'avg + L1.pctloss_1_rcp`r'avg + L2.pctloss_1_rcp`r'avg + L3.pctloss_1_rcp`r'avg + L4.pctloss_1_rcp`r'avg)/5 if fips == 1001 & year > 2014
	gen pctloss2_ma_rcp`r'avg = (pctloss_2_rcp`r'avg + L1.pctloss_2_rcp`r'avg + L2.pctloss_2_rcp`r'avg + L3.pctloss_2_rcp`r'avg + L4.pctloss_2_rcp`r'avg)/5 if fips == 1001 & year > 2014
	gen pctloss3_ma_rcp`r'avg = (pctloss_3_rcp`r'avg + L1.pctloss_3_rcp`r'avg + L2.pctloss_3_rcp`r'avg + L3.pctloss_3_rcp`r'avg + L4.pctloss_3_rcp`r'avg)/5 if fips == 1001 & year > 2014
	gen pctloss4_ma_rcp`r'avg = (pctloss_4_rcp`r'avg + L1.pctloss_4_rcp`r'avg + L2.pctloss_4_rcp`r'avg + L3.pctloss_4_rcp`r'avg + L4.pctloss_4_rcp`r'avg)/5 if fips == 1001 & year > 2014

	// Create prediction intervals for the moving averages
	gen pred1_rcp`r'high = pred1_ma_rcp`r'avg + 1.28*pred1_rcp`r'std
	gen pred1_rcp`r'low = pred1_ma_rcp`r'avg - 1.28*pred1_rcp`r'std
	gen pred2_rcp`r'high = pred2_ma_rcp`r'avg + 1.28*pred2_rcp`r'std
	gen pred2_rcp`r'low = pred2_ma_rcp`r'avg - 1.28*pred2_rcp`r'std
	gen pred3_rcp`r'high = pred3_ma_rcp`r'avg + 1.28*pred3_rcp`r'std
	gen pred3_rcp`r'low = pred3_ma_rcp`r'avg - 1.28*pred3_rcp`r'std
	gen pred4_rcp`r'high = pred4_ma_rcp`r'avg + 1.28*pred4_rcp`r'std
	gen pred4_rcp`r'low = pred4_ma_rcp`r'avg - 1.28*pred4_rcp`r'std

	gen pred1p_rcp`r'high = pred1p_ma_rcp`r'avg + 1.28*pred1p_rcp`r'std
	gen pred1p_rcp`r'low = pred1p_ma_rcp`r'avg - 1.28*pred1p_rcp`r'std
	gen pred2p_rcp`r'high = pred2p_ma_rcp`r'avg + 1.28*pred2p_rcp`r'std
	gen pred2p_rcp`r'low = pred2p_ma_rcp`r'avg - 1.28*pred2p_rcp`r'std
	gen pred3p_rcp`r'high = pred3p_ma_rcp`r'avg + 1.28*pred3p_rcp`r'std
	gen pred3p_rcp`r'low = pred3p_ma_rcp`r'avg - 1.28*pred3p_rcp`r'std
	gen pred4p_rcp`r'high = pred4p_ma_rcp`r'avg + 1.28*pred4p_rcp`r'std
	gen pred4p_rcp`r'low = pred4p_ma_rcp`r'avg - 1.28*pred4p_rcp`r'std
	
	// Average out the full county by county results by model
	gen predfup_1_rcp`r'avg = 0
	gen predfup_2_rcp`r'avg = 0
	gen predfup_3_rcp`r'avg = 0
	gen predfup_4_rcp`r'avg = 0	
	local i = 0
	foreach m in `rcp`r'list' {
		replace predfup_1_rcp`r'avg = predfup_1_rcp`r'avg + predfup_1_`m'
		replace predfup_2_rcp`r'avg = predfup_2_rcp`r'avg + predfup_2_`m'
		replace predfup_3_rcp`r'avg = predfup_3_rcp`r'avg + predfup_3_`m'
		replace predfup_4_rcp`r'avg = predfup_4_rcp`r'avg + predfup_4_`m'

		local i = `i' + 1
	}
	replace predfup_1_rcp`r'avg = predfup_1_rcp`r'avg / `i'
	replace predfup_2_rcp`r'avg = predfup_2_rcp`r'avg / `i'
	replace predfup_3_rcp`r'avg = predfup_3_rcp`r'avg / `i'
	replace predfup_4_rcp`r'avg = predfup_4_rcp`r'avg / `i'
} 
}

// 
// Table 4 Replication
// 

qui {
matrix Table4 = J(12,12,0)
matrix rownames Table4 = RCP85_30 RCP85_50 RCP85_80 RCP85_00 RCP45_30 RCP45_50 RCP45_80 RCP45_00 RCP26_30 RCP26_50 RCP26_80 RCP26_00
matrix colnames Table4 = FENA_Mean FENA_low FENA_high MONA_Mean MONA_low MONA_high FEA_Mean FEA_low FEA_high MOA_Mean MOA_low MOA_high
local i = 0
local rcplist = "85 45 26"
qui foreach r in `rcplist' {

	local yearlist = "2030 2050 2080 2100"
	local b = 1
	foreach y in `yearlist' {
		su pred1p_ma_rcp`r'avg if fips == 1001 & year == `y', meanonly
		matrix Table4[`i'*4 + `b',1] = r(mean)
		su pred1p_rcp`r'low if fips == 1001 & year == `y', meanonly
		matrix Table4[`i'*4 + `b',2] = r(mean)
		su pred1p_rcp`r'high if fips == 1001 & year == `y', meanonly
		matrix Table4[`i'*4 + `b',3] = r(mean)
		
		su pred3p_ma_rcp`r'avg if fips == 1001 & year == `y', meanonly
		matrix Table4[`i'*4 + `b',4] = r(mean)
		su pred3p_rcp`r'low if fips == 1001 & year == `y', meanonly
		matrix Table4[`i'*4 + `b',5] = r(mean)
		su pred3p_rcp`r'high if fips == 1001 & year == `y', meanonly
		matrix Table4[`i'*4 + `b',6] = r(mean)

		su pred2p_ma_rcp`r'avg if fips == 1001 & year == `y', meanonly
		matrix Table4[`i'*4 + `b',7] = r(mean)
		su pred2p_rcp`r'low if fips == 1001 & year == `y', meanonly
		matrix Table4[`i'*4 + `b',8] = r(mean)
		su pred2p_rcp`r'high if fips == 1001 & year == `y', meanonly
		matrix Table4[`i'*4 + `b',9] = r(mean)

		su pred4p_ma_rcp`r'avg if fips == 1001 & year == `y', meanonly
		matrix Table4[`i'*4 + `b',10] = r(mean)
		su pred4p_rcp`r'low if fips == 1001 & year == `y', meanonly
		matrix Table4[`i'*4 + `b',11] = r(mean)
		su pred4p_rcp`r'high if fips == 1001 & year == `y', meanonly
		matrix Table4[`i'*4 + `b',12] = r(mean)

		local b = `b' + 1
	}
	
	local i = `i' + 1
}
noi display "Table 4:"
noi matlist Table4
}

//
// Figure 6 Replication
// 

twoway rarea pred4p_rcp45high pred4p_rcp45low year if fips == 1001 & year > 2014, fcolor("117 112 179") lcolor("117 112 179")  lwidth(none) fintensity(50) ///
|| rarea pred4p_rcp85high pred4p_rcp85low year if fips == 1001 & year > 2014, fcolor("217 95 2") lcolor("217 95 2")  lwidth(none) fintensity(50)  ///
|| rarea pred3p_rcp85high pred3p_rcp85low year if fips == 1001 & year > 2014, fcolor("27 158 119") lcolor("27 158 119")  lwidth(none) fintensity(50)  ///
|| line pred4p_ma_rcp45avg pred4p_ma_rcp85avg pred3p_ma_rcp85avg year if fips == 1001 & year > 2014, lcolor("27 158 119" "217 95 2" "117 112 179") lwidth(medthick medthick medthick)  ///
name(Figure6, replace) graphregion(color(white))  ytitle("Percentage Change in Corn Yield") xtitle("Year") ///
legend(order(4 5 6) label(4 "RCP45 With Adaptation") label(5 "RCP85 With Adaptation") label(6 "RCP45 Without Adaptation") bmargin(none) pos(8) ring(0)) ylabel(0 -0.25 -0.5 -0.75) xlabel(2015[10]2095)

//
// Figure 7 Replication
// 

twoway line pctloss4_ma_rcp26avg pctloss4_ma_rcp45avg pctloss3_ma_rcp85avg pctloss4_ma_rcp85avg year if fips == 1001 & year > 2014,   ///
name(Figure7, replace) graphregion(color(white)) ytitle("Percentage of counties experiencing losses") xtitle("Year") lwidth(medthick medthick medthick medthick medthick medthick) ///
legend(label(1 "RCP26 - Adapt.") label(2 "RCP45 - Adapt.") label(3 "RCP85 - No Adapt.") label(4 "RCP85 - Adapt.") bmargin(none) pos(11) ring(0)) ylabel(0.5 0.6 0.7 0.8 0.9 1) xlabel(2015[10]2095) ///
lcolor(gs10 black gs10 black) lpattern(dash dash solid solid)

//
// Table 5 Replication
// 
qui {
local modelnamelist = "BNUESM CANESM CCSM4 CESM1CAM5 CSIROMK3 ECEARTH FGOALS_G2 FIOESM GFDL2G GFDLCM3 GISSE2R HadGSM IPSL_CM5 IPSL_MR MIROC_CH MIROC MPI_ESM_LR NorESM"

gen temp_mo_ad = .
gen temp_mo_er = .
gen temp_mo_both = .
gen temp_mo_ser = .
gen temp_fe_ad = .
gen temp_fe_both = .

local i = 1
qui foreach m in `modelnamelist' {
	
	// Baseline == that specific models average prediction 2006-2017
	su predm3_`m'_85_r1i1p1 if year > 2005 & year < 2016 & fips == 1001, meanonly
	local ref = r(mean)
	su predm1_`m'_85_r1i1p1 if year > 2005 & year < 2016 & fips == 1001, meanonly
	local ref2 = r(mean)
	gen tempz1 = `ref' - predm3_`m'_85_r1i1p1 if fips == 1001 & year >= 2020
	su tempz1 if fips == 1001
	local sum_mo_loss = r(sum)
	gen tempz2 =  `ref'  - predm4_`m'_85_r1i1p1 if fips == 1001 & year >= 2020
	su tempz2 if fips == 1001
	local sum_mo_adloss = r(sum)
	gen tempz3 =  `ref'  - predm3_`m'_45_r1i1p1 if fips == 1001 & year >= 2020
	su tempz3 if fips == 1001
	local sum_mo_erloss = r(sum)
	gen tempz4 =  `ref'  - predm4_`m'_45_r1i1p1 if fips == 1001 & year >= 2020
	su tempz4 if fips == 1001
	local sum_mo_bothloss = r(sum)
	gen tempz5 =  `ref'  - predm3_`m'_26_r1i1p1 if fips == 1001 & year >= 2020
	su tempz5 if fips == 1001
	local sum_mo_serloss = r(sum)
	
	replace temp_mo_ad = 1- `sum_mo_adloss' / `sum_mo_loss' in `i'
	replace temp_mo_er = 1- `sum_mo_erloss' / `sum_mo_loss' in `i'
	replace temp_mo_both = 1- `sum_mo_bothloss' / `sum_mo_loss' in `i'
	replace temp_mo_ser = 1- `sum_mo_serloss' / `sum_mo_loss' in `i'
	
	gen tempz7 =  `ref2'  - predm1_`m'_85_r1i1p1 if fips == 1001 & year >= 2020
	su tempz7 if fips == 1001
	local sum_fe_loss = r(sum)
	gen tempz8 =  `ref'  - predm2_`m'_85_r1i1p1 if fips == 1001 & year >= 2020
	su tempz8 if fips == 1001
	local sum_fe_adloss = r(sum)
	gen tempz10 =  `ref'  - predm2_`m'_45_r1i1p1 if fips == 1001 & year >= 2020
	su tempz10 if fips == 1001
	local sum_fe_bothloss = r(sum)

	replace temp_fe_ad = 1- `sum_fe_adloss' / `sum_mo_loss' in `i'
	replace temp_fe_both = 1- `sum_fe_bothloss' / `sum_mo_loss' in `i'
	
	local i = `i' + 1
	
	drop tempz*
	
}

matrix table_loss = J(2,12,0)
su temp_mo_ad in 1/`i'
matrix table_loss[1,1] = r(mean)
matrix table_loss[1,2] = r(mean) - 1.28*r(sd)
matrix table_loss[1,3] = r(mean) + 1.28*r(sd)
su temp_fe_ad in 1/`i'
matrix table_loss[2,1] = r(mean)
matrix table_loss[2,2] = r(mean) - 1.28*r(sd)
matrix table_loss[2,3] = r(mean) + 1.28*r(sd)
su temp_mo_er in 1/`i'
matrix table_loss[1,4] = r(mean)
matrix table_loss[1,5] = r(mean) - 1.28*r(sd)
matrix table_loss[1,6] = r(mean) + 1.28*r(sd)
su temp_mo_both in 1/`i'
matrix table_loss[1,7] = r(mean)
matrix table_loss[1,8] = r(mean) - 1.28*r(sd)
matrix table_loss[1,9] = r(mean) + 1.28*r(sd)
su temp_fe_both in 1/`i'
matrix table_loss[2,7] = r(mean)
matrix table_loss[2,8] = r(mean) - 1.28*r(sd)
matrix table_loss[2,9] = r(mean) + 1.28*r(sd)
su temp_mo_ser in 1/`i'
matrix table_loss[1,10] = r(mean)
matrix table_loss[1,11] = r(mean) - 1.28*r(sd)
matrix table_loss[1,12] = r(mean) + 1.28*r(sd)

noi display "Table 5 Replication:" 
matrix rownames table_loss = MOOLS FEOLS
matrix colnames table_loss = RCP85_Ad_Mean RCP85_Ad_Low RCP85_Ad_High RCP45_Mean RCP45_Low RCP45_High RCP45_Ad_Mean RCP45_Ad_Low RCP45_Ad_High RCP26_Mean RCP26_Low RCP26_High 
noi matlist table_loss
drop temp*
}

