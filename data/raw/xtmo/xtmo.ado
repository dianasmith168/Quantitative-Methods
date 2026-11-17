*! Timothy Neal -- 11/01/17
*! This is the first public version of xtmo, used to conduct Mean Observation OLS ('MO-OLS') on large panel data models. 
*! If there are any questions, issues, or comparatibility problems with this procedure, please email tjrneal@gmail.com. 
*!
program define xtmo, eclass prop(xt)
	version 11
	syntax varlist [if] [in] [, MAXITER(integer 1500) TOL(real 0.0001)]

	*! Mark the sample that is usable, identify the panel and time variable, and other panel statistics.
	marksample touse

	qui {
		xtset
		local ivar `r(panelvar)'
		local tvar `r(timevar)'
		levelsof `ivar' if `touse', local(ids)
		levelsof `tvar' if `touse', local(tds)
		sort `ivar' `tvar'
		global tts "`tds'"
		global iis "`ids'"
		local m=wordcount("`varlist'")
		local is = wordcount("`ids'")
		local ts = wordcount("`tds'")
		mata: macro_maxiter = strtoreal(st_local("maxiter"))
		mata: ids = tokens(st_local("ids"))
		mata: tds = tokens(st_local("tds"))

	*! Set up the variables
	tempvar constant
	gen `constant' = 1 if `touse'
	local varlist_cons = "`varlist' `constant'"

	*! Set up empty coefficient variables
	qui forvalues k = 1/`m' {
		tempvar bi`k' bt`k' bpool`k' bit`k'
		gen `bi`k'' = .
		gen `bt`k'' = .
		gen `bpool`k'' = .
		gen `bit`k'' = .
		
		local bilist = "`bilist' `bi`k''"
		local btlist = "`btlist' `bt`k''"
		local bitlist = "`bitlist' `bit`k''"
	}
	
	local num = 1
	foreach x in `varlist' {
		if (`num' == 1) local Yname "`x'"
		
		if `num' > 1 {
			local varnames "`varnames' `x'"
			gen bit_`x' = .
			local bitbclist = "`bitbclist' bit_`x'"
		}
		local num = `num' + 1
	}
	local varnames "`varnames' constant"
	gen bit_cons = .
	local bitbclist = "`bitbclist' bit_cons"
	
	*! Run regressions 1 - pooled
	regress `varlist_cons' if `touse', nocons
	tempvar samplevar
	gen `samplevar' = 0
	replace `samplevar' = 1 if e(sample)
	local obies = e(N)
	matrix beta = e(b)
	forvalues k = 1/`m' {
		replace `bpool`k'' = beta[1,`k'] if `touse'
	}
	*! Run regressions 2 - for each i
	local first = 1
	qui foreach i of global iis {		
		regress `varlist_cons' if `touse' & `ivar' == `i', nocons
		matrix beta = e(b)
		if (`first' == 1) mata: Tlist = st_numscalar("e(N)")
		else mata: Tlist = Tlist \ st_numscalar("e(N)")	
		local first = 0
		forvalues k = 1/`m' {
			replace `bi`k'' = beta[1,`k'] if `touse' & `ivar' == `i'
		}
	}
	*! Run regressions 3 - for each t
	local first = 1
	qui foreach t of global tts {		
		regress `varlist_cons' if `touse' & `tvar' == `t', nocons
		matrix beta = e(b)
		if (`first' == 1) mata: Nlist = st_numscalar("e(N)")
		else mata: Nlist = Nlist \ st_numscalar("e(N)")
		local first = 0
		forvalues k = 1/`m' {
			replace `bt`k'' = beta[1,`k'] if `touse' & `tvar' == `t'
		}
	}
	
	*! Create new it coefficient
	qui forvalues k = 1/`m' {
		replace `bit`k'' = `bi`k'' + `bt`k'' - `bpool`k'' if `touse'
	}
	
	*! Pass these coefficient variables to mata
	sort `ivar' `tvar'
	mata: dataf = st_data(., ("`varlist_cons'"),("`samplevar'"))
	mata: bi = st_data(., ("`bilist'"),("`samplevar'"))
	mata: bt = st_data(., ("`btlist'"),("`samplevar'"))
	mata: bit = st_data(., ("`bitlist'"),("`samplevar'")) 
	mata: ivar = st_data(., ("`ivar'"),("`samplevar'"))
	mata: tvar = st_data(., ("`tvar'"),("`samplevar'"))
	sort `tvar' `ivar'
	mata: dataf_t = st_data(., ("`varlist_cons'"),("`samplevar'"))
	mata: bi_t = st_data(., ("`bilist'"),("`samplevar'"))
	mata: ivar_t = st_data(., ("`ivar'"),("`samplevar'"))
	sort `ivar' `tvar'
	
	
	*! Run bias correction through the mata function
	noi di  in gr "Please Wait: Running Bias Correction"

	mata: bc(`m', `tol',"`samplevar'")
	
	*! Obtain MO-OLS beta vector
	mata: varcalc(`m',"`samplevar'")
	matrix colnames b = `varnames'
	
	*! Obtain asymptotic standard errors of each variable
	local num = 1
	matrix V = J(`m',`m',0)
	
	qui foreach x in `bitbclist' {
		tempvar temp1 temp2 var1 var2
		bysort `ivar': egen `temp1' = mean(`x') if !missing(`x')
		bysort `tvar': egen `temp2' = mean(`x') if !missing(`x')
		gen `var1' = (`x'-`temp1') * (`x' - `temp1')
		gen `var2' = (`x'-`temp2') * (`x' - `temp2')
		su `var1', meanonly
		local comp1 = r(mean)
		su `var2', meanonly
		local comp2 = r(mean)
		matrix V[`num',`num'] = `comp1'/`ts' + `comp2'/`is'		
		
		local num = `num' + 1
	}
	sort `ivar' `tvar'
	matrix rownames V = `varnames'
	matrix colnames V = `varnames'
	}
	
	*! ereturn storing
	ereturn post b V, depname("`Yname'") esample(`samplevar') obs(`obies')
	ereturn local tvar "`tvar'"
	ereturn local ivar "`ivar' `tvar'"
	ereturn local cmd "xtmo"
	ereturn local title "Mean Observation OLS (MO-OLS)"
	mata: poste()
	
	capture test `varnames', min constant
	if (_rc == 0) {
		ereturn scalar chi2 = r(chi2)
		ereturn scalar df_m = r(df)
	}
	else est scalar df_m = 0
	ereturn local chi2type "Wald"

	ereturn scalar iters = iterations
	local iters2 = iterations
	
	// Display Results
	display "`title'" 
		
	// Preliminary stats
	_crcphdr
	
	// Display the main regression results	
	ereturn display
	
	// Postscripts
	display in gr "Iterations used in the Bias Correction: `iters2'"
	display "Individual coefficients have been stored in the variables: `bitbclist'"
end

mata:
void bc(real scalar m, real scalar tol, string scalar touse) {
	external real matrix dataf, dataf_t, bi, bi_t, bt, bit, Tlist, Nlist, ivar, tvar, ivar_t
	external macro_maxiter, macro_tol, ids, tds
		
	data = dataf[|1,2\.,.|]
	data_t = dataf_t[|1,2\.,.|]
	st_view(bitbc, ., st_local("bitbclist"),touse)
	bcsummary = J(1,m+(4*m), 0)
	bcsummary[|1,1\1,m|]= mean(bit)
	macro_is = length(ids)
	macro_ts = length(tds)
	fail = 0
	
	/* Start by calculating the first round of bias correction */
	/* Set up matrices */
	pool_xx = quadcross(data,data) 
	pool_xxinv = invsym(pool_xx) /* M x M */
	pool_xxbi = J(m,1,0)
	pool_xxbt = J(m,1,0)
	i_xxbt = J(m,macro_is,0)
	t_xxbi = J(m,macro_ts,0)
	
	/* Find the correlation between coefficients and variables over the whole sample */
	Tsum = J(macro_is,1,0)
	Nsum = J(macro_ts,1,0)
	Tall = 0
	Nall = 0
	
	for (i=1; i<=macro_is; i++) {
		t_spec = Tlist[i,1]
		if (i > 1) Tsum[i,1] = Tlist[i-1,1] + Tsum[i-1,1]
		Tall = Tall + Tlist[i,1]
		for (t=1; t<=t_spec; t++) {
			VARS = data[Tsum[i,1] + t,.]
			b_i = bi[Tsum[i,1] + t,.]'
			b_t = bt[Tsum[i,1] + t,.]'
			pool_xxbi = pool_xxbi  + (VARS' * VARS * b_i)
			pool_xxbt = pool_xxbt  + (VARS' * VARS * b_t)
			i_xxbt[.,i] = i_xxbt[.,i] + (VARS' * VARS * b_t)
		}
	}
	for (t=1; t<=macro_ts; t++) {
		n_spec = Nlist[t,1]
		if (t > 1) Nsum[t,1] = Nlist[t-1,1] + Nsum[t-1,1]
		Nall = Nall + Nlist[t,1]
		for (i=1; i<=n_spec; i++) {
			VARS = data_t[Nsum[t,1] + i,.]
			b_i = bi_t[Nsum[t,1] + i,.]'
			t_xxbi[.,t] = t_xxbi[.,t] + (VARS' * VARS * b_i)
		}
	}
	/* Calculate the four components of the bias */
	comp1 = J(macro_is,m,0)
	comp3 = J(macro_ts,m,0)
	for (i=1; i<=macro_is; i++) {
		if (i<macro_is) VARS = data[|(Tsum[i,1]+1),1\(Tsum[i+1,1]),m|]
		else VARS = data[|(Tsum[i,1]+1),1\(Tall),m|]
		i_xx = quadcross(VARS,VARS)
		i_xxinv = invsym(i_xx) /* M x M */
		comp1[i,.] = (i_xxinv*i_xxbt[.,i])' 
	}

	for (t=1; t<=macro_ts; t++) {
		if (t<macro_ts) VARS = data_t[|(Nsum[t,1]+1),1\(Nsum[t+1,1]),m|]
		else VARS = data_t[|(Nsum[t,1]+1),1\(Nall),m|]
		t_xx = quadcross(VARS,VARS)
		t_xxinv = invsym(t_xx) /* M x M */
		comp3[t,.] = (t_xxinv*t_xxbi[.,t])'
	}
	comp4 = (pool_xxinv*pool_xxbi)'
	comp2 = (pool_xxinv*pool_xxbt)'

	/* Apply the bias correction to bit */	
	for (i=1; i<=macro_is; i++) {
		tcount = 1
		for (t=1; t<=macro_ts; t++) {
			tdst = strtoreal(tds[t])
			if (tcount  <= Tlist[i,1]) {
				if (tdst == tvar[Tsum[i,1] + tcount]) {
					bitbc[Tsum[i,1] + tcount,.] = bit[Tsum[i,1] + tcount,.] - (comp1[i,.] - comp2) - (comp3[t,.] - comp4) 
					tcount = tcount + 1
				}
			}
		}
	}
	
	bcsummary[|1,m+1\1,(m+4*m)|] = (mean(comp1), mean(comp2), mean(comp3), mean(comp4))

	/* Parameters for the loop */
	sign = 1
	/* Loop for the remaining iterations */	
	for (s=1; s<= macro_maxiter; s++) {
		pool_comp4 = J(m,1,0)
		pool_comp2 = J(m,1,0)
		i_comp1 = J(m,macro_is,0)
		t_comp3 = J(m,macro_ts,0)

		/* Find the correlation between coefficients and variables over the whole sample */
		for (i=1; i<=macro_is; i++) {
			tcount = 1
			for (t=1; t<=macro_ts; t++) {
				tdst = strtoreal(tds[t])
				if (tcount <= Tlist[i,1]) {
					if (tdst == tvar[Tsum[i,1] + tcount]) {
						VARS = data[Tsum[i,1] + tcount,.]
						c3 = comp3[t,.]'
						i_comp1[.,i] = i_comp1[.,i] + (VARS' * VARS * c3)
						pool_comp2 = pool_comp2  + (VARS' * VARS * c3)
						tcount = tcount + 1
					}
				}
			}
		}
		for (t=1; t<=macro_ts; t++) {
			icount = 1
			for (i=1; i<=macro_is; i++) {
				idst = strtoreal(ids[i])
				if (icount <= Nlist[t,1]) {
					if (idst == ivar_t[Nsum[t,1] + icount]) {
						VARS = data_t[Nsum[t,1] + icount,.]
						c1 = comp1[i,.]'
						t_comp3[.,t] = t_comp3[.,t] + (VARS' * VARS * c1)
						pool_comp4 = pool_comp4  + (VARS' * VARS * c1)
						icount = icount + 1
					}
				}
			}
		}

		
		/* Calculate the four components of the bias */
		comp1 = J(macro_is,m,0)
		comp3 = J(macro_ts,m,0)

		for (i=1; i<=macro_is; i++) {
			if (i<macro_is) VARS = data[|(Tsum[i,1]+1),1\(Tsum[i+1,1]),m|]
			else VARS = data[|(Tsum[i,1]+1),1\(Tall),m|]
			i_xx = quadcross(VARS,VARS)
			i_xxinv = invsym(i_xx) /* M x M */
			comp1[i,.] = (i_xxinv*i_comp1[.,i])' 
		}
		for (t=1; t<=macro_ts; t++) {
			if (t<macro_ts) VARS = data_t[|(Nsum[t,1]+1),1\(Nsum[t+1,1]),m|]
			else VARS = data_t[|(Nsum[t,1]+1),1\(Nall),m|]
			t_xx = quadcross(VARS,VARS)
			t_xxinv = invsym(t_xx) /* M x M */
			comp3[t,.] = (t_xxinv*t_comp3[.,t])'
		}
		comp4 = (pool_xxinv*pool_comp4)'
		comp2 = (pool_xxinv*pool_comp2)'
		
		bcsummary = bcsummary \ (mean(bitbc), mean(comp1), mean(comp2), mean(comp3), mean(comp4))
	
		/* Apply the bias correction to bit */	
		bitbcold = bitbc
		for (i=1; i<=macro_is; i++) {
			tcount = 1
			for (t=1; t<=macro_ts; t++) {
				tdst = strtoreal(tds[t])
				if (tcount <= Tlist[i,1]) {
					if (tdst == tvar[Tsum[i,1] + tcount]) {
						bitbc[Tsum[i,1] + tcount,.] = bitbc[Tsum[i,1] + tcount,.] + sign*((comp1[i,.] - comp2) + (comp3[t,.] - comp4)) 
						tcount = tcount + 1
					}
				}
			}
		}
		sign = sign*(-1)
	
		/* Check for convergence */
		change = mean(mean(bitbc)') - mean(mean(bitbcold)')
		if (abs(change) < tol) s_end = s
		if (abs(change) < tol) break
		if (s == macro_maxiter) s_end = s
		if (s == macro_maxiter) fail = 1
	}
	st_matrix("iterresults", bcsummary)
	st_numscalar("iterations", s_end)
	st_numscalar("failstate", fail)
	
}

void varcalc(real scalar m, string scalar touse) {
	st_view(bitbc, ., st_local("bitbclist"),touse)
	b = mean(bitbc)
	st_matrix("b",b)
}

void poste () {
	external real matrix Tlist, Nlist
	NTlist = Tlist \ Nlist
	g1 = min(NTlist)
	g2 = mean(NTlist)
	g3 = max(NTlist)
	g4 = rows(NTlist)
	st_numscalar("e(g_min)",g1)
	st_numscalar("e(g_avg)",g2)
	st_numscalar("e(g_max)",g3)
	st_numscalar("e(N_g)",g4)
}

end
