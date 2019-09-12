subroutine VARprocs_SUR
	'VAR Procedures for Sign Restriction Code
	'Program to collect and construct useful matrices and strings from the VAR and System Objects
	
	'VAR objects
	'============================================
	!coefs = {%varname}.@eqncoef(1) 'number of estimated coefficients per equation
	!nvarx = (!coefs-(!nvar*!nlag)) 'number of exogenous vars
	!nvarXeq = !coefs-!nvarx 'number of endogenous coefficients
	!nobs={%varname}.@eqregobs(1) 'number of observations in the VAR
	!wishdof = !nobs - !coefs 'degrees of freedom for wishart random number
	string varsmpl={%varname}.@smpl 'var sample period
	{%varname}.makeresids'generate residuals
	group resids resid* 'group the residuals
	resids.drop resid 'drop eviews default resid series
	stom(resids,vresiduals) 'convert the VAR residuals to matrix object
	delete resid* 'delete residuals group
	
	'Make VAR objects
	'============================================
	matrix(!nobs,1) mconst=1 'create vector of ones
	'save the endogous variables used in the VAR
	{%varname}.makeendog yvars
	for !abc=1 to !nvar
		string ser{!abc}=yvars.@seriesname(!abc)
	next
	stom(yvars,ymat) 'convert the endog variables into a matrix
	matrix Ctx = {%varname}.@coefmat 'matrix of coefficients from the estimated VAR
	{%varname}.makemodel(mod_{%varname}) 'convert VAR object to a model object
	
	' Assemble VAR data into matrix for impulses/fevd/hdecomp
	'================================================
	'Create Matrix of data in the VAR system
	matrix(!nobs,(!nvar*!nlag+!nvarx)) vardata
	
	'Collect useful VAR data procs
	string s_exvars=mod_{%varname}.@exoglist 'save list of exogenous variables from VAR
	group exvars {s_exvars} 'group exog variables into group
	matrix(!nobs,!nvarXeq) exmat=na 'create exogenous variables matrix
	sym sigma = {%varname}.@residcov 'save residual cov/var matrix from VAR
	matrix invAirf = @cholesky(sigma) 'Cholesky of the residual cov/var matrix
	matrix(!nsteps,!nvar) SE = 0 'For use in the Hdecomp proc
	close {%varname} 'close the VAR object
	
	'Place appropriate exog variables into data matrix
	if !setexog<>0 then
		stom(exvars,exmat_temp) 'stom(exvars,exmat)
		!tempr=@rows(exmat_temp)
		matplace(exmat,@subextract(exmat_temp,!tempr-!nobs+1,1),1,1)

		'Place exogenous data into VAR
		if !setexog=0 then 'no constant
			'Nothing
		else
		if !setexog=1 then 'constnat
			matplace(vardata,mconst,1,1)
		else
		if !setexog=2 or !setexog=3 then 'constant and trends
			'matplace(vardata,mconst,1,1)
			matplace(vardata,@subextract(exmat,!nlag+1,1),1,1)
		else 
		if !setexog=4 then 'includign a dummmy
			'matplace(vardata,mconst,1,1)
			matplace(vardata,@subextract(exmat,!nlag+1,1),1,1)
		endif
		endif
		endif
		endif
	endif
	
	!firstdata=!nlag+1
	!start=0
	'Loop to place endog variables into the matrix ordered by lag
	for !ii=1 to !nvar
		!start=!start+1
		!cc=!ii
		
		for !tt=1 to !nlag	
			if !tt>1 then
				!cc=!cc+!nvar
			endif
			matplace(vardata,@subextract(ymat,!firstdata-!tt,!ii,!nobs+!nlag-!tt,!ii),1,!nvarx+!cc)
			next
	next

	matrix vardatax=@subextract(vardata,1,!nvarx+1)
	
	'Estimate VAR in System Object (provides a nice way of getting vector of coefficients in nicer ordering)
	'==============================================================================	
	c=0
	{surname}.sur
	{surname}.makeresids
	group resids resid*
	resids.drop resid
	stom(resids,sresiduals)
	delete resid*
	matrix Csx = c '(to incorporate the zero coefficients - if use system default procs it skips missing coefficients)
	sym sigma={surname}.@residcov
	matrix invAirf = @cholesky(sigma)
	
	'Setup Coefficient Matrix in different ordering (exog first, ordered by variable, then lag)
	'==================================================================
	matrix(!nvar*!nlag,!nvar) Fs
	
	!a=0
	for !cc = 1 to !nvar
		for !rr = 1 to (!nvar*!nlag)
				Fs(!rr,!cc) = Csx(!rr+!a) 
		next
		!a = !a + (!nvar*!nlag + !nvarx)
	next 
	
	matrix(!coefs,!nvar) Ft 'matrix of coefficients ordered with exogenous variables first
	
	if !setexog=0 then
		matplace(Ft,Fs,1+!nvarx,1)
	else
		matplace(Ft,@subextract(Ctx,!coefs-!nvarx+1,1,!coefs,!nvar))
		matplace(Ft,Fs,1+!nvarx,1)
	endif
	
	'Create Companion Matrix for Impulse Responses
	'===================================================================
	matrix Fst=@transpose(Fs) 'transpose coefficient matrix
	matrix(!nvar*!nlag,!nvar*!nlag) Fcomp = 0 'create the companion matrix
	matplace(Fcomp,Fst,1,1) 'place some coefficients in correct spot in comp. matrix
	
	!r=!nvar
	!c=0
	!zzz=0
	
	'Add 1s to relevant cells of the comp. matrix
	while !zzz < !nvar*!nlag-!nvar
		!zzz=!zzz+1
		!r=!r+1
		!c=!c+1
		Fcomp(!r,!c)=1
	wend
	
	'Store Original Coefficient and Covariance Matrix for VARdraw proc
	matrix Ft_hat = Ft
	matrix Ft_hat_hold = Ft_hat
	matrix Fhat_vec_hold = @vec(Ft_hat_hold)
	matrix(@rows(Fhat_vec_hold),1) betacheck1=0
	matrix(@rows(Fhat_vec_hold),1) betacheck2=@abs(Fhat_vec_hold)
	matrix betacheck = @egt(betacheck2,betacheck1)

	matrix(!nvar,!nvar) blockshocks=1
	if @wcount(block1)>0 then
		for !i=1 to @wcount(block1)
			for !j=@wcount(block1)+1 to !nvar
				blockshocks(!i,!j) = 0
			next
		next
	endif

	sym inv_sigma_hat = @inverse(sigma)			
	sym inv_sigma_hat_hold = inv_sigma_hat

	'Create LU factorization of partition of covariance variance matrix
	matrix(!nvar,!nvar) Vtilde = 0
	matrix(@wcount(block1),@wcount(block1)) V11 = @subextract(sigma,1,1,@wcount(block1),@wcount(block1))
	matrix(@wcount(block1),@wcount(block2)) V12 = @subextract(sigma,1,@wcount(block1)+1,@wcount(block1),!nvar)
	matrix(@wcount(block2),@wcount(block1)) V21 = @subextract(sigma,@wcount(block1)+1,1,!nvar,@wcount(block1))
	matrix(@wcount(block2),@wcount(block2)) V22 = @subextract(sigma,@wcount(block1)+1,@wcount(block1)+1,!nvar,!nvar)
	matplace(Vtilde,V11,1,1)
	matplace(Vtilde,V22-V21*@inverse(V11)*V12,@wcount(block1)+1,@wcount(block1)+1)
	sym symtilde=Vtilde

	'Misc Objects
	'===================================================================	
	string lowerb = @str(!pctilel*100)+"th "+"percentile"
	string upperb = @str(!pctileu*100)+"th "+"percentile"
	string allvars = @wunion(block1,block2)

	svector(!nvar) s_vars
	s_vars=@wsplit(allvars)

endsub


