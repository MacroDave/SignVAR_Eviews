subroutine VARFevd2(matrix Fcomp, matrix invAirf,scalar fsteps)

	string dave=""
	matrix IK = @identity(!nvar)
	matrix(!nvar,!nvar*!nlag) JJ = 0
	matplace(JJ,IK,1,1)
	matrix(!nvar,!nvar) MSE = 0
	matrix(!nvar,!nvar) FE = 0
	matrix(!nvar,!nvar) WFEVD = 0
	for !i=1 to !nvar
		matrix(fsteps,!nvar) FEVD{!i}
	next
	
	'!h=1
	for !i=1 to fsteps
		if !i=1 then
			%s_fcomp="@identity(!nvar*!nlag)"
			matrix tmp = (JJ*({%s_fcomp})*@transpose(JJ))*invAirf
		else
			%s_fcomp=%s_fcomp+"*Fcomp"
			matrix tmp = (JJ*({%s_fcomp})*@transpose(JJ))*invAirf
		endif
		MSE = MSE + tmp*@transpose(tmp)
	
	'	!h=!h+1
	
		for !j=1 to !nvar
			for !k=1 to !nvar
				!FEtemp=(@transpose(@subextract(IK,1,!j,!nvar,!j))*tmp*@subextract(IK,1,!k,!nvar,!k))*(@transpose(@subextract(IK,1,!j,!nvar,!j))*tmp*@subextract(IK,1,!k,!nvar,!k))
				FE(!j,!k) = FE(!j,!k) + !FEtemp
				WFEVD(!j,!k) = FE(!j,!k) / MSE(!j,!j)
			next
		next
		
		for !z=1 to !nvar
			matplace(FEVD{!z},@subextract(WFEVD,!z,1,!z,!nvar),!i,1)
			for !x=1 to !nshocks
				!fevdlocation=success(!x)
				FEVD{!z}(!i,!x) = WFEVD(!z,!fevdlocation)
			next
		next
	next

endsub


