subroutine VARHDecomp2(matrix residuals, matrix invAirf,matrix Fcomp)

	matrix eps = @inverse(invAirf) * @transpose(residuals) 

	'Contribution of each shock
	matrix(!nvarXeq,!nvar) invA_big = 0
	matplace(invA_big,invAirf,1,1)
	
	matrix (!nvar,!nvarXeq) Icomp = 0
	matrix tempi=@identity(!nvar)
	matplace(Icomp,tempi,1,1)
	delete tempi

	for !bb=1 to !nvar
		matrix(!nlag*!nvar,!nobs+1) HDshock_big{!bb} = 0
		matrix(!nvar,!nobs+1) HDshock{!bb} = 0
	next
	
	for !b=1 to !nvar  'for each variable
	
		matrix(!nvar,!nobs+1) eps_big = 0 'matrix of shocks conformable with companion
		matplace(eps_big,@rowextract(eps,!b),!b,2)
	
		for !i = 2 to !nobs+1
	
			matrix temp1
			temp1 = invA_big*@columnextract(eps_big,!i) + Fcomp*@columnextract(HDshock_big{!b},!i-1)
			matplace(HDshock_big{!b},temp1,1,!i)	
		
			matrix temp2
			temp2 = Icomp*@columnextract(HDshock_big{!b},!i)
			matplace(HDshock{!b},temp2,1,!i)
			matrix HDshocke{!b} = @subextract(HDshock{!b},1,2)
		next
	
		delete HDshock_big{!b}
	
	next

	'Initial value
	matrix(!nlag*!nvar,!nobs+1) HDinit_big= 0
	matrix(!nvar,!nobs+1) HDinit = 0
	matplace(HDinit_big,@subextract(@transpose(vardatax),1,1,!nlag*!nvar,1))
    	matplace(HDinit,Icomp*@subextract(HDinit_big,1,1,!nlag*!nvar,1),1,1)

    	for !i = 2 to !nobs+1
        	matplace(HDinit_big,Fcomp*@subextract(HDinit_big,1,!i-1,!nlag*!nvar,!i-1),1,!i)
        	matplace(HDinit,Icomp*@subextract(HDinit_big,1,!i,!nlag*!nvar,!i),1,!i)
    	next

	'Constant
    	matrix(!nlag*!nvar,!nobs+1) HDconst_big = 0
    	matrix(!nvar,!nobs+1) HDconst = 0
	vector(!nlag*!nvar,1) CC = 0
    	if !setexog>0 then
		matplace(CC,@subextract(@transpose(Ft),1,1,!nvar,1),1,1)
        	for !i = 2 to !nobs+1
        		matplace(HDconst_big,CC+Fcomp*@subextract(HDconst_big,1,!i-1,!nlag*!nvar,!i-1),1,!i)
        		matplace(HDconst,Icomp*@subextract(HDconst_big,1,!i,!nlag*!nvar,!i),1,!i)
        	next
    	endif

	'Linear trend
    	matrix(!nlag*!nvar,!nobs+1) HDtrend_big = 0
    	matrix(!nvar,!nobs+1) HDtrend = 0
	vector(!nlag*!nvar,1) TT = 0
    	if !setexog>1 then
		matplace(TT,@subextract(@transpose(Ft),1,2,!nvar,2),1,1)
        	for !i = 2 to !nobs+1
        	matplace(HDtrend_big,TT*(!i-1)+Fcomp*@subextract(HDtrend_big,1,!i-1,!nlag*!nvar,!i-1),1,!i)
        	matplace(HDtrend,Icomp*@subextract(HDtrend_big,1,!i,!nlag*!nvar,!i),1,!i)
        	next
    	endif

	'Reshape Results
	
	matrix(!nlag,1) matnan
	
	for !i=1 to !nvar
	
		matrix(!nobs+!nlag,!nvar) Hshock{!i} = 0 ' [nobs x shock x var]
	
		for !bb=1 to !nvar
			matplace(Hshock{!i},matnan,1,!bb)
			matplace(Hshock{!i},@transpose(@rowextract(hdshocke{!bb},!i)),!nlag+1,!bb)
		next
	
	next

	HDtrend=@transpose(HDtrend)
	HDconst=@transpose(HDconst)
	HDinit=@transpose(HDinit)
	'delete hdshock* hdshocke* eps* *_big

endsub


