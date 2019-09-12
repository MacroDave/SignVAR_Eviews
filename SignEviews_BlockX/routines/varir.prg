subroutine VARir(matrix Fcomp, matrix invAirf, vector cumul)

	%irfname=invAirf.@name
'	show %irfname
	
	''Compute the impulse response
	'===============================================
	
	for !mm=1 to !nvar

		'Initialize the impulse response vector
		matrix(!nvar,!nsteps) {%irfname}{!mm} = 0
		
		'Create the impulse vector
		matrix(!nvar,1) impulse = 0
  			
		'Set the size of the shock
		impulse(!mm,1) = 1 ' one stdev shock
		
		'First period impulse response (=impulse vector)
		matrix temp = invAirf*impulse
		matplace({%irfname}{!mm},temp,1,1)
		
		'Make it comparable with companion
		vector(!nvar*!nlag) impulse_big = 0
		matplace(impulse_big,temp,1,1)
		
		'Recursive computation of impulse response
		matrix Fcomp_eye = @identity(!nvar*!nlag)
		
		for !kk = 2 to !nsteps
  			Fcomp_eye = Fcomp_eye * Fcomp
   			vector(!nvar*!nlag) irf_big{!kk} = Fcomp_eye * impulse_big
  			for !rr = 1 to !nvar
  					{%irfname}{!mm}(!rr,!kk) = irf_big{!kk}(!rr)
				next
		next
		
		delete irf_big*

		for !rr = 1 to !nvar
			if cumul(!rr)=1 then
				for !kk =1 to !nsteps
					if !kk =1 then
						{%irfname}{!mm}(!rr,!kk) = {%irfname}{!mm}(!rr,!kk)
					else
						{%irfname}{!mm}(!rr,!kk) = {%irfname}{!mm}(!rr,!kk) + {%irfname}{!mm}(!rr,!kk-1) 
					endif
				next
			endif
		next

	next
	
	delete temp impulse_big impulse Fcomp_eye 
	
endsub


