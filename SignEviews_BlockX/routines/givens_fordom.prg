subroutine givens(scalar n1, scalar n2, scalar n3)

logmode on

if n3<>0 then

	for !i=1 to n2
		for !j=1 to n2
			if !j>!i then
				matrix(n1,n1) Qfor{!i}{!j}=@identity(n1)
			endif
		next
	next

	for !i=1 to n3
		for !j=1 to n3
			if !j>!i then
				matrix(n1,n1) Qdom{!i}{!j}=@identity(n1)
			endif
		next
	next

	'Draw from the theta vector of uniform disn between 0 and pi
	'Note in Eviews @acos(-1) = pi

	vector(1) givdraws
	rnd(givdraws)
	givdraws=givdraws*@acos(-1)

	for !i=1 to n2
		for !j=!i+1 to n2

			Qfor{!i}{!j}(!i,!i) = @cos(givdraws(1))
			Qfor{!i}{!j}(!i,!j) = -@sin(givdraws(1))
	
			if !j>!i then
				qfor{!i}{!j}(!j,!i) = @sin(givdraws(1))
				qfor{!i}{!j}(!j,!j) = @cos(givdraws(1))
			endif
	
			if !i=1 and !j=2 then
	
			'Generate Q matrix = multiplication of individual Givens Rotation Matrices
			matrix(n1,n1) Qfor= qfor{!i}{!j}
	
			else
	
				Qfor = Qfor * qfor{!i}{!j}
	
			endif		

		 next
	next

	for !i=1 to n3
		for !j=!i+1 to n3

			Qdom{!i}{!j}(!i,!i) = @cos(givdraws(1))
			Qdom{!i}{!j}(!i,!j) = -@sin(givdraws(1))
	
			if !j>!i then
				qdom{!i}{!j}(!j,!i) = @sin(givdraws(1))
				qdom{!i}{!j}(!j,!j) = @cos(givdraws(1))
			endif
	
			if !i=1 and !j=2 then
	
			'Generate Q matrix = multiplication of individual Givens Rotation Matrices
			matrix(n1,n1) Qdom= qdom{!i}{!j}
	
			else	
				Qdom = Qdom * qdom{!i}{!j}
	
			endif		

		 next
	next

	matrix(n1,n1) Qall = @identity(n1)
	matplace(Qall,@subextract(Qfor,1,1,n2,n2),1,1)
	matplace(Qall,@subextract(Qdom,1,1,n3,n3),n2+1,n2+1)
	delete(noerr) qdom* qfor*

else

	for !i=1 to n2
		for !j=1 to n2
			if !j>!i then
				matrix(n1,n1) Qfor{!i}{!j}=@identity(n1)
			endif
		next
	next

	'Draw from the theta vector of uniform disn between 0 and pi
	'Note in Eviews @acos(-1) = pi

	vector(1) givdraws
	rnd(givdraws)
	givdraws=givdraws*@acos(-1)

	for !i=1 to n2
		for !j=!i+1 to n2

			Qfor{!i}{!j}(!i,!i) = @cos(givdraws(1))
			Qfor{!i}{!j}(!i,!j) = -@sin(givdraws(1))
	
			if !j>!i then
				qfor{!i}{!j}(!j,!i) = @sin(givdraws(1))
				qfor{!i}{!j}(!j,!j) = @cos(givdraws(1))
			endif
	
			if !i=1 and !j=2 then
	
			'Generate Q matrix = multiplication of individual Givens Rotation Matrices
			matrix(n1,n1) Qfor= qfor{!i}{!j}
	
			else
	
				Qfor = Qfor * qfor{!i}{!j}
	
			endif		

		 next
	next

	matrix Qall = Qfor
	delete(noerr) qdom* qfor*

endif

endsub


