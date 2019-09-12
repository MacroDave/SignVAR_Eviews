subroutine medtarget
'Routine to calculate the median target of Fry and Pagan (2007,2011)
	
	for !rrr = 1 to !nshocks
		for !ccc = 1 to !nvar
			matrix(!ndraws,!nsteps) target{!rrr}_{!ccc}
			for !nnn = 1 to !ndraws
				for !sss = 1 to !nsteps
					target{!rrr}_{!ccc}(!nnn,!sss) = (irf{!rrr}_{!ccc}(!nnn,!sss) - @median(@rowextract(irf{!rrr}_{!ccc},!nnn)))/@stdev(@rowextract(irf{!rrr}_{!ccc},!nnn))
				next
			next
		next
	next

	'next vectorize all the shocks together (tranpose each target row then stack
	'then calculate distance = @transpose(vectorized)*vectorized

	!min=999999
	scalar distance = 0

	for !nnn = 1 to !ndraws
		matrix(!nsteps,!nshocks*!nvar) targets
	
		!xyz=0
		for !rrr = 1 to !nshocks
			for !ccc=1 to !nvar
				!xyz=!xyz+1
				matplace(targets,@transpose(@rowextract(target{!rrr}_{!ccc},!nnn)),1,!xyz)
			next
		next
		vector vtargets=@vec(targets)
		scalar sdist = @transpose(vtargets)*vtargets
		if sdist<!min then
			distance=!nnn
			!min=sdist
		endif
		
	next

	'delete target*
	'delete vtargets
	show !min

endsub


