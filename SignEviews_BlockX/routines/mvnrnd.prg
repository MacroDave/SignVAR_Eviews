subroutine mvnrnd(matrix mu,sym matsigma)

	'Draw n random d-dimensional vectors from a multivariate Gaussian distribution with mean mu and covariance matrix matsigma.
	'Adapting Code by Iain Murray 2003 
	'http://homepages.inf.ed.ac.uk/imurray2/code/matlab_octave_missing/mvnrnd.m
	
	if @columns(mu)=1 and @rows(matsigma)<>1 and @columns(matsigma)<>1 then
		mu=@transpose(mu)
	endif
	
	!n=@rows(mu)
	!d=@columns(mu)
	!ns=@rows(matsigma)
	!ds=@columns(matsigma)
	
	if !ns<>!d or !ds<>!d then
		logmsg "matsigma must have dimensions dxd where mu is nxd"
		stop
	
	else
	
		!oldcount = @errorcount
		
		matrix uu=@transpose(@cholesky(matsigma))
		
		!newcount = @errorcount
		
		
		if !newcount>!oldcount then
		
			matrix eval = @eigenvalues(matsigma)
			matrix lambda=@eigenvectors(matsigma)
		
			if (min(diag(lambda))<0),error('matsigma must be positive semi-definite.'),end
				
				matrix uu = @sqrt(lambda)*@transpose(eval)
		
			endif
		
		endif
		
		matrix temprand = @mnrnd(!n,!d)
		matrix tempmult = temprand*uu
		vector ss = tempmult + mu 'matrix
		'vector s = @mnrnd(!n,!d)*uu + mu 'matrix 
	
		if @isobject("ftdraw") then
			delete ftdraw
		endif
			
		rename ss Ftdraw
		
	endif

	delete uu

endsub


