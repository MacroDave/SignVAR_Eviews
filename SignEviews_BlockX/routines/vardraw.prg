subroutine VARDraw

	' Draw the VCV matrix (sigma)
	'============================================
	sym inv_sigma_in=inv_sigma_hat/!nobs
	call wish_rnd(inv_sigma_in,!nobs)
	sym sigma_draw = @inverse(inv_sigma_draw)
	
	' Draw the coeffiecient matrix (Ft)
	'============================================
	sym aux1 = @inverse(@transpose(vardata)*vardata)
	sym aux2 = @kronecker(sigma_draw, aux1)
	matrix Fthat_vec = @vec(Ft_hat)
	call mvnrnd(Fthat_vec,aux2)
	Ftdraw = @emult(Ftdraw,@transpose(betacheck)) 'Added to impose block exogeneity
	matrix Ft_draw = @unvec(@transpose(Ftdraw), @rows(Ft))
	'sigma_draw = @emult(sigma_draw,blockshocks) 'Added to impose block exogeneity on B matrix

	'Create Reshaped Coefficient Matrices for Analysis
	'============================================
	matrix(!nvar*!nlag,!nvar) Fs_draw=na
	matrix Fs_draw=@subextract(Ft_draw,1+!nvarx,1)
	matrix Fst_draw=@transpose(Fs_draw)

	matrix(!nvar*!nlag,!nvar*!nlag) Fcomp_draw = 0
	
	matplace(Fcomp_draw,Fst_draw,1,1)
	
	!r=!nvar
	!c=0
	
	for !xxx = 1 to !nvar*!nlag-!nvar
	
		!r=!r+1
		!c=!c+1
	
		Fcomp_draw(!r,!c)=1
	
	next	

'	matrix sigma_draw{!kkk} = sigma_draw
'	delete aux1 aux2 Fthat_vec Fs_draw Fst_draw 
'	rename ft_draw ft_draw_!kkk
'	rename Fcomp_draw Fcomp_draw_!kkk

endsub


