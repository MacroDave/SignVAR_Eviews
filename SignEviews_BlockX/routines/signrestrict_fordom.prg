subroutine SignRestrict
		
	'Create Matrices to Hold Correct Sign Restriction Impulses
	for !i = 1 to !nshocks
		for !j = 1 to !nvar
			matrix(!ndraws,!nsteps) irf{!i}_{!j}    'the !j response to a !i shock
			matrix(1,!nsteps) mIRF{!i}_{!j}    'the median !j response to a !i shock
			matrix(1,!nsteps) mIRFLow{!i}_{!j}    'the lower percentile !j response to a !i shock
			matrix(1,!nsteps) mIRFUpp{!i}_{!j}    'the upper percentile !j response to a !i shock
			if !aux=1 then
				'FEVD
				matrix(!ndraws,!nsteps) fevd{!j}_{!i}    'FEVD !ith variable due to !jth shock
				matrix(1,!nsteps) mFEVD{!j}_{!i}    'FEVD !ith variable due to !jth shock - median
				matrix(1,!nsteps) mFEVDLow{!j}_{!i}    'FEVD !ith variable due to !jth shock - lower percentile
				matrix(1,!nsteps) mFEVDUpp{!j}_{!i}    'FEVD !ith variable due to !jth shock - upper percentile
				'HDECOMP
				matrix(!ndraws,!nobs+!nlag) hd{!j}_{!i}    'HD !j response to a !i shock
				matrix(1,!nobs+!nlag) mHD{!j}_{!i}    'HD !j response to a !i shock
				matrix(1,!nobs+!nlag) mHDlow{!j}_{!i}    'HD !j response to a !i shock
				matrix(1,!nobs+!nlag) mHDupp{!j}_{!i}    'HD !j response to a !i shock
			endif
		next
	next
	
	''Sign restriction routine
	'==========================================================
	!jj = 0 'accepted draws
	!kkk=0 'counter to stop procedure if restrictions are not working
	!ww = 1 ' index for printing on screen; =1 is intervals of 10 loops shown when completed
	
	while !jj < !ndraws 'do until ndraws succesfully found
	!kkk = !kkk+1
	statusline Draw Number: !kkk Success: !jj
	
		if !jj=0 and !kkk>4000000*!ndraws then
			logmsg Error: Max number of iterations reached. Try different sign restrictions
			stop	
		endif
	
		'1. Draw a random orthonormal matrix
		call givens(!nvar,!nshocks_for,!nshocks_dom)

		'2. Use rotation matrix to generate new structural shocks
		matrix invAirf_draw = @cholesky(sigma)*Qall

		'3 Compute IRFs given new rotation
		matrix Fcomp_draw = Fcomp
		call VARir(Fcomp,invAirf_draw,cumul) 'invAirf_draw

		' ... and check whether sign restrictions are satisfied
		matrix(!nvar,1) checkall=0 '!nshocks columns
		matrix(!nvar,1) checkall_flip=0
		vector(!nshocks) success=0

		for !ss1=1 to !nshocks 'loop over the estimated impulses
			for !ss2=1 to !nshocks		'loop over the sign shock restrictions
				checkall=0
				checkall_flip=0
				for !ii = 1 to !nvar
					if R{!ss2}(!ii,1) = 0 then
						checkall(!ii,1) = 1
						checkall_flip(!ii,1) = 1
					else
		
					if R{!ss2}(!ii,3) = 1 then
						!check=0
						for !rr=R{!ss2}(!ii,1) to R{!ss2}(!ii,2)
							if {%irfname}{!ss1}(!ii,!rr)>0 then
								!check=!check+1
							endif
						next
						
						if !check=( R{!ss2}(!ii,2)-R{!ss2}(!ii,1) +1 ) then
							checkall(!ii,1) = 1
						else
							checkall(!ii,1) = 0
						endif

						' Check flipped signs
						!check_flip=0
						for !rr=R{!ss2}(!ii,1) to R{!ss2}(!ii,2)
							if {%irfname}{!ss1}(!ii,!rr)<0 then
								!check_flip=!check_flip+1
							endif
						next
	
						if !check_flip=( R{!ss2}(!ii,2)-R{!ss2}(!ii,1) +1 ) then
							checkall_flip(!ii,1) = 1
						else
							checkall_flip(!ii,1) = 0
						endif 	

					else

					if R{!ss2}(!ii,3) = -1 then
						!check=0
						for !rr=R{!ss2}(!ii,1) to R{!ss2}(!ii,2)
							if {%irfname}{!ss1}(!ii,!rr)<0 then
								!check=!check+1
							endif
						next
						
						if !check=(R{!ss2}(!ii,2)-R{!ss2}(!ii,1)+1) then
							checkall(!ii,1) = 1
						else
							checkall(!ii,1) = 0
						endif			
						' Check flipped signs
						!check_flip=0
						for !rr=R{!ss2}(!ii,1) to R{!ss2}(!ii,2)
							if {%irfname}{!ss1}(!ii,!rr)>0 then
								!check_flip = !check_flip+1
							endif
						next
	
						if !check_flip=(R{!ss2}(!ii,2)-R{!ss2}(!ii,1)+1) then
							checkall_flip(!ii,1) = 1
						endif  
					endif
					endif
					endif
				next

				if @cmean(checkall)=1 then
					success(!ss2) = !ss1				
				else
				if @cmean(checkall_flip)=1 then
					success(!ss2) = !ss1
					{%irfname}{!ss1}=-{%irfname}{!ss1}	
				endif
				endif	
			next
		next

		if @csum(success) = !signcheck then
			!jj = !jj+1
			show success
			stop

			if !aux=1 then 'do FEVD and Hdecomp
	
'				call VARFevd(Fcomp,invAirf_draw,sigma_draw)
				call VARFevd2(Fcomp_draw,invAirf_draw,!nsteps)
				call VARHDecomp2(vresiduals, invAirf_draw,Fcomp_draw)
			
				for !rrr = 1 to !nshocks
					for !ccc = 1 to !nvar
						!location = success(!rrr)
						matrix tempstoreirf = @subextract({%irfname}{!location},!ccc,1,!ccc,!nsteps)
						matplace(irf{!rrr}_{!ccc},tempstoreirf,!jj,1)
						
						matplace(fevd{!ccc}_{!rrr},@subextract(@transpose(fevd{!ccc}),!rrr,1,!rrr,!nsteps),!jj,1)
						matplace(hd{!ccc}_{!rrr},@subextract(@transpose(hshock{!ccc}),!location,1,!location,!nobs+!nlag),!jj,1)
'						matplace(hd{!ccc}_{!rrr},@subextract(@transpose(hshock{!ccc}),!rrr,1,!rrr,!nobs+!nlag),!jj,1)
					next
				next

			else 'Do not do FEVD and HDecomp just impulses
	
				for !rrr = 1 to !nshocks
					for !ccc = 1 to !nvar
						!location = success(!rrr)
						matrix tempstoreirf = @subextract({%irfname}{!location},!ccc,1,!ccc,!nsteps)
						matplace(irf{!rrr}_{!ccc},tempstoreirf,!jj,1)
					next
				next
			endif
		endif
	
		'Reset Ft_draw and inv_sigma matrices, otherwise they seem to get cumulated during the loop for some reason
		Ft_hat=Ft_hat_hold
		inv_sigma_hat=inv_sigma_hat_hold
'		delete Fcomp_draw
'		delete Ft_draw
'		delete sigma_draw
	
	wend 'Finished sign restriction procedure
	
	'Create Summaries of impulses, FEVD etc (median, lower and upper percentiles)
	'================================================================
	
	'use @columnextract to pull out each column 1x1 and take the median then compile into one median and for pctile ranges
	for !shock=1 to !nshocks
		for !response = 1 to !nvar
			for !r = 1 to (!nobs+!nlag)
				if !r<!nsteps+1 then
					'Impulses
					mIRF{!shock}_{!response}(1,!r) = @median(@columnextract(irf{!shock}_{!response},!r))
					mIRFLow{!shock}_{!response}(1,!r) = @quantile(@columnextract(irf{!shock}_{!response},!r),!pctilel)
					mIRFUpp{!shock}_{!response}(1,!r) = @quantile(@columnextract(irf{!shock}_{!response},!r),!pctileu)
	
					if !aux=1 then
						'FEVD
						mFEVD{!response}_{!shock}(1,!r) = @median(@columnextract(FEVD{!response}_{!shock},!r))
						mFEVDLow{!response}_{!shock}(1,!r) = @quantile(@columnextract(FEVD{!response}_{!shock},!r),!pctilel)
						mFEVDUpp{!response}_{!shock}(1,!r) = @quantile(@columnextract(FEVD{!response}_{!shock},!r),!pctileu)
						'HDecomp
						mHD{!response}_{!shock}(1,!r) = @median(@columnextract(Hd{!response}_{!shock},!r))
						mHDLow{!response}_{!shock}(1,!r) = @quantile(@columnextract(Hd{!response}_{!shock},!r),!pctilel)
						mHDUpp{!response}_{!shock}(1,!r) = @quantile(@columnextract(Hd{!response}_{!shock},!r),!pctileu)
					endif
				else
					if !aux=1 then
						'HDecomp
						mHD{!response}_{!shock}(1,!r) = @median(@columnextract(Hd{!response}_{!shock},!r))
						mHDLow{!response}_{!shock}(1,!r) = @quantile(@columnextract(Hd{!response}_{!shock},!r),!pctilel)
						mHDUpp{!response}_{!shock}(1,!r) = @quantile(@columnextract(Hd{!response}_{!shock},!r),!pctileu)
					endif
				endif										
			next
		next
	next
	
endsub


