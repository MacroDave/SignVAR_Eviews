subroutine ChartArrange

	for !i=1 to !nshocks
		string imps=""
		for !j=1 to !nvar
			if !medt=0 then
				matrix(!nsteps,3) cirf{!i}_{!j}

				matplace(cirf{!i}_{!j},@transpose(mirflow{!i}_{!j}),1,1)	
				matplace(cirf{!i}_{!j},@transpose(mirf{!i}_{!j}),1,2)
				matplace(cirf{!i}_{!j},@transpose(mirfupp{!i}_{!j}),1,3)	
	
				freeze(g_irf{!j}) cirf{!i}_{!j}.line
				
				g_irf{!j}.name(1) {lowerb}
				g_irf{!j}.name(2) Median
				g_irf{!j}.name(3) {upperb}
				imps = imps + " " + g_irf{!j}.@name

			else

				matrix(!nsteps,4) cirf{!i}_{!j}

				matplace(cirf{!i}_{!j},@transpose(mirflow{!i}_{!j}),1,1)	
				matplace(cirf{!i}_{!j},@transpose(mirf{!i}_{!j}),1,2)
				matplace(cirf{!i}_{!j},@transpose(mirfupp{!i}_{!j}),1,3)	
				matplace(cirf{!i}_{!j},@transpose(@rowextract(irf{!i}_{!j}),distance),1,4)	
	
				freeze(g_irf{!j}) cirf{!i}_{!j}.line
				
				g_irf{!j}.name(1) {lowerb}
				g_irf{!j}.name(2) Median
				g_irf{!j}.name(3) {upperb}
				g_irf{!j}.name(4) Median_Target
				imps = imps + " " + g_irf{!j}.@name			
	
			endif

		next
		string s{!i}=@replace(s{!i}," ","")
		graph impulse_{s{!i}}.merge {imps}
		delete g_irf*

	next
	delete cirf* g_irf* 

	if !aux=1 then
		if !medaux=0 then		
			for !ii=1 to !nvar
				string fevs=""
				for !jj=1 to !nshocks
					matrix(!nsteps,3) cfevd{!ii}_{!jj}
					
					if !nshocks>1 then
						matrix(!nobs+!nlag,!nshocks) chd{!ii}
						matplace(chd{!ii},@transpose(mhd{!ii}_{!jj}),1,!jj)
					endif	

					matplace(cfevd{!ii}_{!jj},@transpose(mfevdlow{!ii}_{!jj}),1,1)	
					matplace(cfevd{!ii}_{!jj},@transpose(mfevd{!ii}_{!jj}),1,2)
					matplace(cfevd{!ii}_{!jj},@transpose(mfevdupp{!ii}_{!jj}),1,3)	
		
					freeze(g_fevd{!jj}) cfevd{!ii}_{!jj}.line
					g_fevd{!jj}.addtext(t) {s{!jj}}
					g_fevd{!jj}.name(1) {lowerb}
					g_fevd{!jj}.name(2) Median
					g_fevd{!jj}.name(3) {upperb}
					fevs = fevs+ " " + g_fevd{!jj}.@name
		
				next
				graph fevd_{ser{!ii}}.merge {fevs}
				delete g_fevd*
		
				if !nshocks>1 then
					freeze(hd_{ser{!ii}}) chd{!ii}.bar
					hd_{ser{!ii}}.addtext(t) {ser{!ii}}
					for !zz=1 to !nshocks
						hd_{ser{!ii}}.name(!zz) {s{!zz}}			
					next
				endif
			next
			delete cfevd* chd*
			if !nshocks>1 then	
				delete g_fevd*
			endif

		else

			for !ii=1 to !nvar
				string fevs=""
				for !jj=1 to !nshocks
					matrix(!nsteps,1) cfevd{!ii}_{!jj}

					if !nshocks>1 then
						matrix(!nobs+!nlag,!nshocks) chd{!ii}
						matplace(chd{!ii},@subextract(@transpose(hd{!ii}_{!jj}),1,distance,!nobs+!nlag,distance),1,!jj)
					endif

					matplace(cfevd{!ii}_{!jj},@subextract(@transpose(fevd{!ii}_{!jj}),1,distance,!nsteps,distance),1,1)	
		
					freeze(g_fevd{!jj}) cfevd{!ii}_{!jj}.line
					g_fevd{!jj}.addtext(t) {s{!jj}}
					g_fevd{!jj}.name(1) Median_Target
					fevs = fevs+ " " + g_fevd{!jj}.@name
		
				next
				graph fevd_{ser{!ii}}.merge {fevs}
				delete g_fevd*
		
				if !nshocks>1 then
					freeze(hd_{ser{!ii}}) chd{!ii}.bar
					hd_{ser{!ii}}.addtext(t) {ser{!ii}}
					for !zz=1 to !nshocks
						hd_{ser{!ii}}.name(!zz) {s{!zz}}			
					next
				endif
			next
			delete cfevd* chd*
			if !nshocks>1 then
				delete g_fevd*
			endif
		endif
	endif
endsub


