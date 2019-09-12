subroutine sur_append(string surname, string block1, string block2)
	svector(!setexog) exogs
	exogs = @wsplit(g_exog)

if @wcount(block1)=0 then

	!n_coef=0
	for %v1 {block2}
		string sys=%v1+"="
		
		for !l=1 to !nlag
			for %v2 {block2}
				!n_coef=!n_coef+1
				sys = sys + "+" + "c("+@str(!n_coef)+")"+"*"+%v2+"(-"+@str(!l)+")"
			next	
		next

		if !setexog<>0 then		
			for !x=1 to @rows(exogs)
				!n_coef = !n_coef + 1
				sys = sys + "+" + "c("+@str(!n_coef)+")"+"*"+exogs(!x)
			next
		endif
	
		rename sys sys_{%v1}
	
		{surname}.append {sys_{%v1}}
		delete sys_{%v1}
	next

else

	!n_coef=0
	for %v1 {block1} {block2}
		string sys=%v1+"="
		
		for !l=1 to !nlag
			if @wfind(block1,%v1) then
				for %v2 {block1}
					!n_coef=!n_coef+1
					sys = sys + "+" + "c("+@str(!n_coef)+")"+"*"+%v2+"(-"+@str(!l)+")"
				next	
				!n_coef=!n_coef+@wcount(block2)
			else
				for %v2 {block1} {block2}		
					!n_coef=!n_coef+1
					sys = sys + "+" + "c("+@str(!n_coef)+")"+"*"+%v2+"(-"+@str(!l)+")"		
				next
			endif
		next
		
		if !setexog<>0 then
			for !xyx=1 to @rows(exogs)
				!n_coef = !n_coef + 1
				sys = sys + "+" + "c("+@str(!n_coef)+")"+"*"+exogs(!xyx)
			next
		endif
		
		rename sys sys_{%v1}
	
		{surname}.append {sys_{%v1}}
		delete sys_{%v1}
	next
endif

endsub


