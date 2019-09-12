subroutine wish_rnd(sym sigma_in,scalar v)

!n=@rows(sigma_in)
!k=@columns(sigma_in)

if !n <> !k then
	stop
else

if !n < !k then
	logmsg "wish_rnd: n must be >= k+1 for a finite distribution"
endif
endif

matrix t=@transpose(@cholesky(sigma_in))
matrix y = @transpose(t)*@mnrnd(!n,v)
matrix w = y*@transpose(y)

if @isobject("inv_sigma_draw")>0 then
	delete inv_sigma_draw
endif

rename w inv_sigma_draw

endsub


