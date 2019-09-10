'Code adapted from Matlab code: http://statweb.stanford.edu/~candes/acm106a/Handouts/Ur.pdf

subroutine qr(matrix A)

scalar m = @rows(A)
scalar n = @columns(A)

matrix(m,n) MATU	
matrix R = A

for !k = 1 to n
	matrix x = @subextract(R,!k,!k,m,!k)
	vector(@rows(x)) e
	e(1) = 1

	if x(1)> 0 then
		scalar sign=1
	else 

	if x(1)<0 then
		scalar sign=-1
	else

	if x(1)=0 then
		scalar sign=0
	endif
	endif
	endif

	vector MATUU = sign*@norm(x,2)*e + x
	MATUU = MATUU/@norm(MATUU,2)
	vector temp_uu_{!k}=MATUU

	matrix subR = @subextract(R,!k,!k,m,n)
	subR = subR - 2*MATUU*@transpose(MATUU)*subR
	matplace(R,subR,!k,!k)

	matrix(m,n) MATU	
	matplace(MATU,MATUU,!k,!k)

next

'Clean up small values of R matrix to be zero

for !rr = 1 to n
	for !cc= 1 to n
		if abs(R(!rr,!cc))<1e-14 then
			R(!rr,!cc) = 0
		endif
	next
next

matrix S = A*@inverse(R)

'Clean up error whereby the last column of Q is the wrong sign

for !i = 1 to m

	S(!i,n) = S(!i,n)*-1

next
	
endsub


