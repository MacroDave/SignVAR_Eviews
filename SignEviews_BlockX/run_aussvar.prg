'  Description: Sign Restrictions Procedure in Eviews
'  Programmer: David Stephan, World Bank Macroeconomics and Fiscal Global Practice
'  Last Update: 10/12/2015
'___________________________________________________________________________________________________________________________________________________________________________
tic 'start program timer

'rndseed(12345) 'uncomment for reproducibility

close @all
'logmode all
logmode l

'change path to Eviews sample programs filepath
%path = @runpath
cd %path

include .\routines\VARir.prg 'Impulse Response Code
include .\routines\VARHDecomp2.prg 'Historical Decomposition Code
include .\routines\VARFevd2.prg 'Variance Decomp Code
include .\routines\VARDraw.prg 'Code to Draw from the posterior distribution of a VAR model
include .\routines\wish_rnd.prg 'Generates a random wishart matrix
include .\routines\mvnrnd.prg 'Draw n random d-dimensional vectors from a multivariate Gaussian distribution
include .\routines\qr.prg 'QR matrix decomposition
include .\routines\givens_fordom.prg 'Givens Rotation
include .\routines\VARprocs_sur.prg 'Collects and creates useful objects for VAR analysis
include .\routines\signrestrict_fordom.prg 'Runs the sign restriction procedure
include .\routines\medtarget.prg 'Calculates the median target for impulses
include .\routines\ChartArrange.prg 'Rearranges results into nicer WF objects
include .\routines\sur_append.prg 'Writes the system object with block exog variables

' create workfile
'============================================
'Open workfile or create new workfile and read in data in this section
wfopen .\data\aussvar_data.wf1

'Set Parameters
'============================================
!nsteps = 40 'number of periods impulses calculated for
!nlag=2 'number of lags in the VAR
!ndraws=10 'number of successful draws required
%varname = "AUSSVAR" 'name of VAR
!pctilel=0.16 'lower percentile for impulse responses
!pctileu=0.84 'upper percentile for impulse responses
!aux=1 '=1 if you want FEVD and Hdecomp to be computed during the sign restriction procedures
!SR = 1 '=1 if you want SR or 0 if you want standard cholesky decomp
!medt=0 '=1 if you want to calculate the median target
!medaux=0 '=1 if you want FEVD and HD with the median target

'Set exogenous
'============================================
'0 no constant 
'1 constant
'2 constant and trend
'3 constant and trend and trend^2
!setexog=1

if !setexog=0 then
	string g_exog=""
endif 
if !setexog=1 then
	string g_exog="const"
	series const=1
endif
if !setexog=2 then
	string g_exog="const trend"
	series const=1
	series trend=@trend
endif
if !setexog=3 then
	string g_exog="const trend trend2"
	series const=1
	series trend2=@trend^2
endif

if !setexog=4 then
	string g_exog="const dum1"
	series const=1
	series dum1=0
	smpl 1993q1 @last
	dum1=1 
endif

smpl @all

' estimate unrestricted VAR 'Change the variables after !nlag to the ones you want in your model
'===========================================================================
'Don't replace the @ {g_exog} terms
var {%varname}.ls 1 !nlag IRUS GDPUS CPIUS GDP IR CPIEX RER TOT @ {g_exog}

string block1="IRUS GDPUS CPIUS " 'list the block exog variables
string block2="GDP IR CPIEX RER TOT" 'list the non-exog variables
'***NOTE if want all variables in all equations set block1="" and put all vars in block2  ******

!nvar={%varname}.@neqn 'number of variables in VAR system

' set which impulses you want to be accumulated or not (=1 for accumulation)
'==============================================================
vector(!nvar) cumul = 0
'cumul.fill 1,1,1,1,0,0,1             '(set the corresponding row number of the variable to 1 if you want accumulated)

' Collect and create useful objects for VAR analysis
'==============================================================
string surname="sys"+%varname
system {surname}

call sur_append(surname,block1,block2)

call VARprocs_SUR

' Determine the Identification (SR=0 is Cholesky, SR=1 is sign restriction)
'==============================================================

if !SR=0 then

	call VARir(Fcomp,invAirf,cumul)

	if !aux=1 then
	
		call VARFevd(Fcomp,invAirf,sigma)
		call VARHDecomp(sresiduals, invAirf,Fcomp)

	endif

else

	' Set the shocks for the VAR Model
	'==============================================================

	!nshocks_for = 3
	!nshocks_dom = 5
	!nshocks=!nshocks_dom + !nshocks_for

	'Create restriction matrices
	for !r=1 to !nshocks
		matrix(!nvar,3) R{!r}
	next
	
	'Fill the shock matrices according to the following method 
	'  The first column defines the first period from which the restriction is imposed 
	'  The second column defines the last period to which the restriction is imposed 
	'  The last column defines the sign of the restriction: positive (1) or negative (-1)
	'  To leave unrestricted set all three columns to zero

	string s1="Foreign Interest Rate" 'name shock for graphs
	R1.fill(b=r) _
	1,2,1, _
	0,0,0, _ 
	1,2,-1, _
	0,0,0, _
	0,0,0, _
	0,0,0, _
	0,0,0, _
	0,0,0
	
	string s2="Foreign Output"
	R2.fill(b=r) _
	1,2,1, _
	1,2,1, _ 
	0,0,0, _
	0,0,0, _
	0,0,0, _
	0,0,0, _
	0,0,0, _
	0,0,0

	string s3="Foreign Inflation"
	R3.fill(b=r) _
	1,2,1, _
	1,2,-1, _ 
	1,2,1, _
	0,0,0, _
	0,0,0, _
	0,0,0, _
	0,0,0, _
	0,0,0

	string s4="Output" 'name shock for graphs
	R4.fill(b=r) _
	0,0,0, _
	0,0,0, _
	0,0,0, _
	1,2,1, _
	1,2,1, _
	0,0,0, _
	0,0,0, _
	1,2,-1
	

	string s5="Interest Rate"
	R5.fill(b=r) _
	0,0,0, _
	0,0,0, _
	0,0,0, _
	0,0,0, _
	1,2,1, _
	1,2,-1, _
	0,0,0, _
	0,0,0

	string s6="Cost-push"
	R6.fill(b=r) _
	0,0,0, _
	0,0,0, _
	0,0,0, _
	1,2,-1, _
	1,2,1, _
	1,2,1, _
	1,2,1, _
	0,0,0

	string s7="Risk Premium" 'name shock for graphs
	R7.fill(b=r) _
	0,0,0, _
	0,0,0, _
	0,0,0, _
	0,0,0, _
	0,0,0, _
	1,2,-1, _
	1,2,1, _
	0,0,0
	
	string s8="Terms of Trade"
	R8.fill(b=r) _
	0,0,0, _
	0,0,0, _
	0,0,0, _
	0,0,0, _
	1,2,1, _
	1,2,1, _
	1,2,1, _
	1,2,1

	!signcheck = (!nshocks*(!nshocks+1))/2

	' Run the sign restriction procedures
	'==============================================================
	call SignRestrict
	
	' Calc Median Target if wanted
	'==============================================================
	if !medt=1 then
		call medtarget
	endif
	
	' Rearrange Results for Easier Viewing
	'==============================================================
	call ChartArrange

endif

'Cleanup workfile
'============================================
'delete resid* fcomp* eps* csx ctx fs* ft* temp* sigma* invairf* checkall* *_hold temp* exvars exogvars vardata yvars cumul g_exog inv_sigma* matu matuu mconst orthnorm s_exvars sign sresiduals vresiduals ymat varsmpl varmod subr invamult* upperb lowerb s? ser? 'fevd?_? mhd?_? mfevd?_? mirflow?_? mhdlow?_? mfevdlow?_? mirfupp?_? mhdupp?_? mfevdupp?_? hshock* 
if !aux=1 then
	'delete cumulfevd fevd? hshock? FECD* MSE* PSI* mult* icomp inva_big matnan column 
endif
if !setexog<>0 then 
	delete exmat const 
endif	

'Display length of time taken for the program to run (in seconds).
toc
!elapsed = @toc
logmsg Time Taken: !elapsed seconds


