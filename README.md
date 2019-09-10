# SignVAR_Eviews
Sign Restricted VAR in Eviews

Sign Restricted VAR in Eviews
The first code I have uploaded here is a series of eviews programs and subroutines that allows you to estimate a VAR model with identification using sign restrictions. If you unzip the file you should be able to run the codes directly from the folder (at least in theory). There is a fairly general template (run_sr.prg) which you can use to play around with the settings.
The main motivation in writing this code was to get a proper understanding of sign restrictions, which I feel is always best done by coding it yourself and seeing how everything works. 

The two examples include: Uhlig's Results from an agnostic identification procedure (run_uhlig) and 
Stock and Watson's Vector Autoregressions (run_sr). The results are computed using the median response but can also be done using the median target method of Fry and Pagan (this includes the FEVD and Historical Decomposition). 

I used the matlab toolbox from Ambrogio Cesa-Bianchi as the basis for most of the code. This had the benefit of being able to check the codes were working as I went along but may mean that the code is not as efficient as it could be. Some of the sub-routines I adapted from other sources (generally noted in the programs, eg. the random multivariate normal or the QR decomposition). If you find any bugs or problems please email me (email). I hope the code is fairly straightforward to understand (it is heavily commented) but if not let me know.

Multivariate Random Normal
Eviews subroutine to draw from a multivariate Gaussian distribution

QR Decomposition
Eviews subroutine to calculate the QR decomposition of a matrix

Generate Wishart Matrix
Eviews subroutine to generate a random wishart matrix

Some of the subroutines are redundant in the Eviews 11 which has more matrix routines in-built.
