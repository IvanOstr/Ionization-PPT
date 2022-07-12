# Ionization-PPT
Ionization rate in EM field PPT formulation

I've used A.14 formula in the paper ["Ultrashort filaments of light in weakly-ionized, optically-transparent media"][1] to calculate the dependence of ionization rate in oxygen molecules in air on laser intensity.
I'm interested in comparing it with the famous sigma_k*I**k law of MPI ionization. 
However, there is something wrong in the simple simulation I wrote since the ionization rate in PPT (using  is much higher than we expect.
Do you have a caveat for why is it happening?

I think the error is in the Z value in the n_eff
