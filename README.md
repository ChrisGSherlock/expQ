# expQ
Product of non-negative vector and exponential of a rate matrix.
RCppArmadillo package to evaluate nu^T exp(Q) for a non-negative vector nu and a sparse or dense rate matrix Q.
Two method choices: scaling and squaring or uniformisation.
Also includes code to evaluate exp(Q) using scaling and squaring.
## Accompanying paper
Sherlock (2020) Direct statistical inference for finite Markov jump processes via the matrix exponential. Submitted.
## Files
expQ_1.0.tar.gz - library. From the terminal use: R CMD INSTALL expQ_1.0.tar.gz  
rexpQ.cpp - if R version is too early for the library then within your R file use: sourceCpp("rexpQ.cpp")  
Eyam.r - performs the benchmarking reported in the accompanying paper, demonstrating use of the functions in earnest  
genericQ.r - ancillary file that sets up the Q matrices for a few different models and, in particular, the SIR model  

Further, very simple demonstration calls to the functions appear in the help files for the individual functions. If you are unable to load the package then search for the @examples in expQ.cpp

You will need to install the C++ library Boost before using the expQ package.
