# Matlab Code for the test of "structural identifiability of fractional-order impedance spectroscopy models with two constant phase elements and Grunwald-Letnikov approximation"

## T Soleymani-Aghdam, SMM Alavi, M Saif

                                    R1
                              |---/\/\/\----|
                  Rinf        |             |  1/(C2*s^alpha2)
             ---/\/\/\--------|             |-------->>---------
                              |             |
                              |----->>------|
                              1/(C1*s^alpha1)
 
##   Applications
This model is widely used in the study of energy storage systems and biomedical systems
                 

## Model Structure 
By using the Grunwald-Letnikov approximation, a discrete time transfer function of the system is given by
 
             f_(2T+2)*z^(2T+2)+f_(2T+1)*z^(2T+1)+f_(2T)*z^(2T)+..... +f_0
            
H(z,theta)=----------------------------------------------------------------------
 
                  z^(2T+2)+g_(2T+1)*z^(2T+1)+g_(2T)*z^(2T)+..... +g_0
               
 where, 
 
 theta=[Rinf,R1,C1,C2,alpha1,alpha2]
 
 The problem is to determine whether there is one-to-one map between the theta vector and the coefficients of the tranfer function H(z,theta)? 

## Inputs: 
  - Rinf, R1, C1, alpha1, C2, alpha2
  - precision of computations
## Outputs: 
  - Impedance spectra in Nyquits diagram
  - The number of solutions to the structural identifiability equations
  For global identifiability, there must be only one solution 

## Code written by: 
 - Tohid Soleymani Aghdam (Shahid Beheshti University, Tehran, Iran)
 - Seyed Mohammad Mahdi Alavi (Shahid Beheshti University, Tehran, Iran)
 - Mehrdad Saif (University of Windsor, windsor, ON, Canada)
 
## Dates
 - March 15, 2019: First version publihsed on github.
