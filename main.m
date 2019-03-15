%==========================================================================
%  Identifiability test of fractional-order impedance spectroscopy models 
%                with two constant phase elements
%
%                                    R1
%                              |---/\/\/\----|
%                  Rinf        |             |  1/(C2*s^alpha2)
%             ---/\/\/\--------|             |-------->>---------
%                              |             |
%                              |----->>------|
%                              1/(C1*s^alpha1)
% 
%   Applications: This model is widely used in the study of 
%                 - Energy Storage Systems and
%                 - Biomedical systems
%                 
%
% By using the Grunwald-Letnikov approximation, a discrete time transfer
% function of the system is given by
%            f_(2T+2)*z^(2T+2)+f_(2T+1)*z^(2T+1)+f_(2T)*z^(2T)+..... +f_0
%H(z,theta)=---------------------------------------------------------------
%               z^(2T+2)+g_(2T+1)*z^(2T+1)+g_(2T)*z^(2T)+..... +g_0
% where, 
% theta=[Rinf,R1,C1,C2,alpha1,alpha2]
% 
% The problem is to determine whether there is one-to-one map between the
% theta vector and the coefficients of the tranfer function H(z,theta)? 
%
% Code written by: 
%      Tohid Soleymani Aghdam (Shahid Beheshti University, tehran, Iran)
% Under the supervision of:
%      Seyed Mohammad Mahdi Alavi (Shahid Beheshti University, tehran, Iran)
%      Mehrdad Saif (University of Windsor, Canada)
% 
% Inputs: 
%  Rinf, R1, C1, alpha1, C2, alpha2
%  precision of computations
% Outputs: 
%  Impedance spectra in Nyquits diagram
%  The number of solutions to the structural identifiability equations
%  For global identifiability, there must be only one solution 
%  
%
%
% April 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
close all
clear all
tic
XRinf=input('Rinf=');%typically in (0.05 0.3)
XC1=input('C1=');%typically in (0.5 3)
XR1=input('R1=');%typically in (0.05 0.3)
XC2=input('C2=');%typically in (50 400)
Xalpha1=input('alpha1=');%typically in (0.01 0.99)
Xalpha2=input('alpha2=');%typically in (0.01 0.99)

precision=input('precision=');% must be integer, the larger the more accuracy 

tohid=-1.*ones(length(XRinf)*length(XC1)*length(XR1)*length(XC2)*length(Xalpha1)*length(Xalpha2),1);
indice=1;
ij=1;
jk=1;
for t=1:length(Xalpha1)
    for tt=1:length(Xalpha2)
        for ttt=1:length(XRinf)
            for tttt=1:length(XR1)
                for ttttt=1:length(XC1)
                    for tttttt=1:length(XC2)
%                         alfa1=Xalpha1(t);
%                         alfa2=Xalpha1(tt);
%                         rinf=XRinf(ttt);
%                         r1=XR1(tttt);
%                         c1=XC1(ttttt);
%                         c2=XC2(tttttt);
                          tohid(indice)=NumSOL([Xalpha1(t),Xalpha2(tt), XRinf(ttt), XR1(tttt),XC1(ttttt),XC2(tttttt)],precision);
                          if tohid(indice)==0
                              thetaproblem1(jk,:)=[Xalpha1(t),Xalpha2(tt), XRinf(ttt), XR1(tttt),XC1(ttttt),XC2(tttttt)];
                              jk=jk+1;
                          end
                          if tohid(indice)>1
                             thetaproblem2(ij,:)=[Xalpha1(t),Xalpha2(tt), XRinf(ttt), XR1(tttt),XC1(ttttt),XC2(tttttt)];
                             ij=ij+1;
                          end
                          indice=indice+1;
                    end
                end
            end
        end
    end
end
                        
toc

%%
figure
w=logspace(-4,2,10000);
Z=XRinf+XR1./(XR1*XC1*(1i*w).^(Xalpha1)+1)+1./(XC2*(1i*w).^(Xalpha2));
plot(real(Z),-imag(Z))
ylabel('-Imag\{Z\}')
xlabel('Real\{Z\}')

Msg = sprintf('The number of solutions: %d',tohid);
disp(Msg)

Msg = sprintf('If the number of solutions is 1, the model is globally identifiable');
disp(Msg)
