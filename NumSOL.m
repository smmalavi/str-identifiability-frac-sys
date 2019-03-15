% this function determines the number of solutions for  transfer
% function coefficients. 
% its input is the required precision and original parameter vector

function Nanswer=NumSOL(theta,precision)
theta=vpa(theta,precision);
alpha1=theta(1);alpha2=theta(2);Rinf=theta(3);R1=theta(4);C1=theta(5);C2=theta(6);
Ts=vpa(1/2000,precision); % sampling time which is used in f and g calculation.
T=15; % the number of samples which determines the degree of transfer function. 
R2=inf; %for Warburg Term
%% This part calculates the f and g for given input.

d=vpa(Rinf,precision); 
b1=vpa(Ts^alpha1/C1,precision);
b2=vpa(Ts^alpha2/C2,precision);

a10=vpa(alpha1-Ts^alpha1/(R1*C1),precision);
a20=vpa(alpha2-Ts^alpha2/(R2*C2),precision);

for j=1:T
a(1,j)=vpa(-(-1)^(j+1)*gamma(alpha1+1)/(gamma(j+2)*gamma(alpha1-j)),precision);

a(2,j)=vpa(-(-1)^(j+1)*gamma(alpha2+1)/(gamma(j+2)*gamma(alpha2-j)),precision);
end

g2Tp1=vpa(-a10-a20,precision); m=g2Tp1;
g2T=vpa(-a(1,1)-a(2,1)+a10*a20,precision); n=g2T;
g2T_1=vpa(a10*a(2,1)+a(1,1)*a20-a(1,2)-a(2,2),precision); p=g2T_1;
g2T_2=vpa(a10*a(2,2)+a20*a(1,2)+a(1,1)*a(2,1)-a(1,3)-a(2,3)); q=g2T_2;
g2T_3=vpa(a10*a(2,3)+a20*a(1,3)+a(1,1)*a(2,2)+a(1,2)*a(2,1)-a(1,4)-a(2,4),precision); r=g2T_3;
g2T_4=vpa(-a(1,5)-a(2,5)+a10*a(2,4)+a20*a(1,4)+a(1,1)*a(2,3)+a(1,3)*a(2,1)+a(1,2)*a(2,2),precision); s=g2T_4;
g2T_5=vpa(-a(1,6)-a(2,6)+a10*a(2,5)+a20*a(1,5)+a(1,1)*a(2,4)+a(1,4)*a(2,1)+a(1,2)*a(2,3)+a(1,3)*a(2,2),precision);
g2T_6=vpa(-a(1,7)-a(2,7)+a10*a(2,6)+a20*a(1,6)+a(1,1)*a(2,5)+a(1,5)*a(2,1)+a(1,2)*a(2,4)+a(1,4)*a(2,2)+a(1,3)*a(2,3),precision);
g2T_7=vpa(-a(1,8)-a(2,8)+a10*a(2,7)+a20*a(1,7)+a(1,1)*a(2,6)+a(1,6)*a(2,1)+a(1,2)*a(2,5)+a(1,5)*a(2,2)+a(1,3)*a(2,4)+a(1,4)*a(2,3),precision);
g2T_8=vpa(-a(1,9)-a(2,9)+a10*a(2,8)+a20*a(1,8)+a(1,1)*a(2,7)+a(1,7)*a(2,1)+a(1,2)*a(2,6)+a(1,6)*a(2,2)+a(1,3)*a(2,5)+a(1,5)*a(2,3)+a(1,4)*a(2,4),precision);
g2T_9=vpa(-a(1,10)-a(2,10)+a10*a(2,9)+a20*a(1,9)+a(1,1)*a(2,8)+a(1,8)*a(2,1)+a(1,2)*a(2,7)+a(1,7)*a(2,2)+a(1,3)*a(2,6)+a(1,6)*a(2,3)+a(1,4)*a(2,5)+a(1,5)*a(2,4),precision);
g2T_10=vpa(-a(1,11)-a(2,11)+a10*a(2,10)+a20*a(1,10)+a(1,1)*a(2,9)+a(1,9)*a(2,1)+a(1,2)*a(2,8)+a(1,8)*a(2,2)+a(1,3)*a(2,7)+a(1,7)*a(2,3)+a(1,4)*a(2,6)+a(1,6)*a(2,4)+a(1,5)*a(2,5),precision);
g2T_11=vpa(-a(1,12)-a(2,12)+a10*a(2,11)+a20*a(1,11)+a(1,1)*a(2,10)+a(1,10)*a(2,1)+a(1,2)*a(2,9)+a(1,9)*a(2,2)+a(1,3)*a(2,8)+a(1,8)*a(2,3)+a(1,4)*a(2,7)+a(1,7)*a(2,4)+a(1,5)*a(2,6)+a(1,6)*a(2,5),precision);
g2T_12=vpa(-a(1,13)-a(2,13)+a10*a(2,12)+a20*a(1,12)+a(1,1)*a(2,11)+a(1,11)*a(2,1)+a(1,2)*a(2,10)+a(1,10)*a(2,2)+a(1,3)*a(2,9)+a(1,9)*a(2,3)+a(1,4)*a(2,8)+a(1,8)*a(2,4)+a(1,5)*a(2,7)+a(1,7)*a(2,5)+a(1,6)*a(2,6),precision);
g2T_13=vpa(-a(1,14)-a(2,14)+a10*a(2,13)+a20*a(1,13)+a(1,1)*a(2,12)+a(1,12)*a(2,1)+a(1,2)*a(2,11)+a(1,11)*a(2,2)+a(1,3)*a(2,10)+a(1,10)*a(2,3)+a(1,4)*a(2,9)+a(1,9)*a(2,4)+a(1,5)*a(2,8)+a(1,8)*a(2,5)+a(1,6)*a(2,7)+a(1,7)*a(2,6),precision);
g2T_14=vpa(-a(1,15)-a(2,15)+a10*a(2,14)+a20*a(1,14)+a(1,1)*a(2,13)+a(1,13)*a(2,1)+a(1,2)*a(2,12)+a(1,12)*a(2,2)+a(1,3)*a(2,11)+a(1,11)*a(2,3)+a(1,4)*a(2,10)+a(1,10)*a(2,4)+a(1,5)*a(2,9)+a(1,9)*a(2,5)+a(1,6)*a(2,8)+a(1,8)*a(2,6)+a(1,7)*a(2,7),precision);
GT=vpa([g2Tp1 g2T g2T_1 g2T_2 g2T_3 g2T_4 g2T_5 g2T_6 g2T_7 g2T_8 g2T_9 g2T_10 g2T_11 g2T_12 g2T_13 g2T_14],precision);

f2Tp1=vpa(b1+b2+d*g2Tp1,precision); u=f2Tp1;
f2T=vpa(-b1*a20-b2*a10+d*g2T,precision); v=f2T;
f2T_1=vpa(-b1*a(2,1)-b2*a(1,1)+d*g2T_1,precision); w=f2T_1;
f2T_2=vpa(-b1*a(2,2)-b2*a(1,2)+d*g2T_2,precision); w0=f2T_2;
f2T_3=vpa(-b1*a(2,3)-b2*a(1,3)+d*g2T_3,precision);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% This part of code developes the 8th degree equation

coef1(1)=vpa(-12*d^2*m^2+24*d*m*u-12*u^2,precision);
coef1(2)=vpa(-56*d^2*m^3-12*d^2*m^2-8*d^2*m*n+112*d*m^2*u+24*d*m*u+8*d*m*v+8*d*n*u-56*m*u^2-12*u^2-8*u*v,precision);
coef1(3)=vpa(-97*d^2*m^4-34*d^2*m^3-72*d^2*m^2*n+194*d*m^3*u-16*d^2*m*n+24*d^2*m*p+4*d^2*n^2+68*d*m^2*u+48*d*m^2*v+96*d*m*n*u-97*m^2*u^2+16*d*m*v-24*d*m*w+16*d*n*u-8*d*n*v-24*d*p*u-34*m*u^2-48*m*u*v-24*n*u^2-16*u*v+24*u*w+4*v^2);
coef1(4)=vpa(-74*d^2*m^5-26*d^2*m^4-190*d^2*m^3*n+148*d*m^4*u+4*d^2*m^3-84*d^2*m^2*n+106*d^2*m^2*p-8*d^2*m*n^2+52*d*m^3*u+102*d*m^3*v+278*d*m^2*n*u-74*m^3*u^2-8*d^2*m*n+24*d^2*m*p-4*d^2*n^2+16*d^2*n*p-8*d*m^2*u+60*d*m^2*v-88*d*m^2*w+108*d*m*n*u-124*d*m*p*u+16*d*n^2*u-26*m^2*u^2-102*m^2*u*v-88*m*n*u^2+8*d*m*v-24*d*m*w+8*d*n*u+8*d*n*v-16*d*n*w-24*d*p*u-16*d*p*v+4*m*u^2-60*m*u*v+88*m*u*w+8*m*v^2-24*n*u^2-16*n*u*v+18*p*u^2-8*u*v+24*u*w-4*v^2+16*v*w);
coef1(5)=vpa(-21*d^2*m^6+2*d^2*m^5-198*d^2*m^4*n+42*d*m^5*u+8*d^2*m^4-110*d^2*m^3*n+147*d^2*m^3*p-97*d^2*m^2*n^2-4*d*m^4*u+92*d*m^4*v+304*d*m^3*n*u-21*m^4*u^2-12*d^2*m^2*n+44*d^2*m^2*p+12*d^2*m^2*q-58*d^2*m*n^2+108*d^2*m*n*p+8*d^2*n^3-16*d*m^3*u+66*d*m^3*v-114*d*m^3*w+154*d*m^2*n*u+90*d*m^2*n*v-180*d*m^2*p*u+104*d*m*n^2*u+2*m^3*u^2-92*m^3*u*v-106*m^2*n*u^2-8*d^2*n^2+32*d^2*n*p-24*d^2*n*q+12*d*m^2*v-44*d*m^2*w-12*d*m^2*w0+12*d*m*n*u+84*d*m*n*v-84*d*m*n*w-44*d*m*p*u-84*d*m*p*v-12*d*m*q*u+32*d*n^2*u-16*d*n^2*v-48*d*n*p*u+8*m^2*u^2-66*m^2*u*v+114*m^2*u*w-5*m^2*v^2-44*m*n*u^2-80*m*n*u*v+33*m*p*u^2-12*n^2*u^2+16*d*n*v-32*d*n*w+24*d*n*w0-32*d*p*v+24*d*q*v-12*m*u*v+44*m*u*w+12*m*u*w0-26*m*v^2+60*m*v*w-32*n*u*v+24*n*u*w+8*n*v^2+24*p*u*v-8*v^2+32*v*w-24*v*w0,precision);
coef1(6)=vpa(6*d^2*m^6-72*d^2*m^5*n+4*d^2*m^5-24*d^2*m^4*n+68*d^2*m^4*p-174*d^2*m^3*n^2-12*d*m^5*u+30*d*m^5*v+114*d*m^4*n*u+8*d^2*m^3*n+36*d^2*m^3*q-140*d^2*m^2*n^2+228*d^2*m^2*n*p-8*d*m^4*u+16*d*m^4*v-62*d*m^4*w+32*d*m^3*n*u+160*d*m^3*n*v-74*d*m^3*p*u+188*d*m^2*n^2*u+6*m^4*u^2-30*m^4*u*v-42*m^3*n*u^2-8*d^2*m^2*p+12*d^2*m^2*q-32*d^2*m*n^2+132*d^2*m*n*p-72*d^2*m*n*q-36*d^2*m*p^2-8*d^2*n^3+22*d^2*n^2*p-12*d*m^3*w-36*d*m^3*w0-16*d*m^2*n*u+168*d*m^2*n*v-144*d*m^2*n*w+12*d*m^2*p*u-146*d*m^2*p*v-36*d*m^2*q*u+112*d*m*n^2*u-8*d*m*n^2*v-166*d*m*n*p*u+8*d*n^3*u+4*m^3*u^2-16*m^3*u*v+62*m^3*u*w-18*m^3*v^2-8*m^2*n*u^2-124*m^2*n*u*v+6*m^2*p*u^2-32*m*n^2*u^2+16*d^2*n*p-24*d^2*n*q+8*d*m^2*w-12*d*m^2*w0+48*d*m*n*v-108*d*m*n*w+72*d*m*n*w0+8*d*m*p*u-108*d*m*p*v+36*d*m*p*w-12*d*m*q*u+72*d*m*q*v+16*d*n^2*u+16*d*n^2*v-16*d*n^2*w-48*d*n*p*u-28*d*n*p*v+36*d*p^2*u+12*m^2*u*w+36*m^2*u*w0-40*m^2*v^2+80*m^2*v*w+8*m*n*u^2-88*m*n*u*v+64*m*n*u*w+8*m*n*v^2-12*m*p*u^2+66*m*p*u*v-12*n^2*u^2-8*n^2*u*v+18*n*p*u^2-16*d*n*w+24*d*n*w0-16*d*p*v+24*d*q*v-8*m*u*w+12*m*u*w0-16*m*v^2+84*m*v*w-72*m*v*w0-16*n*u*v+24*n*u*w-8*n*v^2+16*n*v*w+24*p*u*v-36*p*u*w+6*p*v^2+16*v*w-24*v*w0,precision);
coef1(7)=vpa(18*d^2*m^5*n+3*d^2*m^5*p-90*d^2*m^4*n^2+12*d^2*m^4*n-26*d^2*m^4*p+33*d^2*m^4*q-72*d^2*m^3*n^2+143*d^2*m^3*n*p-50*d^2*m^2*n^3-6*d*m^5*v-12*d*m^5*w-30*d*m^4*n*u+78*d*m^4*n*v+6*d*m^4*p*u+102*d*m^3*n^2*u-12*d^2*m^3*p+18*d^2*m^3*q-16*d^2*m^2*n^2+102*d^2*m^2*n*p-54*d^2*m^2*n*q-42*d^2*m^2*p^2-60*d^2*m*n^3+101*d^2*m*n^2*p+4*d^2*n^4-4*d*m^4*v+14*d*m^4*w-33*d*m^4*w0-20*d*m^3*n*u+80*d*m^3*n*v-95*d*m^3*n*w+38*d*m^3*p*u-78*d*m^3*p*v-33*d*m^3*q*u+64*d*m^2*n^2*u+68*d*m^2*n^2*v-113*d*m^2*n*p*u+32*d*m*n^3*u+6*m^4*u*v+12*m^4*u*w-9*m^4*v^2+12*m^3*n*u^2-60*m^3*n*u*v-9*m^3*p*u^2-21*m^2*n^2*u^2+24*d^2*m*n*p-36*d^2*m*n*q-16*d^2*n^3+56*d^2*n^2*p-24*d^2*n^2*q-24*d^2*n*p^2+12*d*m^3*w-18*d*m^3*w0+24*d*m^2*n*v-90*d*m^2*n*w+54*d*m^2*n*w0+12*d*m^2*p*u-82*d*m^2*p*v+42*d*m^2*p*w-18*d*m^2*q*u+66*d*m^2*q*v+8*d*m*n^2*u+104*d*m*n^2*v-44*d*m*n^2*w-32*d*m*n*p*u-134*d*m*n*p*v-12*d*m*n*q*u+42*d*m*p^2*u+16*d*n^3*u-8*d*n^3*v-24*d*n^2*p*u+4*m^3*u*v-14*m^3*u*w+33*m^3*u*w0-18*m^3*v^2+45*m^3*v*w+8*m^2*n*u^2-44*m^2*n*u*v+50*m^2*n*u*w-18*m^2*n*v^2-12*m^2*p*u^2+33*m^2*p*u*v-10*m*n^2*u^2-32*m*n^2*u*v+15*m*n*p*u^2-24*d*m*n*w+36*d*m*n*w0-24*d*m*p*v+36*d*m*q*v+32*d*n^2*v-32*d*n^2*w+24*d*n^2*w0-80*d*n*p*v+24*d*n*p*w+24*d*n*q*v+24*d*p^2*v-12*m^2*u*w+18*m^2*u*w0-8*m^2*v^2+70*m^2*v*w-66*m^2*v*w0-8*m*n*u*v+20*m*n*u*w+12*m*n*u*w0-44*m*n*v^2+44*m*n*v*w+12*m*p*u*v-42*m*p*u*w+33*m*p*v^2-16*n^2*u*v+4*n^2*v^2+24*n*p*u*v+24*m*v*w-36*m*v*w0-16*n*v^2+32*n*v*w-24*n*v*w0+24*p*v^2-24*p*v*w,precision);
coef1(8)=vpa(-6*d^2*m^5*p+9*d^2*m^5*q+18*d^2*m^4*n^2+3*d^2*m^4*n*p-48*d^2*m^3*n^3-4*d^2*m^4*p+6*d^2*m^4*q+12*d^2*m^3*n^2-14*d^2*m^3*n*p+6*d^2*m^3*n*q-3*d^2*m^3*p^2-64*d^2*m^2*n^3+105*d^2*m^2*n^2*p+6*d*m^5*w-9*d*m^5*w0-12*d*m^4*n*v-21*d*m^4*n*w+6*d*m^4*p*u-9*d*m^4*q*u-24*d*m^3*n^2*u+66*d*m^3*n^2*v+15*d*m^3*n*p*u+30*d*m^2*n^3*u+12*d^2*m^2*p^2-18*d^2*m^2*p*q-24*d^2*m*n^3+100*d^2*m*n^2*p-48*d^2*m*n^2*q-48*d^2*m*n*p^2-4*d^2*n^4+6*d^2*n^3*p+4*d*m^4*w-6*d*m^4*w0-8*d*m^3*n*v-10*d*m^3*n*w-6*d*m^3*n*w0+4*d*m^3*p*u-6*d*m^3*p*v+12*d*m^3*p*w-6*d*m^3*q*u+18*d*m^3*q*v-16*d*m^2*n^2*u+100*d*m^2*n^2*v-36*d*m^2*n^2*w+44*d*m^2*n*p*u-132*d*m^2*n*p*v-24*d*m^2*n*q*u-6*d*m^2*p^2*u+28*d*m*n^3*u-42*d*m*n^2*p*u-6*m^4*u*w+9*m^4*u*w0+9*m^4*v*w+12*m^3*n*u*v+12*m^3*n*u*w-18*m^3*n*v^2-9*m^3*p*u*v+6*m^2*n^2*u^2-30*m^2*n^2*u*v-9*m^2*n*p*u^2+16*d^2*n^2*p-24*d^2*n^2*q-24*d^2*n*p^2+36*d^2*n*p*q-8*d*m^2*p*v-12*d*m^2*p*w+18*d*m^2*p*w0+12*d*m^2*q*v+40*d*m*n^2*v-52*d*m*n^2*w+48*d*m*n^2*w0+8*d*m*n*p*u-124*d*m*n*p*v+30*d*m*n*p*w-12*d*m*n*q*u+48*d*m*n*q*v-12*d*m*p^2*u+48*d*m*p^2*v+18*d*m*p*q*u+8*d*n^3*u+8*d*n^3*v-24*d*n^2*p*u-12*d*n^2*p*v+18*d*n*p^2*u-4*m^3*u*w+6*m^3*u*w0+18*m^3*v*w-18*m^3*v*w0+8*m^2*n*u*v-8*m^2*n*u*w+24*m^2*n*u*w0-36*m^2*n*v^2+36*m^2*n*v*w-12*m^2*p*u*v-12*m^2*p*u*w+27*m^2*p*v^2+4*m*n^2*u^2-28*m*n^2*u*v-12*m*n*p*u^2+42*m*n*p*u*v+9*m*p^2*u^2-16*d*n^2*w+24*d*n^2*w0-16*d*n*p*v+24*d*n*p*w-36*d*n*p*w0+24*d*n*q*v+24*d*p^2*v-36*d*p*q*v+8*m^2*v*w-12*m^2*v*w0-8*m*n*u*w+12*m*n*u*w0-16*m*n*v^2+52*m*n*v*w-48*m*n*v*w0+12*m*p*u*w-18*m*p*u*w0+24*m*p*v^2-30*m*p*v*w-8*n^2*u*v-4*n^2*v^2+24*n*p*u*v+6*n*p*v^2-18*p^2*u*v+16*n*v*w-24*n*v*w0-24*p*v*w+36*p*v*w0,precision);
coef1(9)=vpa(-6*d^2*m^4*n*p+9*d^2*m^4*n*q+6*d^2*m^3*n^3-9*d^2*m^2*n^4+6*d^2*m^3*n*q+6*d^2*m^3*p^2+4*d^2*m^2*n^3-18*d^2*m^2*n^2*q-18*d^2*m*n^4+6*d*m^4*n*w-9*d*m^4*n*w0-6*d*m^3*n^2*v-9*d*m^3*n^2*w+6*d*m^3*n*p*u-9*d*m^3*n*q*u-6*d*m^2*n^3*u+18*d*m^2*n^3*v+9*d*m^2*n^2*p*u-12*d^2*m*n^2*q+18*d^2*m*n*p*q-8*d^2*n^4+24*d^2*n^3*p-18*d^2*n^2*p^2+4*d*m^3*n*w-6*d*m^3*n*w0-4*d*m^2*n^2*v-18*d*m^2*n^2*w+18*d*m^2*n^2*w0+4*d*m^2*n*p*u-6*d*m^2*n*p*v+9*d*m^2*n*p*w-6*d*m^2*n*q*u+18*d*m^2*n*q*v-6*d*m^2*p^2*u+9*d*m^2*p*q*u-4*d*m*n^3*u+36*d*m*n^3*v+12*d*m*n^2*p*u-54*d*m*n^2*p*v-9*d*m*n*p^2*u-6*m^3*n*u*w+9*m^3*n*u*w0+9*m^3*n*v*w+6*m^2*n^2*u*v-9*m^2*n^2*v^2-9*m^2*n*p*u*v-8*d*m*n^2*w+12*d*m*n^2*w0-8*d*m*n*p*v+12*d*m*n*p*w-18*d*m*n*p*w0+12*d*m*n*q*v-18*d*m*p*q*v+16*d*n^3*v-4*m^2*n*u*w+6*m^2*n*u*w0+18*m^2*n*v*w-18*m^2*n*v*w0+6*m^2*p*u*w-9*m^2*p*u*w0-9*m^2*p*v*w+4*m*n^2*u*v-18*m*n^2*v^2-12*m*n*p*u*v+27*m*n*p*v^2+9*m*p^2*u*v+8*m*n*v*w-12*m*n*v*w0-12*m*p*v*w+18*m*p*v*w0-8*n^2*v^2+24*n*p*v^2-18*p^2*v^2-4*d^2*m^3*n*p-9*d^2*m^3*p*q+6*d^2*m^2*n^2*p+27*d^2*m*n^3*p+8*d^2*m*n^2*p-12*d^2*m*n*p^2-6*d*m^3*p*w+9*d*m^3*p*w0+12*d*m*p^2*v-48*d*n^2*p*v+36*d*n*p^2*v,precision);


coef1=vpa(coef1./coef1(1),precision);
Candid1=[];
Candid1=vpa(roots(coef1),precision); %all possible soultion for alfa2 is stored in Candid.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% This part of code developes the equation for special case of alfa2.


coef2(1)=vpa(-2,precision);
coef2(2)=vpa(m+4,precision);
coef2(3)=vpa(3*m^2-2*n+4*m,precision);
coef2(4)=vpa(3*m*n-6*p+4*n,precision);

Candid2=[];
Candid2=vpa(roots(coef2),precision);

%% This part eliminates the unacceptable solutions. 


SOL1=[]; SOL2=[]; FinalSOL1=[]; FinalSOL2=[]; FinalSOL=[];erer1=[]; erer2=[];
k0=1;k1=1;k2=1;k3=1;

for k=1:length(Candid1)
    if abs(imag(Candid1(k)))<10^-3 && real(Candid1(k))>0 && real(Candid1(k))<1 %the aalfa2 should be real and in interval (0,1).
        y1=vpa(real(Candid1(k)),precision);
        candid1alpha1=vpa(-(6*d*m^3*y1^2+11*d*m^2*y1^3+5*d*m*y1^4-2*d*m^3*y1+9*d*m^2*n*y1-2*d*m^2*y1^2+7*d*m*n*y1^2-d*n*y1^3-6*m^2*u*y1^2-11*m*u*y1^3-5*u*y1^4-2*d*m^2*n+3*d*m*n^2+4*d*m*n*y1-11*d*m*p*y1-d*n^2*y1+4*d*n*y1^2-10*d*p*y1^2+2*m^2*u*y1-3*m^2*v*y1-6*m*n*u*y1+2*m*u*y1^2-2*m*v*y1^2-5*n*u*y1^2+v*y1^3+4*d*n^2-6*d*n*p+3*m^2*w+2*m*n*u-3*m*n*v-3*m*p*u-4*m*v*y1+11*m*w*y1+n*v*y1-4*v*y1^2+10*w*y1^2-4*n*v+6*p*v)/(d*m^3*y1+2*d*m^2*y1^2+d*m*y1^3+d*m^2*n-d*n*y1^2-m^2*u*y1-2*m*u*y1^2-u*y1^3-d*m*p-d*n^2-2*d*p*y1-m*n*u+m*v*y1-n*u*y1+v*y1^2+m*w+n*v+2*w*y1),precision);
        %the candid1alpha1 is the aalfa1
        if candid1alpha1>0 && candid1alpha1<=1 %%aalfa1 should be in range of (0,1)
        SOL1(k0, :)=vpa([candid1alpha1 real(Candid1(k))],precision);
        %SOL1 is the matrix of {aalfa1,aalfa2} sets,all lie in range (0,1)
        %Maxerror function is call to find the error.
            if Maxerror(SOL1(k0, :),GT,d,f2Tp1,f2T,precision)<10^-10
            erer1(k1)=Maxerror(SOL1(k0, :),GT,d,f2Tp1,f2T,precision);
            FinalSOL1(k1,:)=vpa(SOL1(k0, :),precision);
            k1=k1+1;
        end
        k0=k0+1;
        end
    end  
end

%% %% This part eliminates the unacceptable solutions for special case of aalfa1=aalfa2. 

for k=1:length(Candid2)
    if abs(imag(Candid2(k)))<10^-3 && real(Candid2(k))>0 && real(Candid2(k))<1
       SOL2(k2,:)=vpa([real(Candid2(k)) Candid2(k)],precision);
       if Maxerror(SOL2(k2, :),GT,d,f2Tp1,f2T,precision)<10^-10
           erer2(k3)=Maxerror(SOL2(k2, :),GT,d,f2Tp1,f2T,precision);
           FinalSOL2(k3,:)=vpa(SOL2(k2, :),precision);
            k3=k3+1;
      end
        k2=k2+1;
    end
end

%% The solutions are gathered and number of them is assigned as output

FinalSOL=[FinalSOL1;FinalSOL2];
%similar sOL elimination

NanswerNet=length(find(diff(sort(FinalSOL(:,1)))>0.01))+1;

if k1==1 && k3==1
    Nanswer=0;
else
    Nanswer=NanswerNet;
end
