%%This function checks the error between original g and soultions gg
function err=Maxerror(xy,GT,d,f2Tp1,f2T,precision)
%% This function input is the {aalfa1,aalfa2}-the original transfer function coefficiet d,GT={all g} and f.
%% This section checks the {aalfa1,aalfa2}

aalpha1=vpa(xy(1,1),precision);
aalpha2=vpa(xy(1,2),precision);
Rinf=vpa(d,precision);
RR2=inf;
Ts=vpa(1/2000,precision);
T=15; %to reduce the computational task

bb1=vpa(((f2Tp1-d*GT(1))*(GT(1)+aalpha2)-(f2T-d*GT(2)))/(GT(1)+2*aalpha2),precision);
bb2=vpa(f2Tp1-d*GT(1)-bb1,precision);
%bb1 and bb2 should be positive, {aalfa1,aalfa2} with negative bb1 or bb2
%is unacceptable. the error is assigned as 1.
if bb1<0 || bb2<0
    err=1;
else
 %the GT consistency is checked.first the R1,C1,C2 associated with {aalfa1
 %and aalfa2} is obtained and then gg are computed.
CC1=vpa((Ts^aalpha1)/bb1,precision); CC2=vpa((Ts^aalpha2)/bb2,precision);
RR1=vpa(bb1/(aalpha1+aalpha2+GT(1)),precision);

a10=vpa(aalpha1-Ts^aalpha1/(RR1*CC1),precision);
a20=vpa(aalpha2-Ts^aalpha2/(RR2*CC2),precision);

for j=1:T
a(1,j)=vpa(-(-1)^(j+1)*gamma(aalpha1+1)/(gamma(j+2)*gamma(aalpha1-j)),precision);

a(2,j)=vpa(-(-1)^(j+1)*gamma(aalpha2+1)/(gamma(j+2)*gamma(aalpha2-j)),precision);
end
gg2Tp1=vpa(-a10-a20,precision); 
gg2T=vpa(-a(1,1)-a(2,1)+a10*a20,precision); 
gg2T_1=vpa(a10*a(2,1)+a(1,1)*a20-a(1,2)-a(2,2),precision); 
gg2T_2=vpa(a10*a(2,2)+a20*a(1,2)+a(1,1)*a(2,1)-a(1,3)-a(2,3),precision); 
gg2T_3=vpa(a10*a(2,3)+a20*a(1,3)+a(1,1)*a(2,2)+a(1,2)*a(2,1)-a(1,4)-a(2,4),precision); 
gg2T_4=vpa(-a(1,5)-a(2,5)+a10*a(2,4)+a20*a(1,4)+a(1,1)*a(2,3)+a(1,3)*a(2,1)+a(1,2)*a(2,2),precision); 
gg2T_5=vpa(-a(1,6)-a(2,6)+a10*a(2,5)+a20*a(1,5)+a(1,1)*a(2,4)+a(1,4)*a(2,1)+a(1,2)*a(2,3)+a(1,3)*a(2,2),precision);
gg2T_6=vpa(-a(1,7)-a(2,7)+a10*a(2,6)+a20*a(1,6)+a(1,1)*a(2,5)+a(1,5)*a(2,1)+a(1,2)*a(2,4)+a(1,4)*a(2,2)+a(1,3)*a(2,3),precision);
gg2T_7=vpa(-a(1,8)-a(2,8)+a10*a(2,7)+a20*a(1,7)+a(1,1)*a(2,6)+a(1,6)*a(2,1)+a(1,2)*a(2,5)+a(1,5)*a(2,2)+a(1,3)*a(2,4)+a(1,4)*a(2,3),precision);
gg2T_8=vpa(-a(1,9)-a(2,9)+a10*a(2,8)+a20*a(1,8)+a(1,1)*a(2,7)+a(1,7)*a(2,1)+a(1,2)*a(2,6)+a(1,6)*a(2,2)+a(1,3)*a(2,5)+a(1,5)*a(2,3)+a(1,4)*a(2,4),precision);
gg2T_9=vpa(-a(1,10)-a(2,10)+a10*a(2,9)+a20*a(1,9)+a(1,1)*a(2,8)+a(1,8)*a(2,1)+a(1,2)*a(2,7)+a(1,7)*a(2,2)+a(1,3)*a(2,6)+a(1,6)*a(2,3)+a(1,4)*a(2,5)+a(1,5)*a(2,4),precision);
gg2T_10=vpa(-a(1,11)-a(2,11)+a10*a(2,10)+a20*a(1,10)+a(1,1)*a(2,9)+a(1,9)*a(2,1)+a(1,2)*a(2,8)+a(1,8)*a(2,2)+a(1,3)*a(2,7)+a(1,7)*a(2,3)+a(1,4)*a(2,6)+a(1,6)*a(2,4)+a(1,5)*a(2,5),precision);
gg2T_11=vpa(-a(1,12)-a(2,12)+a10*a(2,11)+a20*a(1,11)+a(1,1)*a(2,10)+a(1,10)*a(2,1)+a(1,2)*a(2,9)+a(1,9)*a(2,2)+a(1,3)*a(2,8)+a(1,8)*a(2,3)+a(1,4)*a(2,7)+a(1,7)*a(2,4)+a(1,5)*a(2,6)+a(1,6)*a(2,5),precision);
gg2T_12=vpa(-a(1,13)-a(2,13)+a10*a(2,12)+a20*a(1,12)+a(1,1)*a(2,11)+a(1,11)*a(2,1)+a(1,2)*a(2,10)+a(1,10)*a(2,2)+a(1,3)*a(2,9)+a(1,9)*a(2,3)+a(1,4)*a(2,8)+a(1,8)*a(2,4)+a(1,5)*a(2,7)+a(1,7)*a(2,5)+a(1,6)*a(2,6),precision);
gg2T_13=vpa(-a(1,14)-a(2,14)+a10*a(2,13)+a20*a(1,13)+a(1,1)*a(2,12)+a(1,12)*a(2,1)+a(1,2)*a(2,11)+a(1,11)*a(2,2)+a(1,3)*a(2,10)+a(1,10)*a(2,3)+a(1,4)*a(2,9)+a(1,9)*a(2,4)+a(1,5)*a(2,8)+a(1,8)*a(2,5)+a(1,6)*a(2,7)+a(1,7)*a(2,6),precision);
gg2T_14=vpa(-a(1,15)-a(2,15)+a10*a(2,14)+a20*a(1,14)+a(1,1)*a(2,13)+a(1,13)*a(2,1)+a(1,2)*a(2,12)+a(1,12)*a(2,2)+a(1,3)*a(2,11)+a(1,11)*a(2,3)+a(1,4)*a(2,10)+a(1,10)*a(2,4)+a(1,5)*a(2,9)+a(1,9)*a(2,5)+a(1,6)*a(2,8)+a(1,8)*a(2,6)+a(1,7)*a(2,7),precision);
GTnew=vpa([gg2Tp1 gg2T gg2T_1 gg2T_2 gg2T_3 gg2T_4 gg2T_5 gg2T_6 gg2T_7 gg2T_8 gg2T_9 gg2T_10 gg2T_11 gg2T_12 gg2T_13 gg2T_14],precision);
%the maximum perunit error is set as output.
err=max(abs((GT-GTnew)./GT));



end
