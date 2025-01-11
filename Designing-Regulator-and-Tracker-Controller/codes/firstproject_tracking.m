% Tracking: calculating ua : static pre-compensator
clc;
clear all;
s=tf('s');
num=[0 0 1 0.6 -9.6];
den=[1 0 3.5 4 1.062];
A=[0 1 0 0;0 0 1 0;0 0 0 1;-1.062 -4 -3.5 0];
B=inv([1 0 0 0;0 1 0 0;3.5 0 1 0;4 3.5 0 1])*[0;1;0.6;-9.6];
C=[1 0 0 0];


%state feedback for desierd poles 1
K1=[-2.070746295816953,-3.595041954368031,-2.197204600415871,-1.138424787375386];
%state feedback for desierd poles 2
K2=[-1.088511183014434e+02,-75.738608696318150,-23.532547504161734,-8.691460854871387];
%close loop for K1
Acl_1=A-B*K1;
%close loop for K2
Acl_2=A-B*K2;

%finding ua1 for first close loop system
ua1=inv((-1.*C)*inv(Acl_1)*B);

%finding ua2 for second close loop system
ua2=inv((-1.*C)*inv(Acl_2)*B);


%Tracking : Integral controller
%testing controllability of [B A;0 -C]
phi=[0 0 1 0 0;1 0 0 1 0;0.6 0 0 0 1;-13.1 -1.062 -4 -3.5 0;0 -1 0 0 0];
r=rank(phi);

if(r == 5)
display('[B A;0 -C] is full rank so this system is controllable');
display('so it is ok to designing trackig controller(integral)');
else
display('[B A;0 -C] is not full rank so this system is not controllable');

end
AI=[0 1 0 0 0;0 0 1 0 0;0 0 0 1 0;-1.062 -4 -3.5 0 0;1 0 0 0 0];
BI=[4.520193743116709e-17; 1 ;0.6; -13.1; 0];
CI=[1 0 0 0 0];

dp1 = [-1 -2 -3 -4 -2.5];
dp2 = [-4+0.5i -4-0.5i -8+1.5i -8-1.5i -2.5];
%calculating for both desired poles
KI_1=acker(AI,BI,dp1);
KI1_1=KI_1(1,1:4);
KI1_2=KI_1(1,5);;
KI_2=acker(AI,BI,dp2);
KI2_1=KI_2(1,1:4);
KI2_2=KI_2(1,5);

%A_r for testing robustness
A_r=[0 1 0 0; 0 0 1 0 ;0 0 0 1;-1.062 -4 -3.5 -5 ];
%B_r for testing robustness
B_r=[0;0.8;0.6;-8];
%B_r for testing robustness
C_r=[ 0.9  0.1 0 0];
