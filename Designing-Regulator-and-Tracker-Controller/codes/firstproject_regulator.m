%Regulator



%determining the stability and cotrollability of the system
clc;
clear all;

num=[0 0 1 0.6 -9.6];
den=[1 0 3.5 4 1.062];

A=[0 1 0 0;0 0 1 0;0 0 0 1;-1.062 -4 -3.5 0];
B=inv([1 0 0 0;0 1 0 0;3.5 0 1 0;4 3.5 0 1])*[0;1;0.6;-9.6];
C=[1 0 0 0];
phi_c = ctrb(A,B);

if(rank(phi_c) == length(A))
display('phi_c(A,B) is full rank so this system is controllable');
end

eigen = eig(A);
unstable_poles=eigen(real(eigen)>=0)
stable_poles=eigen(real(eigen)<0)
if isempty(unstable_poles)
    disp('this system is stable');
else 
    disp('this system is not stable');
end


%finding state feedback :
% 1-equivalency 2-Bass and Gura  3-Ackerman 4-canonical controller

%desired poles
desired_poles1 = [-1 -2 -3 -4];
desired_poles2 = [-4+0.5i -4-0.5i -8+1.5i -8-1.5i];

%equivalency
Ke1 = place(A, B, desired_poles1);
Ke2 = place(A, B, desired_poles2);

%Bass and Gura:the function of Bass_Gura is end of the code
Kbg1=Bass_Gura(A,B,desired_poles1);
Kbg2=Bass_Gura(A,B,desired_poles2);

%Ackerman
Ka1 = acker(A, B, desired_poles1);
Ka2 = acker(A, B, desired_poles2);



%canonical controller
Ac=[0 1 0 0;0 0 1 0;0 0 0 1; -1.062 -4 -3.5 0 ];
Bc=[0;0;0;1];
Cc=[-9.6 0.6 1 0];
Dc=0;
m=[0 0 0 1];
a1=m*(-1.*Ac);
s=tf('s');
delta_s_desired1=(s+1)*(s+2)*(s+3)*(s+4);
a2=[24 50 35 10];
delta_s_desired2=(s+4+0.5i)*(s+4-0.5i)*(s+8+1.5i)*(s+8-1.5i);
a3=[1077  790 210.5 24];
Kc1=a2-a1;
Kc2=a3-a1;
phi_c_Ac=ctrb(Ac,Bc);
inv_phi_c=inv(phi_c);
Kcc1=Kc1*phi_c_Ac*inv_phi_c;
Kcc2=Kc2*phi_c_Ac*inv_phi_c;


%in total we have state feedbacks

%state feedback for desierd poles 1
K1=Ka1;
%state feedback for desierd poles 2
K2=Ka2;




function k = Bass_Gura(A,B,pd)
phi_c = ctrb(A,B);
alpha = poly(pd);
alpha = alpha(1,2:end);
n = length(A);
e = eig(A);
a = poly(e);
a = a(1,2:end);
si=eye(n);
for i = 2:n
si = si+diag(a(i-1)*ones(1,n-i+1),i-1);
end
k = (alpha -a)*inv(si)*inv(phi_c);
end


