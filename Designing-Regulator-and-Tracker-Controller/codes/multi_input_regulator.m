% calculating k for multi input_regulator
clc;
clear all;

num=[0 0 1 0.6 -9.6];
den=[1 0 3.5 4 1.062];
A=[0 1 0 0;0 0 1 0;0 0 0 1;-1.062 -4 -3.5 0];
B=inv([1 0 0 0;0 1 0 0;3.5 0 1 0;4 3.5 0 1])*[0;1;0.6;-9.6];
C=[1 0 0 0];

B_B=[4.520193743116709e-17 0;1 1;0.6 0;-13.1 0];
phi_c_2=rank(ctrb(A,B_B));
if(phi_c_2 == size(A))
display('phi_c_2 is full rank so this system is controllable');
else
display('phi_c_2 is not full rank so this system is not controllable');
end
%two group desierd poles
dpoles1=[-3 -4 -5 -6];
dpoles2=[-7 -8 -9 -10];
%calculating special structures and choosing special vectors (for eaech
%group of poles 2 different structure)

%structure one for first group of poles then calculating feed back vector
e=null([A-(-3.*eye(4)) B_B]);
e1_1=e(:,1);
e1_2=e(:,2);
e1=3.*e1_1-2.*e1_2;
v1=e1(1:4,1);
q1=e1(5:6,1);

s=null([A-(-4.*eye(4)) B_B]);
e1_1=s(:,1);
e1_2=s(:,2);
e1=e1_1-2.*e1_2;
v2=e1(1:4,1);
q2=e1(5:6,1);

l=null([A-(-5.*eye(4)) B_B]);
e1_1=l(:,1);
e1_2=l(:,2);
e1=3.*e1_1-e1_2;
v3=e1(1:4,1);
q3=e1(5:6,1);

f=null([A-(-6.*eye(4)) B_B]);
f1_1=f(:,1);
f1_2=f(:,2);
f1=1.5.*f1_1-2.5.*f1_2;
v4=f1(1:4,1);
q4=f1(5:6,1);
%first feedback matrix for first special structure for first group poles
ks1_1=-1.*[q1 q2 q3 q4]*inv([v1  v2 v3 v4]);

%structure two for first group of poles then calculating feed back vector
e=null([A-(-3.*eye(4)) B_B]);
e1_1=e(:,1);
e1_2=e(:,2);
e1=1.*e1_1-3.5*e1_2;
v1=e1(1:4,1);
q1=e1(5:6,1);

s=null([A-(-4.*eye(4)) B_B]);
e1_1=s(:,1);
e1_2=s(:,2);
e1=0.3.*e1_1-2.4.*e1_2;
v2=e1(1:4,1);
q2=e1(5:6,1);

l=null([A-(-5.*eye(4)) B_B]);
e1_1=l(:,1);
e1_2=l(:,2);
e1=4.5.*e1_1-1.25.*e1_2;
v3=e1(1:4,1);
q3=e1(5:6,1);

f=null([A-(-6.*eye(4)) B_B]);
f1_1=f(:,1);
f1_2=f(:,2);
f1=2.25.*f1_1-3.5.*f1_2;
v4=f1(1:4,1);
q4=f1(5:6,1);
%second feedback matrix for second special structure for first group poles
ks1_2=-1.*[q1 q2 q3 q4]*inv([v1  v2 v3 v4]);

%structure one for second group of poles then calculating feed back vector
e=null([A-(-7.*eye(4)) B_B]);
e1_1=e(:,1);
e1_2=e(:,2);
e1=3.*e1_1-2.*e1_2;
v1=e1(1:4,1);
q1=e1(5:6,1);

s=null([A-(-8.*eye(4)) B_B]);
e1_1=s(:,1);
e1_2=s(:,2);
e1=e1_1-2.*e1_2;
v2=e1(1:4,1);
q2=e1(5:6,1);

l=null([A-(-9.*eye(4)) B_B]);
e1_1=l(:,1);
e1_2=l(:,2);
e1=3.*e1_1-e1_2;
v3=e1(1:4,1);
q3=e1(5:6,1);

f=null([A-(-10.*eye(4)) B_B]);
f1_1=f(:,1);
f1_2=f(:,2);
f1=1.5.*f1_1-2.5.*f1_2;
v4=f1(1:4,1);
q4=f1(5:6,1);
%first feedback matrix for first special structure for second group poles
ks2_1=-1.*[q1 q2 q3 q4]*inv([v1  v2 v3 v4]);

%structure two for first group of poles then calculating feed back vector
e=null([A-(-7.*eye(4)) B_B]);
e1_1=e(:,1);
e1_2=e(:,2);
e1=1.*e1_1-3.5*e1_2;
v1=e1(1:4,1);
q1=e1(5:6,1);

s=null([A-(-8.*eye(4)) B_B]);
e1_1=s(:,1);
e1_2=s(:,2);
e1=0.3.*e1_1-2.4.*e1_2;
v2=e1(1:4,1);
q2=e1(5:6,1);

l=null([A-(-9.*eye(4)) B_B]);
e1_1=l(:,1);
e1_2=l(:,2);
e1=4.5.*e1_1-1.25.*e1_2;
v3=e1(1:4,1);
q3=e1(5:6,1);

f=null([A-(-10.*eye(4)) B_B]);
f1_1=f(:,1);
f1_2=f(:,2);
f1=2.25.*f1_1-3.5.*f1_2;
v4=f1(1:4,1);
q4=f1(5:6,1);
%second feedback matrix for second special structure for second group poles
ks2_2=-1.*[q1 q2 q3 q4]*inv([v1  v2 v3 v4]);


%calculating K for these two group of poles with Lyapunov method
K_Lya1=Multi_input(A,B_B,dpoles1);
K_Lya2=Multi_input(A,B_B,dpoles2);

%Lyapunov method function
function k = Multi_input(A,B,pd)
F = diag(pd);
k_bar = [1 1 1 0;0 1 0 1];
if(rank(obsv(F,k_bar)) == length(pd))
T = lyap(A,-F,-B*k_bar);
k = k_bar*inv(T);
else
 
k = NaN;
end
end

