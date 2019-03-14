%This code solves the Overlapping Control Problem in Example IV.B.
%The code requires YALMIP parser for Linear Matrix Inequality, freely avaialbe at https://yalmip.github.io. 
%Any SDP solver can be used. Here we used SDPT3 freely avaialbe at https://github.com/SQLP/SDPT3

clear all;

A=[1 4 0; 1 2 2; 0 -2 3];
B=[1 0; 0 0; 0 1]; 
[n,m]=size(B);
l1=sdpvar(1,1,'symmetric');
l2=sdpvar(1,1,'symmetric');
l3=sdpvar(1,1,'symmetric');
l4=sdpvar(1,1,'symmetric');

R=[l1 l2 0; 0 l3 l4];

Y11=sdpvar(1,1,'symmetric');
Y12=sdpvar(1,1,'symmetric');
Y22=sdpvar(1,1,'symmetric');
Y23=sdpvar(1,1,'symmetric');
Y33=sdpvar(1,1,'symmetric');

Y=[Y11, Y12, 0;
    0, Y22, 0;
    0, Y23, Y33];

P=sdpvar(n,n,'symmetric');
Condition=[zeros(n,n), P; P, zeros(n,n)]+[A*Y+B*R, (A*Y+B*R);-Y, -Y]+([A*Y+B*R, A*Y+B*R;-Y, -Y])';
problem=[Condition<=-1e-6*eye(2*n), P>=1e-6*eye(n)];
options=sdpsettings('solver','sdpt3','verbose',2);
solution=solvesdp(problem,trace(0),options);
P=double(P);
K=double(R)*inv(double(Y));