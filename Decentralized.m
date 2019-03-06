%This code solves the Decentralized Control Problem in Example IV.b.
%The code requires YALMIP parser for Linear Matrix Inequality, freely avaialbe at https://yalmip.github.io. 
%Any SDP solver can be used. Here we used SDPT3 freely avaialbe at https://github.com/SQLP/SDPT3

clear all;

A=[0 0 0; 0 0.1 0; 0 0 0];

S=[1 0 1; 0 1 1];
B=S';
n=max(size(A));
m=min(size(B));
r1=sdpvar(1,1,'symmetric');
r2=sdpvar(1,1,'symmetric');
r3=sdpvar(1,1,'symmetric');
r4=sdpvar(1,1,'symmetric');
r5=sdpvar(1,1,'symmetric');

R=[r1 0 r2; 0 r3 r4];

Y11=sdpvar(1,1,'symmetric');

Y13=sdpvar(1,1,'symmetric');

Y22=sdpvar(1,1,'symmetric');

Y23=sdpvar(1,1,'symmetric');

Y33=sdpvar(1,1,'symmetric');

Y=[Y11 0 Y13; 0 Y22 Y23; 0 0 Y33];

P=sdpvar(3,3,'symmetric');

M=[zeros(n,n), P; P, zeros(n,n)]+[(A*Y+B*R), (A*Y+B*R);-Y, -Y]...
    +([A*Y+B*R, (A*Y+B*R);-Y, -Y])';


problem=[M<=-1e-8*eye(2*n), P>=1e-5*eye(n)];
   
options=sdpsettings('solver','sdpt3','verbose',2);

solution=solvesdp(problem,trace(0),options);

P=double(P);

K=double(R)*inv(double(Y));

