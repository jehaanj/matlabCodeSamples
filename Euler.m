% MACM 316 
% Instructor: Sarah Huber
% Student: Jehaan Joseph
% File name: euler.m

clear

Ti=0; % Start time
Tf=200; % End time
e=0.5; % base of natural logarithm

% Initial conditions
q10=1-e; 
q20=0;
p10=0;
p20=sqrt((1+e)/(1-e));

% Stepsize and mesh
h=0.0005;

tt=Ti:h:Tf; % Mesh
N=length(tt);

qt1=zeros(N,1); % Values for qt1
qt2=zeros(N,1); % Values for qt2
pt1=zeros(N,1); % Values for pt1
pt2=zeros(N,1); % Values for pt2

qtA=zeros(N,1); % Values for qtA
qtB=zeros(N,1); % Values for qtB
ptA=zeros(N,1); % Values for ptA
ptB=zeros(N,1); % Values for ptB

At1=zeros(N,1); % values for angular momentum for normal Euler's
Ht1=zeros(N,1); %   "     "  Hamilton    "     "    "       "

At2=zeros(N,1); % values for angular momentum for symplectic Euler's
Ht2=zeros(N,1); %   "     "  Hamilton    "     "       "       "


qt1(1)=q10; % Initial values
qt2(1)=q20; % Initial values
pt1(1)=p10; % Initial values
pt2(1)=p20; % Initial values

qtA(1)=q10; % Initial values
qtB(1)=q20; % Initial values
ptA(1)=p10; % Initial values
ptB(1)=p20; % Initial values

% Euler steps

for i=1:N-1
    qn1 = qt1(i) + h*pt1(i);
    qn2 = qt2(i) + h*pt2(i);
    pn1 = pt1(i) - h*(1/(( (qt1(i))^2 + (qt2(i))^2 )^(3/2)))*qt1(i);
    pn2 = pt2(i) - h*(1/(( (qt1(i))^2 + (qt2(i))^2 )^(3/2)))*qt2(i);
    At1(i) = (qt1(i))*(pt2(i)) - (qt2(i))*(pt1(i));
    Ht1(i) = (1/2)*((pt1(i))^2 + (pt2(i))^2) - 1/(sqrt((qt1(i))^2 + (qt2(i))^2));
    qt1(i+1)=qn1;
    qt2(i+1)=qn2;
    pt1(i+1)=pn1;
    pt2(i+1)=pn2;
end


% Symplectic Euler steps
tic
for i=1:N-1
    qnA = qtA(i) + h*ptA(i);
    qnB = qtB(i) + h*ptB(i);
    pnA = ptA(i) - h*(1/(( (qnA)^2 + (qnB)^2 )^(3/2)))*qnA;
    pnB = ptB(i) - h*(1/(( (qnA)^2 + (qnB)^2 )^(3/2)))*qnB;
    At2(i) = (qnA)*(pnB) - (qnB)*(pnA);
    Ht2(i) = (1/2)*((pnA)^2 + (pnB)^2) - 1/(sqrt((qnA)^2 + (qnB)^2));
    qtA(i+1)=qnA;
    qtB(i+1)=qnB;
    ptA(i+1)=pnA;
    ptB(i+1)=pnB;
end
toc

%plots the graphs
figure

plot(qt1,qt2,'r')
title('Ordinary Eulers method','fontsize',14)

figure

plot(qtA,qtB,'b')
title('Symplectic Eulers method','fontsize',14)

figure

plot(tt,At1,'r')
title('Angular Momentum 1','fontsize',14)

figure

plot(tt,Ht1,'r')
title('Hamiltonian 1','fontsize',14)

figure

plot(tt,At2,'r')
title('Angular Momentum 2','fontsize',14)

figure

plot(tt,Ht2,'r')
title('Hamiltonian 2','fontsize',14)