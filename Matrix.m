% MACM 316 - HW 2
% Gaussian Elimination for a random tridiagonal matrix
% Description: Computes the mean solution error over Ntr trials for the
% system Ax=b where A is a random sparse tridiagonal NxN matrix and
% x is a vector of ones. Plots the mean solution error for itermax values of N.
% Instructor: Sarah Huber
% File name: TridiagRandomSolve.m

clear all

Nsizes=10; %Number of different sizes of N
Ntr=2000; %Number of trials for each size of N

time=zeros(Nsizes,1);
N_it=zeros(Nsizes,1);
%Kappa=zeros(Nsizes,1);

for iter=1:Nsizes
    N=2^iter; %Matrix size
    times=zeros(Ntr,1); % Vector of errors
    x=ones(N,1); % exact solution vector
    
    for i=1:Ntr
        
        b1=randn(N,1);
b2=rand(N,1);
b3=rand(N,1);
A=spdiags([b1 b2 b3],-1:1,N,N);
        
        tic
        b=A*x; % Compute the right-hand side vector
        z=A\b; % Solve the linear system
        
        times(i)=toc; % Compute the error
    end
    
    % Compute the mean of the error
    time(iter)=mean(times);
    N_it(iter)=N;
    
    disp(['N=' num2str(N_it(iter)) '   ' 'Time=' num2str(time(iter))]);
end

%% Plot the error vs N in a loglog plot
loglog(N_it,time,'r*');
hold on
title(['Time vs N with M=', num2str(Ntr), ' trials'],'fontsize',14)
xlabel(['log_{10}(N)'],'fontsize',12)
ylabel(['log_{10}(Time)'],'fontsize',12)

%% Compute a linear regression of the data

% Polyfit fits a straight line to data of log_10(Nr) vs. log_10(Err)
% This is the data shown in the loglog plot.
p=polyfit(log10(N_it),log10(time),1);

% Extrapolate to find the point where Err(N)=1
% Fill in your work here!

grad = -p(2)/p(1)

N_estimate = 10^(-p(2)/p(1))
x_prime = [N_it(1), N_estimate];
y_prime = 10.^(polyval(p,log10(x_prime)));
loglog(x_prime, y_prime)
