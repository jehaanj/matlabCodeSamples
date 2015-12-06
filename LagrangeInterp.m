% MACM 316- Lagrange Polynomial Interpolation
% Description: computes the interpolating polynomial at equally-spaced
% points using the barycentric formula and plots the result
% Instructor: Sarah Huber
% Student: Jehaan Jacob Joseph
% Name: LagrangeInterp.m

clear all
clc

% Define f(x)
% f=@(x) exp(x);
%f=@(x) 1./(1+25.*x.^2);
f=@(x)((exp(3*x)).*sin(200*x.^2))./(1+20*x.^2);

% Define a grid to evaluate f(x) and P(x) on
%m=10000;
%x_int=linspace(-1,1,m);
%x_int=x_int';
%fg=f(x_int);

% Define a grid with Chebyshev nodes to evaluate f(x) and P(x) on
m=10000;
j=1:m;
x_int=cos(j*pi/m);
x_int=x_int';
fg=f(x_int);

n=100:50:1000; % Max n

error=zeros(1,length(n));

for i = 1:length(n)  
    
    x=linspace(-1,1,n(i)); % Define the set of equally-spaced interpolation points x
    
    y=f(x); % Compute the data
    x=x';
    y=y';

    % Use the barycentric form to compute P(x)
    w=baryweights(x);
    %w=[0.5 [(-1).^[1:n(i)-2]] 0.5*(-1)^n(i)-1]';
    u=baryinterp(x,w,y,x_int);

    error(i)=max(abs(u-fg));
    
    figure(i)
    
    % Plot f(x) and P(x)
    plot(x_int,fg,x_int,u);
    xlabel('x-axis','fontsize',12);
    ylabel('y-axis','fontsize',12);
    title(['Polynomial interpolation with n=' num2str(n(i)) ],'fontsize',14);
    legend({'f(x)','P(x)'},'fontsize',14,'Location','northwest');
    tol=10^-10;
    if error(i) < tol
    disp([num2str(n(i))])
    break
end
end
semilogy(n, error);
xlabel('n','fontsize',12);
ylabel('Error','fontsize',12);
title('Error for f1','fontsize',14);