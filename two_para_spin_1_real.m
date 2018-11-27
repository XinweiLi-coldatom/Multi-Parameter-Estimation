clear
clc
%% definitions
% best rotation angle
theta0 = 0.277*pi; 
% mesh
delta = 0.001*pi;
[x,y] = meshgrid(0:delta:pi);
nsample = 10000; % atomic number of each experiment

nexp = 500; % number of experiments
sigma = 14; % detection noise
var = 0.005; % variance of rotation angle
% real value of parameters
x0 = 0.45*pi;
y0 = 0.35*pi;
% probabilities
p1 = 0.25*sin(theta0)^2*(1+3*cos(theta0)^2-2*cos(x).*cos(theta0)*(1+cos(theta0))...
    +4*cos(y).*cos(theta0)*sin(theta0/2)^2-cos(x-y).*sin(theta0)^2);
p2 = cos(theta0)^4+(cos(x)+cos(y)).*cos(theta0)^2*sin(theta0)^2+0.5*(1+cos(x-y)).*sin(theta0)^4;
p3 = -0.125*sin(theta0)^2*(-5+8*cos(y).*cos(0.5*theta0)^2*cos(theta0)-3*cos(2*theta0)...
    -8*cos(x).*cos(theta0)*sin(theta0/2)^2+2*cos(x-y).*sin(theta0)^2);
%% simulation of experiments
xmvec = []; % estimated x
ymvec = []; % estimated y
for l = 1:nexp
    
    if mod(l,50)==0
        disp(l)
    end
    
    % probability of experiment
    theta = theta0 * (1 + normrnd(0,var));    
    p10 = 0.25*sin(theta)^2*(1+3*cos(theta)^2-2*cos(x0)*cos(theta)*(1+cos(theta))...
        +4*cos(y0)*cos(theta)*sin(theta/2)^2-cos(x0-y0)*sin(theta)^2);
    p20 = cos(theta)^4+(cos(x0)+cos(y0))*cos(theta)^2*sin(theta)^2+0.5*(1+cos(x0-y0))*sin(theta)^4;
    p30 = -0.125*sin(theta)^2*(-5+8*cos(y0)*cos(0.5*theta)^2*cos(theta)-3*cos(2*theta)...
        -8*cos(x0)*cos(theta)*sin(theta/2)^2+2*cos(x0-y0)*sin(theta)^2);
    % generate random number
    alphabet = [1,0,-1];
    prob = [p10,p20,p30];
    data = randsrc(nsample,1,[alphabet; prob]);
    % atomic number detection    
    n1 = numel(find(data==1));
    n0 = numel(find(data==0));
    nm1 = numel(find(data==-1));
    % detection noise
    n1d = round(normrnd(0,sigma));
    n0d = round(normrnd(0,sigma));
    nm1d = round(normrnd(0,sigma));
    % Logarithm of likelihood function
    pf = (n1+n1d)*log(p1)+(n0+n0d)*log(p2)+(nm1+nm1d)*log(p3);

    % maximum likelihood estimation
    [ym,xm]=find(pf==max(max(real(pf))));% estimated x and y
    xmvec = [xmvec,xm*delta];
    ymvec = [ymvec,ym*delta];
end

%%
figure(1)
h=histogram(xmvec,20);
h.FaceColor = [1 0 0];
h.EdgeColor = 'k';
figure(2)
h=histogram(ymvec,20);
h.FaceColor = [0 1 0];
h.EdgeColor = 'k';
%% variance analysis
std(xmvec)
std(ymvec)
sqrt(1.45/nsample)