% Maximal likelihood estimation of two parameters with spin-1 atom
% Total atomic number nsample = 1000
% Repeated experiment nexp = 500
% Detection noise sigma = 6
clear
clc
%% definitions
% best rotation angle
theta = 0.2774*pi; 
% mesh
delta = 0.001*pi;
[x,y] = meshgrid(0:delta:pi,0:-delta:-pi);
% probabilities
p1 = 0.25*sin(theta)^2*(1+3*cos(theta)^2-2*cos(x).*cos(theta)*(1+cos(theta))...
    +4*cos(y).*cos(theta)*sin(theta/2)^2-cos(x-y).*sin(theta)^2);
p2 = cos(theta)^4+(cos(x)+cos(y)).*cos(theta)^2*sin(theta)^2+0.5*(1+cos(x-y)).*sin(theta)^4;
p3 = -0.125*sin(theta)^2*(-5+8*cos(y).*cos(0.5*theta)^2*cos(theta)-3*cos(2*theta)...
    -8*cos(x).*cos(theta)*sin(theta/2)^2+2*cos(x-y).*sin(theta)^2);

nsample = 1000; % atomic number of each experiment
nexp = 500; % number of experiments
sigma = 6; % detection noise

% figure
% mesh(x,y,p1)
% imagesc(p1)
%% Estimation process
% real value of parameters
x0vec = 0.01*pi:0.01*pi:pi;
y0vec = -0.01*pi:-0.01*pi:-pi;
xest = zeros(length(x0vec),length(y0vec)); % estimation of x
yest = zeros(length(x0vec),length(y0vec)); % estimation of y
xstd = zeros(length(x0vec),length(y0vec)); % standard deviation of estimated x
ystd = zeros(length(x0vec),length(y0vec)); % standard deviation of estimated y

for k = 1:length(x0vec)    
    for p = 1:length(y0vec)     
        fprintf('-----k=%d,p=%d-----\n',k,p)
        % parameters to be estimated
        x0 = x0vec(k);
        y0 = y0vec(p);
        % probability of experiment
        p10 = 0.25*sin(theta)^2*(1+3*cos(theta)^2-2*cos(x0)*cos(theta)*(1+cos(theta))...
            +4*cos(y0)*cos(theta)*sin(theta/2)^2-cos(x0-y0)*sin(theta)^2);
        p20 = cos(theta)^4+(cos(x0)+cos(y0))*cos(theta)^2*sin(theta)^2+0.5*(1+cos(x0-y0))*sin(theta)^4;
        p30 = -0.125*sin(theta)^2*(-5+8*cos(y0)*cos(0.5*theta)^2*cos(theta)-3*cos(2*theta)...
            -8*cos(x0)*cos(theta)*sin(theta/2)^2+2*cos(x0-y0)*sin(theta)^2);
        % probability of outcome
        alphabet = [1,0,-1];
        prob = [p10,p20,p30];
        % simulating experiment
        xmvec = zeros(1,nexp);
        ymvec = zeros(1,nexp);
        
        for l = 1:nexp                       
            % generate random number
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
            [ym,xm]=find(real(pf)==max(max(real(pf))));% estimated x and y
            xmvec(l) = xm(1)*delta; % if multiple maximum values, just pick one of them
            ymvec(l) = ym(1)*(-delta);
            if mod(l,500)==0
                fprintf('-----x=%d,y=%d-----\n',xmvec(l),ymvec(l))
            end
        end
        % standard deviation for each phase
        xest(k,p) = mean(xmvec);
        yest(k,p) = mean(ymvec);
        xstd(k,p) = std(xmvec);
        ystd(k,p) = std(ymvec);
             
    end
end
%% save data
save xest_noise_6.dat xest -ascii
save yest_noise_6.dat yest -ascii
save xstd_noise_6.dat xstd -ascii
save ystd_noise_6.dat ystd -ascii