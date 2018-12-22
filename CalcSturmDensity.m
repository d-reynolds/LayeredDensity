% RELEASE NOTES
%   Written by Dylan Reynolds (reyno18@uw.edu), Feb 2018)
%
% SYNTAX
% [density] = CalcSturmDensity(depth,datenums,class)

% INPUTS
%   
% depth = time series of snowpack depth in m
% datenums = matlab datenums of observations
% class = Sturm snowpack classificaiton (see Sturm 2010)


function [density] = CalcSturmDensity(depth,datenums,class)
%This function calculates the snowpack density from the empirical method
%outlined in Sturm et al. 2010


if strcmpi(class,'alpine')
    RhoMax=.5975;
    Rho0=.2237;
    k1=.0012;
    k2=.0038;
elseif strcmpi(class,'maritime')
    RhoMax=.5979;
    Rho0=.2578;
    k1=.001;
    k2=.0038;
end

%Get DOY for bulk density calculation
DOY = GetSturmDOY(datenums);

%Bias Correction constant = 0.027222
density = (RhoMax-Rho0).*(1-min(1,exp(-k1.*depth-k2.*DOY)))+Rho0;

end

