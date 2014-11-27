function [X,H,Zhat]= nmfStandard(Z,K,varargin)

%always same random seed for study
randn('seed', 0);

%Size of signals and number of NTF components
[M,N] = size(Z);


%Parsing input and initialization
p = inputParser;
p.addParamValue('X', [], @(x)isnumeric(x));
p.addParamValue('H', [], @(x)isnumeric(x));
p.addParamValue('fixedX', 0, @(x)isnumeric(x)||islogical(x));
% p.addParamValue('fixedH', 0, @(x)isnumeric(x)||islogical(x));
% p.addParamValue('fixedQ', 0, @(x)isnumeric(x)||islogical(x));
p.addParamValue('nIter', 200, @(x)(isnumeric(x)&&(x==round(x))));
p.addParamValue('display', 0, @(x)isnumeric(x));

p.KeepUnmatched = true;
p.parse(varargin{:})

X         = p.Results.X;
H         = p.Results.H;
fixedX    = p.Results.fixedX;
% fsixedH    = p.Results.fixedH;
% fixedQ    = p.Results.fixedQ;
nIter     = p.Results.nIter;
display   = p.Results.display;

%euclidean distance minimization
% Avoid zero values
%-----------------------------------
Z = max(Z,eps);

% initializing empty parameters
if isempty(X), X=.5*(1.5*abs(randn(M,K))+.5); end;
if isempty(H), H=.5*(1.5*abs(randn(K,N))+.5)*mean(Z(:)); end;

% sources spectrogram model
Zhat = X*H;

for i=1:nIter
    H=H.*(X'*Z)./(X'*X*H); %euclidean distance minimization 
    
    %MSE using frobenious norm
    Zhat=X*H;
    %MSE=norm(Z-Zhat); %diplay the norm after each iteration
    
    if ~fixedX
    X=X.*(Z*H')./(X*H*H');
    
    %MSE using frobenious norm
    Zhat=X*H;
    %MSE=norm(Z-Zhat); %diplay the norm after each iteration
    end
end