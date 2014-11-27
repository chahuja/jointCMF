function [X,H,Xc] = nmfWithPhase4(Z,K,varargin)
%------------------------------------------------
%   Z=Complex input which has to be factorized to X(activation) and
%   H(weights)



%always same random seed for study
randn('seed', 0);

%Size of signals and number of NTF components
[M,N] = size(Z);


%Parsing input and initialization
p = inputParser;

p.addParamValue('nIter', 200, @(x)(isnumeric(x)&&(x==round(x))));
p.addParamValue('nIterNMF', 200, @(x)(isnumeric(x)&&(x==round(x))));
p.addParamValue('fixedX', 0, @(x)isnumeric(x)||islogical(x));
p.addParamValue('Xc', [], @(x)isnumeric(x));


p.KeepUnmatched = true;
p.parse(varargin{:})

nIter     = p.Results.nIter;
nIterNMF  = p.Results.nIterNMF;
fixedX    = p.Results.fixedX;
Xc        = p.Results.Xc;

%%
%separate real and imaginaries
ZRealFixed=real(Z);
ZImagFixed=imag(Z);

%separate positives and negatives
ZRealPos=max(ZRealFixed,0);
ZImagPos=max(ZImagFixed,0);
ZRealNeg=-min(ZRealFixed,0);
ZImagNeg=-min(ZImagFixed,0);

%%
%initialization of factors as random numbers
if isempty(Xc)
    XRealPos=.5*(1.5*abs(randn(M,K))+.5);
    XImagPos=.5*(1.5*abs(randn(M,K))+.5);
    XRealNeg=.5*(1.5*abs(randn(M,K))+.5);
    XImagNeg=.5*(1.5*abs(randn(M,K))+.5);
end

HPos=.5*(1.5*abs(randn(K,N))+.5)*mean(ZRealPos(:))/2; %check this init out
HNeg=.5*(1.5*abs(randn(K,N))+.5)*mean(ZRealNeg(:))/2; %check this init out




%Combine matrices

Z=[ZRealPos,ZRealNeg;ZImagPos,ZImagNeg]; %initialize at the centre

if(isempty(Xc))
    X=[XRealPos,XRealNeg;XImagPos,XImagNeg];
else
    X=Xc;
end

%XNeg=[XRealNeg;XImagNeg];
H=[HPos,HNeg;HNeg,HPos];
%H2=[HNeg,HPos];

%%
%main loop
%initialize 
 %ZRealPos = ZRealFixed;
 %ZImagPos = ZImagFixed;
 for i=1:nIter
     if ~fixedX
         [X,H,Zhat]=nmfStandard(Z,K,'X',X,'H',H,'nIter',nIterNMF);
     else
         [X,H,Zhat]=nmfStandard(Z,K,'X',X,'H',H,'nIter',nIterNMF,'fixedX',1);
     end
    H=[H(K+1:2*K,N+1:2*N)/2+H(1:K,1:N)/2,H(K+1:2*K,1:N)/2+H(1:K,N+1:2*N)/2;H(K+1:2*K,1:N)/2+H(1:K,N+1:2*N)/2,H(K+1:2*K,N+1:2*N)/2+H(1:K,1:N)/2];
    
 end
 
 
%  figure
%  surf(abs(H(K+1:2*K,N+1:2*N)-H(1:K,1:N)))
%  figure
%  surf(abs(H(K+1:2*K,1:N)-H(1:K,N+1:2*N)));
%  figure
%  surf(abs(H(K+1:2*K,N+1:2*N)-H(1:K,N+1:2*N)));
 Xc=X;
 X=X(1:M,1:K)-X(1:M,K+1:2*K)+1j*(X(M+1:2*M,1:K)-X(M+1:2*M,K+1:2*K));
 H=H(K+1:2*K,N+1:2*N)/2+H(1:K,1:N)/2-H(K+1:2*K,1:N)/2-H(1:K,N+1:2*N)/2;


end

