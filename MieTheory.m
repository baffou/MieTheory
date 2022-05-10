%% SCATTERING BY A SPHERICAL GOLD NANOPARTICLE USING MIE THEORY
%% inputs
% n_m optical index of the medium
% lambda0 wavelength in nm
% r0 radius of the particle in nm
% N=5 maximum n-pole

function [Qext,Qsca,Qabs]=MieTheory(lambda0,r0,n_m)
%% parameters
n_Au=indexRead(lambda0,'Au'); %any function that returns the optical index of gold
m=n_Au/n_m;
k=2*pi*n_m/lambda0;
x=k*r0;
z=m*x;
N=round(2+x+4*x^(1/3));
%% computation
j=(1:N);

sqr=sqrt(pi*x/2);
sqrm=sqrt(pi*z/2);

phi=sqr.*besselj(j+0.5,x);
xi=sqr.*(besselj(j+0.5,x)+1i*bessely(j+0.5,x));
phim=sqrm.*besselj(j+0.5,z);
phi1=[sin(x), phi(1:N-1)];
phi1m=[sin(z), phim(1:N-1)];
y=sqr*bessely(j+0.5,x);
y1=[-cos(x), y(1:N-1)];

phip=(phi1-j/x.*phi);
phimp=(phi1m-j/z.*phim);
xip=(phi1+1i*y1)-j/x.*(phi+1i*y);
aj=(m*phim.*phip-phi.*phimp)./(m*phim.*xip-xi.*phimp);
bj=(phim.*phip-m*phi.*phimp)./(phim.*xip-m*xi.*phimp);

Qsca=sum( (2*j+1).*(abs(aj).*abs(aj)+abs(bj).*abs(bj)));
Qext=sum( (2*j+1).*real(aj+bj) );

Qext=Qext*2*pi/(k*k);
Qsca=Qsca*2*pi/(k*k);
Qabs=Qext-Qsca;
