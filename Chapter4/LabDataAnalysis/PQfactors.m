function [P, Q] = PQfactors(k0, k, mu0, mu, asp)
% [P, Q] = PQfactors(k0, k, mu0, mu, asp)
%
% Compute the P and Q factors for spheroidal inclusions in an elastic
% background.  Can have multiple
% Inputs
%     k0, mu0       bulk and shear moduli of background material
%     k, mu           vectors of bulk and shear moduli of inclusion materials
%     asp              vector of aspect ratios of inclusions

nu0=(3*k0-2*mu0)/(2*(3*k0+mu0));

k=k(:); mu=mu(:); asp=asp(:);
indx=find(asp==1); asp(indx)=0.99*ones(size(indx));
theta=zeros(size(asp)); fn=zeros(size(asp));  

obdx=find(asp<1);
theta(obdx)=(asp(obdx)./((1-asp(obdx).^2).^(3/2))).*...
             (acos(asp(obdx)) -asp(obdx).*sqrt(1-asp(obdx).^2));
fn(obdx)=(asp(obdx).^2./(1-asp(obdx).^2)).*(3.*theta(obdx) -2);
prdx=find(asp>1);
theta(prdx)=(asp(prdx)./((asp(prdx).^2-1).^(3/2))).*...
             (asp(prdx).*sqrt(asp(prdx).^2-1)-acosh(asp(prdx)));
fn(prdx)=(asp(prdx).^2./(asp(prdx).^2-1)).*(2-3.*theta(prdx));

a=mu./mu0 -1; 
b=(1/3)*(k./k0 -mu./mu0); 
r=(1-2*nu0)/(2*(1-nu0));

f1=1+a.*((3/2).*(fn+theta)-r.*((3/2).*fn+(5/2).*theta-(4/3)));
f2=1+a.*(1+(3/2).*(fn+theta)-(r/2).*(3.*fn+5.*theta))+b.*(3-4*r);
f2=f2+(a/2).*(a+3.*b).*(3-4.*r).*(fn+theta-r.*(fn-theta+2.*theta.^2));
f3=1+a.*(1-(fn+(3/2).*theta)+r.*(fn+theta));
f4=1+(a./4).*(fn+3.*theta-r.*(fn-theta));
f5=a.*(-fn+r.*(fn+theta-(4/3))) + b.*theta.*(3-4*r);
f6=1+a.*(1+fn-r.*(fn+theta))+b.*(1-theta).*(3-4.*r);
f7=2+(a./4).*(3.*fn+9.*theta-r.*(3.*fn+5.*theta)) + b.*theta.*(3-4.*r);
f8=a.*(1-2.*r+(fn./2).*(r-1)+(theta./2).*(5.*r-3))+b.*(1-theta).*(3-4.*r);
f9=a.*((r-1).*fn-r.*theta) + b.*theta.*(3-4.*r);

P = 3*f1./f2; 
Q = (2./f3) + (1./f4) +((f4.*f5 + f6.*f7 - f8.*f9)./(f2.*f4));
P=P./3;
Q=Q./5 ;

