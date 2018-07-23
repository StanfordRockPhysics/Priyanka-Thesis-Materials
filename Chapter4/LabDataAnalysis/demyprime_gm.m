function yprime=demyprime_gm(t,y)
%function yprime=demyprime_gm(t,y)
%used by DEM

%Written by T. Mukerji
%modified by G. Mavko (moved strain concentration to function PQfactors.m

global DEMINPT;

k1=DEMINPT(1); mu1=DEMINPT(2); k2=DEMINPT(3); mu2=DEMINPT(4);
asp=DEMINPT(5); phic=DEMINPT(6);

ka=k2; mua=mu2;

k=y(1); mu=y(2);
yprime=zeros(2,1);

if asp==1.
   asp=0.99;
end
  
%
[pa, qa] = PQfactors(k, ka, mu, mua, asp);

krhs=(ka-k)*pa;
yprime(1)=krhs/(1-t);

murhs=(mua-mu)*qa;
yprime(2)=murhs/(1-t);
