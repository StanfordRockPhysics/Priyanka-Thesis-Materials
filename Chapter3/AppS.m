function [c33, c44, eps, gam, del, eta] = AppS( c, ar )

% Approximate Reuss summation using Legendre polynomials

e = (c(1,1)-c(3,3))/(2*c(3,3));
g = (c(6,6)-c(4,4))/(2*c(4,4));
d = ((c(1,3)+c(4,4))^2 - (c(3,3)-c(4,4))^2)/(2*c(3,3)*(c(3,3)-c(4,4)));
s = inv(c);

for n = 1:length(ar)
    siso = zeros(6,6);
    
    siso(1,1) = (1/15)*(8*s(1,1)+3*s(3,3)+4*s(1,3)+2*s(4,4));
    siso(2,2) = siso(1,1); siso(3,3) = siso(1,1);
    
    siso(1,2) = (1/15)*(s(1,1)+s(3,3)+5*s(1,2)+8*s(1,3)-s(4,4));
    siso(2,3) = siso(1,2); siso(3,1) = siso(1,2);
    siso(2,1) = siso(1,2);siso(1,3) = siso(3,1);siso(3,2) = siso(2,3);
    
    siso(4,4) = (2/15)*(7*s(1,1)+2*s(3,3)-5*s(1,2)-4*s(1,3)+3*s(4,4));
    siso(5,5) = siso(4,4); siso(6,6) = siso(4,4);
    
    
    
    s1 = s(1,1)+s(3,3)-2*s(1,3)-s(4,4);
    s2 = s(1,1)-3*s(1,2)+2*s(1,3)-s(4,4)/2;
    s3 = 4*s(1,1)-3*s(3,3)-s(1,3)-s(4,4)/2;
    
    
    
    dels = zeros(6,6); a = ar(n);
    
    dels(1,1) = (8*sqrt(10)/105)*(pi.^2)*s3*(w(2,0,0,a) - sqrt(6)*w(2,2,0,a)) + (4*sqrt(2)/35)*(pi.^2)*s1*(w(4,0,0,a) - (2*sqrt(10)/3)*w(4,2,0,a) + (sqrt(70)/3)*w(4,4,0,a)); 
    
    dels(2,2) = (8*sqrt(10)/105)*(pi.^2)*s3*(w(2,0,0,a) + sqrt(6)*w(2,2,0,a)) + (4*sqrt(2)/35)*(pi.^2)*s1*(w(4,0,0,a) + (2*sqrt(10)/3)*w(4,2,0,a) + (sqrt(70)/3)*w(4,4,0,a)); 
    
    dels(3,3) = (16*sqrt(2)/105)*(pi.^2)*(2*s1*w(4,0,0,a) - sqrt(5)*s3*w(2,0,0,a));
    
    dels(1,2) = -(8*sqrt(10)/315)*(pi.^2)*(7*s2 - s3)*w(2,0,0,a) + (4*sqrt(2)/105)*(pi.^2)*s1*(w(4,0,0,a) - sqrt(70)*w(4,4,0,a)); dels(2,1) = dels(1,2);
    
    dels(1,3) = (4*sqrt(10)/315)*(pi.^2)*(7*s2 - s3)*(w(2,0,0,a) + sqrt(6)*w(2,2,0,a)) - (16*sqrt(2)/105)*(pi.^2)*s1*(w(4,0,0,a) - sqrt(5/2)*w(4,2,0,a)); dels(3,1) = dels(1,3);
    
    dels(2,3) = (4*sqrt(10)/315)*(pi.^2)*(7*s2 - s3)*(w(2,0,0,a) - sqrt(6)*w(2,2,0,a)) - (16*sqrt(2)/105)*(pi.^2)*s1*(w(4,0,0,a) + sqrt(5/2)*w(4,2,0,a)); dels(3,2) = dels(2,3);
    
    dels(4,4) = -(8*sqrt(10)/315)*(pi.^2)*(7*s2 + 2*s3)*(w(2,0,0,a) - sqrt(6)*w(2,2,0,a)) - (64*sqrt(2)/105)*(pi.^2)*s1*(w(4,0,0,a) + sqrt(5/2)*w(4,2,0,a));
    
    dels(5,5) = -(8*sqrt(10)/315)*(pi.^2)*(7*s2 + 2*s3)*(w(2,0,0,a) + sqrt(6)*w(2,2,0,a)) - (64*sqrt(2)/105)*(pi.^2)*s1*(w(4,0,0,a) - sqrt(5/2)*w(4,2,0,a));
    
    dels(6,6) = (16*sqrt(10)/315)*(pi.^2)*(7*s2 + 2*s3)*w(2,0,0,a) + (16*sqrt(2)/105)*(pi.^2)*s1*(w(4,0,0,a) - sqrt(70)*w(4,4,0,a));
    
    svti = siso + dels;
    cvti = inv(svti);
    
    c33(n) = cvti(3,3);
    c44(n) = cvti(4,4);
    eps(n) = (cvti(1,1)-cvti(3,3))/(2*cvti(3,3));
    gam(n) = (cvti(6,6)-cvti(4,4))/(2*cvti(4,4));
    del(n) = ((cvti(1,3)+cvti(4,4))^2 - (cvti(3,3)-cvti(4,4))^2)/(2*cvti(3,3)*(cvti(3,3)-cvti(4,4)));
    eta(n) = (eps(n) - del(n))./(1+2*del(n));
end

end

