function [ ksc, gsc ] = IsotropicPolycrystalSCA( cc , ar)

% IsotropicPolycrystalSCA: Moduli of randomly oriented polycrystals with
% varing grain aspect ratios using a self-consistent formulation.

% Inputs
% c: 6X6 matrix: Elastic stiffness tensor of single crystal in Voigt notation
% a: real scalar value: crystal grain aspect ratio, a<1 implies oblate, a =
% 1 implies spherical and a>1 implies prolate grain shapes.

% Outputs
% ksc: real scalar value: effective bulk modulus of randomly oriented polycrystals using a self-consistent formulation.
% gsc: real scalar value: effective shear modulus of randomly oriented polycrystals using a self-consistent formulation.

% References
% Gary's white papers

% Coded by Priyanka Dutta, 2017

% if ar==1
%     ar = .99;
% end

%* initial guess: polycrystal Hill moduli
[kv,kr,kh,gv,gr,gh] = IsotropicPolycrystalVRH(cc); 
a = kh+(4/3)*gh; b = kh-(2/3)*gh; c = gh;
c0 = [a b b 0 0 0;b a b 0 0 0;b b a 0 0 0;...
    0 0 0 c 0 0;0 0 0 0 c 0;0 0 0 0 0 c];
% c02 = [2.5 2.5 2.5 0 0 0;2.5 2.5 2.5 0 0 0;2.5 2.5 2.5 0 0 0;...
%     0 0 0 0 0 0;0 0 0 0 0 0;0 0 0 0 0 0];
% c03 = .3.*c02+.7.*c0;

% Voigt to Kelvin notation
ck = VoigtKelvinSwitch(cc,1);

% Eshelby tensor and stress concentration tensor
err = 1;
n = 0;
%n = 0;
while err>1e-3
%     n = n+1
    nu = (c0(1,1)-2*c0(4,4))/(2*c0(1,1)-2*c0(4,4));
    c0k = VoigtKelvinSwitch(c0,1); 
    mul = [1 1 1 sqrt(2) sqrt(2) sqrt(2);1 1 1 sqrt(2) sqrt(2) sqrt(2);1 1 1 sqrt(2) sqrt(2) sqrt(2);...
    1/sqrt(2) 1/sqrt(2) 1/sqrt(2) 1 1 1;1/sqrt(2) 1/sqrt(2) 1/sqrt(2) 1 1 1;1/sqrt(2) 1/sqrt(2) 1/sqrt(2) 1 1 1];
    sk = makeSpheroidalEshelby(ar,nu).*mul
%     sk = EshelbyTensorIsoSpheroid(nu, ar);
    I = eye(6);
    tk = inv( I + TensorDoubleContraction( TensorDoubleContraction( sk,inv(c0k) ),(ck - c0k) ) );
    tkiso = TensorIsotropicProjection( tk )
    cnew = TensorDoubleContraction(TensorDoubleContraction( ck,tk ),inv(tkiso));
% cnew = TensorDoubleContraction(TensorDoubleContraction( ck,tk ),I);
%     cnewiso = TensorDoubleContraction(TensorDoubleContraction( TensorIsotropicProjection(ck),tkiso ),I);
%     cnew = TensorDoubleContraction(TensorDoubleContraction( ck,tk ),inv(tkiso));
    cnewiso =  TensorDoubleContraction(TensorIsotropicProjection( cnew ),I);
    err = sum(sum((c0k-cnewiso).*(c0k-cnewiso)));
    c0 = VoigtKelvinSwitch(cnewiso,2);
    %n = n+1;
end
tkiso;
%err
%c0 = VoigtKelvinSwitch(c0,1);
ksc = c0(1,1)-(4/3)*c0(4,4);
gsc = c0(4,4);

end

