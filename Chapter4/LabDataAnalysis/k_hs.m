function [k_hsu,k_hsl] = k_hs(f1,k1,g1,f2,k2,g2)
%** Hashin-Shtrikman-Walpole bounds for K, RPH p224 **%

k_hsu = k1 + f2./((k2-k1).^-1 + f1.*(k1+4/3.*max(g1,g2)).^-1);
k_hsl = k1 + f2./((k2-k1).^-1 + f1.*(k1+4/3.*min(g1,g2)).^-1);

end

