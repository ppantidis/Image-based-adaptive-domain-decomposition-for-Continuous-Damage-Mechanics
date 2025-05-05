function [u] = func_deltaU(Khh_inv,CA_inverse,CA_inverseB,K_F,Res_F_F,lenA,lenD)

% -------------------------------------------------------------------------
u2 = ((K_F(lenA+1:lenA+lenD, lenA+1:lenA+lenD) - CA_inverseB)) \ (- Res_F_F(lenA+1:lenA+lenD) + CA_inverse * Res_F_F(1:lenA));
% u2 = sparse((K_F(lenA+1:lenA+lenD, lenA+1:lenA+lenD) - CA_inverseB)) \ (- Res_F_F(lenA+1:lenA+lenD) + CA_inverse * Res_F_F(1:lenA));

% -------------------------------------------------------------------------
% x1 = - P * (Res_F_F(1:lenA) + K_F(1:lenA,lenA+1:lenA+lenD) * u2);
% y1 = L \ x1;
% u1 = U \ y1;

u1 = - Khh_inv * (Res_F_F(1:lenA) + K_F(1:lenA,lenA+1:lenA+lenD) * u2);

u = cat(1,u1,u2);

end

