function [coeff,pred]=p_fit110k0(x,y,z)
z = z-x-y;
A = x.*y;
% coeff=A\z;
% coeff(isnan(coeff)) = 0;
if sum(A) == 0
    coeff = 0;
else
    coeff = (A'*z)/(A'*A);
end
pred=A*coeff+x+y;
end
