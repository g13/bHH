function [coeff,pred]=p_fit2(x,y,z)
z = z-x-y;
A = x.*y;

coeff=A\z;
pred=A*coeff+x+y;

end