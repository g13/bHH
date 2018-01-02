function [coeff,pred]=p_fitb110k0(x,y,z)
z = z-x-y;
A=ones(length(x),2);
A(:,2) = x.*y;
% if sum(A(:,2)) == 0
%     coeff = zeros(2,1);
%     coeff(1) = mean(z);     
coeff=A\z;
% else
%     coeff = (A'*z)/(A'*A);
% end
% coeff = coeff(1);
pred=A*coeff+x+y;

end
