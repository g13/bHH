function [coeff,pred]=p_fitb110k0_beta(x,y,z)
z = z-x-y;
A=ones(length(x),2);
A(:,2) = x.*y;
if sum(A(:,2)) == 0
     coeff = zeros(2,1);
     coeff(1) = mean(z);     
%coeff=A\z;
else
     coeff = A\z;
end

pred=A*coeff+x+y;

end
