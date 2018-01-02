function [coeff,pred]=p_fit(x,y,z,n)
porders = n*(n+3)/2;
A=zeros(length(x),porders);
idx = 0;
for i = 1:n
    for j = 0:i
        idx = idx +1;
        A(:,idx) = x.^j.*y.^(i-j);
    end
end

coeff=A\z;
pred=A*coeff;

end