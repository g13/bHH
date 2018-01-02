function [coeff,pred]=p_fit110k0o(x,y,z,n)
porders = n*(n+3)/2;
A=zeros(length(x),porders-4);
z = z-x-y;
A(:,1) = x.*y;

idx = 1;
for i = 3:n
    for j = 0:i
        idx = idx +1;
        A(:,idx) = x.^j.*y.^(i-j);
    end
end

coeff=A\z;
pred=A*coeff+x+y;

end
