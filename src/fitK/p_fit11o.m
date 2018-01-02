function [coeff,pred]=p_fit11o(x,y,z,n)
z = z-x-y;
porders = n*(n+3)/2;
A=zeros(length(x),porders-2);
idx = 0;
for i = 2:n
    for j = 0:i
        idx = idx +1;
        A(:,idx) = x.^j.*y.^(i-j);
    end
end

coeff=A\z;
pred=A*coeff+x+y;

end
