function [coeff,pred,r2]=p_fit110k0(x,y,z0)
    z = z0-x-y;
    A = x.*y;
    % coeff=A\z;
    % coeff(isnan(coeff)) = 0;
    if sum(A) == 0
        coeff = 0;
    else
        coeff = (A'*z)/(A'*A);
    end
    pred=A*coeff+x+y;
    vpred = z0 - pred;
    SSres = vpred*vpred'
    vtot = z0 - mean(z0);
    SStot = vtot * vtot'
    if (SStot == 0)
        r2 = 1;
    else
        r2 = 1 - SSres/SStot;
    end
end
