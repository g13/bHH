function [point, scale] = autoSpacing(f,n,scale,x)
    if nargin < 4
        assert(length(f) >= n);
        data = f;
        if nargin < 3
            scale = round(length(f)/max(f)); % maximize angle changes in the whole function.
        end
    else
        assert(isa(f,'function_handle'));
        assert(nargin == 4);
        assert(length(x) >= n);
        data = f(x);
    end
    % f should be a column vector
    l = length(data);
    d2 = abs(diff(atan(diff(data)*scale))*scale);
    w0 = [d2(1);d2;d2(l-2)];
    w = w0/sum(w0)*n;
    csw = cumsum(w);
    point = zeros(n,1);
    point(1) = 1;
    ip = 1;
    for i=1:l
        if csw(i) > ip
            point(ip+1) = i;
            ip = ip + 1;
        end
    end
end