function [point, scale, csw] = autoSpacing(f,n,scale,p,x)
    if nargin < 5
        assert(length(f) >= n);
        data = f;
        if nargin < 4
            p = 0.85;       % penalty, 0 penalty non-bias aspect ratio.
            if nargin < 3
                scale = round(length(f)/(max(abs(f))-min(abs(f)))); % maximize angle changes in the whole function.
            end
        end
    else
        assert(isa(f,'function_handle'));
        assert(nargin == 4);
        assert(length(x) >= n);
        data = f(x);
    end
    if scale < 0
        scale = round(length(f)/(max(abs(f))-min(abs(f))));
    end
    % f should be a column vector
    l = length(data);
    d2 = abs(diff(atan(diff(data)*scale))*scale); % for angle turns.
    l2 = abs(diff(data(2:l)));      % for distance traveled.
    w0 = [d2(1)*l2(1);d2.*l2;d2(l-2)*l2(l-2)];
    if p ~= 0
        dev_w0 = w0 - mean(w0);
        pick = dev_w0>0;
        penalty = dev_w0(pick)*p;
        w0(pick) = w0(pick) - penalty;
    end
    w = w0/sum(w0)*n;
    csw = cumsum(w);
    point = zeros(n,1);
    ip = 0;
    for i=1:l
        if csw(i) >= ip
            point(ip+1) = i-1;
            ip = ip + 1;
            if ip == n
                break;
            end
        end
    end
end
