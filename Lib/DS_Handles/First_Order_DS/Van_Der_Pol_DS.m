function xd = Van_Der_Pol_DS(x, alpha)
    if(nargin < 2)
        alpha = 0.9;
    end
    xd = x(2,:);
    xd(2,:) = -x(1,:) + alpha*(1 - x(1,:).^2).*x(2,:);
end
