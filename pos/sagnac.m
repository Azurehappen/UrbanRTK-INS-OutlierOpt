function delta = sagnac(p,x,y)
% Compute Sagnac correction/earth rotate correction
%
% Reference: https://escholarship.org/uc/item/1bf6w7j5
delta = p.omge*(x(1)*y(2)-x(2)*y(1))/p.c;
end