function [Y] = P_T(U,V,X)
%P_T 此处显示有关此函数的摘要
%   此处显示详细说明
Y=U*U'*X+X*V*V'-U*U'*X*V*V';
end

