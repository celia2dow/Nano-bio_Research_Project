function y = w8mean(x,w)
% W8MEAN  Weighted Average or mean value.
x=x(:)';
w=w(:)';
y = sum(w.*x)./sum(w);
end