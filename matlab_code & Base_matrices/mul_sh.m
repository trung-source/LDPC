function y = mul_sh(x,k)
%x: Khoi vao
%k: -1 hoac so lan dich
%y: Khoi ra
if (k==-1)
    y = zeros(1,length(x));
else
    y = [x(k+1:end) x(1:k)];
end
