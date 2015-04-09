function [f_prime_out] = dfunc(x)

%give your own derivative here. it is a vector output in 
% f_prime_out=(del/del(x1) del/del(x2) ... etc)'

%f_prime_out(1) = -sin(x(1)).*cos(x(2));
%f_prime_out(2) = -cos(x(1)).*sin(x(2));

%expf = exp(-((x(1)^2+x(2)^2)^0.5)/4);
%ss = -((x(1)^2+x(2)^2)^0.5)/4;

%see func.m for function definition 
%del func/del x(1)
%f_prime_out(1) = cos(x(2)).*(-sin(x(1)).*expf-cos(x(1)).*expf.*ss.*x(1));
%del func/del x(2)
%f_prime_out(2) = cos(x(1)).*(-sin(x(2)).*expf-cos(x(2)).*expf.*ss.*x(2));

f_prime_out = zeros(100, 1);
for i = 1:100
    f_prime_out(i) = exp(x(i)) - sqrt(i);
end

%here I tried to maximise the function
%f_prime_out = -1.*f_prime_out';
