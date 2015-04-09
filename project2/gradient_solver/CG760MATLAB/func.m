function [f_out] = func(x)

%give your own function here. it is a scalar output

%example function f_out = cos(x(1)).*cos(x(2));

%another example function cos x(1)cos x(2)exp(-sqrt(x(1)^2+x(2)^2)/4) 

%f_out = cos(x(1)).*cos(x(2)).*exp(-(x(1).^2+x(2).^2)^0.5/4);

f_out = 0;
for i = 1:100
    f_out = f_out + exp(x(i));
    f_out = f_out - sqrt(i)*x(i);
end

%here I try to maximise the function
%f_out = -1.*f_out;
