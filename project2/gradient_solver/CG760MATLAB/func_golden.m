function [xmin,golden] = func_golden(ax,bx,cx,x_in,d_in)

%1. Given bracket points ax, bx, cx, this function finds the minimum with tol; 
%   return independent variable as xmin (a scalar!) and function value as golden. 
%   See NR F77, p.394
%2. Uses a user-defined func.m
%3. I added x1t=x_in+x1.*d_in, etc, for the use of multidimensional conjugate 
%   gradient method. This step is performed by f1dim in NR on p.413.
%4. Like in func_mnbrak.m the outputs from this function are again scalars.
%5. You may want to change tol here.

tol=1.e-7;
R=0.61803399;
C=1.-R;

x0=ax;
x3=cx;
if abs(cx-bx) > abs(bx-ax) 
  x1=bx;
  x2=bx+C.*(cx-bx);
else
  x2=bx;
  x1=bx-C.*(bx-ax);
end
x1t=x_in+x1.*d_in;
[f1]=func(x1t);
x2t=x_in+x2.*d_in;
[f2]=func(x2t);

while abs(x3-x0) > tol.*(abs(x1)+abs(x2))
      if f2 < f1
         x0=x1;
         x1=x2;
         x2=R.*x1+C.*x3;
         f1=f2;
         x2t=x_in+x2.*d_in;
         [f2]=func(x2t);
      else
         x3=x2;
         x2=x1;
         x1=R.*x2+C.*x0;
         f2=f1;
         x1t=x_in+x1.*d_in;
         [f1]=func(x1t);
      end
end

if f1 < f2
  golden=f1;
  xmin=x1;
else
  golden=f2;
  xmin=x2;
end
