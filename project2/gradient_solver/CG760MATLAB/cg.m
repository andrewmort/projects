% Martin King, ICTP, 29 March 2005
%1. Conjugate Gradient Method with Flecther-Reeves (or Polak-Ribiere)
%   to find a vector x that gives a MINIMUM of a function (a scalar).
%2. Ideas taken from J.R. Shewchuk and Numerical Recipes.
%3. You must modify your own function to minimise in a Matlab 
%   function file called func.m (scalar output) and the first derivative of that
%   function in a matlab function file called dfunc.m, which has a vector output in 
%   (del/del(x1) del/del(x2) ... etc)'. 
%4. If you want to MAXIMISE a function, multiply -1 to the output of func.m
%   and dfunc.m (be careful here, dfunc.m may use func.m; to be safe, give dfunc.m
%   the original output of func.m and then multiply -1 to dfunc.m at the end). 
%   If you are brave, reverse the search direction r and d for maximisation.
%5. If the method is not converging or is giving you a solution that doesn't make sense, 
%   change the initial guess.
%6. As an example, a simple function is given in func.m and its gradient vector in dfunc.m. 
%   Change the initial guess to x = [1 ; 1] for example, the solution it gives is incorrect. 
%   The reason is obvious if you plot the function (it is the saddle points). 
%7. I have used these scripts to optimise a fairly complicated function. They 
%   seem to work well. If you notice any bug or have any comment, please email me
%   king at ictp at it

%put initial guess here (x is an n-dimensional vector).
%x=[-2 ; -3];
x = ones(100, 1);

[F] = func(x);
[F_prime] = dfunc(x);

%From here in this m-file, I usually follow the notations in the 'pseudocode' of B4 
%in Shewchuk's note.

r=-1.*F_prime;
d=r;
%g is for use in Polak-Ribiere
g=r;
delta_new=r'*r; 
delta_0=delta_new;
fp = F;
%ftol is a convergence tolerance
ftol=1.e-7;
ITERMAX = 10000;
%Don't worry too much about this
EPS=1.e-10;

for iter = 1 : ITERMAX 

%Doing the line search here. First bracket a minimum, then use Golden section to find it.
%Not using Newton-Raphson as in Shewchuk. So you don't need the second derivative.
  [ax,bx,cx,fa,fb,fc] = func_mnbrak(0,1,x,d);
  [xt,golden] = func_golden(ax,bx,cx,x,d);
%To recover vector x, which is along d at xt away from initial x.
   x = x + xt.*d;
%The function value at x is golden as returned by func_golden.
   F = golden;
   [F_prime] = dfunc(x);

   r = -1.*F_prime;
   delta_old = delta_new;
%This is Fletcher-Reeves
   delta_new = r'*r;
%This is Polak-Ribiere
%   delta_new = (F_prime+g)'*F_prime;
   beta = delta_new/delta_old;
   d = r + beta * d ; 
   g = r;
   if r'*d <= 0
      d=r;
   end
%this convergence criterion is taken from NR.
   if 2.*abs(F-fp) < ftol.*(abs(F)+abs(fp)+EPS)
      iter
      x
      F
      break;
   end
   fp = F;

end
