
// Simple example of using the CG_DESCENT (Conjugate Gradient -- CG--
// optimizer code from Hager and Zhang of U Florida.
//
// SETTING UP THE PROBLEM

//  First, you need two included files.  math.h gives you all the basic
// math operations.  Cg_user.h gives you the declarations specific
// to the CG_DESCENT package:

#include <math.h>
#include "cg_user.h"


//
// SUPPLYING THE FUNCTION AND ITS GRADIENT
// To optimize a function, you are required to have a routine that
// computes the function value, and the gradient vector.  The
// function prototypes look like this:


// function myvalue takes the current solution vector x[i] and its
// dimension n, and returns a double which is the value of f(x).

double myvalue
(
    double   *x,    // a pointer to an array of doubles for x[i]
    INT       n     // dimension of the problem
) ;


// function mygrad takes the current solution vector x[i]
// and the dimension n, and computes the gradient, which is
// a vector of doules f length n,  
// g[] = [del_f/del x1, ...  del_f/del_xn]
// Since we supply a pointer to the g[] array, mygrad returns void.

void mygrad
(
    double    *g, // a pointer to an array of doubles for gradient g[i]
    double    *x, // a pointer to an array of doubles for x[i]
    INT        n  // dimension of the problem
) ;

// The authors of CG)_DESCENT note that it is sometimes more efficient
// to compute the function f() and its gradient g=[...del_f/del_xi ...]
// at the SAME TIME, in one function call.  So, there is a prototype
// for this combo evaluation as well.   If you supply it to the
// core optimizer routine cg_descent, it will use it.  If you dont
// build this function, you just pass in a NULL.  As far as I can tell,
// you STILL need to supply 'myvalue' and 'mygrad' anyway. 
// You give it pointer to the x[] vector, to the gradient g[] array, 
// and problem size n, and myvalgrad returns the function value, and
// and also fills the g[] array with the current gradient at x[]


double myvalgrad
(
    double    *g,
    double    *x,
    INT        n
) ;




//  CALLING THE OPTIMIZER -- A SIMPLE DRIVER EXAMPLE

int main (void)
{

    // stuff specifically for cg_descent optimizer call
    double *x ;           // pointer to the x[i] vector
    INT n ;               // dimension of the problem
    double cg_tol;        // tolerance for when to stop optimizing
    cg_stats Stats ;      // cg_descent computes some useful stats while
                          // it runs.  You pass a pointer to this, and
                          // after the optimization, you can poke this for
                          // useful info about the run.
    cg_parameter Parm ;   //  cg_descent takes a bunch of control params,
                          // this struct is how you set them.  Good news is
                          // is you can default pretty much ALL of them
                          // and it should still work
    INT cg_return;        // cg_descent returns an int with 'status' info
                          // about the run.  The status nums and meanings are:
                          //  -2 (function value became nan)
                          //  -1 (starting function value is nan)
                          //  0 (convergence tolerance satisfied)
                          //  1 (change in func <= feps*|f|)
                          //  2 (total iterations exceeded maxit)
                          //  3 (slope always negative in line search)
                          //  4 (number secant iterations exceed nsecant)
                          //  5 (search direction not a descent direction)
                          //  6 (line search fails in initial interval)
                          //  7 (line search fails during bisection)
                          //  8 (line search fails during interval update)
                          //  9 (debugger is on and the func value increases)
                          //  10 (out of memory) 
                          //... so basically you want a ZERO return valyue.


    // other stuff for the driver 
    INT i;                // a simple counter to set up initial soln guess


    // RUNNING THE CG_DESCENT OPTIMIZER

    // allocate space for solution 
    n = 100 ;
    x = (double *) malloc (n*sizeof (double)) ;

    
    // initialize default parameter values for cg_descent
    cg_default (&Parm) ;

    //  cg_descent has LOTS of parameters you can control, but defaults
    //  seem to work well for all of them.  A few worth knowing and perhaps
    //  setting are:
    //
    //    Parm.PrintFinal = FALSE   
    //        The default setting is TRUE, which means it prints to stdout
    //        the result of each call to cg_descent.  This can be annoying,
    //        depending on you are doing output, so you can turn it OFF.
    //
    //    Parm.PrintLevel = integer 0 or 1 or 2 or 3 
    //        Internal debugging/progress for cg_descent.  
    //        Default = 0 = nothing, setting to 1 or 2 or 3 
    //        increases amount of output
    //
    //    Parm.step = size of first step in the direction of the gradient
    //        i.e., we calculate the gradient, then take a first step
    //        of size ALPHA in this descent direction.  cg_descent has
    //        has a default computation for how big this should be, but
    //        but the authors also say that sometimes this is the one
    //        place you might want to use what you know about the problem
    //        and just set it yourself.  If you don't set it, it defaults
    //        to something reasonable (you hope).  But you can set it, like
    //        this:   Parm.step = 1.0;
    //
    //    Parm.feps = stop optimizing when the estimated change in the
    //        value of the function (from myvalue) is smaller than feps*|f|.
    //        In other words, if you want cg_descent to quit when the function
    //        is changing by less the 1%, you do:  Parm.feps = 0.01;
    //        If you don't set anything, cg_descent uses a tolerance value
    //        passed in on the function call to determine when to quit.
    //        This is described later in this code.  
    //
    //    Parm.maxit_fac = (double) num ; abort cg after maxit_fac*n iterations
    //         So this is the other way to tell cg_descent how hard to work.
    //         If your problem is of dimension 1000, and you don't want to
    //         do more than 2000 CG iterations (calls to gradient), then
    //         set maxit_fac = 2.0 

    
    // HOW TO CALL THE OPTIMIZER
    //
    // cg_descent itself is called with 9 parameters:
    //
    //   double       *x,         on input: starting guess, on output: the solution 
    //   INT           n,         problem dimension 
    //   cg_stats     *Stat,      structure with statistics (can be NULL) 
    //   cg_parameter *UParm,     user parameters, if NULL = use default params
    //   double        grad_tol   the stopping tolerance, ie, when are we done?
    //                            remember that we search downhill until we
    //                            cannot make any more progress.  The builtin
    //                            mechanism for stopping wants to find where
    //                            the gradient is approx [0,0,...0] as a vector.
    //                            If we compute the length of the gradient vector
    //                            and it is smaller than grad_tol, we quit.
    //                            NOTE: you can also specify a stopping
    //                            rule in terms of the change in the function val
    //                            |f| using the cg_parameters structure.  But you
    //                            still need to put in a value for grad_tol.
    //   double   (*myvalue)      Pointer the routine to compute f = myvalue
    //   void     (*mygrad)       Pointer to routine to compute gradient(f) 
    //   double   (*myvalgrad)    Pointer to routine to compute both func & grad
    //                            If NULL = compute just using mvalue & mygrad 
    //   double  *Work            either size 4n work array or NULL 
  

    //===================================================================
    // EX1:  set up and solve the problem 
    

    // starting guess for solution = [1 1 1 ... 1]
    for (i = 0; i < n; i++) x [i] = 1.0 ;

    // set the tolerance on the ||gradient vector||
    cg_tol = 1.e-8;

    // parameter settings = all default for this run

    // solve it
    printf("FIRST CALL TO CG_DESCENT, cg_tol=%g\n", cg_tol);
    cg_return = cg_descent(x, n, &Stats, &Parm, cg_tol, myvalue, mygrad, myvalgrad, NULL) ;

    // note that at this point x[] is the final solution, 
    // and cg_descent by default printed a lot of stuff to stdout



    //===================================================================
    // EX2:  set up and solve the problem again

    // reset the starting guess for solution
    for (i = 0; i < n; i++) x [i] = 1.0 ;

    // set the tolerance LOOSER on the ||gradient vector||
    cg_tol = 1.e-5;

    // parameter settings = all default for this run

    // solve it
    printf("SECOND CALL TO CG_DESCENT, cg_tol=%g\n", cg_tol);
    cg_return = cg_descent(x, n, &Stats, &Parm, cg_tol, myvalue, mygrad, myvalgrad, NULL) ;

    // note that at this point x[] is again the final solution,
    //  cg_descent again printed a lot of stuff, and we ran LESS
    //  iterations of the CG solver since we had a LOOSER tolerance on 
    //  the size of the final gradient vector


    //===================================================================
    // EX3:  set up and solve the problem again

    // reset the starting guess for solution
    for (i = 0; i < n; i++) x [i] = 1.0 ;

    // set the tolerance on the ||gradient vector||
    cg_tol = 1.e-8;

    //  Parameter settings: lets NOT print all that stuff to Stdout by default
    Parm.PrintFinal = FALSE ;

    // solve it
    printf("THIRD CALL TO CG_DESCENT, cg_tol=%g\n", cg_tol);
    cg_return = cg_descent(x, n, &Stats, &Parm, cg_tol, myvalue, mygrad, myvalgrad, NULL) ;

    // note that at this point x[] is again the final solution,
    //  but this time, cg_descent did NOT print anything...

    // ..but we can look at the Stats structure and print the info OURSELVES
    printf("Final function value = %g\n", Stats.f ) ;
    printf("||gradient at solution|| = %g\n", Stats.gnorm) ;
    // for some reason, cg-descent uuses INT = long ints for these, we cast so gcc won't yell at us
    printf("Number of calls to function eval = %d\n", (int)Stats.nfunc);
    printf("Number of calls to gradient eval = %d\n", (int)Stats.ngrad);
    printf("Number of CG interations =%d\n", (int)Stats.iter);


    free (x) ; /* free work space */
    return(0);
}


//===================================================================
//
// Example function to optimize.  
//     This is a 100 dimensional function of [x1 x2 ... x100]
//     Function is SUM(i=1 to 100) [ exp(xi) - xi*sqrt(i) ]
//
//===================================================================

double myvalue
(
    double   *x ,
    INT       n
)
{
    double f, t ;
    INT i ;
    f = 0. ;
    for (i = 0; i < n; i++)
    {
        t = i+1 ;
        t = sqrt (t) ;
        f += exp (x [i]) - t*x [i] ;
    }
    return (f) ;
}


//===================================================================
//
// Gradient of our function from myvalue.  
//
//   NOTE:  this one has a ClOSED FORM gradient, ie, a formula.
//          In your placer, you are NOT going to have a formula,
//          you will have to do this numerically, carefully, yourself
//===================================================================


void mygrad
(
    double    *g ,
    double    *x ,
    INT        n
)
{
    double t ;
    INT i ;
    for (i = 0; i < n; i++)
    {
        t = i + 1 ;
        t = sqrt (t) ;
        g [i] = exp (x [i]) -  t ;
    }
    return ;
}


//===================================================================
//
//  Combo eval:  Our function AND its gradient
//
//   NOTE:  it just turns out to be easy to do in this example
//===================================================================


double myvalgrad
(
    double    *g,
    double    *x,
    INT        n
)
{
    double ex, f, t ;
    INT i ;
    f = (double) 0 ;
    for (i = 0; i < n; i++)
    {
        t = i + 1 ;
        t = sqrt (t) ;
        ex = exp (x [i]) ;
        f += ex - t*x [i] ;
        g [i] = ex -  t ;
    }
    return (f) ;
}

