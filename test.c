#include<stdlib.h>
#include<stdio.h>
#include<math.h>  
#include<complex.h>

//include -lm during the gcc compilation to ensure the mathematical functions work

int degree = 8;   //This needs to be replaced by the incoming command line argument which gives the degree of the polynomial

int main()
{

  double complex x;
  double complex xk;
  double complex x_k1;
  
  x = 2 + 4*I;  //This is the main complex number variable being used for the testing the functions written below for newtons method. Maybe has to be replaced with the value of complex number in the +-2 real and imaginary plane (maybe?)

  xk = x; x_k1 = x; // Written just so that xk has some value and can be used for computation of the Newtons Iterations Switch Case part of the program below
  
  printf("%lf\n", creal(xk)); //Printing for reference
  printf("%lf\n", cimag(xk)); //Printing for reference

  //**************************************Checking*******************************************************

  int checks = 1;    //Check variable set to 1 by default to state to assume that all validation checks are true. This will be set to 0 if any of the INITIAL validation cehcks is false
  
  //-----------precomputing of the exact values of the roots of the given polynomial----------------------

  double complex roots[degree-1];

  switch(degree)
    {

    case 1:
      //STATEMENTS FOR DEGREE 1
      roots[0] =  1;
      break;

    case 2:
      //STATEMENTS FOR DEGREE 2
      roots[0] =  1;
      roots[1] =  -1;
      break;


    case 3:
      //STATEMENTS FOR DEGREE 3
      roots[0] =  1;
      roots[1] =  -0.500000000000000 - 0.866025403784439 *I;
      roots[2] =  -0.500000000000000 + 0.866025403784439 *I;

      break;

    case 4:
      //STATEMENTS FOR DEGREE 4
      roots[0] =  1;
      roots[1] =  -1;
      roots[2] =  0+1*I;
      roots[3] =  0-1*I;

      break;

    case 5:
      //STATEMENTS FOR DEGREE 5
      roots[0] =  1;
      roots[1] =  -0.809016994374947 - 0.587785252292473 *I;
      roots[2] =  0.309016994374947 + 0.951056516295154 *I;
      roots[3] =  0.309016994374947 - 0.951056516295154 *I;
      roots[4] =  -0.809016994374947 + 0.587785252292473 *I;

      break;

    case 6:
      //STATEMENTS FOR DEGREE 6
      roots[0] =  1;
      roots[1] =  -1;
      roots[2] =  -0.500000000000000 - 0.866025403784439 *I;
      roots[3] =  0.500000000000000 + 0.866025403784439 *I;
      roots[4] =  0.500000000000000 - 0.866025403784439 *I;
      roots[5] =  -0.500000000000000 + 0.866025403784439 *I;

      break;

    case 7:
      //STATEMENTS FOR DEGREE 7
      roots[0] =  1;
      roots[1] =  -0.900968867902419 - 0.433883739117558 *I;
      roots[2] =  0.623489801858734 + 0.781831482468030 *I;
      roots[3] =  -0.222520933956314 - 0.974927912181824 *I;
      roots[4] =  -0.222520933956314 + 0.974927912181824 *I;
      roots[5] =  0.623489801858734 - 0.781831482468030 *I;
      roots[6] =  -0.900968867902419 + 0.433883739117558 *I;

      break;

    case 8:
      //STATEMENTS FOR DEGREE 8
      roots[0] =  1;
      roots[1] =  -1;
      roots[2] =  0+1*I;
      roots[3] =  0-1*I;
      roots[4] =  -0.707106781186548 - 0.707106781186548 *I;
      roots[5] =  0.707106781186548 + 0.707106781186548 *I;
      roots[6] =  0.707106781186548 - 0.707106781186548 *I;
      roots[7] =  -0.707106781186548 + 0.707106781186548 *I;

      break;

    case 9:
      //STATEMENTS FOR DEGREE 9
      roots[0] =  1;
      roots[1] =  -0.939692620785908 - 0.342020143325669 *I;
      roots[2] =  0.766044443118978 + 0.642787609686539 *I;
      roots[3] =  -0.500000000000000 - 0.866025403784439 *I;
      roots[4] =  0.173648177666930 + 0.984807753012208 *I;
      roots[5] =  0.173648177666930 - 0.984807753012208 *I;
      roots[6] =  -0.500000000000000 + 0.866025403784439 *I;
      roots[7] =  0.766044443118978 - 0.642787609686539 *I;
      roots[8] =  -0.939692620785908 + 0.342020143325669 *I;

      break;


    default:
      fprintf(stderr, "unexpected degree\n");
      //      exit(1);     //Commented out for the purpose of jotting down code - maybe have to remove in the final thing

    }

  //Below 2 print statements only for testing the code out. To be removed in the final one
  printf("Real part of 4th Root : %lf\n",creal(roots[5]));
  printf("Imaginary part of 4th Root : %lf\n", cimag(roots[5]));

  
  //----------upper bound on the absolute value of the real and imaginary part of x-----------------------

  //  if (!((creal(x)*creal(x) <= 100000000000000000000) && (cimag(x)*cimag(x) <= 100000000000000000000)))
  if (!(((exp(log(creal(x))*0.5)) <= 10000000000) && ((exp(log(cimag(x))*0.5)) <= 10000000000)))
    {
      checks = 0;
    }
  
  //------------lower bound on the absolute value of x---------------------------------------------------

  //  if (!((((creal(x)*creal(x)) + (cimag(x)*cimag(x))) * ((creal(x)*creal(x)) + (cimag(x)*cimag(x)))) >= 0.000001))
  if (!((exp(log((creal(x)*creal(x)) + (cimag(x)*cimag(x)))*0.5)) >= 0.001))
    {
      checks = 0;
    }
  

  //------------lower bound on the absolute value of x - x’, where x’ is one of the roots of the given polynomial----------


  // Loops through all the roots to ensure the absolute value for all roots are satusfying the condition. However, just the IF statement can be can be used in isolation or somewhere else if only one particular root needs to be checked.
  for (int i = 0; i < degree; i++)
    {

      if (!((exp(log(((creal(x) - creal(roots[i]))*(creal(x) - creal(roots[i]))) + ((cimag(x) - cimag(roots[i]))*(cimag(x) - cimag(roots[i]))))*0.5)) >= 0.001))
	{
	  checks = 0;
	}
            
    }
  
    
  //---------------------------------------------------------------------------------------------------------------------------------------------------------------

  //**************************Implement details: Iteration*****************************************************************************************

  // These are the checks that need to be performed during the iteration process. Should potentially be included in the loop and validated everytime a new x_k1 value is calculated (maybe?)
  
  //-------------If x_i is closer than 10^-3 to one of the roots of f(x), then abort the iteration.------------------------------------------------
 
for (int i = 0; i < degree; i++)
  {

    
    if (!((exp(log(((creal(x_k1) - creal(roots[i]))*(creal(x_k1) - creal(roots[i]))) + ((cimag(x_k1) - cimag(roots[i]))*(cimag(x_k1) - cimag(roots[i]))))*0.5)) >= 0.001))
      {
	checks = 0;
      }

  }

//------------------------------abort iteration if x_i is closer than 10^-3 to the origin--------------------------------------------------------

 if (!((exp(log(((creal(x_k1) - creal(0 + 0*I))*(creal(x_k1) - creal(0 + 0*I))) + ((cimag(x_k1) - cimag(0 + 0*I))*(cimag(x_k1) - cimag(0 + 0*I))))*0.5)) >= 0.001))
   {
     checks = 0;
   }

 //------------------------------ if the absolute value of its real or imaginary part is bigger than 10^10----------------------------------------

 // if (!((creal(x)*creal(x) <= 100000000000000000000) && (cimag(x)*cimag(x) <= 100000000000000000000)))
 if (!(((exp(log(creal(x_k1))*0.5)) <= 10000000000) && ((exp(log(cimag(x_k1))*0.5)) <= 10000000000)))
   {
     checks = 0;
   }


 //*************************************************************************************************************************************************
  
  
  //************************Computation - finding an expression for the iteration step*****************************************

  /*This switch statement needs to be placed inside the iteration loop. For Example : 

    for (...)
    {

    <switch statement>
     
     xk = x_k1; --> This would update the value of xk for the next iteration with the calculated value of x_k1;

    }

  */
  
  switch ( degree )
    {
    case 1:
      //STATEMENTS FOR DEGREE 1
      x_k1 = xk - (xk-1);
      break;
      
    case 2:
      //STATEMENTS FOR DEGREE 2
      x_k1 = xk - (((xk*xk)-1)/(2*xk));
      break;
      
    case 3:
      //STATEMENTS FOR DEGREE 3
      x_k1 = xk - (((xk*xk*xk)-1)/(3*xk*xk));
      break;
      
    case 4:
      //STATEMENTS FOR DEGREE 4
      x_k1 = xk - (((xk*xk*xk*xk)-1)/(4*xk*xk*xk));
      break;
      
    case 5:
      //STATEMENTS FOR DEGREE 5
      x_k1 = xk - (((xk*xk*xk*xk*xk)-1)/(5*xk*xk*xk*xk));
      break;
      
    case 6:
      //STATEMENTS FOR DEGREE 6
      x_k1 = xk - (((xk*xk*xk*xk*xk*xk)-1)/(6*xk*xk*xk*xk*xk));
      break;
      
    case 7:
      //STATEMENTS FOR DEGREE 7
      x_k1 = xk - (((xk*xk*xk*xk*xk*xk*xk)-1)/(7*xk*xk*xk*xk*xk*xk));
      break;
      
    case 8:
      //STATEMENTS FOR DEGREE 8
      x_k1 = xk - (((xk*xk*xk*xk*xk*xk*xk*xk)-1)/(8*xk*xk*xk*xk*xk*xk*xk));
      break;
      
    case 9:
      //STATEMENTS FOR DEGREE 9
      x_k1 = xk - (((xk*xk*xk*xk*xk*xk*xk*xk*xk)-1)/(9*xk*xk*xk*xk*xk*xk*xk*xk));
      break;
      
    default:
      fprintf(stderr, "unexpected degree\n");
      //      exit(1);     //Commented out for the purpose of jotting down the code - may have to be removed in the final thing
    }

  printf("%lf\n", creal(x_k1));
  printf("%lf\n", cimag(x_k1));

  //*************************************************************************************************************************************

  
  return 0;

}




