// Vanilla Options pricing with Black-Scholes and Monte-Carlo
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <time.h> 
#include "tools.hpp" 

using namespace std;

int get_cube_index(int i,int j,int k,int lat_dim)
{
  return k + lat_dim* j + lat_dim * lat_dim *i;
}

double Inv_Cum_Norm(  double x )
{
  // generates the inverse cumulative normal function
  // using the Moro method
  // provide number x on [0,1]

  double y = x - 0.5;
  double r,t;
  double sm1=0.0;

  if (abs(y) < 0.42 )
    {

      double a[4];
      double b[4];

      a[0] = 2.50662823884;
      a[1] = -18.61500062529;
      a[2] = 41.39119773534;
      a[3] = -25.44106049637;

      b[0] = -8.47351093090;
      b[1] = 23.08336743743;
      b[2] = -21.06224101826;
      b[3] = 3.13082909833;

      r = y*y;

      double sm2=0.0;
      
      for ( int j = 0 ; j < 4; j++ )
	{
	  sm1 = sm1 + y*a[j] * pow(r,j);
	  sm2 = sm2 + b[j] * pow(r,j+1);
	}

      t = sm1/(sm2+1.0);
    }
  else
    {
      double c[9];
      
      c[0] = 0.3374754822726147;
      c[1] = 0.9761690190917186;
      c[2] = 0.1607979714918209;
      c[3] = 0.0276438810333863;
      c[4] = 0.0038405729373609;
      c[5] = 0.0003951896511919;
      c[6] = 0.0000321767881768;
      c[7] = 0.0000002888167364;
      c[8] = 0.0000003960315187;
  
      if ( y < 0.0 )	  
	r = x;
      else
        r = 1-x;

      double s = log(-1*log(r)) ;

      for (int j = 0 ; j < 9 ; j++ )
	sm1 = sm1 + c[j]*pow(s,j);

      if (x > 0.5)
	t = sm1;
      else
	t = -1*sm1;	  
    }

  return t;
}


double Cum_Norm(double x)
{
  // cumulative normal function
  // x is any real number
  double y;
  if (x > 0.0)
    y = x;
  else
    y = -1*x;
  
  double k = 1/(1+0.2316419*y);

  double t = 1.0- exp(-0.5*y*y) / sqrt(2.0 * M_PI ) *
    k * (0.319381530 + k *(-0.356563782 + k * ( 1.781477937 +  k * (-1.821255978 + 1.330274429*k))));

  if (x <=  0.0 )
    t = 1-t;
  return t;
}

void random::init(double mu , double sig ) 
{
  srand(time(NULL));
  mean = mu;
  stdv = sig;
}

void Uniform_Random::initialize(double mu , double sig ) 
{
  init(mu,sig);
}

void Gaussian_Random::initialize(double mu , double sig ) 
{
  init(mu,sig);
}

double Uniform_Random::draw()
{
  int v1 = rand();
  double a1 = v1; 
  double a2 = a1 / RAND_MAX; 
  return a2 * 2  * stdv -  stdv + mean;   
}     
  
double Gaussian_Random::draw()
{
  int v1 = rand(); 
  double a1 = v1;
  double a2 = a1 / RAND_MAX;   
  return Inv_Cum_Norm(a2) * stdv + mean; 
}     
  
