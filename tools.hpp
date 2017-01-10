// basic tools for Monte Carlo
#include <iostream>
#include <math.h>
#include <stdlib.h>

double Inv_Cum_Norm(  double x );
double Cum_Norm(double x);
int get_cube_index(int,int,int,int);
  
class random{
public:
  void init(double,double); 
  virtual double draw(void)=0;    
protected:
  double mean,stdv; 
};  
  
class Gaussian_Random: public random
{
public:
  void initialize(double = 0.0,double = 1.0); 
  double draw(void);
};

class Uniform_Random: public random
{
public:
  void initialize(double = 0.5,double = 0.5); 
  double draw(void);
};
