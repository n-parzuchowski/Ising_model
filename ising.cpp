// c++ ising model code for practice

#include <iostream>
#include <string>
#include <fstream>
#include "tools.hpp"

using namespace std; 

// Ising model Hamiltonian 
double site_energy(int site,int *lattice, int lat_dim, double Jng, double Bex){
  
  // get the indeces of i,j,k dimensions
  int i = site / (lat_dim*lat_dim);
  int j = (site - i * lat_dim *lat_dim) / lat_dim;
  int k = site - i * lat_dim *lat_dim - j * lat_dim;

  // H = -J * s1 * s2 - u B s1
  // neighbors
  double energy=0.0;
  int neighbor;
  int site_val = lattice[site];

  // calculate pair contribution from 6 nearest neighbors
  neighbor = lattice[get_cube_index(i,j,(k+1)%lat_dim,lat_dim)];
  energy += -1* Jng * site_val * neighbor;

  neighbor = lattice[get_cube_index(i,j,(k-1)%lat_dim,lat_dim)];
  energy += -1* Jng * site_val * neighbor;

  neighbor = lattice[get_cube_index(i,(j+1)%lat_dim,k,lat_dim)];
  energy += -1* Jng * site_val * neighbor;

  neighbor = lattice[get_cube_index(i,(j-1)%lat_dim,k,lat_dim)];
  energy += -1* Jng * site_val * neighbor;

  neighbor = lattice[get_cube_index((i+1)%lat_dim,j,k,lat_dim)];
  energy += -1* Jng * site_val * neighbor;  
  
  neighbor = lattice[get_cube_index((i-1)%lat_dim,j,k,lat_dim)];
  energy += -1* Jng * site_val * neighbor;    


  // magnetic field part. 
  energy = energy - Bex *  site_val;
  
  return energy;
}



double magnetization(int *lattice, int lattice_size)
{
  // compute magnetization of the lattice, i.e. average spin up - spin down. 
  
  double sm = 0;
  for (int thing= 0 ; thing < lattice_size ; thing++ ){
    sm += lattice[thing];
  }

  return sm / lattice_size;
}
    

double boltz_fac(double Temp,double Enew, double Eold)
// Statistical Mechanics probability factor
{
  return exp(- (Enew-Eold) / Temp );
}


int main() 
{ 
  
  // obtain input parameters =====================================
  int lattice_dim=20;
  double frac_up=1.0, Jneigh=0.1, Bext=0.0,Temp;
  
  cout << "Input Dimension of Lattice: ";
  cin >> lattice_dim;
  
  cout << "Input approximate fraction of spins oriented up: ";
  cin >> frac_up;
  
  cout << "Input J coupling constant between nearest neighbors: ";
  cin >> Jneigh;
  
  cout << "Input external magnetic field strength: ";
  cin >> Bext;

  //////////////////////////////////////////////////////////////

  int lattice_size = pow(lattice_dim,3);
  int *lattice;

  lattice = new int[lattice_size];

  // Random number generator
  Uniform_Random gen;
  gen.initialize();

  FILE * fx;
  fx = fopen("mag_v_temp.dat","w"); 
  
  Temp = 0.0;
  for (int T=0; T<20; T++ ) {
    Temp+=0.05;
    cout << "Temperature is: "<< Temp << endl;
  
  // fill lattice sites with spins either up or down
  for (int i = 0 ; i < lattice_dim ; i++) {
    for (int j = 0 ; j < lattice_dim ; j++) {
      for (int k = 0 ; k < lattice_dim ; k++) {

        double x = gen.draw();
	int spin = 1;
	if ( x > frac_up ) {
	  spin = -1;
	}

	lattice[get_cube_index(i,j,k,lattice_dim)]=spin;
	
      }
    }
  }

  double mag = magnetization(lattice,lattice_size);

  cout << "Initial Magnetization is: " << mag << endl;
  
  // equilibrate with Metropolis spin-flip

  for (int mc_step = 0 ; mc_step < 1000000; mc_step++ ){
    double x = gen.draw();
    int site = round(x*lattice_size);

    double Eold = site_energy(site,lattice,lattice_dim,Jneigh,Bext);

    // flip spin
    lattice[site] *= -1; 
    double Enew = site_energy(site,lattice,lattice_dim,Jneigh,Bext);

    // if Enew is less than Eold, keep the change (energetically favorable) 
    if ( Enew >= Eold )
      {
	double y = gen.draw();
	if (y  >  boltz_fac(Temp,Enew,Eold) )
	  {
	    // revert
	    lattice[site] *= -1;
	  }
      }
  }


  // measure average magnetism over 10000 steps, continue to use metropolis algorithm

  double mag_sum = 0.0;
  for (int mc_step = 0; mc_step < 10000; mc_step++)
    {
      double x = gen.draw();
      int site = round(x*lattice_size);
      
      double Eold = site_energy(site,lattice,lattice_dim,Jneigh,Bext);
      
      // flip spin
      lattice[site] *= -1; 
      double Enew = site_energy(site,lattice,lattice_dim,Jneigh,Bext);
      
      // if Enew is less than Eold, keep the change (energetically favorable) 
      if ( Enew >= Eold )
	{
	  double y = gen.draw();
	  if (y  >  boltz_fac(Temp,Enew,Eold) )
	    {
	      // revert
	      lattice[site] *= -1;
	    }
	}

      mag_sum += magnetization(lattice,lattice_size);



    }
    
  mag = mag_sum/10000.;
  
  cout << "Equilibrated magnetization is: " << mag << "\n"<< endl;      
  
  // write to file
  fprintf(fx,"%8.3f %8.3f\n",Temp,mag);	  

  
  } // end loop over temperature 
  fclose(fx);
  
  delete [] lattice;
  return 0;
}


