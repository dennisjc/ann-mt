#include "functions.h"

int main(int argc, char *argv[])
{
 
  // test and sample the intensity of blobs at different locations.

  double xi[]={5450,5450,5450,5450,5450,5450,5450,
	      5450,5450,5450,5450,5450,5450,5450,
	      5450,5450,5450,5450,5450,5450,5450,
	      5450};
  vector<double> xIntervals(xi,xi+22);
  /*xIntervals={20000.0, 10000.0, 2500.0, 2500.0, 2500.0, 2500.0, 2500.0, 
	      2500.0, 10000.0, 2500.0, 2500.0, 2500.0, 2500.0, 10000.0, 
	      2500.0, 2500.0, 2500.0, 2500.0, 2500.0, 2500.0, 10000.0, 
	      20000.0}; */
  double yi[]={20000.0, 10000.0, 2500.0, 2500.0, 2500.0, 2500.0, 2500.0, 
	      2500.0, 2500.0, 2500.0, 20000.0, 2500.0, 2500.0, 2500.0, 
	      2500.0, 2500.0, 2500.0, 2500.0, 2500.0, 10000.0, 20000.0};
  vector<double> yIntervals(yi,yi+21);
  double zi[]={500.0, 1000.0, 1500.0, 3000.0, 4000.0, 5000.0, 5000.0, 
	      10000.0, 10000.0, 20000.0, 40000.0};
  vector<double> zIntervals(zi,zi+11);

  vector< vector<double> > params(2);
  params[0]=vector<double>(12);
  params[0][offs_inten]=0;
  params[0][offs_stren]=0.8;
  params[0][offs_atten]=3;
  params[0][offs_xpos]=0.5;
  params[0][offs_ypos]=0.5;
  params[0][offs_zpos]=0.2;
  params[0][offs_xrad]=0.1;
  params[0][offs_yrad]=0.2;
  params[0][offs_zrad]=0.4;
  params[0][offs_xrot]=0.2;
  params[0][offs_yrot]=0.5;
  params[0][offs_zrot]=0.3;

  params[1]=vector<double>(12);
  params[1][offs_inten]=255;
  params[1][offs_stren]=0.5;
  params[1][offs_atten]=20;
  params[1][offs_xpos]=0.5;
  params[1][offs_ypos]=0.5;
  params[1][offs_zpos]=0.2;
  params[1][offs_xrad]=0.3;
  params[1][offs_yrad]=0.1;
  params[1][offs_zrad]=0.2;
  params[1][offs_xrot]=0.2;
  params[1][offs_yrot]=0.8;
  params[1][offs_zrot]=0.3;

 
  cout.precision(2);
  render3DModel(xIntervals,yIntervals,zIntervals,params,0.5);
}

  
