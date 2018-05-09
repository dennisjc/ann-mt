#include "functions.h"
#include <fstream>
#include <string>
int main(int argc, char *argv[])
{

  // usage: ./showParams model_parameters
  // outputs a formatted description of the paramters

  
  if(argc < 2){
    cout << "usage: ./showParams model_parameters" << endl;
  }

  double bgcolour = atof(argv[1]); 
  vector <vector<double > > params;	      
  params = processInput3D(extractParam(argc - 2, 
				       argv + 2));
  // ouput background
  
  cout << cout.precision(2);
  cout << fixed;
  cout << "background: " << bgcolour << endl;
  // output parm headings
  cout << "int\t str\t att\t xpos\t ypos\t zpos\t xrad\t yrad\t zrad\t xrot\t ";
  cout << "yrot\t zrot\t" << endl;
  cout << "------------------------------------------------------------------------------" << endl;
  u_int i=0;
  //cout.flush();
  while(i<params.size()){
    vector<double> blob = params[i];
    cout << blob[offs_inten] << "\t ";
    cout << blob[offs_stren] << "\t ";
    cout << blob[offs_atten] << "\t ";
    cout << blob[offs_xpos] << "\t ";
    cout << blob[offs_ypos] << "\t ";
    cout << blob[offs_zpos] << "\t ";
    cout << blob[offs_xrad] << "\t ";
    cout << blob[offs_yrad] << "\t ";
    cout << blob[offs_zrad] << "\t ";
    cout << blob[offs_xrot] << "\t ";
    cout << blob[offs_yrot] << "\t ";
    cout << blob[offs_zrot] << endl;
    i++;
  }
 cout << "---------------------------------------------------------------------------" << endl;

}
