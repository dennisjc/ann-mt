#include "functions.h"

int main(int argc, char *argv[])
{

  // test and sample the intensity of blobs at different locations.

  double tp[12]={1,1,2,0.5,0.5,0.15,0.1,0.7,0.5,0.2,0.5,0.9};
  vector<double> testparams(tp,tp+12);
  //testparams={1,1,2,0.5,0.5,0.15,0.1,0.7,0.5,0.2,0.5,0.9};
 
  cout.precision(2);

  for (double z=0;z<=1;z+=0.2){
    for(double y=0;y<=1;y+=0.1){
      for(double x=0;x<=1;x+=0.1){
	cout << getBlobInfluence3D(x,y,z,testparams) << " ";
      }
      cout << "\n";
    }
    cout << "\n";
  }

}
	
