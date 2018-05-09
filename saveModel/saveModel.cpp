#include "functions.h"
#include <fstream>
#include <string>
int main(int argc, char *argv[])
{

  // usage: ./saveModel input_model_file [background_colour] [blob parameters] > output_model_file
  // input_model_file is just used for the model intervals in its preamble

  // interval info
  int numXInts;
  int numYInts;
  int numZInts;
  vector<double> xInts;
  vector<double> yInts;
  vector<double> zInts;
  // read in intervals

  string firstLine; // dummy holder for first line of file
  double d;         // holder for integer lengths

  // file
  ifstream myFile;
  myFile.open(argv[1]);
  getline(myFile,firstLine);

  // read in the interval numbers in each dimension
  myFile>>numXInts;
  myFile>>numYInts;
  myFile>>numZInts;
  getline(myFile,firstLine);  // skip possible dummy value at end of line 2
  // now on line 3

  for(int i=0;i<numXInts;i++){
    myFile>>d;
    xInts.push_back(d);
    //cout << d << " ";
  }
  //cout << endl;

  for(int i=0; i<numYInts;i++){
    myFile>>d;
    yInts.push_back(d);
    //cout << d << " ";
  }
  //cout << endl;
  
  for(int i=0; i<numZInts;i++){
    myFile>>d;
    zInts.push_back(d);
    //cout << d << " ";
  }
  //cout << endl;


    if(argc < 3){
        cout << "usage: ./saveModel input_model_file [background_colour] [blob parameters] > output_model_file" << endl;
        exit(1);
    }

    double bgcolour = atof(argv[2]);

    renderAndPrint3DModel(xInts,
		  yInts,
		  zInts,
		  processInput3D( extractParam(argc - 3, argv + 3)),
		  bgcolour);

}


