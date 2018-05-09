/*
    evaluates the mismatch between a target file and a model made by the parameters that follow.

    compile: g++ -O3 -o evalFitness evalFitness.cpp functions.cpp `pkg-config --cflags --libs opencv`
*/

#include "functions.h"
#include <fstream>
#include <string>

int main(int argc, char *argv[])
{
    if(argc < 3){
        cout << "usage: ./evalFitness3D [target_model_name] [background_colour] [optional blob parameters]" << endl;
        exit(1);
    } 

    double bgcolour = atof(argv[2]);
    
    // declare both files
    ifstream targFile;

    // placeholder variables
    string dummy;
    double d;
    // parameters for target
    vector<double> targParams;
    // open up and read the target file into a vector
    targFile.open(argv[1]);
    // skip the first line
    getline(targFile,dummy);
    
    // number of intervals (derived from target file)
    int numXInts;
    int numYInts;
    int numZInts;
    // actual intervals (derived from target file)
    vector<double> xInts;
    vector<double> yInts;
    vector<double> zInts;

    // int model dims in x, y and z dims.
    double  xDim=0;
    double  yDim=0;
    double  zDim=0;

    // read in the interval numbers in each dimension
    targFile>>numXInts;
    targFile>>numYInts;
    targFile>>numZInts;

    //cout << " number of samples: " << numXInts << " " << numYInts << " " << numZInts << " " << numXInts*numYInts*numZInts << endl;

    getline(targFile,dummy); // skip to end of the line

    // read in the x, y and z intervals from target file
    // we will use these later to make the values in an model
    // to compare to the target

    for(int i=0;i<numXInts;i++){
      targFile>>d;
      xInts.push_back(d);
      xDim+=d;
    }

    for(int i=0; i<numYInts;i++){
      targFile>>d;
      yInts.push_back(d);
      yDim+=d;
    }
    
    for(int i=0; i<numZInts;i++){
      targFile>>d;
      zInts.push_back(d);
      //cout << d << " ";
      zDim+=d;
    }

    // calculate model volume and number of cells
    // and mean cell size.
    double modelVolume=xDim*yDim*zDim;
    int  numCells=numXInts*numYInts*numZInts;
    double meanCellSize=modelVolume/numCells;
    //cout << meanCellSize << " " << modelVolume << " " << xDim << " " << yDim << " " << zDim << endl;
    // hop sizes for converting an index i into an index in the x, y and z dims
    int  bigHopSize=numXInts*numYInts;
    int smallHopSize=numXInts;
    
    // now read the rest of the model into a vector for comparison
    // will push an extra copy of the last sample onto targParams.
    // ok as long as we don't use this in our comparisons!
    double t;
    while(targFile){
      targFile>>t;
      targParams.push_back(t);
      //cout << t << endl;
    }
    targFile.close();


    vector<double> modelParams 
      = render3DModel(xInts,yInts,zInts,
		      processInput3D( extractParam(argc - 3, argv + 3)),
		      bgcolour);

    // now read both sets of model parameters and compare 
    double o;
    double sumSquared=0.0;
    //cout << "sizes " << modelParams.size() << " " << targParams.size() << endl;
    //cout << "last parm of targParams " << targParams[targParams.size()-1] << endl;
    for(u_int i=0;i<modelParams.size();i++){
      t=targParams[i];
      o=modelParams[i];
      int zIndex=i/bigHopSize;
      int yIndex=(i%bigHopSize)/smallHopSize;
      int xIndex=(i%bigHopSize)%smallHopSize;
      double xSize=xInts[xIndex];
      double ySize=yInts[yIndex];
      double zSize=zInts[zIndex];
      double cellSize=xSize*ySize*zSize;

      cellSize=1; // comment out to activate cell size weighting
      meanCellSize=1; // comment out to activate cell size weighting
      //cout << "cell size " << cellSize <<  " " << zIndex << " " << yIndex << " " << xIndex <<endl;
      sumSquared+=pow(log10(t)-log10(o),2.0)*(cellSize/meanCellSize);
      //cout<< o << " " << t <<  " sumsq: " << sumSquared << endl;
    }

    // print out the RMS error
    cout.precision(12);
    cout << fixed;
    cout << pow(sumSquared/modelParams.size(),0.5) << endl;
    
    
}
