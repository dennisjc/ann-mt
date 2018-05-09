#include "functions.h"

//extracts blob parameters from the command line
vector<double> extractParam(int argc, char *argv[])
{
    vector<double> input;
    for(int i = 0; i < argc; i++) {
        char* token = strtok(argv[i],",");
        while(token != NULL) {
            input.push_back(atof(token));
            token = strtok(NULL,", ");
        }
    }
    return input;
}

// splits and normalises the blob parameters 
// returns a list of normalised blob parameters
vector< vector<double> > processInput(vector<double> input)
{
    vector< vector<double> > params = splitlist(input);
    vector< vector<double> > newparams;
    double max_strength = 0.0000001; // this variable is used to normalise strength values

    // normalise the parameters
    for(int i = 0; i < (int)params.size(); i++) {

        double *param = &params[i][0];

        //check for illegal parameters
        for(int j = 0; j < GENE_SIZE; j++) {
            if(param[j] < 0.0 || param[j] > 1.0) {
                cerr << "invalid parameter at blob" << i << "[" << j << "] = " << param[j] << endl;
                param[j] = (param[j] < 0) ? 0 : param[j];
                param[j] = (param[j] > 1) ? 1 : param[j];
                cerr << "readjusted parameter to " << param[j] << endl;
                cerr << endl;
            }
        }
        if(max_strength < param[offs_stren]) max_strength = param[offs_stren];

        vector<double> newparam;
        newparam.push_back(param[offs_inten]*MAX_INTENSITY);
        if(param[offs_stren] > 0.0000001) newparam.push_back(param[offs_stren]);
        else newparam.push_back(0.0000001);
        newparam.push_back(pow(param[offs_atten], 2.0)*ATTEN_MULTIPLIER);
        // modify xpos and ypos from the initial range of 0.0 ~ 1.0 to a range of -1.0 ~ 1.0
        double xpos = 2 * (param[offs_xpos] - 0.5);
        double ypos = 2 * (param[offs_ypos] - 0.5);
        newparam.push_back(xpos*(MAX_POS/2.0) + MAX_POS/2.0);
        newparam.push_back(ypos*(MAX_POS/2.0) + MAX_POS/2.0);
        newparam.push_back(param[offs_xrad]*MAX_RAD);
        newparam.push_back(param[offs_yrad]*MAX_RAD);
        double theta = param[offs_theta]*M_PI;
        newparam.push_back(sin(theta));
        newparam.push_back(cos(theta));

        newparams.push_back(newparam);
    }

    /*// debug
    for(int i = 0; i < (int)newparams.size(); i++) {
        for(int j = 0; j < (int)newparams[i].size(); j++)
            cout << newparams[i][j] << " ";
        cout << endl;
    }//*/

    return newparams;
}


// splits and normalises the blob parameters for 3D
// because this is for MT models some of the normalisation factors
// will be different
// returns a list of normalised blob parameters
vector< vector<double> > processInput3D(vector<double> input)
{
    vector< vector<double> > params = splitlist3D(input);
    vector< vector<double> > newparams;
    double max_strength = 0.0000001; // this variable is used to normalise strength values

    // normalise the parameters
    for(int i = 0; i < (int)params.size(); i++) {

        double *param = &params[i][0];

        //check for illegal parameters
        for(int j = 0; j < GENE_SIZE; j++) {
            if(param[j] < 0.0 || param[j] > 1.0) {
                cerr << "invalid parameter at blob" << i << "[" << j << "] = " << param[j] << endl;
                param[j] = (param[j] < 0) ? 0 : param[j];
                param[j] = (param[j] > 1) ? 1 : param[j];
                cerr << "readjusted parameter to " << param[j] << endl;
                cerr << endl;
            }
        }
        if(max_strength < param[offs_stren]) max_strength = param[offs_stren];

        vector<double> newparam;
	//cout << "intensity raw" << param[offs_inten] << endl;
        newparam.push_back(param[offs_inten]*MAX_INTENSITY);
        if(param[offs_stren] > 0.0000001) newparam.push_back(param[offs_stren]);
        else newparam.push_back(0.0000001);
        newparam.push_back(pow(param[offs_atten], 2.0)*ATTEN_MULTIPLIER_3D);
        newparam.push_back(param[offs_xpos]);
        newparam.push_back(param[offs_ypos]);
        newparam.push_back(param[offs_xrad]);
        newparam.push_back(param[offs_yrad]);
        newparam.push_back(param[offs_zrot]*M_PI);
	newparam.push_back(param[offs_zpos]);
	newparam.push_back(param[offs_zrad]);
	newparam.push_back(param[offs_xrot]*M_PI);
	newparam.push_back(param[offs_yrot]*M_PI);
	//cout << " processed intensity " << newparam[offs_inten] << endl;
        newparams.push_back(newparam);
    }

    /*// debug
    for(int i = 0; i < (int)newparams.size(); i++) {
        for(int j = 0; j < (int)newparams[i].size(); j++)
            cout << newparams[i][j] << " ";
        cout << endl;
	}// */

    return newparams; 
}

// splits the blob parameters into a list for easier access
vector< vector<double> > splitlistGeneric(vector<double> input,int geneSize)
{
    if(input.size()%(int)geneSize != 0) {
        cerr << "invalid blob parameters!" << endl;
        exit(1);
    }
    vector< vector<double> > params;

    for(int i = 0; i < (int)input.size(); i += geneSize) {
        // insert blob in the order of strength, strongest blobs are placed last in the list
        vector<double> param(&input[i], &input[i + geneSize]);
        int index = 0;
        while(index < (int)params.size() && params[index][offs_stren] < param[offs_stren]) index++;
        params.insert(params.begin() + index, param);
    }

    // debug
    /*for(int i = 0; i < (int)params.size(); i++) {
        for(int j = 0; j < (int)params[i].size(); j++)
            cout << params[i][j] << " ";
        cout << endl;
    }//*/

    return params;
}

// splits the blob parameters into a list for easier access to 2D genomes
vector< vector<double> > splitlist(vector<double> input)
{
  return splitlistGeneric(input,GENE_SIZE);
}

// splits the blob parameters into a list for easier access to 3D genomes
vector< vector<double> > splitlist3D(vector<double> input)
{
  return splitlistGeneric(input,GENE_SIZE3D);
}

// range of background_inten is [0.0, 1.0]
void render(IplImage* img, double background_inten, vector< vector<double> > params)
{
   
    if(background_inten > 1.0 || background_inten < 0.0) {
        cerr << "Invalid background intensity!" << endl;
        exit(1);
    }

    vector<double> res;
    double bgcol = background_inten * MAX_INTENSITY;

    // get image data    
    int height,width,step;
    uchar *data;
    height = img->height;
    width = img->width;
    step = img->widthStep;
    data = (uchar *)img->imageData;

    // collect data for the combine function
    vector<double> strengths, intensities;
    for(int i = 0; i < (int)params.size(); i++) {
        strengths.push_back(params[i][offs_stren]);
        intensities.push_back(params[i][offs_inten]);
    }

    for(int row = 0; row < height; row++) { // height
        for(int col = 0; col < width; col++) { //width
            vector<double> influences;
            for(int i = 0; i < (int)params.size(); i++) {
                influences.push_back(getBlobInfluence(row, col, params[i]));
            }
            data[col*step + row] = combine2(influences, strengths, intensities, bgcol);
        }
    }
}

vector<double> render3DModel(vector<double> xIntervals,
                   vector<double> yIntervals,
		   vector<double> zIntervals,
		   vector< vector<double> > params,
		   double backgroundIntensity){

  double bgColour = backgroundIntensity * MAX_INTENSITY;
  vector<double> res;

  //cout.precision(2); // sets precision to two digits
  // cout<< fixed;      // makes cout in fixed decimal format rather than the 
                     // default scientific notation

  // renders blobs as expressed by params into a 3D MT model as 
  // specified by intervals
   
  // set up Xsamples, Ysamples and zSamples given respective intervals
  double xStartInterval=0;
  double xEndInterval;
  vector<double> xSamples;
  for(unsigned int i=0; i<xIntervals.size();i++){
    xEndInterval=xIntervals[i]+xStartInterval;
    xSamples.push_back((xStartInterval+xEndInterval)/2.0);
    xStartInterval=xEndInterval;
  }

  // Y samples
  double yStartInterval=0;
  double yEndInterval;
  vector<double> ySamples;
  for(unsigned int i=0; i<yIntervals.size();i++){
    yEndInterval=yIntervals[i]+yStartInterval;
    ySamples.push_back((yStartInterval+yEndInterval)/2.0);
    yStartInterval=yEndInterval;
  }
  
  // Z samples
  double zStartInterval=0;
  double zEndInterval;
  vector<double> zSamples;
  for(unsigned int i=0; i<zIntervals.size();i++){
    zEndInterval=zIntervals[i]+zStartInterval;
    zSamples.push_back((zStartInterval+zEndInterval)/2.0);
    zStartInterval=zEndInterval;
  }

  // find max in all dimensions

  double maxExtent = max(xEndInterval,max(yEndInterval,zEndInterval));

 
  //cout << "xSamples" << xSamples << endl;


  // find the maximum extent in all dimensions

  // print out the file preamble
  // first line - just a placeholder
  //cout << "#Artificial Data" << endl;

  // second line - number of intervals in each dimension
  //cout << " " << xSamples.size() << " " << ySamples.size() << " " << 
  //  zSamples.size() << " " << "0" << endl;
  /*
  cout.precision(8);
  cout << fixed;
  // print sample Intervals
  printSampleIntervals(xIntervals);
  printSampleIntervals(yIntervals);
  printSampleIntervals(zIntervals);
  */
  // collect data for the combine function
  vector<double> strengths, intensities;
  for(int i = 0; i < (int)params.size(); i++) {
    strengths.push_back(params[i][offs_stren]);
    intensities.push_back(params[i][offs_inten]);
  }
  
  // iterate through the intervals sampling and printing out file
  /*for( auto z: zSamples){
    for(auto y: ySamples){
    for(auto x: xSamples){    // c++11 only I'm afraid.. */
  for(u_int i=0; i<zSamples.size();i++){
    double z=zSamples[i];
    for(u_int j=0; j<ySamples.size(); j++){
      double y=ySamples[j];
      for(u_int k=0; k<xSamples.size(); k++){
	double x=xSamples[k];
	//cout << "(" << x/maxExtent << " " << 
	//  y/maxExtent << " " << 
	//  z/maxExtent << ")" << endl;
	vector<double> influences;
	for(int i = 0; i < (int)params.size(); i++) {
	  //cout << params[i][0] << " " << params[i][1] << endl;
	  // cout << x << " " << y << " " << z << endl;

	  influences.push_back(getBlobInfluence3D(x/maxExtent,
						  y/maxExtent,
						  z/maxExtent,params[i]));
	  //cout << "influence " << influences[i] << endl;
	  //cout << "influence " << getBlobInfluence3D(x/maxExtent,
	  // 					  y/maxExtent,
	  //			       	     z/maxExtent,params[i]) << " ";
	}
        // color in range 0.0 to MAX_INTENSITY [0.0,255.0]
	double colour = combine2(influences, strengths, intensities, bgColour);
        
	// convert to resitivity
	double resist = pow(10.0,
			    ((colour/INTENSE_RANGE)*RESIST_RANGE)+MIN_RESIST);
	res.push_back(resist);
	//cout << resist << endl;
      }
	
    }
  }
  return res;
}

void renderAndPrint3DModel(vector<double>xIntervals,
                   vector<double> yIntervals,
		   vector<double> zIntervals,
		   vector< vector<double> > params,
				   double bgColor){
  // calls render 3D model and prints the result
  cout << "#Artificial Data" << endl;

  // second line - number of intervals in each dimension
  cout.precision(2);
  cout<<fixed;
  cout << " " << xIntervals.size() << " " << yIntervals.size() << " " << 
    zIntervals.size() << " " << "0" << endl;
  // print sample Intervals
  printSampleIntervals(xIntervals);
  printSampleIntervals(yIntervals);
  printSampleIntervals(zIntervals);
  vector<double> res=render3DModel(xIntervals,
				   yIntervals,
				   zIntervals,
				   params,
				   bgColor);
  cout.precision(12);
  cout << fixed;
  for(u_int i=0; i<res.size(); i++){
    cout << res[i] << endl;
  }


}

// prints the intervals in samples, seven to a line
void printSampleIntervals(vector<double> samples){
  const int SAMPLES_PER_LINE=7;
  int numIntervalsRemaining = samples.size();
  while(numIntervalsRemaining>0){
    for(int i=0;i<SAMPLES_PER_LINE;i++){
      cout << samples[samples.size()-numIntervalsRemaining] << " ";
      numIntervalsRemaining--;
      if (numIntervalsRemaining<=0) break;
    }
    cout << endl;
  }
}


// combines the colors of all blobs according to influence and strength
// this function assumes that the list of values have been sorted according to strength
// in which the largest blob values are placed at the end of the list
double combine(vector<double> influences, vector<double> strengths, vector<double> intensities, double bgcol)
{
    // strength determines solely the ordering in which the blobs are painted
    // blobs with the same strength will be painted at the same time
    // the color of each blob will be determined by the influence of the blob
    double color = bgcol;
    cerr << "colour " << bgcol << endl;
    
    for(int i = 0; i < (int)intensities.size(); ) {
        vector<double> newcolors, newinfluences;
        double norminfluences_sum = 0.0;
        int j = i;
        // first find all the blobs with the same strength in this current iteration
        while(j < (int)intensities.size() && strengths[i] == strengths[j]) j++;
        while(i < j) {
	  cerr << "intensitiy: " << intensities[i] << endl;
            double difference = intensities[i] - color;
            double newcolor = color + difference*influences[i];
            newcolors.push_back(newcolor);
            newinfluences.push_back(influences[i]);
            // normalise influences such that the sum of the influences of all blobs with the same strength = 1.0;
	    cerr << "influence: " << influences[i] << endl;
            norminfluences_sum += influences[i];
            i++;
        }
        double newcolor = 0;
        for(int k = 0; k < (int)newcolors.size(); k++) {
            newinfluences[k] /= norminfluences_sum;
            newcolor += newcolors[k]*newinfluences[k];
        }
        color = newcolor;
    }
    if(color > 255.0) return 255.0;
    else if(color < 0.0) return 0.0;
    return color;
}

// combine function where order doesn't matter.
double combine2(vector<double> influences, vector<double> strengths, vector<double> intensities, double bgcol){

  // append background to influences, strengths, and intensities
  influences.push_back(1.0);  // intensity of foreground is uniform bgcol
  intensities.push_back(bgcol); // intensity of background is normalised bgcol
  strengths.push_back(0.0);   // strength of background is 0.0
  // increment each strength element by tiny amount
  // zeros are no fun.. replace with epsilon
  double epsilon=1.0e-80;
  double epsilon2=1.0e-2;
  for(int i = 0; i < (int)intensities.size(); i++) {
    intensities[i]/=255.0;
    intensities[i]=(intensities[i]<epsilon2)?epsilon2:intensities[i];
    intensities[i]=(intensities[i]>1.0)?1.0:intensities[i];
    
    // normalise influence with the max intensity of a blob
    // note that influence in this case is different from my old understanding
    // of influence which is blob intensity at that pixel. Influence in this case
    // is always between 0.0 and 1.0 and must be multiplied by max-intensity in
    // order to get the influence I want
    influences[i]*=intensities[i];
    //cout << " intensity " << intensities[i] << " norm_infl " << influences[i];
  }

  // add tiny epsilon value to weights.
  // Now to normalise strengths to range epslion  to 1.0
  // exaggerate differences between weights by raising to a high power

  // find max strength
  double max_strength = epsilon;
  for(int i=0;i< (int)strengths.size();i++)
    if (strengths[i]>max_strength)
      max_strength=strengths[i];
  // normalise each strength and exaggerate within range 0 to 1.0
  for(int i=0;i<(int)strengths.size();i++){
    strengths[i]+=epsilon;
    strengths[i]/=max_strength;
    strengths[i]=pow(strengths[i],6.0);
    //cout << strengths[i] << ' ';
    //cerr << "pow strengths[i]" << strengths[i] << endl;
  }
  
  // make it so all strengths add to one
  // find sum of strengths
  double sum_strength=0.0;
  for(int i=0;i<(int)strengths.size();i++){
    sum_strength+=strengths[i];
  }
  // now divide each by sum
  for(int i=0;i<(int)strengths.size();i++){
    strengths[i]/=sum_strength;
  }

  // further normalise strengths make it so that only blobs with 
  // significant influence at current point get full strength
  vector<double> norm_weights;
  double sum_normweights=0.0;
  for(int i=0;i<(int)strengths.size();i++){
    double norm_weight = strengths[i]*pow(influences[i]/intensities[i],6.0);
    //cerr << "strength,influence,intensities,norm_weight " <<  strengths[i] << " " << influences[i] << " " << intensities[i] << " " << norm_weight << endl;
    sum_normweights+=norm_weight;
    norm_weights.push_back(norm_weight);
  }
  double background_term=0.0;
  double blob_factor=0.0;
  double max_intensity=0.0;
  for(int i=0;i<(int)strengths.size();i++){
    // now make all normalised weights sum to one
    norm_weights[i]/=sum_normweights;
    // and calculate terms
    background_term+=(1.0-norm_weights[i])*influences[i];
    blob_factor+=norm_weights[i]*influences[i];
    max_intensity+=norm_weights[i]*intensities[i];
  }
  // calculate combined color
  double combined_colour = background_term+blob_factor*(1-(background_term/max_intensity));
  return combined_colour*255.0;
}


// returns the influence of the specified blob at (row, col)
double getBlobInfluence(int row, int col, vector<double> param)
{
    double atten = param[offs_atten];
    double xpos = param[offs_xpos];
    double ypos = param[offs_ypos];
    double xrad = param[offs_xrad];
    double yrad = param[offs_yrad];
    double sin = param[offs_sin];
    double cos = param[offs_cos];

    if(xrad == 0 || yrad == 0) return 0.0;

    double dx = row - xpos;
    double dy = col - ypos;

    double dxprime = ( dx*cos - dy*sin )/xrad;
    double dyprime = ( dx*sin + dy*cos )/yrad;
    double distance = pow(dxprime*dxprime + dyprime*dyprime, 0.5);

    double result = 1.0/(pow(distance,atten) + 1);

    return result;
}


double getBlobInfluence3D(double x, double y, double z, vector<double> param)
{
    double atten = param[offs_atten];
    double xpos = param[offs_xpos];
    double ypos = param[offs_ypos];
    double zpos = param[offs_zpos];
    double xrad = param[offs_xrad];
    double yrad = param[offs_yrad];
    double zrad = param[offs_zrad];
    
    double xrot = param[offs_xrot];
    double yrot = param[offs_yrot];
    double zrot = param[offs_zrot];

    if(xrad == 0 || yrad == 0 || zrad == 0) return 0.0;

    double dx = x - xpos;
    double dy = y - ypos;
    double dz = z - zpos;
 
    // rotate about z axis

    double dxRz = ( dx*cos(zrot) - dy*sin(zrot));
    double dyRz = ( dx*sin(zrot) + dy*cos(zrot));
    double dzRz = dz;

    // rotate about the y axis
    double dxRzRy = (dzRz*sin(yrot) + dxRz*cos(yrot));
    double dyRzRy = dyRz;
    double dzRzRy = (dzRz*cos(yrot) - dxRz*sin(yrot));

    // rotate about the x axis
    double dxRzRyRx = dxRzRy;
    double dyRzRyRx = (dyRzRy*cos(xrot) - dzRzRy*sin(xrot));
    double dzRzRyRx = (dyRzRy*sin(xrot) + dzRzRy*cos(xrot));

    // scale in x y and z directions
    double dxRzRyRxScaled = dxRzRyRx/xrad;
    double dyRzRyRxScaled = dyRzRyRx/yrad;
    double dzRzRyRxScaled = dzRzRyRx/zrad;


    double distance = pow(pow(dxRzRyRxScaled,2.0)+
                          pow(dyRzRyRxScaled,2.0)+
                          pow(dzRzRyRxScaled,2.0), 0.5);

    double result = 1.0/(pow(distance,atten) + 1);

    return result;
}


double fitness(IplImage* img1, IplImage* img2)
{
    int height = img1->height;
    int width = img1->width;
    int step = img1->widthStep;

    /*// check for potential errors
    if(height != img2->height || width != img2->width || step != img2->widthStep) {
        if(step != img2->widthStep)
            cerr << "unable to compare images of different types (inten/greyscale)" << endl;
        if(height != img2->height || width != img2->width)
            cerr << "unable to compare images of different sizes" << endl;
        exit(1);
    }//*/

    uchar *data1, *data2;
    data1 = (uchar *)img1->imageData;
    data2 = (uchar *)img2->imageData;
    
    // calculate error squared
    double err_squared = 0;
    for(int row=0; row < height; row++) {
        for(int col=0; col < width; col++) {
            err_squared += pow((double)data1[col*step + row] - (double)data2[col*step + row], 2);
        }
    }
    // calculate mean
    double pixelCount = height * width;
    double mean = err_squared/pixelCount;

    // calculate root
    double root = pow(mean, 0.5);

    // return rms error
    return root;
}

