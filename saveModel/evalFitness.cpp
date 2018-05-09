/*
    evaluates the fitness between the given blob parameters and the test.png image

    compile: g++ -O3 -o evalFitness evalFitness.cpp functions.cpp `pkg-config --cflags --libs opencv`
*/

#include "functions.h"

int main(int argc, char *argv[])
{
    if(argc < 3){
        cout << "usage: ./evalFitness [target_name.image_format] [background_colour] [optional blob parameters]" << endl;
        exit(1);
    }

    double bgcolour = atof(argv[2]);

    // load an image in greyscale
    //IplImage* comparison = cvLoadImage("test.png", 0);
    IplImage* comparison = cvLoadImage(argv[1], 0);
    if(!comparison) {
        cerr << "unable to load image!" << endl;
        exit(1);
    }
    //cvSaveImage("img.png", comparison); return 0;

    // render an image
    IplImage* img = cvCreateImage(cvSize(MAX_POS,MAX_POS),IPL_DEPTH_8U,1);
    render(img, bgcolour, processInput( extractParam(argc - 3, argv + 3) ));
    
    // use the fitness function to compare the rendered image and the loaded image
    cout << fitness(img,comparison) << endl;

    // release the image
    cvReleaseImage(&img);
    cvReleaseImage(&comparison);

    return 0;
}
