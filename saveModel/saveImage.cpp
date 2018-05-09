/*
    renders, displays, and saves the image specified by the blob parameter

    compile: g++ -O3 -o saveImage saveImage.cpp functions.cpp `pkg-config --cflags --libs opencv`
*/

#include "functions.h"

int main(int argc, char *argv[])
{
    if(argc < 3){
        cout << "usage: ./saveImage [image_name.image_format] [background_colour] [optional blob parameters]" << endl;
        exit(1);
    }

    double bgcolour = atof(argv[2]);

    // render the image
    IplImage* img = cvCreateImage(cvSize(MAX_POS,MAX_POS),IPL_DEPTH_8U,1);
    render(img, bgcolour, processInput( extractParam(argc - 3, argv + 3) ));

    // save/export image 
    cvSaveImage(argv[1], img);
    
    /*cout << "Image saved. Loading preview..." << endl;

    // create a window, show the image, and wait for a keypress
    cvNamedWindow("mainWin", CV_WINDOW_AUTOSIZE); 
    cvMoveWindow("mainWin", 100, 100);
    cvShowImage("mainWin", img);
    cvWaitKey(0); //*/

    // release the image
    cvReleaseImage(&img);

    return 0;
}
