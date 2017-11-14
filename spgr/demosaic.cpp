/* --------------------------------------------------------------------------
 * File:    demosaic.cpp
 * Created: 2015-10-01
 * --------------------------------------------------------------------------
 *
 *
 *
 * ------------------------------------------------------------------------*/


#include "demosaic.h"
#include <cmath>

using namespace std;


Image basicGreen(const Image &raw, int offset) {
    // --------- HANDOUT  PS03 ------------------------------
    // Takes as input a raw image and returns a single-channel
    // 2D image corresponding to the green channel using simple interpolation
    Image output(raw.width(), raw.height(), 1);

    for (int i=0; i < raw.width(); i++) {
        for (int j=0; j < raw.height(); j++) {
            if (i == 0 || j == 0 || i == (raw.width()-1) || j == (raw.height()-1)) {
                // copy over pixel values from the raw image for the first and last rows and columns
                output(i,j) = raw(i,j); 
            } else {
                if ((i + j%2) % 2 == offset) {
                    output(i,j) = raw(i,j);
                } else {
                    float neighbors = raw(i-1, j) + raw(i, j+1) + raw(i+1, j) + raw(i, j-1);
                    output(i,j) = neighbors/4.0;
                }
            }

        }
    }
    return output;
}


Image basicRorB(const Image &raw, int offsetX, int offsetY) {
    // --------- HANDOUT  PS03 ------------------------------
    //  Takes as input a raw image and returns a single-channel
    //  2D image corresponding to the red or blue channel using simple interpolation
    
    Image output = Image(raw.width(), raw.height(), 1);

    for (int i=0; i < raw.width(); i++) {
        for (int j=0; j < raw.height(); j++) {
            if (i == 0 || j == 0 || i == (raw.width()-1) || j == (raw.height()-1)) {
                // copy over pixel values from the raw image for the first and last rows and columns
                output(i,j) = raw(i,j); 
            } else {
                if (i%2 != offsetX && j%2 != offsetY) {
                    float neighbors = raw(i-1, j-1) + raw(i-1, j+1) + raw(i+1, j+1) + raw(i+1, j-1);
                    output(i,j) = neighbors/4.0;
                } else if (i%2 == offsetX && j%2 != offsetY) {
                    output(i,j) = (raw(i,j-1) + raw(i,j+1))/2.0;
                } else if (i%2 != offsetX && j%2 == offsetY) {
                    output(i,j) = (raw(i-1,j) + raw(i+1,j))/2.0;
                } else {
                    output(i,j) = raw(i,j);
                }
            }
        }
    }

    return output;

}

Image basicDemosaic(const Image &raw, int offsetGreen,
                    int offsetRedX, int offsetRedY, int offsetBlueX, int offsetBlueY) {
    // --------- HANDOUT  PS03 ------------------------------
    // takes as input a raw image and returns an rgb image
    // using simple interpolation to demosaic each of the channels
    Image green = basicGreen(raw, offsetGreen);
    Image red = basicRorB(raw, offsetRedX, offsetRedY);
    Image blue = basicRorB(raw, offsetBlueX, offsetBlueY); 

    Image output = Image(green.width(), green.height(), 3);

    for (int i=0; i < output.width(); i++) {
        for (int j=0; j < output.height(); j++) {
            output(i,j,0) = red(i,j);
            output(i,j,1) = green(i,j);
            output(i,j,2) = blue(i,j);
        }
    }

    return output;
}

Image edgeBasedGreen(const Image &raw, int offset){
    // --------- HANDOUT  PS03 ------------------------------
    // Takes a raw image and outputs a single-channel
    // image corresponding to the green channel taking into account edges
    Image output(raw.width(), raw.height(), 1);

    for (int i=0; i < raw.width(); i++) {
        for (int j=0; j < raw.height(); j++) {
            if ( i == 0 || j == 0 || i == (raw.width()-1) || j == (raw.height()-1) ) {
                output(i,j) = raw(i,j);
            } else {
                if ( (i + j%2) % 2 == offset ) {
                    output(i,j) = raw(i,j);
                } else {
                    float horizontal = abs(raw(i-1,j) - raw(i+1,j));
                    float vertical = abs(raw(i,j-1) - raw(i,j+1));

                    if (horizontal < vertical) {
                        output(i,j) = (raw(i-1,j) + raw(i+1,j)) / 2.0; 
                    } else {
                        output(i,j) = (raw(i,j-1) + raw(i,j+1)) / 2.0;
                    }
                }
            }
        }
    }

    return output;

}

Image edgeBasedGreenDemosaic(const Image &raw, int offsetGreen,
                             int offsetRedX, int offsetRedY,
                             int offsetBlueX, int offsetBlueY) {
    // --------- HANDOUT  PS03 ------------------------------
    // Takes as input a raw image and returns an rgb image
    // using edge-based green demosaicing for the green channel and
    // simple interpolation to demosaic the red and blue channels
    Image green = edgeBasedGreen(raw, offsetGreen);
    Image red = basicRorB(raw, offsetRedX, offsetRedY);
    Image blue = basicRorB(raw, offsetBlueX, offsetBlueY);

    Image output = Image(green.width(), green.height(), 3);

    for (int i=0; i < output.width(); i++) {
        for (int j=0; j < output.height(); j++) {
            output(i,j,0) = red(i,j);
            output(i,j,1) = green(i,j);
            output(i,j,2) = blue(i,j);
        }
    }

    return output;

}


Image greenBasedRorB(const Image &raw, Image &green, int offsetX, int offsetY){
    // --------- HANDOUT  PS03 ------------------------------
    // Takes as input a raw image and returns a single-channel
    // 2D image corresponding to the red or blue channel using green based interpolation
    
    Image output = Image(raw.width(), raw.height(), 1);

    for (int i=0; i < raw.width(); i++) {
        for (int j=0; j < raw.height(); j++) {
            if ( i == 0 || j == 0 || i == (raw.width()-1) || j == (raw.height()-1) ) {
                output(i,j) = raw(i,j);
            } else {
                if (i%2 != offsetX && j%2 != offsetY) {
                    float neighbors = raw(i-1,j-1) - green(i-1,j-1) 
                                        + raw(i-1,j+1) - green(i-1,j+1) 
                                        + raw(i+1,j+1) - green(i+1,j+1)
                                        + raw(i+1,j-1) - green(i+1,j-1);
                    output(i,j) = neighbors/4.0 + green(i,j);
                } else if (i%2 == offsetX && j%2 != offsetY) {
                    float up_down = raw(i,j-1) - green(i,j-1) + raw(i,j+1) - green(i,j+1);
                    output(i,j) = up_down/2.0 + green(i,j);
                } else if (i%2 != offsetX && j%2 == offsetY) {
                    float left_right = raw(i-1,j) - green(i-1,j) + raw(i+1,j) - green(i+1,j);
                    output(i,j) = left_right/2.0 + green(i,j);
                } else {
                    output(i,j) = raw(i,j);
                }
            }
        }
    }

    return output;
}

Image improvedDemosaic(const Image &raw, int offsetGreen,
                       int offsetRedX, int offsetRedY,
                       int offsetBlueX, int offsetBlueY) {
    // // --------- HANDOUT  PS03 ------------------------------
    // Takes as input a raw image and returns an rgb image
    // using edge-based green demosaicing for the green channel and
    // simple green based demosaicing of the red and blue channels
    Image green = edgeBasedGreen(raw, offsetGreen);
    Image red = greenBasedRorB(raw, green, offsetRedX, offsetRedY);
    Image blue = greenBasedRorB(raw, green, offsetBlueX, offsetBlueY);

    Image output = Image(green.width(), green.height(), 3);

    for (int i=0; i < output.width(); i++) {
        for (int j=0; j < output.height(); j++) {
            output(i,j,0) = red(i,j);
            output(i,j,1) = green(i,j);
            output(i,j,2) = blue(i,j);
        }
    }

    return output;    

}
