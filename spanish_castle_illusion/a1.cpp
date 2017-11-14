/* -----------------------------------------------------------------
 * File:    a1.cpp
 * Created: 2015-09-15
 * -----------------------------------------------------------------
 *
 * Assignment 01
 *
 * ---------------------------------------------------------------*/


#include "a1.h"
#include <cmath>
#include <numeric>
using namespace std;

// Create a surprise image
Image create_special() {
    // // --------- HANDOUT  PS01 ------------------------------
    // create the image outlined in the handout
    Image im = Image(290, 150, 3);
    im.set_color(1.0, 1.0, 1.0);

    im.create_rectangle(0, 0, 31, 149, 0.64, 0.12, 0.20);
    im.create_rectangle(52, 0, 83, 102, 0.64, 0.12, 0.20);
    im.create_rectangle(103, 0, 134, 149, 0.64, 0.12, 0.20);
    im.create_rectangle(155, 0, 186, 30, 0.64, 0.12, 0.20);
    im.create_rectangle(155, 48, 186, 149, 0.55, 0.55, 0.55);
    im.create_rectangle(207, 0, 289, 30, 0.64, 0.12, 0.20);
    im.create_rectangle(207, 48, 238, 149, 0.64, 0.12, 0.20);

    return im;
}

// Change the brightness of the image
// const Image & means a reference to im will get passed to the function,
// but the compiler won't let you modify it within the function.
// So you will return a new image
Image brightness(const Image &im, float factor) {
    // --------- HANDOUT  PS01 ------------------------------
	// Image output(im.width(), im.height(), im.channels());
	// Modify image brightness
	// return output;

    Image output = im * factor;
    return output;
}


Image contrast(const Image &im, float factor, float midpoint) {
    // --------- HANDOUT  PS01 ------------------------------
    // Image output(im.width(), im.height(), im.channels());
    // Modify image contrast
    // return output;

    Image output = (im - midpoint) * factor + midpoint;   
    return output; 
}


Image color2gray(const Image &im, const std::vector<float> &weights) {
    // --------- HANDOUT  PS01 ------------------------------
    // Image output(im.width(), im.height(), 1);
    // Convert to grayscale
    Image output = Image(im.width(), im.height(), 1);

    for (int i=0; i < im.width(); i++) {
        for (int j=0; j < im.height(); j++) {
            output(i,j,0) = (weights[0]*im(i,j,0) + weights[1]*im(i,j,1) + weights[2]*im(i,j,2));
        }
    }

    return output;
}


// For this function, we want two outputs, a single channel luminance image
// and a three channel chrominance image. Return them in a vector with luminance first
std::vector<Image> lumiChromi(const Image &im) {
    // --------- HANDOUT  PS01 ------------------------------
    // Create the luminance image
    // Create the chrominance image
    // Create the output vector as (luminance, chrominance)

    // create the luminance image
	Image lumi = color2gray(im);

    // initialize the chrominance image
    Image chromi = im;

    for (int i=0; i < im.width(); i++) {
        for (int j=0; j < im.height(); j++) {
            for (int k=0; k < im.channels(); k++) {
                if (lumi(i,j) == 0.0f) {
                    chromi(i,j,k) = 0.0f;
                } else {
                    chromi(i,j,k) = chromi(i,j,k)/lumi(i,j);
                }
            }
        }
    }

    std::vector<Image> output;
    output.push_back(lumi);
    output.push_back(chromi);

    return output;

}


// Modify brightness then contrast
Image brightnessContrastLumi(const Image &im, float brightF, float contrastF, float midpoint) {
    // --------- HANDOUT  PS01 ------------------------------
    // Modify brightness, then contrast of luminance image

    std::vector<Image> decomposed = lumiChromi(im);
    Image lumi = decomposed[0];
    Image chromi = decomposed[1];

    lumi = brightness(lumi, brightF);
    chromi = contrast(chromi, contrastF, midpoint);

    Image output = Image(im.width(), im.height(), im.channels());

    for (int i=0; i < im.width(); i++) {
        for (int j=0; j < im.height(); j++) {
            for (int k=0; k < im.channels(); k++) {
                output(i,j,k) = lumi(i,j) * chromi(i,j,k);
            }
        }
    }

    return output;
}


Image rgb2yuv(const Image &im) {
    // --------- HANDOUT  PS01 ------------------------------
    // Create output image of appropriate size
    // Change colorspace
    Image output = Image(im.width(), im.height(), im.channels());

    for (int i=0; i < im.width(); i++) {
        for (int j=0; j < im.height(); j++) {
            float r = im(i,j,0); 
            float g = im(i,j,1); 
            float b = im(i,j,2);

            output(i,j,0) = 0.299*r + 0.587*g + 0.114*b; // Y
            output(i,j,1) = -0.147*r - 0.289*g + 0.436*b; // U
            output(i,j,2) = 0.615*r - 0.515*g - 0.1*b; // V
        }
    }

    return output;
}


Image yuv2rgb(const Image &im) {
    // --------- HANDOUT  PS01 ------------------------------
    // Create output image of appropriate size
    // Change colorspace
    Image output = Image(im.width(), im.height(), 3);
    
    for (int i=0; i < im.width(); i++) {
        for (int j=0; j < im.height(); j++) {
            float y = im(i,j,0); 
            float u = im(i,j,1); 
            float v = im(i,j,2);

            output(i,j,0) = y + 1.14*v; // R
            output(i,j,1) = y - 0.395*u - 0.581*v; // G
            output(i,j,2) = y + 2.032*u; // B
        }
    }

    return output;
}


Image saturate(const Image &im, float factor) {
    // --------- HANDOUT  PS01 ------------------------------
    // Create output image of appropriate size
    // Saturate image
    // return output;
    Image output = rgb2yuv(im);

    for (int i=0; i < im.width(); i++) {
        for (int j=0; j < im.height(); j++) {
            output(i,j,1) = output(i,j,1) * factor;
            output(i,j,2) = output(i,j,2) * factor;
        }
    }

    return yuv2rgb(output);
}

// Gamma encodes the image
Image gamma_code(const Image &im, float gamma) {
    // // --------- HANDOUT  PS01 ------------------------------
    // Image output(im.width(), im.height(), im.channels());
    // Gamma encodes the image
    // return output;
    Image output = Image(im.width(), im.height(), im.channels());

    for (int i=0; i < im.width(); i++) {
        for (int j=0; j < im.height(); j++) {
            for (int k=0; k < im.channels(); k++) {
                output(i,j,k) = pow(im(i,j,k), 1.0/gamma);
            }
        }
    }

    return output;
}

// Quantizes the image to 2^bits levels and scales back to 0~1
Image quantize(const Image &im, int bits) {
    // // --------- HANDOUT  PS01 ------------------------------
    // Image output(im.width(), im.height(), im.channels());
    // Quantizes the image to 2^bits levels
    // return output;
    const int Q = pow(2, bits);
    Image output = im * Q;

    for (int i=0; i < im.width(); i++) {
        for (int j=0; j < im.height(); j++) {
            for (int k=0; k < im.channels(); k++) {
                output(i,j,k) = (float)round(output(i,j,k));        
            }
        }
    }

    return output/Q;
}

// Compare between first quantize then gamma_encode and first gamma_encode then quantize
std::vector<Image> gamma_test(const Image &im, int bits, float gamma) {
    // // --------- HANDOUT  PS01 ------------------------------
    // // im1 = quantize then gamma_encode the image
    // // im2 = gamma_encode then quantize the image
    // // Remember to create the output images and the output vector
    // // Push the images onto the vector
    // // Do all the required processing
    // // Return the vector, color image first
    
    std::vector<Image> output;

    Image qim1 = quantize(im, bits);
    Image im1 = gamma_code(qim1, gamma);

    Image gim2 = gamma_code(im, gamma);
    Image im2 = quantize(gim2, bits);

    output.push_back(im1);
    output.push_back(im2);

    return output;
}


// Return two images in a C++ vector
std::vector<Image> spanish(const Image &im) {
    // --------- HANDOUT  PS01 ------------------------------
    // Remember to create the output images and the output vector
    // Push the images onto the vector
    // Do all the required processing
    // Return the vector, color image first

    Image color_image = rgb2yuv(im);
    for (int i=0; i < im.width(); i++) {
        for (int j=0; j < im.height(); j++) {
            color_image(i,j,0) = 0.5;
            color_image(i,j,1) *= -1;
            color_image(i,j,2) *= -1;
        }
    }

    color_image = yuv2rgb(color_image);

    Image gray_image = color2gray(im);

    for (int c=0; c < im.channels(); c++) {
        color_image(floor(im.width()/2.0), floor(im.height()/2.0), c) = 0.0f;
    }
    gray_image(floor(im.width()/2.0), floor(im.height()/2.0)) = 0.0f;

    std::vector<Image> output;

    output.push_back(color_image);
    output.push_back(gray_image);

    return output;
}


// White balances an image using the gray world assumption
Image grayworld(const Image & im) {
    // --------- HANDOUT  PS01 ------------------------------
    // Implement automatic white balance by multiplying each channel
    // of the input by a factor such that the three channel of the output image
    // have the same mean value. The mean value of the green channel
    // is taken as reference.
    
    std::vector<float> r_values;
    std::vector<float> g_values;
    std::vector<float> b_values;

    for (int i=0; i < im.width(); i++) {
        for (int j=0; j < im.height(); j++) {
            r_values.push_back(im(i,j,0));
            g_values.push_back(im(i,j,1));
            b_values.push_back(im(i,j,2));
        }
    }

    float r_mean = std::accumulate(r_values.begin(), r_values.end(), 0.0f)/r_values.size();
    float g_mean = std::accumulate(g_values.begin(), g_values.end(), 0.0f)/g_values.size();
    float b_mean = std::accumulate(b_values.begin(), b_values.end(), 0.0f)/b_values.size();

    float r_factor = g_mean/r_mean;
    float b_factor = g_mean/b_mean;

    Image output = Image(im.width(), im.height(), im.channels());

    for (int i=0; i < im.width(); i++) {
        for (int j=0; j < im.height(); j++) {
            output(i,j,0) = im(i,j,0) * r_factor;
            output(i,j,1) = im(i,j,1);
            output(i,j,2) = im(i,j,2) * b_factor;
        }
    }

    return output;
}
