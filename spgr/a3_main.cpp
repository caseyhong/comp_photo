/* -----------------------------------------------------------------
 * File:    a3_main.cpp
 * Author:  Michael Gharbi <gharbi@mit.edu>
 * Created: 2015-09-30
 * -----------------------------------------------------------------
 *
 *
 *
 * ---------------------------------------------------------------*/


#include "Image.h"
#include "basicImageManipulation.h"
#include "demosaic.h"
#include "align.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>


using namespace std;


// This is a way for you to test your functions.
// We will only grade the contents of demosaic.cpp and align.cpp

void testDenoiseSeq() {
    cout << "starting denoiseSeq..." << endl;
    vector<Image> seq;
    for (int i=1; i < 19; i++) {
        // Image im = Image("./Input/aligned-ISO3200/1D2N-iso3200-" + to_string(i) + ".png");
        Image im = Image("./Input/green/noise-small-" + to_string(i) + ".png");
        seq.push_back(im);
    }
    Image output = denoiseSeq(seq);
    output.write("./Output/denoiseSeq_green.png");
}

void testLogSNR() {
    cout << "starting log SNR..." << endl;
    cout << "testing on ISO 3200..." << endl;
    vector<Image> seq_3200;
    for (int i=1; i < 25; i++) {
        Image im = Image("./Input/aligned-ISO3200/1D2N-iso3200-" + to_string(i) + ".png");
        seq_3200.push_back(im);
    }

    Image snr_3200 = logSNR(seq_3200);
    snr_3200.write("./Output/logsnr_3200.png");

    cout << "testing on ISO 400..." << endl;
    vector<Image> seq_400;
    for (int i=1; i < 25; i++) {
        Image im = Image("./Input/aligned-ISO400/1D2N-iso400-under-" + to_string(i) + ".png");
        seq_400.push_back(im);
    }

    Image snr_400 = logSNR(seq_400);
    snr_400.write("./Output/logsnr_400.png");
}

void testAlignment() {
    cout << "starting alignment..." << endl;
    vector<Image> green;
    for (int i=1; i < 19; i++) {
        Image im = Image("./Input/green/noise-small-" + to_string(i) + ".png");
        green.push_back(im);
    }
    Image aligned = alignAndDenoise(green, 10);
    aligned.write("./Output/aligned_offset10.png");
}

void testBasicGreen() {
    cout << "starting basic green..." << endl;
    Image raw = Image("./Input/raw/signs-small.png");
    Image green = basicGreen(raw, 1);
    green.write("./Output/basic_green.png");

}

void testBasicRorB() {
    cout << "starting basic red or green..." << endl;
    Image raw = Image("./Input/raw/signs-small.png");
    Image red = basicRorB(raw, 1, 1);
    Image blue = basicRorB(raw, 0, 0);

    red.write("./Output/basic_red_v2.png");
    red.write("./Output/basic_green_v2.png");
}

void testBasicDemosaic() {
    cout << "starting basic demosaic..." << endl;
    Image raw = Image("./Input/raw/signs-small.png");
    Image rgb = basicDemosaic(raw, 1, 1, 1, 0, 0);
    rgb.write("./Output/basic_demosaic_v2.png");
}

void testEdgeGreen() {
    cout << "starting edge based green..." << endl;
    Image raw = Image("./Input/raw/signs-small.png");
    Image green = edgeBasedGreen(raw, 1);
    green.write("./Output/edge_green.png");
}

void testEdgeGreenDemosaic() {
    cout << "starting edge based green demosaic..." << endl;
    Image raw = Image("./Input/raw/signs-small.png");
    Image rgb = edgeBasedGreenDemosaic(raw, 1, 1, 1, 0, 0);
    rgb.write("./Output/edge_demosaic_v2.png");
}

void testImprovedDemosaic() {
    cout << "starting improved demosaic..." << endl;
    Image raw = Image("./Input/raw/signs-small.png");
    Image rgb = improvedDemosaic(raw, 1, 1, 1, 0, 0);
    rgb.write("./Output/improved_demosaic_v3.png");
}

void testSergey(string file_path, string pid) {
    Image sergeyImg(file_path);
    Image naive_rgb = split(sergeyImg);
    naive_rgb.write("./Output/sergey_naive_" + pid + ".png");
    Image aligned_rgb = sergeyRGB(sergeyImg, 10);
    aligned_rgb.write("./Output/sergey_aligned_" + pid + ".png");
}

int main() {
    testDenoiseSeq();

    // testLogSNR();

    // testAlignment();

    // testBasicGreen();

    // testBasicRorB();

    // testBasicDemosaic();

    // testEdgeGreen();

    // testEdgeGreenDemosaic();

    // testImprovedDemosaic();

    // testSergey("./Input/Sergey/00106v_third.png", "00106");
    // testSergey("./Input/Sergey/00137v_third.png", "00137");
    // testSergey("./Input/Sergey/00757v_third.png", "00757");
    // testSergey("./Input/Sergey/00888v_third.png", "00888");
    // testSergey("./Input/Sergey/00889v_third.png", "00889");
    // testSergey("./Input/Sergey/00907v_third.png", "00907");
    // testSergey("./Input/Sergey/00911v_third.png", "00911");
    // testSergey("./Input/Sergey/01031v_third.png", "01031");
    // testSergey("./Input/Sergey/01880v_third.png", "01880");

    // cout << "nothing done in a3_main.cpp, debug me !" << endl;

    // Image im1 = Image(100,100,1);
    // Image im2 = Image(100,100,1);
    // im1.create_rectangle(40,40,60,50, 1.0f,1.0f,1.0f);
    // im2.create_rectangle(50,60,70,70, 1.0f,1.0f,1.0f);
    // im1.write("./Input/im1_rectangle.png");
    // im2.write("./Input/im2_rectangle.png");

    // std::vector<int> offset = align(im1, im2, 20);
    // cout << "horizontal offset: " << offset[0] << endl;
    // cout << "vertical offset: " << offset[1] << endl;
    
    

    // // Denoise ---------------------------
    // // Load sequence
    // vector<Image> seq;
    // int n_images = 4;
    // for (int i = 1; i <= n_images; ++i) {
    //     ostringstream fname;
    //     // fname << "./Input/aligned-ISO400/1D2N-iso400-under-";
    //     fname << "./Input/aligned-ISO3200/1D2N-iso3200-";
    //     fname << i;
    //     fname << ".png";
    //     seq.push_back(Image(fname.str()));
    // }
    //
    // // Denoise
    // Image out = denoiseSeq(seq);
    // out.write("./Output/denoised.png");
    //
    // Image SNRIm = logSNR(seq,1/30.0);
    // SNRIm.write("./Output/snr_map.png");
    //
    // // Demosaic ---------------------------
    // Image raw("./Input/raw/signs-small.png");
    // Image green = basicGreen(raw, 1);
    // green.write("./Output/demosaic_green.png");
    // Image red = basicRorB(raw, 1, 1);
    // red.write("./Output/demosaic_red.png");
    // Image blue = basicRorB(raw, 0, 0);
    // blue.write("./Output/demosaic_blue.png");
    // Image rgb = basicDemosaic(raw, 1, 1,1,0,0);
    // rgb.write("./Output/demosaiced.png");
    //
    //
    // // Sergey ---------------------------
    // Image sergeyImg("./Input/Sergey/00911v_third.png");
    // Image rgb2 = split(sergeyImg);
    // rgb2.write("./Output/Sergey_split.png");
    // Image rgbAlign = sergeyRGB(sergeyImg,10);
    // rgbAlign.write("./Output/Sergey_aligned.png");

    return 0;
}

