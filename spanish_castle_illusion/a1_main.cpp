#include "a1.h"
#include <iostream>

using namespace std;

// This is a way for you to test your functions. 
// We will only grade the contents of a1.cpp and Image.cpp
int main() {
 
    cout << "Testing create_special..." << endl;
    Image im1 = create_special();
    im1.write("./Output/create_special.png");

    cout << "Testing brightness..." << endl;
    Image boston = Image("./Input/Boston_underexposed.png");
    Image boston_brighter = brightness(boston, 2.0);
    Image boston_dimmer = brightness(boston, 0.2);
    boston_brighter.write("./Output/boston_underexposed_brighter.png");
    boston_dimmer.write("./Output/boston_underexposed_dimmer.png");

    cout << "Testing contrast..." << endl;
    Image boston_locontrast = contrast(boston, 0.5, 0.5);
    Image boston_hicontrast = contrast(boston, 2.0, 0.5);
    boston_locontrast.write("./Output/boston_locontrast.png");
    boston_hicontrast.write("./Output/boston_hicontrast.png");

    cout << "Testing color2gray..." << endl;
    Image castle = Image("./Input/castle_small.png");
    Image gray_castle = color2gray(castle);
    gray_castle.write("./Output/castle_gray.png");

    std::vector<float> uneven;
    uneven.push_back(0.0f);
    uneven.push_back(0.0f);
    uneven.push_back(1.0f);
    Image unevenly_gray_castle = color2gray(castle, uneven);
    unevenly_gray_castle.write("./Output/castle_unevenly_gray.png");

    cout << "Testing lumiChromi..." << endl;
    Image nature = Image("./Input/skies_and_trees.png");
    std::vector<Image> out = lumiChromi(nature);
    out[0].write("./Output/skies_and_trees_lumi.png");
    out[1].write("./Output/skies_and_trees_chromi.png");
    
    cout << "Testing brightnessContrastLumi..." << endl;
    Image bcl_boston = brightnessContrastLumi(boston, 2.0f, 2.0f, 0.5);
    Image bcl_nature = brightnessContrastLumi(nature, 2.0f, 2.0f, 0.5);
    bcl_boston.write("./Output/boston_bcl.png");
    bcl_nature.write("./Output/skies_and_trees_bcl.png");

    cout << "Testing rgb2yuv..." << endl;
    Image castle_yuv = rgb2yuv(castle);
    Image nature_yuv = rgb2yuv(nature);
    castle_yuv.write("./Output/castle_yuv.png");
    nature_yuv.write("./Output/skies_and_trees_yuv.png");

    cout << "Testing yuv2rgb..." << endl;
    Image castle_rgb = yuv2rgb(castle_yuv);
    Image nature_rgb = yuv2rgb(nature_yuv);
    castle_rgb.write("./Output/castle_rgb.png");
    nature_rgb.write("./Output/skies_and_trees_rgb.png");

    cout << "Testing saturate..." << endl;
    Image saturated_neg = saturate(castle, -1.0f);
    Image saturated_pos = saturate(castle, 2.0f);
    saturated_neg.write("./Output/castle_saturated_neg.png");
    saturated_pos.write("./Output/castle_saturated_pos.png");

    cout << "Testing gamma..." << endl;
    Image gamma_nature = gamma_code(nature, 1.8);
    gamma_nature.write("./Output/skies_and_trees_gamma.png");

    cout << "Testing quantize..." << endl;
    Image quantized_nature_3 = quantize(nature, 3);
    Image quantized_nature_4 = quantize(nature, 4);
    quantized_nature_3.write("./Output/skies_and_trees_q3.png");
    quantized_nature_4.write("./Output/skies_and_trees_q4.png");

    cout << "Testing spanish illusion..." << endl;
    Image zebra = Image("./Input/zebra.png");
    std::vector<Image> zebra_output = spanish(zebra);
    std::vector<Image> castle_output = spanish(castle);
    zebra_output[0].write("./Output/zebra_spanish1.png");
    zebra_output[1].write("./Output/zebra_spanish2.png");
    castle_output[0].write("./Output/castle_spanish1.png");
    castle_output[1].write("./Output/castle_spanish2.png");

    cout << "Testing grayworld white balancing..." << endl;
    Image flower = Image("./Input/flower.png");
    Image flower_wb = grayworld(flower);
    flower_wb.write("./Output/flower_white_balance.png");

    cout << "Testing gamma_test..." << endl;
    std::vector<Image> gamma_test_output = gamma_test(nature, 4, 1.8);
    Image quantize_first = gamma_test_output[0];
    Image gamma_first = gamma_test_output[1];
    quantize_first.write("./Output/nature_quantize_first.png");
    gamma_first.write("./Output/nature_gamma_first.png");
}
