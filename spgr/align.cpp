/* --------------------------------------------------------------------------
 * File:    align.cpp
 * Created: 2015-10-01
 * --------------------------------------------------------------------------
 *
 *
 *
 * ------------------------------------------------------------------------*/



#include "align.h"

using namespace std;

Image denoiseSeq(const vector<Image> &imSeq){
    // // --------- HANDOUT  PS03 ------------------------------
    // Basic denoising by computing the average of a sequence of images
    int w = imSeq[0].width();
    int h = imSeq[0].height();
    int c = imSeq[0].channels();
    Image output = Image(w, h, c);

    for (int x=0; x < imSeq.size(); x++) {
        for (int i=0; i < w; i++) {
            for (int j=0; j < h; j++) {
                for (int k=0; k < c; k++) {
                    output(i,j,k) += imSeq[x](i,j,k);
                }
            }
        }
    }
    
    float norm_const = 1.0f/imSeq.size();
    return output * norm_const;

}


Image logSNR(const vector<Image> &imSeq, float scale){
    // // --------- HANDOUT  PS03 ------------------------------
    // returns an image visualizing the per-pixel and
    // per-channel log of the signal-to-noise ratio scaled by scale.
    int w = imSeq[0].width();
    int h = imSeq[0].height();
    int c = imSeq[0].channels();
    Image output = Image(w, h, c);

        
    for (int i=0; i < w; i++) {
        for (int j=0; j < h; j++) {
            for (int k=0; k < c; k++) {
                float EI = 0.0f;
                float EI2 = 0.0f;
                for (int x=0; x < imSeq.size(); x++) {
                    float val = imSeq[x](i,j,k);
                    EI += val;
                    EI2 += pow(val, 2);
                }
                EI /= imSeq.size();
                EI2 /= imSeq.size();

                float sigma = 0.0f;
                for (int x=0; x < imSeq.size(); x++) {
                    sigma += pow(EI - imSeq[x](i,j,k), 2);
                }
                sigma /= (imSeq.size()-1); // bessels correction
                sigma += 0.000001; // numerical lubrication

                output(i,j,k) = scale*10*log10(EI2/sigma);
            }
        }
    }

    return output;
    

}


vector<int> align(const Image &im1, const Image &im2, int maxOffset){
    // // --------- HANDOUT  PS03 ------------------------------
    // returns the (x,y) offset that best aligns im2 to match im1.
    
    int x_offset = 0;
    int y_offset = 0;

    float min_loss = -1;

    for (int dx = -maxOffset; dx < maxOffset + 1; dx++) {
        for (int dy = -maxOffset; dy < maxOffset + 1; dy++) {
            float sse = 0.0f;
            for (int i = maxOffset; i < im1.width() - maxOffset; i++) {
                for (int j = maxOffset; j < im1.height() - maxOffset; j++) {
                    for (int k = 0; k < im1.channels(); k++) {
                        sse += pow(im1(i,j,k) - im2(i+dx,j+dy,k), 2);
                    }
                }
            }
            if (sse < min_loss || min_loss < 0) {
                min_loss = sse;
                x_offset = dx;
                y_offset = dy;
        }
      }
    }

    std::vector<int> output;
    output.push_back(-x_offset); // negative because we want to push im2 to match im1
    output.push_back(-y_offset);

    return output;

}

Image alignAndDenoise(const vector<Image> &imSeq, int maxOffset){
    // // --------- HANDOUT  PS03 ------------------------------
    // Registers all images to the first one in a sequence and outputs
    // a denoised image even when the input sequence is not perfectly aligned.
    Image ref_im = imSeq[0];
    std::vector<Image> aligned_images;
    aligned_images.push_back(ref_im);

    for (int i=1; i < imSeq.size(); i++) {
        Image im = imSeq[i];
        std::vector<int> offset = align(ref_im, im, maxOffset);
        int dx = offset[0];
        int dy = offset[1];
        Image rolled = roll(im,dx,dy);
        aligned_images.push_back(rolled);
    }

    return denoiseSeq(aligned_images);

}

// extra credit
Image split(const Image &sergeyImg){
    // --------- HANDOUT  PS03 ------------------------------
    // 6.865 only:
    // split a Sergey images to turn it into one 3-channel image.
    int height = floor(sergeyImg.height() / 3);
    Image output = Image(sergeyImg.width(), height, 3);

    for (int i=0; i < sergeyImg.width(); i++) {
        for (int j=0; j < height; j++) {
            output(i,j,2) = sergeyImg(i,j);
            output(i,j,1) = sergeyImg(i,j+height);
            output(i,j,0) = sergeyImg(i,j+2*height);
        }
    }

    return output;

}

// extra credit
Image sergeyRGB(const Image &sergeyImg, int maxOffset){
    // // --------- HANDOUT  PS03 ------------------------------
    // 6.865 only:
    // aligns the green and blue channels of your rgb channel of a sergey
    // image to the red channel. This should return an aligned RGB image
    Image rgb = split(sergeyImg);

    Image output = Image(rgb.width(), rgb.height(), 3);
    Image red = Image(rgb.width(), rgb.height(), 1);
    Image green = Image(rgb.width(), rgb.height(), 1);
    Image blue = Image(rgb.width(), rgb.height(), 1);

    for (int i=0; i < rgb.width(); i++) {
        for (int j=0; j < rgb.height(); j++) {
            red(i,j) = rgb(i,j,0);
            green(i,j) = rgb(i,j,1);
            blue(i,j) = rgb(i,j,2);
        }
    }

    // align green and blue channels to the RED channel
    std::vector<int> shift_green = align(red, green, maxOffset);
    std::vector<int> shift_blue = align(red, blue, maxOffset);

    Image rolled_green = roll(green, shift_green[0], shift_green[1]);
    Image rolled_blue = roll(blue, shift_blue[0], shift_blue[1]);

    for (int i=0; i < rgb.width(); i++) {
        for (int j=0; j < rgb.height(); j++) {
            output(i,j,0) = red(i,j);
            output(i,j,1) = rolled_green(i,j);
            output(i,j,2) = rolled_blue(i,j);
        }
    }

    return output;

}


/**************************************************************
 //               DON'T EDIT BELOW THIS LINE                //
 *************************************************************/

// This circularly shifts an image by xRoll in the x direction and
// yRoll in the y direction. xRoll and yRoll can be negative or postive
Image roll(const Image &im, int xRoll, int yRoll){

    int xNew, yNew;
    Image imRoll(im.width(), im.height(), im.channels());

    // for each pixel in the original image find where its corresponding
    // location is in the rolled image
    for (int x=0; x<im.width(); x++){
        for (int y=0; y<im.height(); y++){

            // use modulo to figure out where the new location is in the
            // rolled image. Then take care of when this returns a negative number
            xNew = (x + xRoll) % im.width();
            yNew = (y + yRoll) % im.height();
            xNew = (xNew<0)*(imRoll.width() + xNew) + (xNew>=0)*xNew;
            yNew = (yNew<0)*(imRoll.height() + yNew) + (yNew>=0)*yNew;

            // assign the rgb values for each pixel in the original image to
            // the location in the new image
            for (int z=0; z<im.channels(); z++){
                imRoll(xNew, yNew, z) = im(x,y,z);
            }
        }
    }

    // return the rolled image
    return imRoll;
}
