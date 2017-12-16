#include "blending.h"
#include "matrix.h"
#include <ctime>

using namespace std;

Image blendingweight(int imwidth, int imheight) {
	// // --------- HANDOUT  PS07 ------------------------------
	
	Image output = Image(imwidth, imheight, 1);
	for (int i=0; i < imwidth; i++) {
		for (int j=0; j < imheight; j++) {
			float x = 1.0f - fabs(2.0f * (0.5f - (i/float(imwidth))));
			float y = 1.0f - fabs(2.0f * (0.5f - (j/float(imheight))));
			output(i,j) = x*y;
		}
	}

	return output;

}

//  ****************************************************************************
//  * blending related functions re-written from previous assignments
//  ****************************************************************************

// instead of writing source in out, *add* the source to out based on the weight
// so out(x,y) = out(x, y) + weight * image
void applyhomographyBlend(const Image &source,
						  const Image &weight,
						  Image &out,
						  const Matrix &H,
						  bool bilinear) {
	// // --------- HANDOUT  PS07 ------------------------------
	
	Matrix H_inv = H.inverse();
    for (int i=0; i < out.width(); i++) {
        for (int j=0; j < out.height(); j++) {
            Vec3f loc = H_inv * Vec3f(i,j,1);
            float w_prime = loc[2];
            float x_prime = loc[0] / w_prime;
            float y_prime = loc[1] / w_prime;

            for (int k=0; k < out.channels(); k++) {
                if (x_prime < source.width() && y_prime < source.height() && x_prime >= 0 && y_prime >= 0) {
                    if (bilinear) {
                        out(i,j,k) += interpolateLin(source, x_prime, y_prime, k, true) * interpolateLin(weight, x_prime, y_prime, 0, true);
                    } else {
                        out(i,j,k) += source(x_prime, y_prime, k) * weight(x_prime, y_prime, 0);
                    }
                }
            }
        }
    }

    return;

}

Image stitchLinearBlending(const Image &im1,
						   const Image &im2,
						   const Image &we1,
						   const Image &we2,
						   const Matrix &H) {
	// // --------- HANDOUT  PS07 ------------------------------
	// stitch using image weights.
	// note there is no weight normalization.

	BoundingBox bbox = bboxUnion(computeTransformedBBox(im1.width(), im1.height(), H), BoundingBox(0, im2.width()-1, 0, im2.height()-1));
	Matrix T = makeTranslation(bbox);

	Image output = Image(bbox.x2-bbox.x1, bbox.y2-bbox.y1, im1.channels());
	applyhomographyBlend(im1, we1, output, T*H, true);
	applyhomographyBlend(im2, we2, output, T, true);

	return output;
}


/*****************************************************************************
 * blending functions Pset 08
 *****************************************************************************/


// low freq and high freq (2-scale decomposition)
vector<Image> scaledecomp(const Image &im, float sigma) {
	vector <Image> ims;
	ims.push_back(gaussianBlur_separable(im, sigma));
	ims.push_back(im - ims[0]);
	return ims;
}

// stitch using different blending models
// blend can be 0 (none), 1 (linear) or 2 (2-layer)
Image stitchBlending(const Image &im1,
					 const Image &im2,
					 const Matrix &H,
					 BlendType blend) {
	// // --------- HANDOUT  PS07 ------------------------------

	switch (blend) {
		case BlendType::BLEND_LINEAR: {
			Image we1 = blendingweight(im1.width(), im1.height());
			Image we2 = blendingweight(im2.width(), im2.height());

			Image slb = stitchLinearBlending(im1, im2, we1, we2, H);

			Image weights = Image(slb.width(), slb.height(), 1);
			BoundingBox bbox = bboxUnion(computeTransformedBBox(we1.width(), we1.height(), H), BoundingBox(0, we2.width()-1, 0, we2.height()-1));
			Matrix T = makeTranslation(bbox);
			Image w1 = Image(we1.width(), we1.height(), 1);
			Image w2 = Image(we2.width(), we2.height(), 1);
			for (int i=0; i < w1.number_of_elements(); i++) {
				w1(i) = 1.0f;
			} 
			for (int i=0; i < w2.number_of_elements(); i++) {
				w2(i) = 1.0f;
			}
			applyhomographyBlend(w1, we1, weights, T*H, true);
			applyhomographyBlend(w2, we2, weights, T, true);

			for (int i=0; i < slb.width(); i++) {
				for (int j=0; j < slb.height(); j++) {
					for (int k=0; k < slb.channels(); k++) {

						if (weights(i,j) != 0) {
							slb(i,j,k) /= weights(i,j);
						}

					}
				}
			}

			return slb;
			break;

		} 
		case BlendType::BLEND_2LAYER: {
			Image we1 = blendingweight(im1.width(), im1.height());
			Image we2 = blendingweight(im2.width(), im2.height());

			vector<Image> sd1 = scaledecomp(im1, 2.0f);
			vector<Image> sd2 = scaledecomp(im2, 2.0f);
			Image lo_freq_1 = sd1[0];
			Image lo_freq_2 = sd2[0];
			Image hi_freq_1 = sd1[1];
			Image hi_freq_2 = sd2[1];

			Image lf = stitchBlending(lo_freq_1, lo_freq_2, H, BlendType::BLEND_LINEAR);

			// only keep the high frequency of the image with the highest weight
			BoundingBox bbox = bboxUnion(computeTransformedBBox(we1.width(), we1.height(), H), BoundingBox(0, we2.width()-1, 0, we2.height()-1));
			Matrix T = makeTranslation(bbox);
			Image out1 = Image(lf.width(), lf.height(), 1);
			Image out2 = Image(lf.width(), lf.height(), 1);
			Image hf1 = Image(lf.width(), lf.height(), 3);
			Image hf2 = Image(lf.width(), lf.height(), 3);
			
			applyHomography(we1, T*H, out1, true);
			applyHomography(we2, T, out2, true);
			applyHomography(hi_freq_1, T*H, hf1, true);
			applyHomography(hi_freq_2, T, hf2, true);
			
			Image hf = hf1;
			for (int i=0; i < lf.width(); i++) {
				for (int j=0; j < lf.height(); j++) {
					if (out1(i,j) < out2(i,j)) {
						for (int k=0; k < hf1.channels(); k++) {
							hf(i,j,k) = hf2(i,j,k);
						}
					}
				}
			}

			return lf + hf;
			break;
		}

		default: {
			Image we1 = blendingweight(im1.width(), im1.height());
			Image we2 = blendingweight(im2.width(), im2.height());

			BoundingBox bbox = bboxUnion(computeTransformedBBox(im1.width(), im1.height(), H), BoundingBox(0, im2.width()-1, 0, im2.height()-1));
			Matrix T = makeTranslation(bbox);
			Image output = Image(bbox.x2-bbox.x1, bbox.y2-bbox.y1, im1.channels());
			
			applyHomography(im1, T*H, output, true);
			applyHomography(im2, T, output, true);
			
			return output;
			break;
		}

	}
}

// auto stitch
Image autostitch(const Image &im1,
				 const Image &im2,
				 BlendType blend,
				 float blurDescriptor,
				 float radiusDescriptor) {
	// // --------- HANDOUT  PS07 ------------------------------
	vector<Point> corners_1 = HarrisCorners(im1);
	vector<Point> corners_2 = HarrisCorners(im2);
	vector<Feature> features_1 = computeFeatures(im1, corners_1, blurDescriptor, radiusDescriptor);
	vector<Feature> features_2 = computeFeatures(im2, corners_2, blurDescriptor, radiusDescriptor);
	vector <FeatureCorrespondence> correspondences = findCorrespondences(features_1, features_2);
	Matrix H = RANSAC(correspondences);

	return stitchBlending(im1, im2, H, blend);
}

/************************************************************************
 * Tricks: mini planets.
 ************************************************************************/

Image pano2planet(const Image &pano, int newImSize, bool clamp) {
	// // --------- HANDOUT  PS07 ------------------------------
	Image planet = Image(newImSize, newImSize, pano.channels());

	float center = newImSize / 2.0f;
	float angle_factor = pano.width() / (2*M_PI);
	float rad_factor = pano.height() / center;

	for (int i=0; i < newImSize; i++) {
		for (int j=0; j < newImSize; j++) {
			float theta = -atan2(j-center, i-center);
	    	float r = sqrt(pow(i-center, 2) + pow(j-center, 2));
	    	if (theta > 0) {
	    		for (int k=0; k < pano.channels(); k++) {
	    			planet(i,j,k) = interpolateLin(pano, angle_factor*theta, pano.height()-rad_factor*r, k, true);
	    		}
	    	} else {
	    		for (int k=0; k < pano.channels(); k++) {
	    			planet(i,j,k) = interpolateLin(pano, pano.width()+angle_factor*theta, pano.height()-rad_factor*r, k, true);
	    		}
	    	}
	    	
	  	}
	}

	return planet;
}


/************************************************************************
 * 6.865: Stitch N images into a panorama
 ************************************************************************/

// Pset07-865. Compute sequence of N-1 homographies going from Im_i to Im_{i+1}
// Implement me!
vector<Matrix> sequenceHs(const vector<Image> &ims,
						  float blurDescriptor,
						  float radiusDescriptor) {
	// // --------- HANDOUT  PS07 ------------------------------
	
	vector<Matrix> Hs;
	for (int i=0; i < ims.size()-1; i++) {
	 	vector<Point> corners_im1 = HarrisCorners(ims[i]);
	 	vector<Point> corners_im2 = HarrisCorners(ims[i+1]);

	  	vector<Feature> features_im1 = computeFeatures(ims[i], corners_im1, blurDescriptor, radiusDescriptor);
	  	vector<Feature> features_im2 = computeFeatures(ims[i+1], corners_im2, blurDescriptor, radiusDescriptor);

	 	vector<FeatureCorrespondence> correspondences = findCorrespondences(features_im1, features_im2);
	 	Matrix H = RANSAC(correspondences);
	 	Hs.push_back(H);
	}

	return Hs;
}

// stack homographies:
//   transform a list of (N-1) homographies that go from I_i to I_i+1
//   to a list of N homographies going from I_i to I_refIndex.
vector<Matrix> stackHomographies(const vector<Matrix> &Hs, int refIndex) {
	// // --------- HANDOUT  PS07 ------------------------------
	
	vector<Matrix> homographies;
	for (int i = 0; i <= Hs.size(); i++) {
	 	Matrix H = Matrix::Identity(3,3);

	 	if (i < refIndex) {
	 		H = Hs[i];
	 		for (int j=i+1; j < refIndex; j++) {
	 			H = Hs[j] * H;
	 		}
	 	} else if (i > refIndex) {
	 		H = Hs[i-1].inverse();
	 		for (int j=i-2; j >= refIndex; j--) {
	 			H = Hs[j].inverse() * H;
	 		}
	 	}
	 	homographies.push_back(H);
	}

	return homographies;
}


// Pset07-865: compute bbox around N images given one main reference.
BoundingBox bboxN(const vector<Matrix> &Hs, const vector<Image> &ims) {
	// // --------- HANDOUT  PS07 ------------------------------
	
	BoundingBox bbox = BoundingBox(0,0,0,0);
	for (int i = 0; i < ims.size(); i++) {
		bbox = bboxUnion(computeTransformedBBox(ims[i].width(), ims[i].height(), Hs[i]), bbox);
	}

	return bbox;
}

// Pset07-865.
// Implement me!
Image autostitchN(const vector<Image> &ims,
				  int refIndex,
				  float blurDescriptor,
				  float radiusDescriptor) {
	// // --------- HANDOUT  PS07 ------------------------------
	vector<Matrix> Hs = sequenceHs(ims, blurDescriptor, radiusDescriptor);
	vector<Matrix> stacked_Hs = stackHomographies(Hs, refIndex);
	BoundingBox bbox = bboxN(stacked_Hs, ims);
	Matrix T = makeTranslation(bbox);

	Image output = Image(bbox.x2-bbox.x1, bbox.y2-bbox.y1, ims[0].channels());
	Image weights = Image(bbox.x2-bbox.x1, bbox.y2-bbox.y1, 1);

	for (int i=0; i<ims.size(); i++) {
		Image weight = blendingweight(ims[i].width(), ims[i].height());
		Image white = Image(weight.width(), weight.height(), 1);
		for (int j = 0; j < white.number_of_elements(); j++) {
			white(j) = 1.0f;
		}
	  	applyhomographyBlend(ims[i], weight, output, T*stacked_Hs[i], true);
	  	applyhomographyBlend(white, weight, weights, T*stacked_Hs[i], true);
	}

	for (int i=0; i < output.width(); i++) {
		for (int j=0; j < output.height(); j++) {
			for (int k=0; k < output.channels(); k++) {
				if (weights(i,j) != 0) {
					output(i,j,k) = output(i,j,k)/weights(i,j);
				}
			}
	  	}
	}

	return output;
}


/******************************************************************************
 * Helper functions
 *****************************************************************************/

// copy a single-channeled image to several channels
Image copychannels(const Image &im, int nChannels) {
	// image must have one channel
	assert(im.channels() == 1);
	Image oim(im.width(), im.height(), nChannels);

	for (int i = 0; i < im.width(); i++) {
		for (int j = 0; j < im.height(); j++) {
			for (int c = 0; c < nChannels; c++) {
				oim(i, j, c) = im(i, j);
			}
		}
	}
	return oim;
}

