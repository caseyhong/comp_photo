#include "homography.h"
#include "matrix.h"

using namespace std;


void applyHomography(const Image &source, const Matrix &H, Image &out, bool bilinear) {
    // // --------- HANDOUT  PS06 ------------------------------
    // Transform image source using the homography H, and composite in onto out.
    // if bilinear == true, using bilinear interpolation. Use nearest neighbor
    // otherwise.
    // 
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
                        out(i,j,k) = interpolateLin(source, x_prime, y_prime, k, true);
                    } else {
                        out(i,j,k) = source(x_prime, y_prime, k);
                    }
                }
            }
        }
    }

    return;
}




Matrix computeHomography(const CorrespondencePair correspondences[4]) {
    // --------- HANDOUT  PS06 ------------------------------
    // Compute a homography from 4 point correspondences.
    
    Matrix A = Matrix::Zero(8,8);
    Matrix B = Matrix::Zero(8,1);

    for (int pair=0; pair < 4; pair++) {
        Vec3f p1 = correspondences[pair].point1;
        Vec3f p2 = correspondences[pair].point2;

        int a_row = 2*pair;
        int d_row = 2*pair + 1;

        A(a_row,0) = p1[0];
        A(a_row,1) = p1[1];
        A(a_row,2) = 1;
        A(a_row,6) = -p1[0]*p2[0];
        A(a_row,7) = -p1[1]*p2[0];

        A(d_row,3) = p1[0];
        A(d_row,4) = p1[1];
        A(d_row,5) = 1;
        A(d_row,6) = -p1[0]*p2[1];
        A(d_row,7) = -p1[1]*p2[1];

        B(a_row,0) = p2[0];
        B(d_row,0) = p2[1];
    }

    Matrix x = A.fullPivLu().solve(B);

    Matrix H = Matrix::Zero(3,3);

    H(0,0) = x(0,0);
    H(1,0) = x(3,0);
    H(2,0) = x(6,0);

    H(0,1) = x(1,0);
    H(1,1) = x(4,0);
    H(2,1) = x(7,0);

    H(0,2) = x(2,0);
    H(1,2) = x(5,0);
    H(2,2) = 1;

    return H;
}


BoundingBox computeTransformedBBox(int imwidth, int imheight, Matrix H) {
    // --------- HANDOUT  PS06 ------------------------------
    // Predict the bounding boxes that encompasses all the transformed
    // coordinates for pixels frow and Image with size (imwidth, imheight)
    
    Vec3f ll = H * Vec3f(0,0,1);
    Vec3f lr = H * Vec3f(imwidth,0,1);
    Vec3f ur = H * Vec3f(imwidth, imheight, 1);
    Vec3f ul = H * Vec3f(0,imheight,1);

    ll /= ll[2];
    lr /= lr[2];
    ur /= ur[2];
    ul /= ul[2];

    int x0 = min(min(min(ll[0], ul[0]), lr[0]), ur[0]);
    int x1 = max(max(max(ll[0], ul[0]), lr[0]), ur[0]);

    int y0 = min(min(min(ll[1], ul[1]), lr[1]), ur[1]);
    int y1 = max(max(max(ll[1], ul[1]), lr[1]), ur[1]);

    return BoundingBox(x0,x1,y0,y1);
}


BoundingBox bboxUnion(BoundingBox B1, BoundingBox B2) {
    // --------- HANDOUT  PS06 ------------------------------
    // Compute the bounding box that tightly bounds the union of B1 an B2.
    
    int x0 = min(B1.x1, B2.x1);
    int x1 = max(B1.x2, B2.x2);

    int y0 = min(B1.y1, B2.y1);
    int y1 = max(B1.y2, B2.y2);

    return BoundingBox(x0,x1,y0,y1);

}


Matrix makeTranslation(BoundingBox B) {
    // --------- HANDOUT  PS06 ------------------------------
    // Compute a translation matrix (as a homography matrix) that translates the
    // top-left corner of B to (0,0).

    Matrix T = Matrix::Identity(3,3);
    T(0,2) = -B.x1;
    T(1,2) = -B.y1;
    return T;
}


Image stitch(const Image &im1, const Image &im2, const CorrespondencePair correspondences[4]) {
    // --------- HANDOUT  PS06 ------------------------------
    // Transform im1 to align with im2 according to the set of correspondences.
    // make sure the union of the bounding boxes for im2 and transformed_im1 is
    // translated properly (use makeTranslation)
    
    Matrix H = computeHomography(correspondences);
    BoundingBox bbox = bboxUnion(computeTransformedBBox(im1.width(), im1.height(), H), BoundingBox(0, im2.width()-1, 0, im2.height()-1));
    Matrix T = makeTranslation(bbox);

    Image output = Image(bbox.x2-bbox.x1, bbox.y2-bbox.y1, im1.channels());
   
    applyHomography(im1, T*H, output, true);
    applyHomography(im2, T, output, true);

    return output;
}

// debug-useful
Image drawBoundingBox(const Image &im, BoundingBox bbox) {
    // // --------- HANDOUT  PS06 ------------------------------
    /*
      ________________________________________
     / Draw me a bounding box!                \
     |                                        |
     | "I jumped to my                        |
     | feet, completely thunderstruck. I      |
     | blinked my eyes hard. I looked         |
     | carefully all around me. And I saw a   |
     | most extraordinary small person, who   |
     | stood there examining me with great    |
     | seriousness."                          |
     \              Antoine de Saint-Exupery  /
      ----------------------------------------
             \   ^__^
              \  (oo)\_______
                 (__)\       )\/\
                     ||----w |
                     ||     ||
    */
    
    Image newImage(im.width(), im.height(), im.channels());
    for (int x = 0; x < newImage.width(); x++) {
      for (int y = 0; y < newImage.height(); y++) {
        for (int c = 0; c < newImage.channels(); c++) {
          newImage(x,y,c) = im(x,y,c);
        }
        if (x == bbox.x1 & x < bbox.x2 & y > bbox.y1 & y < bbox.y2) {
          newImage(x, y, 0) = 0.0;
        } else if (x > bbox.x1 & x == bbox.x2 & y > bbox.y1 & y < bbox.y2) {
          newImage(x, y, 0) = 0.0;
        } else if (x > bbox.x1 & x < bbox.x2 & y == bbox.y1 & y < bbox.y2) {
          newImage(x, y, 0) = 0.0;
        } else if (x > bbox.x1 & x < bbox.x2 & y > bbox.y1 & y == bbox.y2) {
          newImage(x, y, 0) = 0.0;
        }
      }
    }
    return newImage;
}

void applyHomographyFast(const Image &source, const Matrix &H, Image &out, bool bilinear) {
    // // --------- HANDOUT  PS06 ------------------------------
    // Same as apply but change only the pixels of out that are within the
    // predicted bounding box (when H maps source to its new position).
     
    Matrix H_inv = H.inverse();
    BoundingBox bbox = computeTransformedBBox(source.width(), source.height(), H);
    int x_start = max(0,bbox.x1);
    int x_end = min(out.width(), bbox.x2+1);
    int y_start = max(0,bbox.y1);
    int y_end = min(out.width(), bbox.y2+1);

    for (int i=x_start; i < x_end; i++) {
        for (int j=y_start; j < y_end; j++) {
            Vec3f loc = H_inv * Vec3f(i,j,1);
            float w_prime = loc[2];
            float x_prime = loc[0] / w_prime;
            float y_prime = loc[1] / w_prime;

            for (int k=0; k < out.channels(); k++) {
                if (x_prime < source.width() && y_prime < source.height() && x_prime >= 0 && y_prime >= 0) {
                    if (bilinear) {
                        out(i,j,k) = interpolateLin(source, x_prime, y_prime, k, true);
                    } else {
                        out(i,j,k) = source(x_prime, y_prime, k);
                    }
                }
            }
        }
    }

    return;
}
