#ifndef A10_H_PHUDVTKB
#define A10_H_PHUDVTKB

// Write your declarations here, or extend the Makefile if you add source
// files
#include "Image.h"
#include "basicImageManipulation.h"
#include "blending.h"

std::string request(std::string query);
std::vector<std::string> getUrls(std::string query, int count, int purpose);
Image getImage(std::string &url);
Image singleColorMosaic(int w, int h, int tile_size, std::string query);
Image flatColorWheel(int w, int h, int tile_size, std::vector<std::string> allQueries);
Image whitePadding(Image im, int px, int py);
Image mosaic(Image target, int tile_size);
Image amalgamation(int w, int h, std::vector<std::string> queries, int count);

#endif /* end of include guard: A10_H_PHUDVTKB */

