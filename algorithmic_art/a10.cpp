#include <iostream>
#include "a10.h"

#include "curl_easy.h"
#include "curl_exception.h"
#include "curl_header.h"
#include "curl_ios.h"

// #include <curl/curl.h>
#include <Magick++.h>
#include <json/json.h>

using namespace std;

using curl::curl_easy;
using curl::curl_ios;
using curl::curl_easy_exception;
using curl::curlcpp_traceback;
using curl::curl_header;

// Write your implementations here, or extend the Makefile if you add source
// files
string request(string query) {
	ostringstream url;
	url << "https://api.flickr.com/services/rest/?method=flickr.photos.search&api_key=dd61e5495c5eae3c68c9ec5503904db5&sort=relevance&safe_search=true&text=" 
		<< query 
		<< "&format=json&nojsoncallback=1";
	return url.str();
}

string getPhotoSource(Json::Value data, int purpose) {
	string size_suffix;
	if (purpose == 0) {
		size_suffix = "_s.jpg"; // s --> small size (75 x 75)
	} else {
		size_suffix = "_n.jpg"; // n --> small, 320 on longest size
	}
	
	ostringstream srcUrl;
	srcUrl << "https://farm" << data["farm"] 
			<< ".staticflickr.com/" 
			<< data["server"].asString() 
			<< "/" << data["id"].asString()
			<< "_" << data["secret"].asString()
			<< size_suffix; 
	return srcUrl.str();
}

vector<string> getUrls(string query, int count, int purpose) {
	vector<string> urls;
	string url = request(query);

	ostringstream str;
	curl_ios<ostringstream> writer(str);
	
	curl_easy easy(writer);
	curl_header header;

	header.add("Accept: application/json");
	header.add("Content-type: application/json");
	header.add("charsets: utf-8");

	easy.add<CURLOPT_HTTPHEADER>(header.get());
	easy.add<CURLOPT_URL>(url.c_str());

	try {
		easy.perform();
	} catch (curl_easy_exception error) {
		error.print_traceback();
	}

	Json::Value jsonData; 
	Json::Reader jsonReader;
	if (jsonReader.parse(str.str(), jsonData)) {
		Json::Value img = jsonData["photos"]["photo"];
		for (int i=0; i < count; i++) {
			string srcUrl = getPhotoSource(img[i], purpose);
			urls.push_back(srcUrl.c_str());
		}
		return urls;
	} else {
		cout << "failed to read json" << endl;
	}

	return urls;

}

Image getImage(string &url) {
	ostringstream str;
	curl_ios<ostringstream> writer(str);

	curl_easy easy(writer);
	curl_header header;

	header.add("Accept: image/jpeg");
	header.add("Content-type: image/jpeg");

	easy.add<CURLOPT_HTTPHEADER>(header.get());
	easy.add<CURLOPT_URL>(url.c_str());

	try {
		easy.perform();

		string data = str.str();

		Magick::Image img;
		Magick::Blob blob(static_cast<const void *>(data.c_str()), data.size());
		img.read(blob);
		
		img.magick("png");
		img.write("./Output/temp.png");
		Image png("./Output/temp.png");
		return png;

	} catch (...) {
		return Image(75,75,3);
	}

}

bool isAllZero(Image im) {
	for (int i=0; i < im.width(); i++) {
		for (int j=0; j < im.height(); j++) {
			for (int k=0; k < im.channels(); k++) {
				if (im(i,k,k) != 0) {
					return false;
				}
			}
		}
	}

	return true;
}

Image singleColorMosaic(int x_chunk, int y_chunk, int tile_size, string query) {
	Image output(x_chunk*tile_size, y_chunk*tile_size, 3);

	vector<string> urls = getUrls(query, x_chunk * y_chunk, 0);
	int idx = 0;
	for (int x=0; x < x_chunk; x++) {
		for (int y=0; y < y_chunk; y++) {
			Image im = getImage(urls[idx]);
			idx += 1;

			for (int i=0; i < tile_size; i++) {
				for (int j=0; j < tile_size; j++) {
					for (int k=0; k < im.channels(); k++) {
						output(i + x*tile_size, j + y*tile_size, k) = im(i,j,k);
					}
				}
			} 
		}
	}

	return output;
}

Image flatColorWheel(int w, int h, int tile_size, vector<string> allQueries) {
	Image colorWheel(w*tile_size*allQueries.size(), h*tile_size, 3);

	for (int q=0; q < allQueries.size(); q++) {
		Image im = singleColorMosaic(w, h, tile_size, allQueries[q]);
		ostringstream filename;
		filename << "./Output/" << allQueries[q] << "Aggregate.png";
		im.write(filename.str());
		for (int x=0; x < im.width(); x++) {
			for (int y=0; y < im.height(); y++) {
				for (int c=0; c < im.channels(); c++) {
					colorWheel(x + q*im.width(), y, c) = im(x,y,c);
				}
			}
		}
	}

	return colorWheel;
}

Image whitePadding(Image im, int px, int py) {
	Image out(im.width() + 2*px, im.height() + 2*py, im.channels());
	out.set_color(255.0f, 255.0f, 255.0f);

	for (int i=0; i < im.width(); i++) {
		for (int j=0; j < im.height(); j++) {
			for (int k=0; k < im.channels(); k++) {
				out(i + px, j + py, k) = im(i,j,k);
			}
		}
	}
	return out;
}

Image mosaic(Image target, int tile_size) { 
	// for turtle tile size should be 97, yields 20x25 tiles perfectly
	int num_x_tiles = target.width() / tile_size;
	int num_y_tiles = target.height() / tile_size;

	vector<string> redUrls = getUrls("orange", 5000, 0);
	vector<string> greenUrls = getUrls("green", 1000, 0);
	vector<string> blueUrls = getUrls("blue", 10000, 0);

	int red_idx = 0;
	int green_idx = 0;
	int blue_idx = 0;

	Image output(target.width(), target.height(), 3);

	for (int x=0; x < num_x_tiles; x++) {
		for (int y=0; y < num_y_tiles; y++) {
			// get dominant rgb signal
			float tile_r = 0.0f;
			float tile_g = 0.0f;
			float tile_b = 0.0f;
			for (int i=0; i < tile_size; i++) {
				for (int j=0; j < tile_size; j++) {
					tile_r += target(i + x*tile_size, j + y*tile_size, 0);
					tile_g += target(i + x*tile_size, j + y*tile_size, 1);
					tile_b += target(i + x*tile_size, j + y*tile_size, 2);
					
				}
			}

			Image tile(tile_size, tile_size, 3);

			if (tile_r > tile_g && tile_r > tile_b) {
 				tile = getImage(redUrls[red_idx]);
				red_idx += 1;
 				
 				while (isAllZero(tile)) {
	 				tile = getImage(redUrls[red_idx]);
					red_idx += 1;

					if (red_idx > redUrls.size() - 1) {
						vector<string> moreUrls = getUrls("red,orange", 1000, 0);
						redUrls.insert(redUrls.end(), moreUrls.begin(), moreUrls.end());
					}

				}

			} else if (tile_b > tile_r && tile_b > tile_g) {
				tile = getImage(blueUrls[blue_idx]);
				blue_idx += 1;
				while (isAllZero(tile)) {
					tile = getImage(blueUrls[blue_idx]);
					blue_idx += 1;

					if (blue_idx > blueUrls.size() - 1) {
						vector<string> moreUrls = getUrls("blue", 1000, 0);
						blueUrls.insert(blueUrls.end(), moreUrls.begin(), moreUrls.end());
					}

				}

			} else {
				tile = getImage(greenUrls[green_idx]);
				green_idx += 1;
				while (isAllZero(tile)) {
					tile = getImage(greenUrls[green_idx]);
					green_idx += 1;

					if (green_idx > greenUrls.size() - 1) {
						vector<string> moreUrls = getUrls("green", 1000, 0);
						greenUrls.insert(greenUrls.end(), moreUrls.begin(), moreUrls.end());
					}

 				}		
			}

			for (int i=0; i < tile_size; i++) {
				for (int j=0; j < tile_size; j++) {
					for (int k=0; k < 3; k++) {
						output(i + x*tile_size, j + y*tile_size, k) = tile(i,j,k);
					}
				}
			}

		}
	}

	return output;
} 

Image amalgamation(int w, int h, vector<string> queries, int count) {
	Image output(w, h, 3);
	vector<string> urls;
	for (int i=0; i < queries.size(); i++) {
		vector<string> tmp = getUrls(queries[i], count, 1);
		urls.insert(urls.end(), tmp.begin(), tmp.end());
	}

	vector<int> idxs;
	for (int q=0; q < urls.size(); q++) {
		Image im = getImage(urls[q]);
		if (im.width() == w && im.height() == h) {
			idxs.push_back(q);
		}
	} 

	vector<string> refinedUrls;
	for (int i=0; i < idxs.size(); i++) {
		refinedUrls.push_back(urls[idxs[i]]);
	}

	for (int u=0; u < refinedUrls.size(); u++) {
		Image im = getImage(refinedUrls[u]);

		for (int i=0; i < im.width(); i++) {
			for (int j=0; j < im.height(); j++) {
				for (int k=0; k < im.channels(); k++) {
					output(i,j,k) += im.smartAccessor(i,j,k,true) / float(refinedUrls.size());
				}
			}
		}
		
	}

	return output;
}