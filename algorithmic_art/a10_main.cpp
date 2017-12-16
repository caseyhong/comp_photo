#include <iostream>

#include "a10.h"

#include "matrix.h"
#include "panorama.h"
#include "blending.h"

using namespace std;

void testGetUrls() {
	cout << "Testing getUrls()" << endl;
	int count = 10;
	vector<string> urls = getUrls("red", count, 0);
    cout << "Retrieved " << count << " urls for tag RED" << endl;
    for (int i=0; i < urls.size(); i++) {
		cout << urls[i] << endl;
    }
};

void testDownload() {
	cout << "Testing getImage()" << endl;
	vector<string> urls_s = getUrls("red", 1, 0);
	cout << "url : " << urls_s[0] << endl;
	Image im = getImage(urls_s[0]);
	im.write("./Output/sampleDownload_s.png");

	vector<string> urls_n = getUrls("red", 1, 1);
	Image im2 = getImage(urls_n[0]);
	im2.write("./Output/sampleDownload_n.png");
}

void testRed() {
	cout << "Testing singleColorMosaic()" << endl;
	Image im = singleColorMosaic(10, 5, 75, "red");
	cout << "im size : " << im.width() << " x " << im.height() << endl;
	im.write("./Output/redAggregate.png");
}

void testFlatColorWheel() {
	cout << "Testing flatColorWheel()" << endl;
	vector<string> colorTags = {"red", "red,orange", "orange", "orange,yellow", "yellow", "yellow,green", "green", "green,blue", "blue", "blue,purple", "purple"};
	Image im = flatColorWheel(15, 5, 75, colorTags);
	cout << "im size : " << im.width() << " x " << im.height() << endl;
	im.write("./Output/colorWheel.png");
}

void testWheelify() {
	cout << "Wheelifying the flat color wheel..." << endl;
	Image im("./Output/colorWheel.png");
	Image im2 = whitePadding(im, 0, 200);
	Image wheel = pano2planet(im2, 1200);
	wheel.write("./Output/colorWheel-2.png");
}

void testTurtle() {
	cout << "Creating turtle mosaic..." << endl;
	Image turtle("./Input/turtle.png");
	turtle = scaleLin(turtle, (float)1500.0f / 1940.0f);
	Image out = mosaic(turtle, 75);
	out.write("./Output/turtleMosaic.png");
}

void rainieramalgamation() {
	cout << "Rainier" << endl;
	vector<string> mtns = {"mt,rainier", "mount,rainier", "rainier"};
	Image out = amalgamation(320, 213, mtns, 50);
	out.write("./Output/rainier-amalgamation.png");
}

void denaliamalgamation() {
	cout << "Denali" << endl;
	vector<string> mtns = {"denali"};
	Image out = amalgamation(320, 213, mtns, 50);
	out.write("./Output/denali-amalgamation.png");
}

void matterhornamalgamation() {
	cout << "Matterhorn" << endl;
	vector<string> mtns = {"matterhorn"};
	Image out = amalgamation(320, 213, mtns, 50);
	out.write("./Output/matterhorn-amalgamation.png");
}

void everestamalgamation() {
	cout << "Everest" << endl;
	vector<string> mtns = {"mount,everest", "mt,everest"};
	Image out = amalgamation(320, 213, mtns, 50);
	out.write("./Output/everest-amalgamation.png");
}

void chicagoamalgamation() {
	cout << "Chicago" << endl;
	vector<string> chicago = {"chicago", "Chicago"};
	Image out = amalgamation(320, 213, chicago, 100);
	out.write("./Output/chicago-amalgamation.png");
}

void seattleamalgamation() {
	cout << "Seattle" << endl;
	vector<string> seattle = {"seattle"};
	Image out = amalgamation(320, 213, seattle, 100);
	out.write("./Output/seattle-amalgamation.png");
}

void nycamalgamation() {
	cout << "NYC" << endl;
	vector<string> nyc = {"new,york,city", "nyc"};
	Image out = amalgamation(320, 213, nyc, 100);
	out.write("./Output/nyc-amalgamation.png");
}

void yosemiteamalgamation() {
	cout << "Yosemite" << endl;
	vector<string> yosemite = {"yosemite"};
	Image out = amalgamation(320, 213, yosemite, 100);
	out.write("./Output/yosemite-amalgamation.png");
}

void halfdomeamalgamation() {
	cout << "Half Dome" << endl;
	vector<string> hd = {"half,dome", "halfdome"};
	Image out = amalgamation(320, 213, hd, 100);
	out.write("./Output/hd-amalgamation.png");
}

void bostonamalgamation() {
	cout << "Boston" << endl;
	vector<string> boston = {"boston"};
	Image out = amalgamation(320, 213, boston, 100);
	out.write("./Output/boston-amalgamation.png");
}

void seoulamalgamation() {
	cout << "Seoul" << endl;
	vector<string> seoul = {"seoul"};
	Image out = amalgamation(320, 213, seoul, 100);
	out.write("./Output/seoul-amalgamation.png");
}

void parisamalgamation() {
	cout << "Paris" << endl;
	vector<string> paris = {"paris"};
	Image out = amalgamation(320, 213, paris, 100);
	out.write("./Output/paris-amalgamation.png");
}

int main() {
	testGetUrls();
	testDownload();
	testRed();	
	testFlatColorWheel();
	testWheelify();
	testTurtle();
	
	rainieramalgamation();
	denaliamalgamation();
	matterhornamalgamation();
	everestamalgamation();
	
	bostonamalgamation();
	chicagoamalgamation();
	seoulamalgamation();
	parisamalgamation();
	nycamalgamation();
	seattleamalgamation();
	yosemiteamalgamation();
	halfdomeamalgamation();
}