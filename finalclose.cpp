/*
* Simulation code for project:
* PMMA Additive-Induced Active Layer Self-Assembly in Polymer Solar Cells: Toward Low-Cost, Ultra-Portable, Mass-Produced Solar Energy
* Version: close columns
*/

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>

#define PI 3.14159265
#define E 2.718281828459
#define ALPHA 3.5 * (pow(10,-3)) //absorption coefficient of P3HT:PCBM blend
//all length units in nanometers

using namespace std;

double thickness;
double index_pol_one, index_pol_two;
double totalIntensity = 0.0; /absorbed light in the presence of nanostructures
double initialIntensity = 1.0; 

double distance(double x1, double x2, double y1, double y2){
	return sqrt(pow((x1-x2),2)+pow((y1-y2),2));
}

int main() {
	int pol_one, pol_two, pol_three;
	cout << "Enter cell thickness (in nm)" << endl;
	cin >> thickness;
	cout << "What is the type of polymer additive?" << endl;   //0=PMMA
	cin >> pol_one;
	cout << "What is the type of polymer blend?" << endl;   //0=P3HT:PCBM
	cin >> pol_two;

	switch(pol_one) {
		case 0:
			index_pol_one = 1.4893;
            break;
		case 1:
			index_pol_one = 1;
            break;
		default:
			cout << "Polymer not supported" << endl;
	}

	switch(pol_two) {
		case 0:
			index_pol_two = 1.83;
            break;
		case 2:
			index_pol_two = 1;
            break;
		default:
			cout << "Polymer not supported" << endl;
	}

	double crit_angle = asin(index_pol_one/index_pol_two)*(180/PI);
	cout << "Critical angle (degrees): " << crit_angle << endl;

	for (double i=156.0; i<656.0; i++) { //test from edge to edge between columns, incrementing by 1 nanometer at a time
		for (double j=1.0; j<90.0; j++) { //test incident angles from 1 to 90 degrees
			cout << i-156.0 << " nanometers from left edge, " << j << " degrees" << endl;
			double slope = tan (j * PI / 180.0);
			double xCoord = i;
			double yCoord = 158.0;
			double xStart = i; double yStart = 158.0;
			double xIntersect, yIntersect, tangent, incidentAngle;
			double parabolaTangent[4]; //x1 x2 y1 y2
			double reflectance = 1.0;
			initialIntensity = 1.0;

			while (true) {
				if (slope > 0){                     
					xCoord -= 0.001;
					yCoord -= 0.001 * slope;
				}
				else {
					xCoord += 0.001;
					yCoord += 0.001 * slope;
				}
				if (yCoord < 0 or yCoord > 158) { //end if light ray exits solar cell 
					double dist2 = distance(xStart, xCoord, yStart, yCoord);
					//totalDist += dist2;
					totalIntensity += initialIntensity * (1 - pow (E, (-1) * ALPHA * dist2)); //absorbed light intensity
					break;
				}
				if ((abs(pow(yCoord - 79.0, 2)/40.0 - xCoord) <= 0.001) 
				or (abs(810.0 - pow(yCoord - 79.0, 2)/40.0 - xCoord) <= 0.001)) { //extend ray incrementally till close enough to "hit" parabola
					xIntersect = xCoord; 
					yIntersect = yCoord;
					if (xIntersect <= 0.001) { //at the middle point of left parabola
						parabolaTangent[0] = 0.002; parabolaTangent[1] = 0.002;
						parabolaTangent[2] = sqrt(40.0 * 0.002) + 79.0; 
						parabolaTangent[3] = -1 * sqrt(40.0 * 0.002) + 79.0;
					}
					else if (xIntersect >= 809.999) { //at the middle point of right parabola
						parabolaTangent[0] = 809.998; parabolaTangent[1] = 809.998;
						parabolaTangent[2] = sqrt(40.0 * (810 - 809.998)) + 79.0;
						parabolaTangent[3] = - sqrt(40.0 * (810 - 809.998)) + 79.0;
					} 
					else { 
						if (xCoord < 500) { //if intersect at left parabola
							parabolaTangent[0] = xIntersect - 0.001; parabolaTangent[1] = xIntersect + 0.001;
							parabolaTangent[2] = sqrt(40.0 * parabolaTangent[0]) + 79.0;
							parabolaTangent[3] = sqrt(40.0 * parabolaTangent[1]) + 79.0; 
						}
						else { //if intersect at right parabola
							parabolaTangent[0] = xIntersect - 0.001; parabolaTangent[1] = xIntersect + 0.001;
							parabolaTangent[2] = sqrt(40.0 * (810.0 - parabolaTangent[0])) + 79.0;
							parabolaTangent[3] = sqrt(40.0 * (810.0 - parabolaTangent[1])) + 79.0; 
						}                                                                                                                             
					}
					tangent = (parabolaTangent[3]-parabolaTangent[2])/(parabolaTangent[1]-parabolaTangent[0]); //slope of tangent at intersection
					double normal = (-1.0) / (tangent); //slope of normal

					incidentAngle = abs((atan(slope) - atan(normal))*(180 / PI));
					if (incidentAngle <= 54.5)
						reflectance = 0.08;
					else 
						reflectance = 1.00;

					double dist2 = distance(xStart, xIntersect, yStart, yIntersect);
					totalIntensity += initialIntensity * (1 - pow (E, (-1) * ALPHA * dist2)); //absorbed light intensity
					initialIntensity = reflectance * initialIntensity * (pow(E, (-1) * ALPHA * dist2)); //new initial light intensity
					xStart = xIntersect;
					yStart = yIntersect;
					slope = tan(incidentAngle * PI / 180 + atan(normal)); 
					cout << "Reflection" << endl;
					cout << "Distance: " << dist2 << endl;
					cout << "Total intensity: " << totalIntensity << endl;
					cout << "New initial intensity:" << initialIntensity << dist2 << endl;
					if (initialIntensity < 0.001) //if new light intensity is negligible, end to shorten simulation runtime
						break;
				} 
			}
		} 
	}

	initialIntensity = 1.0;
	double originalIntensity = 0.0; //absorbed light without nanostructures
	for(int a=156.0; a<656.0; a++) {
		for (int b=1; b<90; b++) {
			double xCoord = a;
			double yCoord = 158.0;
			double slope = tan (b * PI / 180.0);
			while (true) {
				xCoord -= 0.001;
				yCoord -= 0.001 * slope;
				if (yCoord < 0) { //end if light ray exits solar cell 
					double dist3 = distance(a, xCoord, 158.0, yCoord); 
					originalIntensity += initialIntensity * (1 - pow (E, (-1) * ALPHA * dist3)); //absorbed light intensity
					break;
				}
				if (xCoord <= 78 or xCoord >= 740) { //limit to the same scope tested above
					double dist3 = distance(a, xCoord, 158.0, yCoord);
					originalIntensity += initialIntensity * (1 - pow(E, ((-1) * ALPHA * dist3)));
					break;
				}
			}
			cout << a << " nanometers from edge," << b << " degrees." << endl;
		}
	}
	
	cout << "Intensity with columns: " << totalIntensity << endl;
	cout << "Intensity without columns: " << originalIntensity << endl;
	cout << "Percentage increase: " << ((totalIntensity / originalIntensity) - 1.0) * 100;

	return 0;
}

