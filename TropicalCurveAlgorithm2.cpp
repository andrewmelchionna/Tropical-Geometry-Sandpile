#include <iostream>
#include<cmath>
#include<cstdlib>
#include <vector>
#include<fstream>
#include<ctime>

using namespace std;

void firstgenvertexfinder(vector<vector<int> > &vertices, int numberofpoints, vector<vector<int> > points, int size);
bool pointsearch(int numberofpoints, vector<vector<int> > points, vector<int> vertextocheck);
vector<int> rectdirconverter(int dir);
vector<int> diagdirconverter(int dir);
	
int main() {

	int size = 60;
	int numberofpoints = size/2;
	int i,j;
	int maxgen = 6;

	vector<vector<int> > vertices(2);
	vector<vector<int> > points(2);
	points[0].resize(numberofpoints);
	points[1].resize(numberofpoints);
	vertices[0].resize(numberofpoints);
	vertices[1].resize(numberofpoints);
	

	vertices[0][0] = 28;
	vertices[1][0] = 31;

	points[0] = {28, 34, 16, 18, 24, 36, 36, 20, 21, 44, 35, 40, 39, 24, 22, 39, 27, 20, 42, 34, 36, 18, 37, 31, 40, 21, 25, 29, 28, 25};
	points[1] = {23, 19, 32, 37, 37, 20, 31, 24, 24, 33, 19, 29, 40, 43, 28, 26, 45, 23, 26, 42, 19, 40, 27, 42, 32, 39, 22, 17, 31, 28};



	firstgenvertexfinder(vertices, numberofpoints, points, size);


	for(i = 1; i < 5; i++){
		cout << vertices[0][i] << ", " << vertices[1][i] << endl; 
	}

	return 0;
}

void firstgenvertexfinder(vector<vector<int> > &vertices, int numberofpoints, vector<vector<int> > points, int size) {

	vector<int> internallimits(4);
	int i, j, dir;
	bool hit;
	vector<int> vertextocheck(2);

	//FIND INTERNAL LIMITS
	for(dir = 0; dir < 4; dir++){
		i = 1;
		hit = false;
		do{
			vertextocheck[0] = vertices[0][0] + i*rectdirconverter(dir)[0];
			vertextocheck[1] = vertices[1][0] + i*rectdirconverter(dir)[1];
			if(pointsearch(numberofpoints, points, vertextocheck) == true ){
				hit = true;
				internallimits[dir] = i;
			}
			i++;
			if(vertextocheck[0] < 0 || vertextocheck[0]  >= size || vertextocheck[1] < 0 || vertextocheck[1] >= size){
				hit = true;
				internallimits[dir] = i;
			}
		}while(hit == false);
	}
	
	

	//CHECK FOR EXTERNAL HITS
	
	bool externalhit;
	for(dir = 0; dir < 4; dir++){

		externalhit = false;
		i = 1;
		do{

			j=1;
			do{
				if(i + j < internallimits[dir]) {
					vertextocheck[0] = vertices[0][0] + i*diagdirconverter(dir)[0] + j*rectdirconverter(dir)[0];
					vertextocheck[1] = vertices[1][0] + i*diagdirconverter(dir)[1] + j*rectdirconverter(dir)[1];
					if(pointsearch(numberofpoints, points, vertextocheck) == true){
						externalhit = true;
						vertices[0][dir + 1] = vertices[0][0] + i*diagdirconverter(dir)[0];
						vertices[1][dir + 1] = vertices[1][0] + i*diagdirconverter(dir)[1];
					}
				}			

				if(externalhit == false && i + j < internallimits[(dir+1)%4]){
					vertextocheck[0] = vertices[0][0] + i*diagdirconverter(dir)[0] + j*rectdirconverter((dir+1)%4)[0];
					vertextocheck[1] = vertices[1][0] + i*diagdirconverter(dir)[1] + j*rectdirconverter((dir+1)%4)[1];
					if(pointsearch(numberofpoints, points, vertextocheck) == true){
						externalhit = true;
						vertices[0][dir + 1] = vertices[0][0] + i*diagdirconverter(dir)[0];
						vertices[1][dir + 1] = vertices[1][0] + i*diagdirconverter(dir)[1];
					}
				}
				j++;	
			} while((i + j < internallimits[dir] || i + j < internallimits[(dir+1)%4]) && externalhit == false);
			i++;
		} while(i < internallimits[dir] && i < internallimits[(dir+1)%4] && externalhit == false);

		if(externalhit == false) {
			vertices[0][dir + 1] = vertices[0][0];
			vertices[1][dir + 1] = vertices[1][0];
		}
//cout<< "hello"<< endl;	
	}


}

bool pointsearch(int numberofpoints, vector<vector<int> > points, vector<int> vertextocheck){
	int k;
	bool hit = false;
	for(k = 0; k < numberofpoints; k++){
		if(points[0][k] == vertextocheck[0]){
			if(points[1][k] == vertextocheck[1]){
				hit = true;
				break;
			}
		}
	}
	return hit;
}

vector<int> rectdirconverter(int dir){
	//0 POSITIVE J, 1 POSITIVE I, 2 NEGATIVE J, 3 NEGATIVE I
	vector<int> directionvector(2,0);
	if(dir == 3){
		directionvector[0] = -1;
	} else if(dir == 0){
		directionvector[1] = 1;
	} else if(dir == 1){
		directionvector[0] = 1;
	} else if(dir == 2){
		directionvector[1] = -1;
	}
	
	return directionvector;	

}

vector<int> diagdirconverter(int dir){
	//0 POS I/POS J, 1 POS I/NEG J, 2 NEG I/NEG J, 3 NEG I/ POS J
	vector<int> directionvector(2,0);
	if(dir == 0){
		directionvector[0] = 1;
		directionvector[1] = 1;
	} else if(dir == 1){
		directionvector[1] = -1;
		directionvector[0] = 1;
	} else if(dir == 2){
		directionvector[0] = -1;
		directionvector[1] = -1;
	} else if(dir == 3){
		directionvector[1] = 1;
		directionvector[0] = -1;
	}
	return directionvector;
}