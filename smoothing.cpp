#include <iostream>
#include<cmath>
#include<cstdlib>
#include <vector>
#include<fstream>
#include<ctime>

using namespace std;

vector<vector<int> > Fgenerator(int slope, int n);
vector<vector<int> > OneSmoothing(vector<vector<int> > F, int n, int extra) ;
vector<vector<int> > ConfigGenerator(int k, int numberofconfigs, int n);
bool SuperHarmonicCheck(vector<vector<int> > Config, int n, int k);
int minimum (int i, int j, int k, int l);
void WriteMathematicaFile1 (vector<vector<int> > Function, int x_m, int y_m, string name, int n, int k);


int main () {

	int n, slope, k;

	/*cout << "Input n: " << endl;
	cin >> n;
	cout << "Input slope: " << endl;
	cin >> slope;*/

	//ASSUME K-1 IS DIVISIBLE BY 4
	n=6;
	k = 0;
	slope=2; 

	int i,j;

	vector<vector<int> > F(n+1);
	vector<vector<int> > LoweringList(n+1);
	for (i = 0 ; i < n+1 ; i++ ) {
		F[i].resize(n+1);
   		LoweringList[i].resize(n+1);
	}

	F = Fgenerator(slope, n);
	//F[n/2][n/2] = F[n/2][n/2] + (k-1)/4 + 1 - slope;

	
	

	
	bool terminate;
	int testtally = 1;
	do{
		WriteMathematicaFile1 (F, n, n, "smooth", n, 0);
		terminate = true;

		LoweringList = OneSmoothing(F,n, k);
		cout << LoweringList[1][1];
		//cout << testtally << endl;
		testtally++;
	
		for(i = 0; i < n+1; i++){
			for(j = 0; j < n+1; j++){
				if(LoweringList[i][j] != 0){
					terminate = false;
					break;
				}	
			}
			if (terminate == false){
				break;
			}
		}
		

		if(terminate == false){
			for(i = 0; i < n+1; i++){
				for(j = 0; j < n+1; j++){
					F[i][j] = F[i][j] - LoweringList[i][j];
					//cout << F[1][1] << endl;
				}
			}
		}
		//cout << F[1][1] << endl;

	} while(terminate == false);

	if(SuperHarmonicCheck(F,n,k) == false){
		cout << "NO SOLUTION, SLOPE NOT BIG ENOUGH" << endl;
	}
	
	WriteMathematicaFile1 (F, n, n, "smooth", n, 0);
	

	return 0;
}

vector<vector<int> > Fgenerator(int slope, int n){

	int i,j;

	vector<vector<int> > F(n+1);
	
	for (i = 0 ; i < n+1 ; i++) {
		F[i].resize(n+1);
	}
	
	for(i = 0; i < n+1; i++){
		for(j = 0; j < n+1; j++){
			F[i][j] = slope*minimum(i,j,n-i,n-j);
		}
	}
	
	return F;
	
}

vector<vector<int> > OneSmoothing(vector<vector<int> > F, int n, int extra) {
	//cout << F[1][1] << endl;
	int i,j,k;

	vector<vector<int> > LoweringList(n+1);
	vector<vector<int> > Config(n+1);
	vector<vector<int> > tempvector(n+1);

	for(i = 0; i < n+1; i++){
		LoweringList[i].resize(n+1);
		Config[i].resize(n+1);
		tempvector[i].resize(n+1);
	}	

	for(i = 0; i < n+1; i++){
		Config[0][i] = 0;
		Config[n][i] = 0;
		Config[i][0] = 0;
		Config[i][n] = 0;
	}

	for(i = 0; i < n+1; i++){
		for(j = 0; j < n+1; j++){
			LoweringList[i][j] = 0;
		}
	}

	int numberofconfigs = pow(2,(n-1)*(n-1));
	//int numberofconfigs = pow(2,32);
	cout << numberofconfigs << endl;
	//cout << "336 total" << endl;
	for(k = 0; k < numberofconfigs; k++){
if (k%100000 == 0){cout << "made it here" << endl;}
		Config = ConfigGenerator(k, numberofconfigs,n);
if (k%100000 == 0){cout << "made it here" << endl;}
		for(i = 0; i < n+1; i++){
				for(j = 0; j < n+1; j++){
					tempvector[i][j] = F[i][j] + Config[i][j];
				}
		}
//cout << tempvector[1][1] << "  " << tempvector[0][1] << endl;
//if (k%100000 == 0){cout << "made it here" << endl;}
		if(SuperHarmonicCheck(tempvector, n, extra) == true){
//if (k%100000 == 0){cout << "made it here" << endl;}
			//cout << "made it here" << endl;
			//cout << F[1][1] << "   " << F[0][1] << endl;
			for(i = 0; i < n+1; i++){
				for(j = 0; j < n+1; j++){
					if(LoweringList[i][j] == 0 && Config[i][j] == -1){
						//cout << "got here" << endl;
						LoweringList[i][j] = 1;
					}	
				}
			}
		}
	}
//cout << "made it here" << endl;

	return LoweringList;

	//WriteMathematicaFile1(ConfigGenerator(pow(2,0) +pow(2,6) +pow(2,8) +pow(2,12) +pow(2,16) +pow(2,18) +pow(2,30) +pow(2,32) +pow(2,36) +pow(2,40) +pow(2,42) +pow(2,48), numberofconfigs, n),n,n,"test",n,extra);

}

vector<vector<int> > ConfigGenerator(int k, int numberofconfigs, int n){

	int i,j;
	int level = numberofconfigs/2;

	vector<vector<int> > Config(n+1);
	for(i = 0; i < n+1; i++){
		Config[i].resize(n+1);
	}	

	for(i = 0; i < n+1; i++){
		Config[0][i] = 0;
		Config[n][i] = 0;
		Config[i][0] = 0;
		Config[i][n] = 0;
	}

	for(i = 0; i < n-1; i++){
		for(j = 0; j < n-1; j++){
			Config[i+1][j+1] = (-1)*(k/(level));
			k = k + Config[i+1][j+1]*level;
			level = level/2;
		}
	}

	return Config;

}

bool SuperHarmonicCheck(vector<vector<int> > Config, int n, int k){
	
	bool superharmonic = true;
	int i,j;
	for(i = 0; i < n-1; i++){
		for(j = 0; j < n-1; j++){
			if(-4*Config[i+1][j+1] + Config[i+2][j+1] + Config[i][j+1]+ Config[i+1][j+2] + Config[i+1][j] > 0 ){
				//cout << "not superharmonic!!! " << endl;
				superharmonic = false;
				break;
			}
		}
		if (superharmonic == false){
			break;
		}
	}

	if(-4*Config[n/2][n/2] + Config[n/2+1][n/2] + Config[n/2][n/2+1]+ Config[n/2][n/2-1] + Config[n/2-1][n/2] > -k){
		superharmonic = false;
	}

	//cout << superharmonic << endl;
	return superharmonic;
	
}

int minimum (int i, int j, int k, int l){

	int minimum;
	int winner1, winner2;

	if(i <= j){
		winner1 = i;
	} else{
		winner1 = j;
	}

	if(k <= l){
		winner2 = k;
	} else{
		winner2 = l;
	}

	if(winner1 <= winner2){
		minimum = winner1;
	} else{
		minimum = winner2;
	}

	return minimum;

}

void WriteMathematicaFile1 (vector<vector<int> > Function, int x_m, int y_m, string name, int n, int k){
	int i, j;
	string title;
	title = name + "_" + to_string(k) + ".nb";	
	//title = "sandpile.nb";

	ofstream outfile;
	
	outfile.open(title, ios::out);
	outfile << "CustomColorFunc[z_] := RGBColor[KroneckerDelta[z, 2] + KroneckerDelta[z, 3] + 0.75*KroneckerDelta[z, -1], KroneckerDelta[z, 1] + KroneckerDelta[z, 2] + 0.75*KroneckerDelta[z, -1], KroneckerDelta[z, 0] + KroneckerDelta[z, 1] + 0.75*KroneckerDelta[z, -1]]";

	outfile << endl;

	outfile << name + to_string(k) << " = ";

	for (j = n; j >= 0; j--) {
		if (j == n) {
			outfile << "{{";
		} else {
			outfile << "{";
		}

		for (i = 0; i < n + 1; i++) {
			outfile << Function[i][j];
			if (i != n){
				outfile << ",";
			} else if (j != 0){
				outfile << "},";
			} else if (j == 0){
				outfile << "}};";
			}
		}
	}

	outfile << endl;

	outfile << "ArrayPlot[" << name + to_string(k) << ", ColorFunctionScaling -> False,  ColorFunction -> CustomColorFunc]";
	
	outfile.close();
}
