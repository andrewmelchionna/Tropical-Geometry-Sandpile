#include <iostream>
#include<cmath>
#include<cstdlib>
#include <vector>
#include<fstream>
#include<ctime>

using namespace std;

int routine(int n, double a,double b,double c,double r,vector<vector<int> > v,double vnorm,vector<vector<int> >tile,vector<int> sandpilelegend,vector<vector<int> >tilecodeassignment,int tilesize,int tilewidth, double cptx, double cpty, double alpha, double con, int k, int trialnumber);
bool tileconfirmation(vector<int> tileassignment, vector<vector<int> > tile, int tilesize);
vector<int> goodpoints (vector<vector<int> > s, vector<vector<int> > tilecode, vector<vector<bool> > Domain, vector<int> sandpilelegend, int x_m, int y_m, vector<vector<int> > &goodpointsdomain);
void WriteMathematicaFile1 (vector<vector<int> > Function, int x_m, int y_m, string name, int n, int k);
void WriteMathematicaFile1point1 (vector<vector<int> > Function1, vector<vector<int> > Function2, int x_m, int y_m, string name, int n, int k, int totalchipsadded);
void WriteMathematicaFile1Point25 (vector<vector<double> > Function, int x_m, int y_m, string name, int n, double cptx, double cpty);
void WriteMathematicaFile1Point5 (vector<vector<double> > Function, int x_m, int y_m, string name, int n);
void WriteMathematicaFile2 (vector<vector<bool> > Function, int x_m, int y_m, string name);
void WriteMathematicaFile3 (int trials, vector<vector<int> > data, double alpha, bool overflow);
void WriteMathematicaFile4 (int x_m, int y_m, vector<vector<bool> > bf, vector<vector<int> > doesthisyvaluetouchfrombelow);
void Toppling(vector<vector<int> > &s, vector<vector<int> > &o, vector<vector<bool> > Domain, int n);
void TopplingOriginSink(vector<vector<int> > &s, vector<vector<int> > &o, vector<vector<bool> > Domain, int n);
void BoundaryTopple(vector<vector<int> > &s, vector<vector<int> > &o, vector<vector<bool> > Domain, vector<vector<bool> > bf, int n, int m);
vector<int> rgoodpointsfinder(int tilesize, vector<vector<int> > s, vector<int> sandpilelegend, vector<vector<int> > tile, int x_m, int y_m, int tilewidth, int r, vector<vector<int> > &goodpointsdomain, vector<vector<bool> > Domain, vector<vector<bool> > RrDomain, vector<vector<int> > v, vector<vector<int> > tilecodeassignment);
bool imagematching(int i, int j, int tilesize, vector<vector<int> > s, vector<int> sandpilelegend, vector<vector<int> > tile, int x_m, int y_m, int tilewidth, double r, vector<vector<int> > v, vector<vector<int> > tilecodeassignment);
int configkiller(int i, int j, int tilesize, vector<vector<int> > s, vector<int> sandpilelegend,  vector<vector<int> > tile, int tilewidth, vector<vector<int> > v, vector<vector<int> > tilecodeassignment);
void sandpilecheck(int i, int j, vector<vector<int> > s, int tilesize);
int tilecoordsassignment(int xdisp, int ydisp, int centerpointconfig, vector<vector<int> > tile, vector<vector<int> > v, int tilesize, vector<vector<int> > tilecodeassignment, bool debug);
void timeteller(int t, int n);
void tilecodeassignmentgenerator(int tilesize, int tilewidth, vector<vector<int> > tile, vector<vector<int> > &tilecodeassignment);
void LemmaTest(double vnorm, int tilesize, vector<vector<int> > s, vector<int> sandpilelegend, vector<vector<int> > tile, int x_m, int y_m, int tilewidth, vector<vector<int> > v, vector<vector<int> > tilecodeassignment, vector<vector<bool> > Domain, double bigR, vector<vector<int> > periodicityvectors, vector<vector<int> > &doesthisyvaluegivetouchfrombelow, int n, double bigeval, double a, double b, double c);
bool TouchFromBelowTest(int x_m, int y_m, double bigR, vector<vector<int> > v, vector<vector<double> > phiy, vector<vector<double> > o);
void odometergenerator(vector<vector<double> > &odometer, int L, int M, double a, double b, double c, bool altlinearterms);
int x(int i, int half);
int y(int j, int half);
void Laplacian(vector<vector<double> > &L, vector<vector<double> > v, int x_m, int y_m);
void touchingmapgenerator(vector<vector<vector<int> > > &touchingmap, int i, int j, int x_m, int y_m, vector<vector<double> > phiy, vector<vector<int> > v, vector<vector<bool> > Domain, int n, double bigeval, vector<vector<int> > &rangeoftouchingmap);
void touchingmapyfinder(vector<vector<vector<int> > >  touchingmap, int x_m, int y_m, int xyx, int xyy);
void randomsprinkle (int m, int k, vector<vector<int> > &s);
int littleboxlegend (int number, int m, bool x);
vector<bool> alphahighorlow(int trials, vector<vector<double> > data);


//WARNING: IN ALL BUT LemmaTest AND TouchFromBelowTest, o IS THE IDENTITY ELEMENT AND v ARE THE PERIODICITY VECTORS. IN THE AFOREMENTIONED FUNCTIONS, o AND v ARE AS THEY ARE IN THE LITERATURE

int main (){
	

	srand(time(NULL));
	
	double alpha;
	double con = 1.0;	
	int k;

	//cout << "k" << "  " << "o" << endl;

	int trials = 10;
	vector<int> trialparameters(trials);

	/*trialparameters[0] =200;
	trialparameters[1] = 300;
	trialparameters[2] = 400;
	trialparameters[3] = 500;
	trialparameters[4] = 600;
	trialparameters[5] = 700;
	
	trialparameters[6] = 800;
	trialparameters[7] = 900;
	trialparameters[8] = 1000;
	trialparameters[9] = 620;

	trialparameters[0] = 1100;
	trialparameters[1] = 1200;
	trialparameters[2] = 1300;
	trialparameters[3] = 1400;
	trialparameters[4] = 1500;
	trialparameters[5] = 1600;
	trialparameters[6] = 1700;*/
	
	/*int militarytime;
	cout << "Input current time, in military: ";
	cin >> militarytime;
	cout << endl << endl;*/

	int n;
	int i, j;
	double a, b, c, r, cptx, cpty;
	int tilesize, tilewidth;
	double vnorm;
	vector<vector<int> > v(2);

	for (i = 0 ; i < 2 ; i++ ) {
   		v[i].resize(2);
	}

	v[1][0] = 0;

	

	
	
	
	//AUTOMATIC INPUT

	
	
	/*cout << "Input k: " << endl;
	cin >> k;*/

	//cout << "Input n: " << endl;
	//cin >> trialparameters[0];

	/*cout << "Would you like to compute the fraction of good points?" << endl;
	cin >> yesno;

	cout << "Input Sqrt k:" << endl;
	cin >> trialparameters[0];
	cout << "Input cptx: " << endl;
	cin >> cptx;
	cout << "Input cpty: " << endl;
	cin >> cpty;*/


	/*cout << "Input matrix A" << endl
	<< "a:" << endl;
	cin >> a;
	cout << "b:" << endl;
	cin >> b;
	cout << "c:" << endl;
	cin >> c;*/
	/*cout << "input the x component of the first basis vector: " << endl;
	cin >> v[0][0];
	cout << "input the y component of the first basis vector: " << endl;
	cin >> v[0][1];
	cout << "input the y component of the second basis vector: " << endl;
	cin >> v[1][1];
	cout << "input the basepoint's x component (in true coordinates, origin in middle): " << endl;
	cin >> p[0];
	cout << "input the basepoint's y component (in true coordinates, origin in middle): " << endl;
	cin >> p[1];
	cout << "input vnorm: " << endl;
	cin >> vnorm;	
	cout << "input tile width (in x):" << endl;
	cin >> tilewidth;*/
	
	//MANUAL INPUT AND NECESSARY CALCULATIONS

	cptx = 0;
	cpty = 0;
	//n = 50;
	/*a = 1.14285714286;
	b = 0.357142857;
	c = .8928571429;*/
	v[0][0] = 2;
	v[0][1] = 1;
	v[1][1] = 2;
	vnorm = 2.828427;
	tilewidth = 2;
	r = 3*vnorm;

	

	//tilesize = floor(1/(a + c - 2.0) + 0.5);
	tilesize = 4;

	vector<vector<int> > tile(tilesize);
	vector<vector<int> > tilecodeassignment(tilewidth);

	
	for (i = 0 ; i < tilesize ; i++ ) {
   		tile[i].resize(2);
	}

	for (i = 0 ; i < tilewidth ; i++ ) {
   		tilecodeassignment[i].resize(tilewidth);
	}

	vector<int> sandpilelegend(tilesize);

		

	//n = 9;
	//a = 1.0156;
	//b = 0.125;
	//c = 1.0;

	
		tile[0][0] = 0;
		tile[0][1] = 0;
		tile[1][0] = 1;
		tile[1][1] = 0;
		tile[2][0] = 0;
		tile[2][1] = 1;
		tile[3][0] = 1;
		tile[3][1] = 1;
	
	
		sandpilelegend[0] = 3;
		sandpilelegend[1] = 3;
		sandpilelegend[2] = 3;
		sandpilelegend[3] = 0;

		
		tilecodeassignmentgenerator(tilesize, tilewidth, tile, tilecodeassignment);

	
		
	//cout<< "hello" << endl;
	
	vector<vector<int> > overflowdata(9);
	//vector<vector<double> > halfwayoutdata(trials);
	for (i = 0 ; i < 9 ; i++ ) {
   		overflowdata[i].resize(trials);
		//halfwayoutdata[i].resize(2);
	}

	/*for (i = 0 ; i < trials ; i++ ) {
   		overflowdata[i][0] = trialparameters[i];
		halfwayoutdata[i][0] = trialparameters[i];
	}*/

	int tempdata;

	
	clock_t time;
	time = clock();

	vector<double> alphavector(6);
	alphavector[0] = 1.85;
	alphavector[1] = 1.95;
	alphavector[2] = 2.1;
	alphavector[3] = 1.83;
	alphavector[4] = 1.87;
	alphavector[5] = 1.75;

	int alphaindex = 0;
	
	
	//cout<< "hello" << endl;
	do{
		n = 200;
		j=0;
		do{	
			for(i = 0; i < trials; i++){
				//overflowdata[j][i] = routine(n, a,b,c,r,v,vnorm,tile,sandpilelegend,tilecodeassignment,tilesize,tilewidth,cptx,cpty,alpha,con,k,i);
				overflowdata[j][i] = routine(n, a,b,c,r,v,vnorm,tile,sandpilelegend,tilecodeassignment,tilesize,tilewidth,cptx,cpty,alphavector[alphaindex],con,k,i);
				//cout << "(n,alpha,trial) = (" << n << ", " << alpha << ", " << i << ")" << endl;
				cout << "(n,alpha,trial) = (" << n << ", " << alphavector[alphaindex] << ", " << i << ")" << endl;
				time = clock() - time;
				timeteller(time, trialparameters[i]);
			}
			n = n + 100;
			j++;
		} while (n < 601);
		//} while(n < 1001);
		WriteMathematicaFile3(trials,overflowdata,alphavector[alphaindex], true);
		//alpha = alpha + 0.1;
	//} while(alpha <= 1.91);
		alphaindex++;
	} while(alphaindex < 6);
	
	
		

	
	/*WriteMathematicaFile3(trials,halfwayoutdata,alpha, false);*/
	
	return 0;
}
int routine(int n, double a,double b,double c,double r,vector<vector<int> > v,double vnorm,vector<vector<int> >tile,vector<int> sandpilelegend,vector<vector<int> >tilecodeassignment,int tilesize,int tilewidth, double cptx, double cpty, double alpha, double con, int k, int trialnumber) {

	
	//cout << "n equals: " << n << endl;
	double r_to_use = 8.5;

	//COMPUTE BOUNDS ON ELLIPSE AND EVALS
	int i, j;
	int x_m, y_m;
	double det;
	double x_max, y_max;
	double bigeval, smalleval;

	smalleval = (a+c-sqrt((a-c)*(a-c) + 4*b*b))/2;
	bigeval = (a+c+sqrt((a-c)*(a-c) + 4*b*b))/2;
	

	det = a*c - b*b;
	x_max = n*sqrt(2*c/det);
	y_max = n*sqrt(2*a/det); 
	

	//x_m = ceil(x_max);
	//y_m = ceil(y_max);


	x_m = ceil(x_max) + 3;
	y_m = ceil(y_max) + 3;

	//cout << "x_m is " << x_m << "   y_m is " << y_m << endl;
	
	//INITIALIZE DOMAIN VECTOR



	vector<vector<bool> > Domain((n+1));
	vector<vector<double> > odometer((n+1));

	
	for (i = 0 ; i < n+1 ; i++ ) {
   		Domain[i].resize(n+1);
		odometer[i].resize(n+1);
	}
		
	//odometergenerator(odometer, n+1, n+1, a, b, c, false);

	//BUILD DOMAIN VECTOR
	int x, y;
	for (i = 0; i < n+1; i++) {
			//x = i - x_m;
		for (j = 0; j < n+1; j++) {
			//y = j - y_m;
		
			
			if (i>0 && i < n && j > 0 && j < n){
			//if(odometer[i][j] < n*n){
				Domain[i][j] = true;
	
			} else {
				Domain[i][j] = false;
			}
	
		}
	}

	//INITIALIZE SANDPILE
	
	bool allthrees = true;
	
	

	vector<vector<int> > s(n+1);
	//vector<vector<int> > p(2*x_m + 1);
	//vector<vector<double> > ptest(2*x_m + 1);
	
	for (i = 0 ; i < n+1 ; i++ ) {
   		s[i].resize(n+1);
		//p[i].resize(2*y_m+1);
		//ptest[i].resize(2*y_m+1);
	}

	for(i = 0; i < n+1; i++) {
		for (j = 0; j < n + 1; j++){
			if(Domain[i][j] == true && allthrees == true) {
				s[i][j] = 3;
			} else {
				s[i][j] = 0;
			}
		}
	}

	//SPRINKLE
	
	//cout << "alpha equals: " << alpha << endl;
	//int k = round(con*pow(n, alpha));
	//k = pow(n/4,1.8);

	//s[3*n/4][3*n/4] = s[3*n/4][3*n/4] + k;
	/*for(i=0; i < 5; i++){
		for(j=0; j < 5; j++){
			if(i == 1 && j == 2){
				s[n/4+(i*n)/8-5][n/4+(j*n)/8+5] = 4;
			} else{
				s[n/4+(i*n)/8][n/4+(j*n)/8] = 4;
			}
		}
	}*/
	
	
	
	int totalchipsadded = 0;
	vector<vector<int> > pointsadded (2);
	for(i=0; i < 2; i++){
		pointsadded[i].resize(0);
	}

	//RANDOM ON EACH SQUARE. DENSITY IS FORCED TO BE UNIFORM
	
/*	
	i = 0;
	do{
	k=0;
		do{
			j = (rand()%(n/2 - 2*i)) - (n/4 - i);
			if(i%4 == 0){
				s[n/2+j][i+n/4]++;
				pointsadded[0].push_back(n/2+j);
				pointsadded[1].push_back(i + n/4);
			}else if(i%4 == 2){
				s[n/2+j][3*n/4 - i]++;
				pointsadded[0].push_back(n/2+j);
				pointsadded[1].push_back(3*n/4 - i);
			}else if(i%4==1){
				s[i+n/4][n/2+j]++;
				pointsadded[0].push_back(i + n/4);
				pointsadded[1].push_back(n/2+j);
			} else if(i%4 == 3){
				s[3*n/4 - i][n/2+j]++;
				pointsadded[0].push_back(3*n/4 - i);
				pointsadded[1].push_back(n/2+j);
			}

			k++;
			totalchipsadded++;
		//} while(k <  2*(n/4 - i)/sqrt(((double) (n/4))));
		} while(k < 1);
		i++;
	} while(i < n/4);*/
	

	//TRUE RANDOM

	
	k=0;
	do{
		i = rand()%(n/2);
		j = rand()%(n/2);	
		s[n/4 + i][n/4 + j]++;
		pointsadded[0].push_back(n/4 + i);
		pointsadded[1].push_back(n/4 + j);
		k++;
		totalchipsadded++;
	}while(k < pow(n,alpha));
	//}while(k < 5);

	/*for(i=n/4; i < n/2; i++){
		s[3*n/8][i] = floor(sqrt(n));
	}
	
	s[n/4][3*n/8] = floor(sqrt(n));*/

	/*s[round((9*((double) n))/16)][round((7*((double) n))/16)] = 4;
	s[round((7*((double) n))/16)][round((7*((double) n))/16)] = 4;
	s[round((7*((double) n))/16)][round((9*((double) n))/16)] = 4;
	s[round((9*((double) n))/16)][round((9*((double) n))/16)] = 4; 
	s[round((10*((double) n))/16)][round((10*((double) n))/16)] = 4;
	s[round((6*((double) n))/16)][round((6*((double) n))/16)] = 4;
	s[round((6*((double) n))/16)][round((10*((double) n))/16)] = 4;
	s[round((10*((double) n))/16)][round((6*((double) n))/16)] = 4; 
	s[n/2][n/2] = 4;*/

	/*s[round((17*((double) n))/32)][round((15*((double) n))/32)] = 4;
	s[round((15*((double) n))/32)][round((15*((double) n))/32)] = 4;
	s[round((15*((double) n))/32)][round((17*((double) n))/32)] = 4;
	s[round((17*((double) n))/32)][round((17*((double) n))/32)] = 4; 
	s[round((18*((double) n))/32)][round((18*((double) n))/32)] = 4;
	s[round((14*((double) n))/32)][round((14*((double) n))/32)] = 4;
	s[round((14*((double) n))/32)][round((18*((double) n))/32)] = 4;
	s[round((18*((double) n))/32)][round((14*((double) n))/32)] = 4; 
	s[n/2][n/2] = 4;*/

	/*cout << round((511*((double) n))/1024) << endl;
	cout << round((513*((double) n))/1024) << endl;
	cout << round((17*((double) n))/32) << endl;
	cout << round((18*((double) n))/32) << endl;*/
	
	/*s[n/2+1][n/2] = 4;
	s[n/2-1][n/2] = 4;
	s[n/2][n/2-1] = 4;
	s[n/2][n/2+1] = 4;
	s[n/2+2][n/2] = 4;
	s[n/2-2][n/2] = 4;
	s[n/2][n/2-2] = 4;
	s[n/2][n/2+2] = 4;

	s[n/2-1][n/2+1] = 4;
	s[n/2+1][n/2-1] = 4;
	s[n/2-1][n/2-1] = 4;
	s[n/2+1][n/2+1] = 4;

	s[n/4][5*n/8] = s[n/2][n/2] + k;
	s[n/4][3*n/8] = s[n/2][n/2] + k;
	s[5*n/8][n/4] = s[n/2][n/2] + k;
	s[3*n/8][n/4] = s[n/2][n/2] + k;
	s[3*n/4][5*n/8] = s[n/2][n/2] + k;
	s[3*n/4][3*n/8] = s[n/2][n/2] + k;
	s[5*n/8][3*n/4] = s[n/2][n/2] + k;
	s[3*n/8][3*n/4] = s[n/2][n/2] + k;*/
	

	/*s[n/2][n/2] = -1;
	s[n/2+1][n/2] = 4;
	s[n/2-1][n/2] = 4;
	s[n/2][n/2+1] = 4;
	s[n/2][n/2-1] = 4;*/
	
	
	
	/*s[101][101] = 4;
	s[99][101] = 4;
	s[101][99] = 4;
	s[99][99] = 4;
	s[102][102] = 4;
	s[102][98] = 4;
	s[98][98] = 4;
	s[98][102] = 4;
	s[100][100] = 8;*/

	/*s[round((33*((double) n))/64)][round((31*((double) n))/64)] = 4;
	s[round((31*((double) n))/64)][round((31*((double) n))/64)] = 4;
	s[round((31*((double) n))/64)][round((33*((double) n))/64)] = 4;
	s[round((33*((double) n))/64)][round((33*((double) n))/64)] = 4; 
	s[round((34*((double) n))/64)][round((34*((double) n))/64)] = 4;
	s[round((30*((double) n))/64)][round((30*((double) n))/64)] = 4;
	s[round((30*((double) n))/64)][round((34*((double) n))/64)] = 4;
	s[round((34*((double) n))/64)][round((30*((double) n))/64)] = 4; 
	s[n/2][n/2] = 4;*/


	/*s[round((65*((double) n))/128)][round((63*((double) n))/128)] = 4;
	s[round((63*((double) n))/128)][round((63*((double) n))/128)] = 4;
	s[round((63*((double) n))/128)][round((65*((double) n))/128)] = 4;
	s[round((65*((double) n))/128)][round((65*((double) n))/128)] = 4; 
	s[round((66*((double) n))/128)][round((66*((double) n))/128)] = 4;
	s[round((62*((double) n))/128)][round((62*((double) n))/128)] = 4;
	s[round((62*((double) n))/128)][round((66*((double) n))/128)] = 4;
	s[round((66*((double) n))/128)][round((62*((double) n))/128)] = 4; 
	s[n/2][n/2] = 4;*/

	/*s[round((513*((double) n))/1024)][n/2] = 4;
	s[round((511*((double) n))/1024)][n/2] = 4;
	s[n/2][round((513*((double) n))/1024)] = 4;
	s[n/2][round((511*((double) n))/1024)] = 4; 
	s[round((514*((double) n))/1024)][n/2] = 4;
	s[round((510*((double) n))/1024)][n/2] = 4;
	s[n/2][round((514*((double) n))/1024)] = 4;
	s[n/2][round((510*((double) n))/1024)] = 4;
	s[n/2][n/2] = 4;*/

	
	//randomsprinkle(n,k,s);
	
	
	

	//INITIALIZE ODOMETER
	
	vector<vector<int> > o(n+1);
	vector<vector<double> > rescaledo(n+1);
	for (i = 0 ; i < n+1 ; i++ ) {
   		o[i].resize(n+1);
		rescaledo[i].resize(n+1);
	}
	
	//INITIALIZE BOUNDARY
	
	vector<vector<bool> > bf(n+1);
	for (i = 0 ; i < n+1 ; i++ ) {
   		bf[i].resize(n+1);
	}

	for(i = 0; i < n + 1; i++) {
		for (j = 0; j < n + 1; j++){
			bf[i][j] = false;
			if (Domain[i][j] == false) {
				if (i < n && Domain[i+1][j] == true){
					bf[i][j] = true;
				} else if (i > 0 && Domain[i-1][j] == true){
					bf[i][j] = true;
				} else if (j < n && Domain[i][j+1] == true){
					bf[i][j] = true;
				} else if (j > 0 && Domain[i][j-1] == true){
					bf[i][j] = true;
				}
			}
		}
	}
	


	//MAIN ROUTINE

	int m = 10*n*n + n;
	
	//TOPPLE THE BOUNDARY M TIMES
	//BoundaryTopple(s, o, Domain, bf, n, m);
	
	
	//TOPPLE INTERIOR UNTIL STABILITY
	Toppling(s, o, Domain, n);
	//TopplingOriginSink(s, o, Domain, n);

	/*for(i = 0; i < n + 1; i++){
		for(j = 0; j < n + 1; j++){
			if(Domain[i][j] == 0){
				//o[i][j] = m;
				s[i][j] = 4;
			}
		}
	}*/
	
	/*//GET MAXIMAL ODOMETER VALUE
	
	int max = 0;
	double rescaledmax;
	int ispec, jspec;

	for(i = 0; i < n + 1; i++){
		for(j = 0; j < n + 1; j++){
			if(o[i][j] > max){
				max = o[i][j];
				ispec = i;
				jspec = j;
			}
		}
	}

	rescaledmax = ((double) max)/(n*n);

	cout << "Rescaled maximum of rescaled odometer is: " << rescaledmax << endl;
	cout << "location is: " << ispec << " " << jspec << endl;*/

	//RESCALE ODOMETER

	/*for(i = 0; i < n + 1; i++){
		for(j = 0; j < n + 1; j++){
			rescaledo[i][j] = ((double) o[i][j])/(n*n);
		}
	}*/

	//GET NUMBER OF CHIPS THAT FELL OFF

	int overflow = 0;

	for(i = 0; i < n+1; i++){
		overflow = overflow + s[0][i] + s[n][i] + s[i][0] + s[i][n];
	}

	/*int interior = 0;
	for (i = 0 ; i < n+1 ; i++ ) {
		for (j = 0 ; j < n+1 ; j++ ) {
			if(Domain[i][j] == 1){
			interior = interior + s[i][j];
			}
		}
	}*/

	//cout << "overflow is: " << overflow  << endl;
	

	
	//cout << "log overflow over log n is: " << log((double) overflow )/log(n) << endl;

	//GET ODOMETER HALFWAY OUT
	//cout << "Rescaled odometer halfway out is: " << rescaledo[round(((double) n) * 0.75)][round(((double) n)/2.0)] << endl;
	//cout << "Log odometer halfway out over log n is: " << log((double) o[round(((double) n) * 0.75)][round(((double) n)/2.0)] )/log(n) << endl;
	



	
	//CHECK IF RECURRENT. IF NOT, TOPPLE THE BOUNDARY ONCE. CHECK AGAIN. CONTINUE UNTIL RECURRENT
	bool RecurrencyCheck = false;

/*
	vector<vector<int> > test(2*x_m+1);
	
	for (i = 0 ; i < 2*x_m+1 ; i++ ) {
   		test[i].resize(2*y_m+1);
	}


	
	while(RecurrencyCheck == false) {
		RecurrencyCheck = RecurrencyTest(s, o, test, Domain, bf, x_m, y_m);
	} */





	//cout << k << "  " << o[n/2][n/2] << endl;



	//WriteMathematicaFile1(s, n, n, "sandpile1", n, k);
	//WriteMathematicaFile1Point25(rescaledo,x_m,y_m,"IdentityElement",n,cptx,cpty);
	//WriteMathematicaFile1point1(s, pointsadded, n, n, "sandpile", n, trialnumber, totalchipsadded);
	//WriteMathematicaFile1(o, n, n, "IdentityElement", n,k);
	
	//TOPPLE THE BOUNDARY M TIMES
	//BoundaryTopple(s, o, Domain, bf, n, 1);
	
	
	//TOPPLE INTERIOR UNTIL STABILITY
	//Toppling(s, o, Domain, n);

	//WriteMathematicaFile1(s, n, n, "sandpile2", n,cptx,cpty);

	//return 0;

	
	//resultsvec[0] = (log((double) overflow )/log(n));

	//resultsvec[1] = log((double) o[round(((double) n) * 0.75)][round(((double) n)/2.0)] )/log(n);

	
	//cout << "n equals: " << n << endl << /*"n/4^1.8 " << ( (double) n)*sqrt((double) n)/8.0 << endl <<*/ "total chips added equals: " << k << endl << "overflow over n squared is: " << resultsvec[0] << endl;

	return overflow;

}



bool tileconfirmation(int xdisp, int ydisp, vector<vector<int> > tile, int tilesize){
	bool tileconfirmed = false;
	int i;
	
	for (i = 0; i < tilesize; i++){
		if(xdisp == tile[i][0] && ydisp == tile[i][1]){
			tileconfirmed = true;
			break;
		}
	}
	return tileconfirmed;
}


vector<int> goodpoints (vector<vector<int> > s, vector<vector<int> > tilecode, vector<vector<bool> > Domain, vector<int> sandpilelegend, int x_m, int y_m, vector<vector<int> > &goodpointsdomain){
	int num = 0;
	int den = 0;
	int i,j,k;

	
	for(i = 0; i < 2*x_m+1; i++){
		for(j = 0; j < 2*y_m+1; j++){
			goodpointsdomain[i][j] = 1;
			if(Domain[i][j] == 1){
				den++;
				if(s[i][j] == sandpilelegend[tilecode[i][j]]){
					num++;
					goodpointsdomain[i][j] = 0;
				}
			}
		}
	}

	vector<int> numdengoodpoints(2);
	numdengoodpoints[0] = num;
	numdengoodpoints[1] = den;
	
	return numdengoodpoints;
}
			

void WriteMathematicaFile1 (vector<vector<int> > Function, int x_m, int y_m, string name, int n, int k){
	int i, j;
	string title;
	title = name + ".nb";	
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

	outfile << "ListPlot3D[" << name + to_string(k) << "]";
	
	outfile.close();
}

void WriteMathematicaFile1point1 (vector<vector<int> > Function1, vector<vector<int> > Function2, int x_m, int y_m, string name, int n, int k, int totalchipsadded){
	int i, j;
	string title;
	title = name + "_" + to_string(k+1) + ".nb";	
	//title = "sandpile.nb";

	ofstream outfile;
	
	outfile.open(title, ios::out);
	outfile << "CustomColorFunc[z_] := RGBColor[KroneckerDelta[z, 2] + KroneckerDelta[z, 3] + 0.75*KroneckerDelta[z, -1], KroneckerDelta[z, 1] + KroneckerDelta[z, 2] + 0.75*KroneckerDelta[z, -1], KroneckerDelta[z, 0] + KroneckerDelta[z, 1] + 0.75*KroneckerDelta[z, -1]]";

	outfile << endl;

	outfile << name + to_string(k+1) << " = ";

	for (j = n; j >= 0; j--) {
		if (j == n) {
			outfile << "{{";
		} else {
			outfile << "{";
		}

		for (i = 0; i < n + 1; i++) {
			outfile << Function1[i][j];
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

	//outfile << "points" + to_string(k+1) << " = ";

	for (j = 0; j < totalchipsadded; j++) {
		if (j == 0) {
			outfile << "{{";
		} else {
			outfile << "{";
		}

		outfile << Function2[0][j] + 1 << ", " << Function2[1][j] + 1 << "}";

		if(j != totalchipsadded - 1){
			outfile << ", ";
		} else{
			outfile << "};";
		}
	}

	outfile << endl;

	//outfile << "ArrayPlot[" << name + to_string(k) << ", ColorFunctionScaling -> False,  ColorFunction -> CustomColorFunc]";
	outfile << "ArrayPlot[" << name + to_string(k+1) << ", ColorFunctionScaling -> False, ColorFunction -> CustomColorFunc]";
	
	outfile.close();
}

void WriteMathematicaFile1Point25 (vector<vector<double> > Function, int x_m, int y_m, string name, int n, double cptx, double cpty){
	int i, j;
	string title;
	title = name + "_" + to_string(n) + "_" + to_string(cptx) +"_" + to_string(cpty) + ".nb";	
	//title = "sandpile.nb";

	ofstream outfile;
	
	outfile.open(title, ios::out);

	outfile << name << " = ";

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

	outfile << "ListPlot3D[IdentityElement]";
	outfile << endl;
	
	outfile.close();
}




void WriteMathematicaFile1Point5 (vector<vector<double> > Function, int x_m, int y_m, string name, int n){
	int i, j;
	string title;
	title = "patternmatching_firstmatrix_biggerconstant" + name + "_" + to_string(n) + ".nb";	
	//title = "sandpile.nb";

	ofstream outfile;
	
	outfile.open(title, ios::out);
	outfile << "CustomColorFunc[z_] := RGBColor[KroneckerDelta[z, 2] + KroneckerDelta[z, 3] + 0.75*KroneckerDelta[z, -1], KroneckerDelta[z, 1] + KroneckerDelta[z, 2] + 0.75*KroneckerDelta[z, -1], KroneckerDelta[z, 0] + KroneckerDelta[z, 1] + 0.75*KroneckerDelta[z, -1]]";

	outfile << endl;

	outfile << name << " = ";

	for (j = 2*y_m-1; j >= 1; j--) {
		if (j == 2*y_m-1) {
			outfile << "{{";
		} else {
			outfile << "{";
		}

		for (i = 0; i < 2*x_m; i++) {
			outfile << Function[i][j];
			if (i != 2*x_m - 1){
				outfile << ",";
			} else if (j != 1){
				outfile << "},";
			} else if (j == 1){
				outfile << "}};";
			}
		}
	}

	outfile << endl;

	outfile << "ArrayPlot[" << name << ", ColorFunctionScaling -> False,  ColorFunction -> CustomColorFunc]";
	
	outfile.close();
}


void WriteMathematicaFile2 (vector<vector<bool> > Function, int x_m, int y_m, string name){
	int i, j;
	string title;
	title = "patternmatching_firstmatrix_biggerconstant" + name + ".nb";	

	ofstream outfile;
	
	outfile.open(title, ios::out);
	outfile << "CustomColorFunc[z_] := RGBColor[KroneckerDelta[z, 2] + KroneckerDelta[z, 3] + 0.75*KroneckerDelta[z, -1], KroneckerDelta[z, 1] + KroneckerDelta[z, 2] + 0.75*KroneckerDelta[z, -1], KroneckerDelta[z, 0] + KroneckerDelta[z, 1] + 0.75*KroneckerDelta[z, -1]]";

	outfile << endl;

	outfile << name << " = ";

	for (j = 2*y_m; j >= 0; j--) {
		if (j == 2*y_m) {
			outfile << "{{";
		} else {
			outfile << "{";
		}

		for (i = 0; i < 2*x_m + 1; i++) {
			outfile << Function[i][j];
			if (i != 2*x_m){
				outfile << ",";
			} else if (j != 0){
				outfile << "},";
			} else if (j == 0){
				outfile << "}};";
			}
		}
	}

	outfile << endl;

	outfile << "ArrayPlot[" << name << ", ColorFunctionScaling -> False,  ColorFunction -> CustomColorFunc]";
	
	outfile.close();
}

void WriteMathematicaFile3 (int trials, vector<vector<int> > data, double alpha, bool overflow){
		
	int i,j;
	string title;
	
	if(overflow == true){
		title = "overflow" + to_string(alpha) + ".nb";
	} else{
		title = "explosion_halfwayout_alpha_" + to_string(alpha) + ".nb";
	}

	ofstream outfile;
	
	outfile.open(title, ios::out);
	
	outfile << "(* this data is overflow *)" << endl;
	
	outfile << "data  = {" ;
	for(j=0; j < 5; j++){
		outfile << "{";
		for(i = 0; i < trials; i++){
			outfile << data[j][i];
			if(i != trials -1){
				outfile << ", ";
			}else{
				outfile << "}";
			}
		}
		if(j != 4){
			outfile << ", ";
		}else{
			outfile << "};" << endl;
		}
	}


	outfile << "averagedata = {";
	for(j=0; j < 5; j++){
		outfile << "{" << to_string(200 + j*100) << ", 0}";
		if(j!=4){
			outfile << ", ";
		}
	}
	outfile << "};" << endl;
	outfile << "trials = " << to_string(trials) << ";" << endl;
	outfile << "Do[total = 0; Do[total = total + data[[j,i]], {i,1,trials}]; averagedata[[j,2]] = N[total/trials], {j,1,5}];" << endl;	

	outfile << "logaveragedata = ConstantArray[0,{5,2}];" << endl;
	outfile << "Do[Do[logaveragedata[[j,i]] = Log[averagedata[[j,i]]],{j,1,5}],{i,1,2}];" << endl;
	outfile << "ListPlot[logaveragedata]" << endl;
	outfile << "FindFit[logaveragedata, {a*x + b}, {a, b}, x]";

	
	
	outfile.close();
}

void WriteMathematicaFile4 (int x_m, int y_m, vector<vector<bool> > bf, vector<vector<int> > doesthisyvaluetouchfrombelow){
		
	int i,j;
	
	vector<vector<int> > composite (2*x_m+1);
	for(i = 0; i < 2*x_m+1; i++){
		composite[i].resize(2*y_m+1);
	}
	

	for(i = 0; i < 2*x_m+1; i++){
		for(j = 0; j < 2*y_m +1; j++){
			if(bf[i][j] == true){
				composite[i][j] = 1 + doesthisyvaluetouchfrombelow[i][j];
			} else{
				composite[i][j] = doesthisyvaluetouchfrombelow[i][j];
			}
		}
	}

	ofstream outfile;

	
	outfile.open("doesthisyvaluetouchfrombelow.nb", ios::out);

	outfile << "CustomColorFunc[z_] := RGBColor[KroneckerDelta[z, 2] + KroneckerDelta[z, 3] + 0.75*KroneckerDelta[z, -1], KroneckerDelta[z, 1] + KroneckerDelta[z, 2] + 0.75*KroneckerDelta[z, -1], KroneckerDelta[z, 0] + KroneckerDelta[z, 1] + 0.75*KroneckerDelta[z, -1]]";

	outfile << endl;

	outfile << "touchfrombelowdata=";

	for (j = 2*y_m; j >= 0; j--) {
		if (j == 2*y_m) {
			outfile << "{{";
		} else {
			outfile << "{";
		}

		for (i = 0; i < 2*x_m + 1; i++) {
			outfile << composite[i][j];
			if (i != 2*x_m){
				outfile << ",";
			} else if (j != 0){
				outfile << "},";
			} else if (j == 0){
				outfile << "}};";
			}
		}
	}


	outfile << "ArrayPlot[touchfrombelowdata, ColorFunctionScaling -> False,  ColorFunction -> CustomColorFunc]";
	
	outfile.close();
}



void Toppling(vector<vector<int> > &s, vector<vector<int> > &o, vector<vector<bool> > Domain, int n){
	int i, j;
	bool clear = false;
	while (clear == false) {
		for(i = 0; i < n + 1; i++) {
			if (i == 0) {
				clear = true;
			}

			for (j = 0; j < n + 1; j++){
				if (s[i][j] >= 4 && Domain[i][j] == 1){
					o[i][j]++;
					clear = false;
					s[i][j] = s[i][j] - 4;
					s[i+1][j]++;
					s[i-1][j]++;
					s[i][j+1]++;		
					s[i][j-1]++;
				}
			}
		}
	}
}	

void TopplingOriginSink(vector<vector<int> > &s, vector<vector<int> > &o, vector<vector<bool> > Domain, int n){
	int i, j;
	bool clear = false;
	while (clear == false) {
		for(i = 0; i < n + 1; i++) {
			if (i == 0) {
				clear = true;
			}

			for (j = 0; j < n + 1; j++){
				//if (s[i][j] >= 4 && Domain[i][j] == 1 && (i != n/2 || j != n/2)){
				if (s[i][j] >= 4 && Domain[i][j] == 1 && !(i <= 9*n/16 && i >= 7*n/16 && j <= 9*n/16 && j >= 7*n/16)){
					o[i][j]++;
					clear = false;
					s[i][j] = s[i][j] - 4;
					s[i+1][j]++;
					s[i-1][j]++;
					s[i][j+1]++;		
					s[i][j-1]++;
				}
			}
		}
	}
}	

void BoundaryTopple(vector<vector<int> > &s, vector<vector<int> > &o, vector<vector<bool> > Domain, vector<vector<bool> > bf,  int n, int m){
	
	int i,j;
	for (i = 0; i < n + 1; i++) {
		for (j = 0; j < n + 1; j++) {
			if (bf[i][j] == true){
					
				o[i][j] = o[i][j] + m;				

				if (i < n  && Domain[i+1][j] == true) {
					s[i+1][j] = s[i+1][j] + m;
				}
				if (i > 0 && Domain[i-1][j] == true) {
					s[i-1][j] = s[i-1][j] + m;
				}
				if (j < n  && Domain[i][j+1] == true) {
					s[i][j+1] = s[i][j+1] + m;
				}
				if (j > 0 && Domain[i][j-1] == true) {
					s[i][j-1] = s[i][j-1] + m;
				}
			}
		}
	}
}

vector<int> rgoodpointsfinder(int tilesize, vector<vector<int> > s, vector<int> sandpilelegend, vector<vector<int> > tile, int x_m, int y_m, int tilewidth, int r, vector<vector<int> > &goodpointsdomain, vector<vector<bool> > Domain, vector<vector<bool> > RrDomain, vector<vector<int> > v, vector<vector<int> > tilecodeassignment){

	int i,j;
	int circletotalpoints = 0;
	int circlergoodpoints = 0;
	int ellipsetotalpoints = 0;
	int ellipsergoodpoints = 0;

	for(i = 0; i < 2*x_m + 1; i++){
		for(j=0; j<2*y_m + 1; j++){

			goodpointsdomain[i][j] = 0;

			if(Domain[i][j] == 1){
				ellipsetotalpoints++;
				if(RrDomain[i][j] == 1){
					circletotalpoints++;
					if(imagematching(i, j, tilesize, s, sandpilelegend, tile, x_m, y_m, tilewidth, double(r), v, tilecodeassignment) == true){
						goodpointsdomain[i][j] = 1;
						circlergoodpoints++;
						ellipsergoodpoints++;
					} 
				}else{
					if(imagematching(i, j, tilesize, s, sandpilelegend, tile, x_m, y_m, tilewidth, double(r), v, tilecodeassignment) == true){
						//goodpointsdomain[i][j] = 1;
						ellipsergoodpoints++;
					}
				}
			}

		}
	}

	vector<int> points(4);
	points[0] = circlergoodpoints;
	points[1] = circletotalpoints;
	points[2] = ellipsergoodpoints;
	points[3] = ellipsetotalpoints;
	
	return points;
	
}

bool imagematching(int i, int j, int tilesize, vector<vector<int> > s, vector<int> sandpilelegend, vector<vector<int> > tile, int x_m, int y_m, int tilewidth, double r, vector<vector<int> > v, vector<vector<int> > tilecodeassignment){

	int k,m;
	bool thispointisrgood = true;

	//BUILD AN r DOMAIN
	vector<vector<bool> > rDomain(2*x_m + 1);
	for(k = 0; k < 2*x_m + 1; k++){
		rDomain[k].resize(2*y_m+1);
	}
	
	for(k = 0; k < 2*x_m + 1; k++){
		for(m = 0; m < 2*y_m + 1; m++){
			rDomain[k][m] = false;
		}
	}

	for(k = 0; k < 2*x_m + 1; k++){
		for(m = 0; m < 2*y_m + 1; m++){
			if(pow(i - k,2)+pow(j - m,2) < pow(r,2)){
				rDomain[k][m] = true;
			}
		}
	}


	//FIRST RUN CONFIG KILLER TO GET WHICH CONFIG WE ARE WORKING WITH
	int config;

	config = configkiller(i, j, tilesize, s, sandpilelegend, tile, tilewidth, v, tilecodeassignment);

	if(config != tilesize){
		//THAT IS, THERE ACTUALLY EXISTS A CONFIG THAT WE CAN CHECK
		for(k = 0; k < 2*x_m + 1; k++){
			for(m = 0; m < 2*y_m + 1; m++){
				if(rDomain[k][m] == true){
					if(s[k][m] != sandpilelegend[tilecoordsassignment(k-i,m-j,config,tile,v,tilesize,tilecodeassignment,false)]){
						thispointisrgood = false;
						break;
					}
				}
			}
			if(thispointisrgood == false){
				break;
			}
		}
	} else{
		//IN THIS CASE, WE SEE R-BADNESS IMMEDIATELY!
		thispointisrgood = false;
	}	

	return thispointisrgood;
	
	
}

int configkiller(int i, int j, int tilesize, vector<vector<int> > s, vector<int> sandpilelegend, vector<vector<int> > tile, int tilewidth, vector<vector<int> > v, vector<vector<int> > tilecodeassignment){
	int winningconfig;
	vector<bool> configs(tilesize, true);
	int k;
	bool atleastoneremaining, morethanoneremaining;
	
	

	//KILL ON THE BASIS OF SANDPILE HEIGHT
	for(k=0; k < tilesize; k++){
		if(sandpilelegend[k] != s[i][j]){
			configs[k] = false;
		}
	}

	//CHECK
	atleastoneremaining = false;
	morethanoneremaining = false;
	for(k = 0; k < tilesize; k++){
		if(configs[k] == true){
			if(atleastoneremaining == true){
				morethanoneremaining = true;
				break;
			} else{
				atleastoneremaining = true;
				winningconfig = k;
			}
		}
	}
	if(morethanoneremaining == true){
		//NEXT ELIMINIATION PHASE: KILL ON THE BASIS OF CARDINAL DIRECTIONS
		for(k = 0; k<tilesize; k++){
			if(configs[k] == true){
				if(s[i][j+1] != sandpilelegend[tilecoordsassignment(0,1,k,tile,v,tilesize,tilecodeassignment,false)] || s[i+1][j] != sandpilelegend[tilecoordsassignment(1,0,k,tile,v,tilesize,tilecodeassignment,false)] ||s[i][j-1] != sandpilelegend[tilecoordsassignment(0,-1,k,tile,v,tilesize,tilecodeassignment,false)] ||s[i-1][j] != sandpilelegend[tilecoordsassignment(-1,0,k,tile,v,tilesize,tilecodeassignment,false)] ){
					configs[k] = false;
				}
			}
		}


		//CHECK
		atleastoneremaining = false;
		morethanoneremaining = false;
		for(k = 0; k < tilesize; k++){
			if(configs[k] == true){
				if(atleastoneremaining == true){
					morethanoneremaining = true;
					cout << "ERROR: MORE THAN ONE CONFIG REMAINS: INPUT MORE CHECKS!" << endl;
					cout << i << "  " << j << endl;

					break;
				} else{
					atleastoneremaining = true;
					winningconfig = k;
				}
			}
		}
	}else if(atleastoneremaining == false){
		//ALL CONFIGURATIONS KILLED -> NOT AN R-GOOD POINT
		winningconfig = tilesize;
	}
	
	return winningconfig;	
}
	


void sandpilecheck(int i, int j, vector<vector<int> > s, int tilesize){
	cout << "sandpile for point " << i << " " << j << ", and its NESW neighbors, is: " << s[i][j] << s[i][j+1] << s[i+1][j] << s[i][j-1] << s[i-1][j] << endl;
	cout << "remaining configs are: ";
	int k;
	for (k = 0; k < tilesize; k++){
	//	cout << config[k];
	}
	cout << endl;
}

int tilecoordsassignment(int xdisp, int ydisp, int centerpointconfig, vector<vector<int> > tile, vector<vector<int> > v, int tilesize, vector<vector<int> > tilecodeassignment, bool debug){
	int tilespot;
	//CONSIDER A TILE WITH 'BASEPOINT' AT ORIGIN. BELOW ARE THE COORDS OF THE CENTERPOINT IN THIS FUNDAMENTAL TILE
	int centerpointfundamentalx, centerpointfundamentaly;
	centerpointfundamentalx = tile[centerpointconfig][0];
	centerpointfundamentaly = tile[centerpointconfig][1];
	xdisp = xdisp + centerpointfundamentalx;
	ydisp = ydisp + centerpointfundamentaly;
	
	

	bool belongsintile = false;
	if(xdisp < 0){
		ydisp = ydisp + (((-1*xdisp)/v[0][0])+1)*v[0][1];
		xdisp = xdisp + (((-1*xdisp)/v[0][0])+1)*v[0][0];
	}
if(debug == true){cout << "xdisp and ydisp " << xdisp << " " << ydisp << endl;}
	ydisp = ydisp - (v[0][1]*(xdisp/(v[0][0])));
	xdisp = xdisp%(v[0][0]);
	

	if(ydisp < 0){
		ydisp = ydisp + (((-1*ydisp)/v[1][1])+1)*v[1][1];
	}
	ydisp = ydisp % v[1][1];

	belongsintile = tileconfirmation(xdisp, ydisp, tile, tilesize);	

	if (belongsintile == true){
		tilespot = tilecodeassignment[xdisp][ydisp];
	}else{				
		xdisp = xdisp + v[0][0];
		ydisp = ydisp + v[0][1];
		belongsintile = tileconfirmation(xdisp, ydisp, tile, tilesize);
		if(belongsintile == true){
			tilespot = tilecodeassignment[xdisp][ydisp];
		} else{
			xdisp = xdisp - v[0][0];
			ydisp = ydisp - v[0][1] - v[1][1];
			belongsintile = tileconfirmation(xdisp, ydisp, tile, tilesize);
			if (belongsintile == true){
				tilespot = tilecodeassignment[xdisp][ydisp];
			} else{
				cout << "COULD NOT FIND A TILE ASSIGNMENT" << endl;	
			}
		}
	}
	//if(debug == true){cout << "xdisp and ydisp: " << xdisp << " " << ydisp << endl;}
	return tilespot;
	
}


void timeteller(int t, int n){
	int hours, minutes, seconds;
	seconds = t/1000;
	minutes = seconds/60;
	hours = minutes/60;
	minutes = minutes % 60;
	seconds = seconds % 60;

	cout << "Time elapsed: " << hours;
	if (minutes < 10){
		cout << ":0" << minutes;
	}else{
		cout << ":" << minutes;
	}

	if (seconds < 10){
		cout << ":0" << seconds;
	}else{
		cout << ":" << seconds;
	}

	cout << endl << endl;
}
	

void tilecodeassignmentgenerator(int tilesize, int tilewidth, vector<vector<int> > tile, vector<vector<int> > &tilecodeassignment){
	int i,j,k;
	
	for (i = 0; i < tilewidth; i++){
		for (j=0; j < tilewidth; j++){
			for (k=0; k < tilesize; k++){
				if (tile[k][0] == i && tile[k][1] == j){
					tilecodeassignment[i][j] = k;
					break;
				} else if (k == tilesize - 1){
					tilecodeassignment[i][j] = tilesize;
				}
			}
		}
	}
}


void LemmaTest(double vnorm, int tilesize, vector<vector<int> > s, vector<int> sandpilelegend, vector<vector<int> > tile, int x_m, int y_m, int tilewidth, vector<vector<int> > v, vector<vector<int> > tilecodeassignment, vector<vector<bool> > Domain, double bigR, vector<vector<int> > periodicityvectors, vector<vector<int> > &doesthisyvaluegivetouchfrombelow, int n, double bigeval, double a, double b, double c){
	//TOUCH FROM BELOW INDICATES PHI FROM HYPOTHESIS OF LEMMA 14 EXISTS, PERFECT MATCH TESTS WHETHER THE CONCLUSION OF LEMMA 14 HOLDS

	bool touchfrombelow, perfectmatch, isorigingood;
	int i,j,k,l;
	//double bigR = 3*vnorm*vnorm*vnorm;
	

	vector<vector<double> > phiy(2*x_m+1);
	for (i = 0 ; i < 2*x_m+1 ; i++ ) {
   		phiy[i].resize(2*y_m+1);
	}

	
	vector<vector<double> > o(2*x_m+1);
	for (i = 0 ; i < 2*x_m+1 ; i++ ) {
   		o[i].resize(2*y_m+1);
	}


	isorigingood = imagematching(x_m, y_m, tilesize, s, sandpilelegend, tile, x_m, y_m, tilewidth, 2*bigR, periodicityvectors, tilecodeassignment);


	vector<vector<vector<int> > > touchingmap (2*x_m+1, vector<vector<int> >(2*y_m+1,vector <int>(2)));
	for(k = 0; k < 2*x_m + 1; k++){
		for(l = 0; l < 2*y_m + 1; l++){
			touchingmap[k][l][0] = 2*x_m+1;
			touchingmap[k][l][1] = 2*y_m + 1;
		}
	}

	vector<vector<int> > rangeoftouchingmap(2*x_m+1);
	for (i = 0 ; i < 2*x_m+1 ; i++ ) {
   		rangeoftouchingmap[i].resize(2*y_m+1);
	}

	for(i = 0; i < 2*x_m + 1; i++) {
		for (j = 0; j < 2*y_m + 1; j++){
			rangeoftouchingmap[i][j] = 0;
		}
	}

	odometergenerator(o, 2*x_m+1, 2*y_m+1,a,b,c,true);
	

	for(i = 0; i < 2*x_m + 1; i++) {
		for (j = 0; j < 2*y_m + 1; j++){
			//if(Domain[i][j] == 1){
			//if(i == x_m  && j == y_m ){
				for(k = 0; k < 2*x_m + 1; k++) {
					for (l = 0; l < 2*y_m + 1; l++){
						phiy[k][l] = o[k][l] - 0.5*vnorm*vnorm*((i-k)*(i-k)+(j-l)*(j-l))/(bigR*bigR) + (v[x_m][y_m] - o[x_m][y_m] +  0.5*vnorm*vnorm*((i-x_m)*(i-x_m)+(j-y_m)*(j-y_m))/(bigR*bigR));
						//phiy[k][l] = o[k][l] + v[x_m][y_m] - o[x_m][y_m];
						//phiy[k][l] = o[k][l];
						//phiy[k][l] = -2.0*(i - x_m) + o[k][l] - 0.5*vnorm*vnorm*((i-k)*(i-k)+(j-l)*(j-l))/(bigR*bigR) + (v[x_m][y_m] - o[x_m][y_m] +  0.5*vnorm*vnorm*((i-x_m)*(i-x_m)+(j-y_m)*(j-y_m))/(bigR*bigR));
					}
				}
				//cout << "phi at origin is " << phiy[x_m][y_m] << endl;
				//cout << "v at origin is " << v[x_m][y_m] << endl;

				touchingmapgenerator(touchingmap, i, j, x_m, y_m, phiy, v, Domain, n, bigeval, rangeoftouchingmap);
				
			/*	if(TouchFromBelowTest(x_m, y_m, bigR, v, phiy, o) == true){
					cout << "TOUCH FROM BELOW" << "      (" << i << ", " << j << ")" << endl;
					doesthisyvaluegivetouchfrombelow[i][j] = 2;
					//imagematching(x_m, y_m, tilesize, s, sandpilelegend, tile, x_m, y_m, tilewidth, vnorm*vnorm*vnorm, v, tilecodeassignment);
					if (isorigingood == false){
						cout << "CONTRADICTION: origin is not good but we have found a function which touches from below" << endl;
					}*/
				/*} else{cout << "DOES NOT TOUCH FROM BELOW" << endl;}*/
			//}
			//}
		}
	}

	//WriteMathematicaFile1Point5 (o,x_m, y_m, "odometer", n);
	WriteMathematicaFile1Point5(phiy,x_m,y_m,"phiy",n);


	for(i = x_m - 1; i < x_m + 2; i++) {
		for (j = y_m - 1; j < y_m + 2; j++){
			if((i-x_m)*(i-x_m) + (j-y_m)*(j-y_m) < 1.22*n*n){
				cout << i << ", " << j << "   maps to    " << touchingmap[i][j][0] << ", " << touchingmap[i][j][1] << endl;
			}
		}
	}



	touchingmapyfinder(touchingmap, x_m, y_m, x_m, y_m);

	//WriteMathematicaFile1(rangeoftouchingmap, x_m, y_m, "rangeoftouchingmap", n);	

	vector<vector<double> > diff(2*x_m + 1);
	for(i = 0; i < 2*x_m + 1; i++){
		diff[i].resize(2*y_m+1);
	}

	for(i = 0; i < 2*x_m + 1; i++){
		for(j = 0; j < 2*y_m + 1; j++){
			diff[i][j] = v[i][j] - phiy[i][j];
		}
	}

	//cout << "diff is " << diff[x_m][y_m + 1] << endl;
	

	WriteMathematicaFile1Point5 (diff, x_m, y_m, "diff", n);
	
	if(isorigingood == true){
		cout << "origin is good" << endl;
	}else {
		cout << "origin is not good" << endl;
	}
}

bool TouchFromBelowTest(int x_m, int y_m, double bigR, vector<vector<int> > v, vector<vector<double> > phiy, vector<vector<double> > o){
	int i,j;
	bool touchfrombelow = true;
	for(i = 0; i < 2*x_m + 1; i++) {
		for (j = 0; j < 2*y_m + 1; j++){
			//MAKE SURE bigR IS SUCH THAT B_R (0) IS CONTAINED IN THE ELLIPSE
			if (bigR*bigR > (i-x_m)*(i-x_m) + (j-y_m)*(j-y_m)){ 
				//cout << "probing the following i, j: " << i << ", " << j << "     ";
				//cout << "phiy: " << phiy[i][j] << "   v: " << v[i][j] << "   o: " << o[i][j] << endl;
				if (phiy[i][j] - v[i][j] > 0.00001){
					touchfrombelow = false;
					//cout << "problematic i,j: " << i << " , " << j << endl; 
					//cout << "phi - v at problem " << phiy[i][j] - v[i][j] << endl;
					break;
				}
			}
		}
		if (touchfrombelow == false){
			break;
		}
	}
	return touchfrombelow;		
}


void odometergenerator(vector<vector<double> > &odometer, int L, int M, double a, double b, double c, bool altlinearterms){
	
//cout << endl << endl << endl;
//cout << "testing uncorrected odometer" << endl;


	//HALF IS x_m, DEMI IS y_m
	//L IS 2*x_m + 1 AND M IS 2*y_m + 1
	int half = (L - 1)/2;
	int demi = (M - 1)/2;

	int i, j;

	vector<vector<int> > Uncorrected_O(L);

	for (i = 0 ; i < L ; i++) {
   		Uncorrected_O[i].resize(M);
	}


	


	Uncorrected_O[half][demi] = 0;
	Uncorrected_O[half][demi+1] = 0;	
	Uncorrected_O[half][demi+2] = 0;	
	Uncorrected_O[half+1][demi] = 0;
	Uncorrected_O[half+1][demi+1] = 1;	
	Uncorrected_O[half+1][demi+2] = 1;
	Uncorrected_O[half+2][demi] = 0;	
	Uncorrected_O[half+2][demi+1] = 1;	
	Uncorrected_O[half+2][demi+2] = 2;

	for (i = half; i < L-2; i++) {
		Uncorrected_O[i+2][y(i, half)+1+demi] = Uncorrected_O[i][y(i, half)+demi] + x(i, half) + y(i, half) + 1;
		Uncorrected_O[i+2][y(i, half)+2+demi] = Uncorrected_O[i][y(i,half)+1+demi] + x(i, half) + (y(i, half)+1) + 1;
	}
	
//cout << Uncorrected_O[half + 1][demi + 1] << endl;

	for (i = half+1; i > 1; i--){
		Uncorrected_O[i-2][y(i-2, half)+1+demi] = Uncorrected_O[i][y(i, half)+1+demi] - (x(i-2,half) + (y(i-2, half)+1)+1);
		Uncorrected_O[i-2][y(i-2, half)+demi] = Uncorrected_O[i][y(i, half)+demi] - (x(i-2,half) + y(i-2, half) + 1);
	}
		
//cout << Uncorrected_O[half + 1][demi + 1] << endl;
	

	for (i = 0; i < L; i++) {
		j = y(i, half) + demi;
		while (j + 2 < M) {
//if(j == demi - 1 && i == half + 1){cout << "WTF?" << endl;}
			Uncorrected_O[i][j+2] = Uncorrected_O[i][j] + x(i, half);
			j++;
		}
//if(i==half + 1){cout << Uncorrected_O[half + 1][demi + 1] << endl;}
		j = y(i, half) + 1 + demi;
		while (j - 2 >= 0){
			Uncorrected_O[i][j-2] = Uncorrected_O[i][j] - x(i, half);
			j--;
		}
	}
//cout << Uncorrected_O[half + 1][demi + 1] << endl;

//cout << "done testing uncorrected odometer" << endl << endl << endl;

	double f, g;

	
			//odometer[i][j] = (Uncorrected_O[i][j] - 0.5*((a-1)*f*f + 2*b*f*g + (c-1)*g*g)+f/4) + 0.5*((a)*f*f + 2*b*f*g + (c)*g*g) + .387*f + .5918*g;
	if(altlinearterms == false){
		for(i = 0; i < L; i++){
			f = i - half;
			for(j=0; j< M; j++){
				g = j - demi;
				odometer[i][j] = (Uncorrected_O[i][j] - 0.5*((a-1)*f*f + 2*b*f*g + (c-1)*g*g)+f/4)+ 0.5*((a)*f*f + 2*b*f*g + (c)*g*g) + 0.25*f + 0.5*g;
			}
		}
	} else{
		for(i = 1; i < L; i++){
			f = i - half;
			for(j=0; j< M; j++){
				g = j - demi;
				odometer[i][j] = (Uncorrected_O[i-1][j] - 0.5*((a-1)*(f-1)*(f-1) + 2*b*(f-1)*g + (c-1)*g*g)+(f-1)/4)+ 0.5*((a)*(f-1)*(f-1) + 2*b*(f-1)*g + (c)*g*g) + 1.25*f + 0.5*g;
			}
		}
		for(j=0; j< M; j++){
			odometer[0][j] = odometer[1][0];
		}
	}
			
	

}

int x(int i, int half) {
	int x = i - half;
	return x;
}

int y(int i, int half){
	int temp, y;
	temp = x(i, half);
	if (temp >= 0) {
		y = (temp/2);
	} else {
		y = (temp - 1)/2;
	}
	return y;
}


void Laplacian(vector<vector<double> > &L, vector<vector<int> > v, int x_m, int y_m) {
	
	int i,j;
	
	for (i = 1; i < 2*x_m; i++) {
		for (j = 1; j < 2*y_m; j++) {
				L[i][j] = -4*v[i][j] + v[i+1][j] + v[i-1][j] + v[i][j+1] + v[i][j-1];
		}
	}
}


void touchingmapgenerator(vector<vector<vector<int> > > &touchingmap, int i, int j, int x_m, int y_m, vector<vector<double> > phiy, vector<vector<int> > v, vector<vector<bool> > Domain, int n, double bigeval, vector<vector<int> > &rangeoftouchingmap){

	double runningmin;
	int runningminlocx, runningminlocy;
	bool firstpointchecked = true;
	int k,l;

	for(k = 0; k < 2*x_m+1; k++){
		for(l = 0; l < 2*y_m+1; l++){
			if ((k-x_m)*(k-x_m) + (l-y_m)*(l-y_m) < 2*n*n/bigeval){
				if(firstpointchecked == true){
					runningmin = v[k][l] - phiy[k][l];
					runningminlocx = k;
					runningminlocy = l;
					firstpointchecked = false;
				}else if(v[k][l] - phiy[k][l] < runningmin){
					runningmin = v[k][l] - phiy[k][l];
					runningminlocx = k;
					runningminlocy = l;
				}
			}
		}
	}

	touchingmap[i][j][0] = runningminlocx;
	touchingmap[i][j][1] = runningminlocy;

	rangeoftouchingmap[runningminlocx][runningminlocy] = 1;

if(i == x_m && j ==y_m){cout << "min of diff corresponding to origin: " << runningmin << "and it occurs at" << runningminlocx << ", " << runningminlocy <<  "and v is " << v[runningminlocx][runningminlocy] << "and o is " << phiy[runningminlocx][runningminlocy] << endl;}
}


void touchingmapyfinder(vector<vector<vector<int> > > touchingmap, int x_m, int y_m, int xyx, int xyy){
	
	int i,j;
	bool foundapoint = false;
	
	for(i = 0; i < 2*x_m+1; i++){
		for(j = 0; j < 2*y_m+1; j++){
			if(touchingmap[i][j][0] == xyx){
				if(touchingmap[i][j][1] == xyy){
					cout << "Found a point in the preimage of x_y: (" << i << ", " << j << ")" << endl;
					foundapoint = true;
				}
			}
		}
	}	

	cout << "a test of touching map: touching map of origin is: " << touchingmap[92][103][0] << ", " << touchingmap[92][103][1] << endl;

	if(foundapoint == false){
		cout << "could not find a point (within the rectangle) in the preimage of x_y" << endl;
	}

}

void randomsprinkle (int n, int k, vector<vector<int> > &s){
	int m = n/2;
	int i, number;

	for(i = 0; i < k; i++){
		number = rand ()%(m*m);
		s[littleboxlegend(number,m,true)][littleboxlegend(number,m,false)]++;
	}

}

int littleboxlegend (int number, int m, bool x){

	int coord;

	if(x == true){
		coord = number/m + m - m/2;
	} else {
		coord = number%m + m - m/2;
	}

	return coord;
}

