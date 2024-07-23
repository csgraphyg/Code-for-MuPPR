#include "prototypes.h"
#include <string>
#include <cstring>
#include <cstdio>
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <math.h>
#include <queue>
#include <stdlib.h>
//============================================================================//
using namespace std;
//============================================================================//
//               Global variables that can be accessed anywhere               //
//============================================================================//
uint N=0; 			// N is the highest node index.
uint realN=0;		// Real number of nodes
uint M=0; 			// M is the total number of edges.
uint n_countInt=0;				// The counter of countInt[]
uint maxDegree=0, index=0;
uint **adjList=NULL; 	// The adjacent list
uint *degree=NULL;		// The degree list
uint *countInt=NULL;	// countInt[] array to find intersection, etc
uint *m_Map=NULL;
uint abc=0;
uint *finall=NULL;
uint *why=NULL;uint *XX=NULL;double *R=NULL;double *X=NULL;
uint* mapp; double max_triangle;
bool lwt; double afal;
//============================================================================//
//                         Details of global prototypes                       //
//============================================================================//
void serr(const char str[]) {
	cout << str << endl;
	exit(0);
}
//============================================================================//
FILE *initOpenFile(const char filename[], const char mode[]) {
	FILE *f = NULL;
	if ( strcmp(mode, "read")==0 || strcmp(mode, "rt")==0 ) {// Open the file to read
		fopen_s(&f, filename,"rt");
	} else if ( strcmp(mode, "write")==0 || strcmp(mode, "wt")==0 ) { // Open the file to write
		fopen_s(&f, filename,"wt");
	} else {
		serr("Error: <mode> is incorrect in initOpenFile()");
	}
	if ( !f ) {
		serr("Error: Can not open input file in initOpenFile()");
	}
	return f;
}
//============================================================================//
void initDegree() {
	degree = new uint [index];
	if (!degree ) {
		serr("Error: Not enough memory for degree[] in initDegree()");
	}
	memset(degree, 0, UINT_SIZE*index); // Set up the degree arrary
}
//============================================================================//
void initm_Map() {
	m_Map = new uint [2*index];
	mapp = new uint [2*index];
	if (!m_Map||!mapp ) {
		serr("Error: Not enough memory for m_Map[] in initm_Map()");
	}
}
//============================================================================//
void initAdjList() {	
	adjList = new uint* [index]; // Initialize the adjacent list
	if ( !adjList ) {
		serr("Error: Not enough memory for adjList[] in readFromFile()");
	}
	maxDegree = 0;
	for(uint a=0; a<index; a++) 
	{ 
		adjList[a] = new uint [ degree[a] ];
		if ( !adjList[a] ) {
			serr("Error: !adjList[a] in initAdjList()");
		}
	}
}
//============================================================================//
void initcountInt() {
	countInt = new uint [index];
	memset(countInt, 0, UINT_SIZE*(index)); // Set up the countInt array
}
//----------------------------------------------------------------------------//
void initR(){
	R= new double[index];
	X= new double[index];
	XX= new uint[index];
	why= new uint[index];
}
//----------------------------------------------------------------------------//
void updateCounter() { // update the counter
	if (n_countInt >= SAFE_ULONGMAX) {
		memset(countInt, 0, index*UINT_SIZE); // clear the counter
		n_countInt = 0; // clear the counting variable
	}
	n_countInt++; // increase the counter
}
//============================================================================//
// Function to read the network from file
void readFromFile(const char fileName[]) {
	std::map<uint, uint>* mapa = new std::map<uint, uint>();
	if (!mapa) {
		serr("\t Graph: Error allocating mapa\n");
	}
	uint a, b;
	FILE* f = initOpenFile(fileName, "read"); // open file to read
	N = MAX_N;
	while (!feof(f)) { // read file to find the degree and weight for each node
		fscanf_s(f, "%u %u", &a, &b);
		if (a == b || a < 0 || b < 0) {
			continue;
		}
		if (a > b)
			continue;
		M++;
		if (a > realN) {
			realN = a;
		}
		if (b > realN) {
			realN = b;
		}
		if (!mapa->count(a)) {
			mapa->insert(std::pair<uint, uint>(a, index));
			index++;
		}
		if (!mapa->count(b)) {
			mapa->insert(std::pair<uint, uint>(b, index));
			index++;
		}
	}
	initm_Map();
	for (std::map<uint, uint>::iterator it = mapa->begin(); it != mapa->end(); it++) {
		m_Map[(*it).second] = (*it).first;
		mapp[(*it).first] = (*it).second;
	}
	fseek(f, 0, SEEK_SET);
	initDegree(); // initialize the degree and weight arrays
	while (!feof(f)) {
		fscanf_s(f, "%u %u", &a, &b);
		if (a == b || a < 0 || b < 0) {
			continue;
		}
		if (a > b)
			continue;
		degree[(*mapa->find(a)).second]++;
		degree[(*mapa->find(b)).second]++;
	}
	cout << "realN=" << realN << '\n';
	cout << "index=" << index << '\n'; cout << "M=" << M << '\n';

	initAdjList(); // initialize the adjList[]
	memset(degree, 0, UINT_SIZE * index); // memset 0 for degree
	fseek(f, 0, SEEK_SET); // go back and read the adjacent list
	while (!feof(f)) { // note that we assume the input is friendly, that is there are no duplicate edges.
		fscanf_s(f, "%u %u", &a, &b);
		if (a == b || a < 0 || b < 0) {
			continue;
		}
		if (a > b)
			continue;
		adjList[(*mapa->find(a)).second][degree[(*mapa->find(a)).second]] = (*mapa->find(b)).second;
		degree[(*mapa->find(a)).second]++;
		adjList[(*mapa->find(b)).second][degree[(*mapa->find(b)).second]] = (*mapa->find(a)).second;
		degree[(*mapa->find(b)).second]++;
	}
	fclose(f);
	delete mapa;
}
//============================================================================//
int muppr(double lest,uint seed){
	std::queue<uint> q; std::queue<uint> faker; 
	uint depth=0;
	updateCounter();
	if (lwt == 1)
	{
		faker.push(seed); countInt[seed] = n_countInt;
		R[seed] = 1;
		if (R[seed] >= degree[seed] * lest)
		{
			q.push(seed);
			why[seed] = 1;
		}
		if (q.empty())
		{
			R[seed] = 0;
			return 0;
		}
		while (!q.empty())
		{
			X[q.front()] += R[q.front()] * afal;
			double R_value = R[q.front()];
			R[q.front()] = (1 - afal) * R[q.front()] / 2.0;
			int front_num = q.front();
			if (R[q.front()] < degree[q.front()] * lest)
			{
				why[q.front()] = 0;
				q.pop();
			}
			for (uint k = 0; k < degree[front_num]; k++)
			{
				R[adjList[front_num][k]] += (1 - afal) * R_value  / (2.0 * degree[front_num]);
				if (countInt[adjList[front_num][k]] != n_countInt)
				{
					faker.push(adjList[front_num][k]);
					countInt[adjList[front_num][k]] = n_countInt;
				}
				if (R[adjList[front_num][k]] >= degree[adjList[front_num][k]] * lest && why[adjList[front_num][k]] == 0)
				{
					q.push(adjList[front_num][k]);
					why[adjList[front_num][k]] = 1;
				}
			}
			depth++;
		}
		cout << "iterations=" << depth << endl;
	}
//----------------------------------------------------------------------------//
	else
	{
		faker.push(seed); countInt[seed] = n_countInt;
		R[seed] = afal;
		if (R[seed] >= degree[seed] * lest* afal)
		{
			q.push(seed);
			why[seed] = 1;
		}
		if (q.empty())
		{
			R[seed] = 0;
			return 0;
		}
		while (!q.empty())
		{
			X[q.front()] += R[q.front()];
			if (countInt[q.front()] != n_countInt)
			{
				faker.push(q.front());//向量p(x)中的节点
				countInt[q.front()] = n_countInt;
			}//无法清空其它数据
			double R_value = R[q.front()];
			R[q.front()] = (1 - afal) * R[q.front()] / pow(degree[q.front()] + 1, 1.5);
			int front_num = q.front();
			if (R[q.front()] < degree[q.front()] * lest * afal)
			{
				why[q.front()] = 0;
				q.pop();
			}
			for (uint k = 0; k < degree[front_num]; k++)
			{
				R[adjList[front_num][k]] += (1 - afal) * R_value / ((degree[front_num] + 1) * sqrt(degree[adjList[front_num][k]] + 1));
				if (R[adjList[front_num][k]] >= degree[adjList[front_num][k]] * lest * afal && why[adjList[front_num][k]] == 0)
				{
					q.push(adjList[front_num][k]);
					why[adjList[front_num][k]] = 1;
				}
			}
			depth++;
		}
		cout << "iterations=" << depth << endl;
	}


	while(!q.empty())
	{
		cout<<"can q be empty?"<<endl;
		q.pop();			
	}
	uint position=0; uint postion1=0; uint p=0; double ave_triangle=0.0;max_triangle=0.0; double maxcond=1.0;
	std::multimap<double, uint, std::greater<double> > * map = new std::multimap<double, uint, std::greater<double> >();
	while(!faker.empty())
	{
		if(X[faker.front()] > 0.0&&degree[faker.front()]!=0)
		{
			map->insert(std::pair<double, uint>(1.0*X[faker.front()]/degree[faker.front()], faker.front()));
			R[faker.front()]=0.0;
			X[faker.front()]=0.0;
			why[faker.front()]=0;
			faker.pop();
		}
		else
		{
			R[faker.front()]=0.0;
			X[faker.front()]=0.0;
			why[faker.front()]=0;
			faker.pop();
		}
	}
	for(std::multimap<double, uint, std::greater<double> >::iterator it = map->begin();it != map->end();it++)
	{
		XX[position]=it->second;
		position++;
		//cout << "it->first="<<it->first << endl;
	}
	cout << "length=" << position << endl;

	updateCounter();
	double vol = 0, cut = 0, conductance=1;
	for (uint i = 0; i < position; i++)
	{
		uint el = XX[i];
		vol += degree[el];
		for (uint j = 0; j < degree[el]; j++)
		{
			uint ele = adjList[el][j];
			if (countInt[ele] != n_countInt)
				cut += 1;
			else
				cut -= 1;
		}
		conductance = cut / vol;
		//cout << "conductance=" << conductance << endl;
		countInt[el] = n_countInt;
		if (conductance <= maxcond)
		{
			p = i;
			maxcond = conductance;
		}
	}
	finall=new uint[p+1];
	abc = p+1;	
	for(std::multimap<double, uint, std::greater<double> >::iterator it = map->begin();it != map->end()&&postion1<=p; it++)
	{
		finall[postion1] = it->second;
      	postion1++;
	}
	delete map;
	return 1;
}
//============================================================================//
void tryFormingNewCommunity(){
	initR();
	initcountInt();
	for(uint k=0;k<index;k++)
	{
			R[k]=0.0;
			X[k]=0.0;
			why[k]=0;
	}
	double lest=0.0001; lwt = 1; afal = 0.1;
	int query = 6037;
		//154846;
	for(uint i=0;i<1;i++)
	{
		muppr(lest, mapp[query]);
	}
	if (R != NULL) 
	{
		delete [] R;
		R = NULL;
	}
	if (X != NULL) 
	{
		delete [] X;
		X = NULL;
	}
	if (XX != NULL) 
	{
		delete [] XX;
		XX = NULL;
	}
	if (why != NULL) 
	{
		delete [] why;
		why = NULL;
	}
}
//============================================================================//
void printData(const char FILENAME[]) {
 // Print everything to a file specified by the FILENAME
		double part_name=0.75;
		char filename[256], tmp[5];
		strcpy_s(filename,FILENAME);
		cout<<"FILENAME="<<FILENAME<<'\n';
		size_t len = strlen(FILENAME);
		cout<<"strlen="<<len<<'\n';
		filename[len-4] = '\0';
		strcat_s(filename, "_CID_AMU_");
		_itoa_s(int(part_name*100), tmp, 10);
		strcat_s(filename,tmp);
		strcat_s(filename,".txt");
		cout<<"Results = "<<filename<<endl;
		FILE *f = initOpenFile(filename, "write");
		for (uint j = 0; j < abc; j++)
		{
			fprintf(f, "%u ", m_Map[finall[j]]);
		}
		fprintf(f, "\n");
		fclose(f);
}
//============================================================================//
void MU(const char filename[]) {
	cout<<"Reading base file: "<<filename;
	readFromFile(filename);
	cout<<". Done."<<endl;
	clock_t tstart = clock(); // Start the timer
	cout<<"+++++++++++++++++++++++++++++++++++++++++++++"<<endl;
	cout<<"+++++++++++++++ Starting MU +++++++++++++++"<<endl;
	cout<<"+++++++++++++++++++++++++++++++++++++++++++++"<<endl;
	tryFormingNewCommunity();
	clock_t tend = clock(); // Stop the timer	
	cout<<"----- Printing MU results -----"<<endl;
	cout<<"Time taken = "<<(double)(tend - tstart)/CLOCKS_PER_SEC<<"s"<<endl;
	printData(filename);
	cout<<"DONE MU"<<endl;
}
//============================================================================//
//========================== END ADAPTIVE FUNCTION ===========================//
//============================================================================//
void doAMU(const char base_file[], const char file_path[], const uint numFile){
	MU(base_file);
	cout << "*********** ALL DONE ***********" << endl;
}
