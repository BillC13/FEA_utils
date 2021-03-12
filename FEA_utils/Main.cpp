#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <vector>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <chrono>

// Global variables
double E;
double Pn;
double ro;
double h;
int NNod;
int NEle;
int NDOF;
int NDOFN;
double** NDat;
double** EDat;
double Dee[3][3];
double** KG;
double** MMat;
double** PVres;
double* BC;
double** Lods;
double** Acc0;
double s0[4], t0[4];


class Timer
{
private:
	// Type aliases to make accessing nested type easier
	using clock_t = std::chrono::high_resolution_clock;
	using second_t = std::chrono::duration<double, std::ratio<1> >;

	std::chrono::time_point<clock_t> m_beg;

public:
	Timer() : m_beg(clock_t::now())
	{
	}

	void reset()
	{
		m_beg = clock_t::now();
	}

	double elapsed() const
	{
		return std::chrono::duration_cast<second_t>(clock_t::now() - m_beg).count();
	}
};

void split(std::string strToSplit, char delimeter, double temp[]) // Inlet file processor
{
	std::stringstream ss(strToSplit);
	std::string item;
	int i = 0;
	while (std::getline(ss, item, delimeter))
	{
		temp[i] = (stod(item));
		i++;
	}
}
void getDee(double type, double Dee[3][3]) // D matrix
{
	if (type == 1)
	{
		Dee[0][0] = E / (1 - pow(Pn, 2));
		Dee[0][1] = Pn * E / (1 - pow(Pn, 2));
		Dee[0][2] = 0;
		Dee[1][0] = Pn * E / (1 - pow(Pn, 2));
		Dee[1][1] = E / (1 - pow(Pn, 2));
		Dee[1][2] = 0;
		Dee[2][0] = 0;
		Dee[2][1] = 0;
		Dee[2][2] = ((1 - Pn) / 2 * E) / (1 - pow(Pn, 2));
	}
	else
	{
		Dee[0][0] = E * (1 - Pn) / ((1 + Pn) * (1 - 2 * Pn));
		Dee[0][1] = Pn * E / ((1 + Pn) * (1 - 2 * Pn));
		Dee[0][2] = 0;
		Dee[1][0] = Pn * E / ((1 + Pn) * (1 - 2 * Pn));
		Dee[1][1] = E * (1 - Pn) / ((1 + Pn) * (1 - 2 * Pn));
		Dee[1][2] = 0;
		Dee[2][0] = 0;
		Dee[2][1] = 0;
		Dee[2][2] = E * ((1 - 2 * Pn) / 2) / ((1 + Pn) * (1 - 2 * Pn));
	}
}

int readFile(std::string filename) // Reads the inlet file
{
	std::ifstream myfile(filename);
	if (!myfile.is_open()) return 1;
	std::string dummy, spam;
	getline(myfile, dummy);
	getline(myfile, dummy);
	double type = stod(dummy);
	getline(myfile, dummy);
	E = stod(dummy);
	getline(myfile, dummy);
	Pn = stod(dummy);
	getline(myfile, dummy);
	ro = stod(dummy);
	getline(myfile, dummy);
	h = stod(dummy);
	getline(myfile, dummy);
	double AccX = stod(dummy);
	getline(myfile, dummy);
	double AccY = stod(dummy);
	getline(myfile, dummy);
	NNod = stoi(dummy);
	getline(myfile, dummy);
	NEle = stoi(dummy);
	getline(myfile, dummy);
	int NBC = stoi(dummy);
	getline(myfile, dummy);
	int NLod = stoi(dummy);
	NDOF = 2 * NNod;
	NDOFN = 2;
	int BCn = 0, Lodn = 0;
	getline(myfile, dummy);
	std::string data;
	std::stringstream datas(data);

	NDat = (double**)malloc(NNod * sizeof(double*));
	for (int i = 0; i < NNod; i++)
	{
		NDat[i] = (double*)malloc(3 * sizeof(double));
	}

	for (int i = 0; i < NNod; i++) // Node data
	{
		double temp[3];
		getline(myfile, dummy);
		split(dummy, ',', temp);
		for (int j = 0; j < 3; j++)
		{
			NDat[i][j] = temp[j];
		}
	}

	for (int i = 0; i < NNod; i++) // Node data
	{
		for (int j = 0; j < 3; j++)
		{
			std::cout << NDat[i][j] << " ";
		}
		std::cout << std::endl;
	}

	getline(myfile, dummy);
	EDat = (double**)malloc(NEle * sizeof(double*));
	for (int i = 0; i < NEle; i++)
	{
		EDat[i] = (double*)malloc(5 * sizeof(double));
	}

	for (int i = 0; i < NEle; i++) // Element data
	{
		double temp[5];
		getline(myfile, dummy);
		split(dummy, ',', temp);
		for (int j = 0; j < 5; j++)
		{
			EDat[i][j] = temp[j];
		}
	}

	for (int i = 0; i < NEle; i++) // Node data
	{
		for (int j = 0; j < 5; j++)
		{
			std::cout << EDat[i][j] << " ";
		}
		std::cout << std::endl;
	}

	BC = (double*)malloc(NDOF * sizeof(double));

	for (int i = 0; i < NDOF; i++)
	{
		BC[i] = 0;
	}

	(getline(myfile, dummy));
	for (int i = 0; i < NBC; i++) // Boundary conditions (1 = fixed)
	{
		double temp[3]; // LOOKY HERE
		getline(myfile, dummy);
		split(dummy, ',', temp);
		BC[int(2 * temp[0] - 2)] = temp[1];
		BC[int(2 * temp[0] - 1)] = temp[2];
	}

	for (int i = 0; i < NDOF; i++) // Node data
	{
		std::cout << BC[i] << std::endl;
	}

	Lods = (double**)malloc(NDOF * sizeof(double*));
	for (int i = 0; i < NDOF; i++)
	{
		Lods[i] = (double*)malloc(1 * sizeof(double));
	}
	for (int i = 0; i < NDOF; i++)
	{
		Lods[i][0] = 0;
	}

	(getline(myfile, dummy));
	for (int i = 0; i < NLod; i++) // Loads
	{
		double temp[3]; // LOOKY HERE
		getline(myfile, dummy);
		split(dummy, ',', temp);
		Lods[int(2 * temp[0] - 2)][0] = temp[1];
		Lods[int(2 * temp[0] - 1)][0] = temp[2];
	}

	for (int i = 0; i < NDOF; i++) // Node data
	{
		std::cout << Lods[i][0] << std::endl;
	}


	Acc0 = (double**)malloc(NDOF * sizeof(double*));
	for (int i = 0; i < NDOF; i++)
	{
		Acc0[i] = (double*)malloc(1 * sizeof(double));
	}

	for (int i = 0; i < NDOF; i++) // Acceleration
	{
		for (int j = 0; j < 1; j++)
		{
			if (i % 2 == 0)
			{
				Acc0[i][0] = AccX;
			}
			else
			{
				Acc0[i][0] = AccY;
			}
		}
	}

	for (int i = 0; i < NDOF; i++) // Node data
	{
		std::cout << Acc0[i][0] << std::endl;
	}
	double s1 = 1 / sqrt(3);
	double s0[2][2];
	double t0[2][2];
	s0[0][0] = -s1; // Shape functions
	s0[0][1] = s1;
	s0[1][0] = s1;
	s0[1][1] = -s1;
	t0[0][0] = -s1;
	t0[0][1] = -s1;
	t0[1][0] = s1;
	t0[1][1] = s1;

	double Dee[3][3];
	getDee(type, Dee);

	for (int i = 0; i < 3; i++) // Node data
	{
		for (int j = 0; j < 3; j++)
		{
			std::cout << Dee[i][j] << " ";
		}
		std::cout << std::endl;
	}

	KG = (double**)malloc(NDOF * sizeof(double*));
	for (int i = 0; i < NDOF; i++)
	{
		KG[i] = (double*)malloc(NDOF * sizeof(double));
	}
	MMat = (double**)malloc(NDOF * sizeof(double*));
	for (int i = 0; i < NDOF; i++)
	{
		MMat[i] = (double*)malloc(NDOF * sizeof(double));
	}
	for (int i = 0; i < NDOF; i++) // Node data
	{
		for (int j = 0; j < NDOF; j++)
		{
			KG[i][j] = 0;
			MMat[i][j] = 0;
		}
	}
	for (int i = 0; i < NDOF; i++) // Node data
	{
		for (int j = 0; j < NDOF; j++)
		{
			std::cout << KG[i][j] << " ";
		}
		std::cout << std::endl;
	}
	for (int i = 0; i < NDOF; i++) // Node data
	{
		for (int j = 0; j < NDOF; j++)
		{
			std::cout << MMat[i][j] << " ";
		}
		std::cout << std::endl;
	}
	double Wee[3][3];
	//matadd(Dee[3][3], Dee[3][3], Wee[3][3]);
	return 0;
}

int main(int argc, char* argv[])
{
	Timer t;
	//std::string inFile = argv[1];
	std::string inFile = "input4.txt";
	readFile(inFile);
	/*stiffness();
	KG2 = (double**)malloc(NDOF * sizeof(double*));
	for (int i = 0; i < NDOF; i++)
	{
		KG2[i] = (double*)malloc(NDOF * sizeof(double));
	}
	LodsTot = (double*)malloc(NDOF * sizeof(double));

	BFS(LodsTot);
	BCS(KG2);

	results(KG2, LodsTot);*/
	std::cout << "Time elapsed: " << t.elapsed() << " seconds\n";
	return 0;
}
