// This program computes a lower bound for \zeta_q as defined in the paper "On tournament inversion"

/* Prerequisites:
*
* 1) The api of the package lp_solve for linear programming see: https://lpsolve.sourceforge.net/5.5/.
*
* 2) A file containing the databse D_q of all tournaments on q vertices. This is taken from https://users.cecs.anu.edu.au/~bdm/data/digraphs.html for q <= 9.
*
* 3) A text file containing all q! permutations of [q] (e.g., if q=9, the first line of the file is 123456789  and the last line is 987654321.
*    The program loads this file into a two-dimensional array of dimensions q! by q. This is done to improve performance, as we can choose a random permutation by selecting a random row index into this
*	 array. Alternatively, one can avoid using this file and write a function which generates a random permutation of order q (this is less efficient as this function
*    must be called many times throughout the program run).
* 
*/

#include <ctime>           
#include "lp_lib.h"    // inluce file with the function definitons of the api of lp_solve.

const int q = 9;      // order of tournaments to check.
int numTournaments[] = { 1,1,1,2,4,12,56,456,6880,191536 };       // number of tournaments on q vertices for q = 0,...,9.  
int numPermutations[] = { 1,1,2,6,24,120,720,5040,40320,362880 }; // number of permutations of q for q =0,...,9.

int tournament[q][q];			// adjacency matrix of current tournament Q \in D_q (see the pseudocode in the paper for notation).
int graph[q][q];				// adjacency matrix of currnet permutation graph Q_L(\sigma)  (see the pseudocode in the paper for notation).
int** permutations;				// two-dimensional array containing all permutations of {0,..,q-1}.
int permutationSet[q * (q - 1)][q]; // the current orthogonal family of permutations to scan. This is called P^* in the paper (see the pseudocode in the paper for notation).
int rowNum[q][q];				// labeling the rows of the LP. The labels start from 1 according to lp_solve specification. The LP contains one inequality for each edge of Q_L(\sigma).
int colNum[q][q][q];			// labeling the variables of the LP. The labels start from 1 according to lp_solve specification. The LP contains one variable for each triangle of Q_L(\sigma).


FILE* datafile;

#define DATAFILE "d:\\research\\general\\combinatorial data\\tour"			// prefix of file name containing the tournament database file D_q.
#define PERMFILE "d:\\research\\general\\combinatorial data\\permutations"	// prefix of file name containing all permutations of [q].

// various definitions of function pointers to lp_solve api functions.
HINSTANCE lpsolve;
lprec* lp;
make_lp_func* _make_lp;
delete_lp_func* _delete_lp;
set_row_func* _set_row;
set_obj_fn_func* _set_obj_fn;
solve_func* _solve;
set_maxim_func* _set_maxim;
get_objective_func* _get_objective;
print_lp_func* _print_lp;
set_rh_vec_func* _set_rh_vec;
set_verbose_func* _set_verbose;

/* The seed orthogonal family of order 9 given in Table 1 of the paper. This is denoted P in the paper (see the pseudocode in the paper for notation). Comment out if using q != 9 */
int P[72][9] = {
{0,1,2,3,4,5,6,7,8},
{1,2,0,4,5,3,7,8,6},
{2,0,1,5,3,4,8,6,7},
{3,4,5,6,7,8,0,1,2},
{4,5,3,7,8,6,1,2,0},
{5,3,4,8,6,7,2,0,1},
{6,7,8,0,1,2,3,4,5},
{7,8,6,1,2,0,4,5,3},
{8,6,7,2,0,1,5,3,4},
{0,2,1,6,8,7,3,5,4},
{1,0,2,7,6,8,4,3,5},
{2,1,0,8,7,6,5,4,3},
{3,5,4,0,2,1,6,8,7},
{4,3,5,1,0,2,7,6,8},
{5,4,3,2,1,0,8,7,6},
{6,8,7,3,5,4,0,2,1},
{7,6,8,4,3,5,1,0,2},
{8,7,6,5,4,3,2,1,0},
{0,3,6,4,7,1,8,2,5},
{1,4,7,5,8,2,6,0,3},
{2,5,8,3,6,0,7,1,4},
{3,6,0,7,1,4,2,5,8},
{4,7,1,8,2,5,0,3,6},
{5,8,2,6,0,3,1,4,7},
{6,0,3,1,4,7,5,8,2},
{7,1,4,2,5,8,3,6,0},
{8,2,5,0,3,6,4,7,1},
{0,4,8,7,2,3,5,6,1},
{1,5,6,8,0,4,3,7,2},
{2,3,7,6,1,5,4,8,0},
{3,7,2,1,5,6,8,0,4},
{4,8,0,2,3,7,6,1,5},
{5,6,1,0,4,8,7,2,3},
{6,1,5,4,8,0,2,3,7},
{7,2,3,5,6,1,0,4,8},
{8,0,4,3,7,2,1,5,6},
{0,5,7,1,3,8,2,4,6},
{1,3,8,2,4,6,0,5,7},
{2,4,6,0,5,7,1,3,8},
{3,8,1,4,6,2,5,7,0},
{4,6,2,5,7,0,3,8,1},
{5,7,0,3,8,1,4,6,2},
{6,2,4,7,0,5,8,1,3},
{7,0,5,8,1,3,6,2,4},
{8,1,3,6,2,4,7,0,5},
{0,6,3,8,5,2,4,1,7},
{1,7,4,6,3,0,5,2,8},
{2,8,5,7,4,1,3,0,6},
{3,0,6,2,8,5,7,4,1},
{4,1,7,0,6,3,8,5,2},
{5,2,8,1,7,4,6,3,0},
{6,3,0,5,2,8,1,7,4},
{7,4,1,3,0,6,2,8,5},
{8,5,2,4,1,7,0,6,3},
{0,7,5,2,6,4,1,8,3},
{1,8,3,0,7,5,2,6,4},
{2,6,4,1,8,3,0,7,5},
{3,1,8,5,0,7,4,2,6},
{4,2,6,3,1,8,5,0,7},
{5,0,7,4,2,6,3,1,8},
{6,4,2,8,3,1,7,5,0},
{7,5,0,6,4,2,8,3,1},
{8,3,1,7,5,0,6,4,2},
{0,8,4,5,1,6,7,3,2},
{1,6,5,3,2,7,8,4,0},
{2,7,3,4,0,8,6,5,1},
{3,2,7,8,4,0,1,6,5},
{4,0,8,6,5,1,2,7,3},
{5,1,6,7,3,2,0,8,4},
{6,5,1,2,7,3,4,0,8},
{7,3,2,0,8,4,5,1,6},
{8,4,0,1,6,5,3,2,7} };

/* A seed orthogonal family of order 8. This is denoted P in the paper (see the pseudocode in the paper for notation). Comment out if using q != 8 
int P[56][8] = {
	{0,1,2,3,4,5,6,7},
	{1,0,3,2,5,4,7,6},
	{2,3,0,1,6,7,4,5},
	{3,2,1,0,7,6,5,4},
	{4,5,6,7,0,1,2,3},
	{5,4,7,6,1,0,3,2},
	{6,7,4,5,2,3,0,1},
	{7,6,5,4,3,2,1,0},
	{0,2,4,6,5,7,1,3},
	{1,3,5,7,4,6,0,2},
	{2,0,6,4,7,5,3,1},
	{3,1,7,5,6,4,2,0},
	{4,6,0,2,1,3,5,7},
	{5,7,1,3,0,2,4,6},
	{6,4,2,0,3,1,7,5},
	{7,5,3,1,2,0,6,4},
	{0,3,6,5,1,2,7,4},
	{1,2,7,4,0,3,6,5},
	{2,1,4,7,3,0,5,6},
	{3,0,5,6,2,1,4,7},
	{4,7,2,1,5,6,3,0},
	{5,6,3,0,4,7,2,1},
	{6,5,0,3,7,4,1,2},
	{7,4,1,2,6,5,0,3},
	{0,4,5,1,7,3,2,6},
	{1,5,4,0,6,2,3,7},
	{2,6,7,3,5,1,0,4},
	{3,7,6,2,4,0,1,5},
	{4,0,1,5,3,7,6,2},
	{5,1,0,4,2,6,7,3},
	{6,2,3,7,1,5,4,0},
	{7,3,2,6,0,4,5,1},
	{0,5,7,2,3,6,4,1},
	{1,4,6,3,2,7,5,0},
	{2,7,5,0,1,4,6,3},
	{3,6,4,1,0,5,7,2},
	{4,1,3,6,7,2,0,5},
	{5,0,2,7,6,3,1,4},
	{6,3,1,4,5,0,2,7},
	{7,2,0,5,4,1,3,6},
	{0,6,1,7,2,4,3,5},
	{1,7,0,6,3,5,2,4},
	{2,4,3,5,0,6,1,7},
	{3,5,2,4,1,7,0,6},
	{4,2,5,3,6,0,7,1},
	{5,3,4,2,7,1,6,0},
	{6,0,7,1,4,2,5,3},
	{7,1,6,0,5,3,4,2},
	{0,7,3,4,6,1,5,2},
	{1,6,2,5,7,0,4,3},
	{2,5,1,6,4,3,7,0},
	{3,4,0,7,5,2,6,1},
	{4,3,7,0,2,5,1,6},
	{5,2,6,1,3,4,0,7},
	{6,1,5,2,0,7,3,4},
	{7,0,4,3,1,6,2,5}
}; */


/* A seed orthogonal family of order 7. This is denoted P in the paper (see the pseudocode in the paper for notation). Comment out if using q != 7 
int P[42][7] = {
	{0,1,2,3,4,5,6},
	{1,2,3,4,5,6,0},
	{2,3,4,5,6,0,1},
	{3,4,5,6,0,1,2},
	{4,5,6,0,1,2,3},
	{5,6,0,1,2,3,4},
	{6,0,1,2,3,4,5},
	{0,2,4,6,1,3,5},
	{1,3,5,0,2,4,6},
	{2,4,6,1,3,5,0},
	{3,5,0,2,4,6,1},
	{4,6,1,3,5,0,2},
	{5,0,2,4,6,1,3},
	{6,1,3,5,0,2,4},
	{0,3,6,2,5,1,4},
	{1,4,0,3,6,2,5},
	{2,5,1,4,0,3,6},
	{3,6,2,5,1,4,0},
	{4,0,3,6,2,5,1},
	{5,1,4,0,3,6,2},
	{6,2,5,1,4,0,3},
	{0,4,1,5,2,6,3},
	{1,5,2,6,3,0,4},
	{2,6,3,0,4,1,5},
	{3,0,4,1,5,2,6},
	{4,1,5,2,6,3,0},
	{5,2,6,3,0,4,1},
	{6,3,0,4.1,5,2},
	{0,5,3,1,6,4,2},
	{1,6,4,2,0,5,3},
	{2,0,5,3,1,6,4},
	{3,1,6,4,2,0,5},
	{4,2,0,5,3,1,6},
	{5,3,1,6,4,2,0},
	{6,4,2,0,5,3,1},
	{0,6,5,4,3,2,1},
	{1,0,6,5,4,3,2},
	{2,1,0,6,5,4,3},
	{3,2,1,0,6,5,4},
	{4,3,2,1,0,6,5},
	{5,4,3,2,1,0,6},
	{6,5,4,3,2,1,0}
}; */

/* A seed orthogonal family of order 5. This is denoted P in the paper (see the pseudocode in the paper for notation). Comment out if using q != 5 
int P[20][5] = {
	{0,1,2,3,4},
	{1,2,3,4,0},
	{2,3,4,0,1},
	{3,4,0,1,2},
	{4,0,1,2,3},
	{0,2,4,1,3},
	{1,3,0,2,4},
	{2,4,1,3,0},
	{3,0,2,4,1},
	{4,1,3,0,2},
	{0,3,1,4,2},
	{1,4,2,0,3},
	{2,0,3,1,4},
	{3,1,4,2,0},
	{4,2,0,3,1},
	{0,4,3,2,1},
	{1,0,4,3,2},
	{2,1,0,4,3},
	{3,2,1,0,4},
	{4,3,2,1,0}
}; 
*/

/* A seed orthogonal family of order 4. This is denoted P in the paper (see the pseudocode in the paper for notation). Comment out if using q != 4 
int P[12][4] = {
	{0,3,1,2},
	{1,2,0,3},
	{2,1,3,0},
	{3,0,2,1},
	{0,2,3,1},
	{1,3,2,0},
	{2,0,1,3},
	{3,1,0,2},
	{0,1,2,3},
	{1,0,3,2},
	{2,3,0,1},
	{3,2,1,0}
};
*/

// Load the next tournament from the tournament database file into its adjacency matrix.
void nextTournamentFromFile()
{
	char c;
	for (int i = 0; i < q; i++)
	{
		tournament[i][i] = 0;
		for (int j = i + 1; j < q; j++)
		{
			c = fgetc(datafile);
			tournament[i][j] = c - '0';
			tournament[j][i] = 1 - tournament[i][j];
		}
	}
	fgetc(datafile);
}

// Load the permutations of [q] from text file into the two-dimensional array "permutations" (the permutations in the array are shifted to be in {0,...q-1}.
void loadPermutations() {
	permutations = new int* [numPermutations[q]];
	for (int i = 0; i < numPermutations[q]; i++)
		permutations[i] = new int[q];
	FILE* permfile;

	char permFileName[200];
	sprintf_s(permFileName, "%s%d.txt", PERMFILE, q);
	fopen_s(&permfile, permFileName, "r");
	for (int i = 0; i < numPermutations[q]; i++) {
		for (int j = 0; j < q; j++)
		{
			char c;
			fscanf_s(permfile, "%c", &c);
			permutations[i][j] = c - '1';
		}
		fgetc(permfile);
	}
	fclose(permfile);
}

// construct P^* = {\pi \sigma : \sigma \in P}.
void setPermutationSet() {
	int pi = rand() % numPermutations[q]; // choose a random permutation pi.

	for (int i = 0; i < q * (q - 1); i++)
		for (int k = 0; k < q; k++)
			permutationSet[i][k] = permutations[pi][P[i][k]];
}


// create adjacency matrix of Q_L(\sigma)
void createPermutationGraph(int permNum) {
	for (int i = 0; i < q; i++)
		for (int j = 0; j < q; j++)
			graph[i][j] = 0;
	for (int i = 0; i < q; i++)
		for (int j = 0; j < q; j++)
			if (tournament[i][j] == 1 && permutationSet[permNum][i] < permutationSet[permNum][j]) {
				graph[i][j] = 1;
				graph[j][i] = 1;
			}
}

// set the row numbers of the LP. One row for each possible edge.
void setRowNum()
{
	int count = 1;
	for (int i = 0; i < q - 1; i++)
		for (int j = i + 1; j < q; j++)
			rowNum[i][j] = count++;
}

// set the column numbers of the LP. One column for each possible triangle.
void setColNum()
{
	int count = 1;
	for (int i = 0; i < q - 2; i++)
		for (int j = i + 1; j < q - 1; j++)
			for (int k = j + 1; k < q; k++)
				colNum[i][j][k] = count++;
}

// write the target of the LP according to the lpSolve api specification.
void writeTargetLP()
{
	int i, j, k;
	static REAL row[q * (q - 1) * (q - 2) / 6];

	for (i = 0; i < q - 2; i++)
		for (j = i + 1; j < q - 1; j++)
			for (k = j + 1; k < q; k++)
				if (graph[i][j] == 1 && graph[i][k] == 1 && graph[j][k] == 1)
					row[colNum[i][j][k]] = 1;
				else
					row[colNum[i][j][k]] = 0;
	_set_obj_fn(lp, row);
	_set_maxim(lp);  // tell the LP that we want to maximize the total value on all traingles of Q_L(\sigma).
}

// write the constraints of the LP according to the lpSolve api specification.
void writeConstraintsLP()
{
	int i, j, k;
	static REAL row[q * (q - 1) * (q - 2) / 6];
	static REAL column[q * (q - 1) / 2];

	for (i = 0; i < q - 1; i++)
		for (j = i + 1; j < q; j++)
		{
			row[0] = 1;
			for (k = 1; k < 1 + q * (q - 1) * (q - 2) / 6; k++)
				row[k] = 0;
			for (k = 0; k < q; k++)
				if (k < i)
					row[colNum[k][i][j]] = 1;
				else if (k < j && k > i)
					row[colNum[i][k][j]] = 1;
				else if (k > j)
					row[colNum[i][j][k]] = 1;
			column[rowNum[i][j]] = 1;
			_set_row(lp, rowNum[i][j], row);
		}
	_set_rh_vec(lp, column);
}

int main()
{
	srand((unsigned)time(0));
	loadPermutations();

	lpsolve = LoadLibrary(TEXT("d:\\Projects\\Legacy\\LPSolveIDE\\lpsolve55.dll"));
	_make_lp = (make_lp_func*)GetProcAddress(lpsolve, "make_lp");
	_delete_lp = (delete_lp_func*)GetProcAddress(lpsolve, "delete_lp");
	_set_row = (set_row_func*)GetProcAddress(lpsolve, "set_row");
	_set_obj_fn = (set_obj_fn_func*)GetProcAddress(lpsolve, "set_obj_fn");
	_solve = (solve_func*)GetProcAddress(lpsolve, "solve");
	_set_maxim = (set_maxim_func*)GetProcAddress(lpsolve, "set_maxim");
	_get_objective = (get_objective_func*)GetProcAddress(lpsolve, "get_objective");
	_print_lp = (print_lp_func*)GetProcAddress(lpsolve, "print_lp");
	_set_rh_vec = (set_rh_vec_func*)GetProcAddress(lpsolve, "set_rh_vec");
	_set_verbose = (set_verbose_func*)GetProcAddress(lpsolve, "set_verbose");

	lp = _make_lp(q * (q - 1) / 2, q * (q - 1) * (q - 2) / 6);   // set the dimensions of the LP.
	_set_verbose(lp, SEVERE);
	setRowNum();
	setColNum();

	char dataFileName[200];
	sprintf_s(dataFileName, "%s%d.txt", DATAFILE, q);    //suffix the file name according to the value of q.
	printf("Checking %d instances on %d vertices.\n", numTournaments[q], q);
	fopen_s(&datafile, dataFileName, "r");

	double zeta = 1000000;  // set zeta to "infinity". See pseudocode in the paper.
	int permutationSetSize = q * (q - 1);   // the size of the orthogonal family.
	for (int counter = 0; counter < numTournaments[q]; counter++)  // loop on all tournaments in the tournamnt database.
	{
		if (counter % 1000 == 0 && counter > 0)
			printf("Up until now we have listed %d; the beest is %.10f\n", counter, zeta);  // printing information on program's progression.
		nextTournamentFromFile();  // read the next tournament Q from the tournament database file.
		double best = -1;          // initialize best to a negative value. See pseudocode in the paper.
		int loop = 0;		       // initialize loop counter. See pseudocode in the paper.
		while ( best < 3.8 &&  loop < 1000 && best < zeta) {
			loop++;
			setPermutationSet(); // set P^*
			double sum = 0;      // initialize sum. See pseudocode in the paper. 
			for (int j = 0; j < permutationSetSize; j++) { // for all \sigma \in P^*
				createPermutationGraph(j);  // consdtruct adjacency matrix of Q_L(\sigma)
				writeTargetLP();	        // arrange LP target
				writeConstraintsLP();       // arrange LP constraints
				_solve(lp);					// compute \nu^*(Q_L(\\sigma))
				sum += _get_objective(lp);; // add \nu^*(Q_L(\\sigma)) to currnet sum
			}
			double avg = sum /= permutationSetSize;   // compute avg_{P^*}(Q)
			if (avg > best)  // if the current loop trial improves the present best then improve best
				best = avg;
		}
		printf("%.10f\n", best);  // print the best that has been found for the current tournament
		if (best < zeta)          // if currnet tournament is the worst thus found, then update \zeta
			zeta = best;
	}
	printf("The lower bound for zeta_q computed by the program is %.10f\n", zeta); 
	printf("This implies a value of beta = %.10f\n", zeta/(q*(q-1)));
	printf("This implies a value of %.10f for the lower bound constant in Theorem 1.3.\n", 1.0/8.0 - zeta / (2*q * (q - 1)));
	_delete_lp(lp);
	FreeLibrary(lpsolve);
	return 0;
}