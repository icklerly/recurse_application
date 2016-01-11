//============================================================================
// Name        : LEGO.cpp
// Author      : Lydia Ickler
// Version     : HMM-like
//============================================================================
#include <iostream>
#include <array>
using namespace std;

void maxR(double foo [], int N) {
	int maxIndex = -1;
    double maxVal(-1.0);
    for (int k=0; k<N; ++k)
      if (foo[k] > maxVal) {
    	  maxVal = foo[k];
    	  maxIndex = k;
      }
	std::cout << "Best guess: "<< maxIndex << " white brick(s) with probability: "<< maxVal << ".\n";
  }

int main() {

	//number of LEGO mushrooms
	int l = 4;
	//number of white bricks
	int k = 2;
	int N = k*l+1;
	// allocate result array
	double * results = new double[N];
	//fill results array with 1 to multiply with initially
	std::fill_n(results, N, double(1));
	//help array to save previous maximums
	double * help = new double[N];

	//create weight-matrix [l][k+1]
	double WB [l][k+1];
	for(size_t j = 0; j < l; ++j){
		int sum=0;
		for(int i=0; i<k+1; i++){
			WB[j][i] = (rand()%100)+1;
		    sum = sum + WB[j][i];
		 }
		for(int i=0; i<k+1; i++){
			WB[j][i] = WB[j][i]/sum;
		}
	}


	// INITIALIZATION: LEGO 0
	for(int i = 0; i<=k; i++){
		results[i] = double(WB[0][i]);
		help[i]= results[i];
	}

	// ITERATE over LEGO 1 to l-1
	for(int a = 1; a<l; a++){
		// ITERATE over Levels
		for(int b = 0; b<a*k+1; b++){

			//update all values in the help-array (except in first run)
			if(b==0 && a !=1){
				for(int x =0; x<= a*k; x++){
					help[x] = results[x];
				}
			}

			// ITERATE over k+1 weights
			for(int c = 0; c<=k; c++){
				//calculating at position?
				int pos = b+c;
				//save new value in variable current - calc with the help array
				double current = help[b] * double(WB[a][c]);

				//if in level 0 set initial result (for later comparison)
				if(b==0){
					results[pos] = current;
					continue;
				}
				//if weight k set also initial result (for later comparison)
				if(c==k){
					results[pos] = current;
					continue;
				}
				//if current is bigger than the val in results
				if(results[pos] < current){
					//save current as better result
					results[pos] = current;
				}
			}
		}
	}

	for(int n =0;n<N; n++){
		std::cout << n <<": " << results[n] << "\n";
	}

	//get index of max
	maxR(results,N);
	return 0;
}
