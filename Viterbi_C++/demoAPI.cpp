#include <iostream>
#include <fstream>
#include <array>
#include "LogDouble.hpp"
#include <sstream>

// for a model with K states
template <std::size_t K>
class HMMParameters {
private:
  std::array<double, K> priorProbabilities;
  std::array<std::array<double, K>, K> transitionProbabilities;
  std::array<std::array<double, 256>, K> emissionProbabilities;

public:

  static const std::array<char, 4> emissionCharacters;

  // Prior: Pr(H_0 = state)
  double getPriorProbability(int state) const {
    return priorProbabilities[state];
  }

  // Transition: Pr(H_{i+1} = toState | H_i = fromState)
  double getTransitionProbability(int fromState, int toState) const {
    return transitionProbabilities[fromState][toState];
  }

  // Emission: Pr(D_i = output | H_i = state)
  double getEmissionProbability(int state, char output) const {
    double result = emissionProbabilities[state][(int) output];
    // Given the state, the table of emission probabilities for
    // different characters is essentially a map<char, double>.
    //
    // It's faster to use vector<double>, and use the char as an
    // integer key (i.e. the index).
    if (result >= 0.0)
      return result;

    // A negative probability indicates that an emission for that
    // character has never been initialized, and is thus not
    // supported:
    throw std::exception();
  }

  // Edit: Serang
  // Added update in priors to setParameters:
  void setParameters(const std::array<double, K> & newPriors, const std::array<std::array<double, K>, K> & newTransitionProbabilities, const std::array<std::array<double, 256>, K> & newEmissionProbabilities) {
    for (int state=0; state<K; ++state) {
      priorProbabilities[state] = newPriors[state];
      for (int toState=0; toState<K; ++toState)
	transitionProbabilities[state][toState] = newTransitionProbabilities[state][toState];
      for (int emissionChar=0; emissionChar<256; ++emissionChar)
	emissionProbabilities[state][emissionChar] = newEmissionProbabilities[state][emissionChar];
    }
  }

  // Edit: Serang
  // Nice for debugging:
  void print() {
    std::cerr << "Priors:" << std::endl;
    for (int state=0; state<K; ++state) {
      std::cerr << getPriorProbability(state) << "\t";
    }
    std::cerr << std::endl;

    std::cerr << "Transitions:" << std::endl;
    for (int state=0; state<K; ++state) {
      for (int toState=0; toState<K; ++toState)
	std::cerr << getTransitionProbability(state, toState) << "\t";
      std::cerr << std::endl;
    }

    std::cerr << "Emissions:" << std::endl;
    for (int state=0; state<K; ++state) {
      for (char emissionChar : HMMParameters<K>::emissionCharacters)
	std::cerr << getEmissionProbability(state, emissionChar) << "\t";
      std::cerr << std::endl;
    }

    std::cerr << std::endl;
  }

  HMMParameters(char*path) {
    // Initialize all probabilities with -1 so that we can easily
    // detect a query for a disallowed state/emission:
	#pragma omp parallel for
    for (int state=0; state<K; ++state) {
      priorProbabilities[state] = -1.0;
	#pragma omp parallel for
      for (int toState=0; toState<K; ++toState)
	transitionProbabilities[state][toState] = -1.0;
	#pragma omp parallel for
      for (int emissionCharIndex=0; emissionCharIndex<256; ++emissionCharIndex)
	emissionProbabilities[state][emissionCharIndex] = -1.0;
    }
    std::ifstream fin(path);

    double probability, remainingProbability;

    // read the prior probabilities (do not read the final one; it
    // will be 1 - the sum of the others)
    remainingProbability = 1.0;

    for (int state=0; state<K-1; ++state) {
      fin >> probability;
      remainingProbability -= probability;
      priorProbabilities[state] = probability;
    }
    priorProbabilities[K-1] = remainingProbability;

    // read the transition probabilities

    for (int fromState=0; fromState<K; ++fromState) {
      remainingProbability = 1.0;

      for (int toState=0; toState<K-1; ++toState) {
	fin >> probability;
	remainingProbability -= probability;
	transitionProbabilities[fromState][toState] = probability;
      }
      transitionProbabilities[fromState][K-1] = remainingProbability;
    }

    // read the emission characters and their probabilities

    for (int state=0; state<K; ++state) {
      remainingProbability = 1.0;

      for (int charIndex=0; charIndex<emissionCharacters.size()-1; ++charIndex) {
	fin >> probability;
	remainingProbability -= probability;
	emissionProbabilities[state][ emissionCharacters[charIndex] ] = probability;
      }
      emissionProbabilities[state][ emissionCharacters[emissionCharacters.size()-1] ] = remainingProbability;
    }
  }
};

template <std::size_t K>
const std::array<char, 4> HMMParameters<K>::emissionCharacters = {'G', 'A', 'T', 'C'};

// For a model with K states
template <std::size_t K>
class HMM {
private:
  const std::string * emissions;
  HMMParameters<K>*hmmParams;
public:
  HMM(const std::string & genome, HMMParameters<K> & parameters) {
    emissions = & genome;
    hmmParams = & parameters;
  }

  int argmax(const std::array<LogDouble, K> & arr) {
    int maxIndex = -1;
    LogDouble maxVal(-1.0);
    for (int k=0; k<K; ++k)
      if (arr[k] > maxVal) {
	maxVal = arr[k];
	maxIndex = k;
      }
    return maxIndex;
  }

  LogDouble max(const std::array<LogDouble, K> & arr) {
    return arr[argmax(arr)];
  }

  void getViterbiStates() {
	std::array<LogDouble, K> T1;
	std::array<LogDouble, K> comp;
	int L = emissions->size();
	// Allocate memory
	int ** T2 = new int*[2];
	#pragma omp parallel for
	for (int i = 0; i < 2; ++i)
	T2[i] = new int[L+1];

    //Initialization
	#pragma omp parallel for
    for(int i=0;i<=K-1;i++){
    	T1[i] = LogDouble(hmmParams->getPriorProbability(i)) * LogDouble(hmmParams->getEmissionProbability(i,emissions->at(0)) );
    	T2[i][0] = 0;
    }
    //Left to right
    for( int i = 1; i < L; i++ ){
    	for( int j = 0; j <= K-1; j++ ){
    		comp[0] = T1[0]*LogDouble(hmmParams->getTransitionProbability(0,j))*LogDouble(hmmParams->getEmissionProbability(j,emissions->at(i)));
    		comp[1] = T1[1]*LogDouble(hmmParams->getTransitionProbability(1,j))*LogDouble(hmmParams->getEmissionProbability(j,emissions->at(i)));
    		T1[j] = max(comp);
    		T2[j][i] = argmax(comp);
    	}
    }
    //get h_n*
    int z = argmax(T1);
    int x =T2[z][L-1];
    T2[x][L] = z;
    //likelihood
    LogDouble likelihood = LogDouble(hmmParams->getEmissionProbability(z,emissions->at(L-1)));
    //Right to left + calculate Likelihood
    for(int i = (L-1); i >= 1; i-- )
      {
    	z = T2[z][i];
    	T2[x][i] = z;
    	likelihood *= LogDouble(hmmParams->getEmissionProbability(z,emissions->at(i-1)));
      }
    for (int k=1; k<L+1; ++k)
    std::cout << T2[x][k];
    std::cout << std::endl;
    std::cout << "Likelihood: " << likelihood;
    return;
  }
};


std::string loadFile(char*path) {
  std::ifstream fin(path);
  return std::string(std::istreambuf_iterator<char>(fin), std::istreambuf_iterator<char>());
}

int main(int argc, char**argv) {

  if (argc == 3) {
    std::string genome = loadFile(argv[1]);
	std::cout << genome << std::endl;
    HMMParameters<2> hmmParams(argv[2]);
    HMM<2> hmm(genome, hmmParams);
    hmm.getViterbiStates();
  }
  else
    std::cerr << "usage: GCRichGenomeViterbi <genome> <hmm params>" << std::endl;

 	    return 0;
}
