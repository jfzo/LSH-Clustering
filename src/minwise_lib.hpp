/*
 * minwise_lib.hpp
 *
 *  Created on: 07-11-2014
 *      Author: juan
 */

#ifndef MINWISE_LIB_HPP_
#define MINWISE_LIB_HPP_
#include <random>
#include <fstream>
#include <algorithm>
#include <cstring>
#include <limits>
#include <unordered_set>
#include <cstdlib>
#include <set>
//#include "data_definitions.hpp"
#include <cstdlib>
//#include "linear_hashing.hpp"
#include <cmath>
#include <chrono>
#include <functional>

#include <time.h>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/uniform_int_distribution.hpp>

#define UINT_T uint32_t
#define SGN(X) copysign(1.0, X) // signum function macro.

/*
 class BinarySet {
 private:
 bool * data;

 public:
 size_t size;

 BinarySet(size_t _size) {
 size = _size;
 data = new bool[size];
 memset(data, 0, size * sizeof(bool));
 }

 ~BinarySet() {
 delete[] data;
 }

 bool operator [](int i) const {
 assert(i < size);
 return data[i];
 }
 bool & operator [](int i) {
 assert(i < size);
 return data[i];
 }

 // TO IMPLEMENT
 unsigned int s_intersect(const BinarySet & other) const {
 unsigned int int_s = 0;

 return int_s;
 }

 unsigned int s_difference(const BinarySet & other) const {
 unsigned int int_s = 0;
 // A - B

 return int_s;
 }
 };
 */


UINT_T hash(const UINT_T* v, size_t size) {
	uint64_t z, z2, s, zi, xi;
	int bigprime = 0xfffffffb;//2^32 - 5
	z = 0x64b6055aL;
	z2 = 0x5067d19d;
	s = 0;
	zi = 1;

	//std::cout << "applying the right function for int array of size "<< size <<"\n";

	for (size_t i = 0; i < size; i++) {
		xi = (v[i] * z2) >> 1;
		s = (s + zi * xi) % bigprime;
		zi = (zi * z) % bigprime;
	}

	s = (s + zi * (bigprime - 1)) % bigprime;
	UINT_T hval = (int) s;
	//if (hval < 0)
	//	hval += bigprime;

	return hval;
}

/*
Not the best implementation, but allows to
avoid using Murmurhash library.
*/
UINT_T hash(int key, int a, int b){	
	int P = 0xfffffffb;//2^32 - 5	
	return (unsigned)( a * key + b) % P;
}


class BinarySet {
private:
	//bool * data;
	std::unordered_set<int> data;
	//std::set<int> data;

public:
	size_t size;

	BinarySet(size_t _size) {
		size = _size;
		//data = new bool[size];
		//memset(data, 0, size * sizeof(bool));
	}

	~BinarySet() {
		//delete[] data;
	}

	void clear(){
		data.clear();
	}

	bool operator [](int i) const {
		//assert(i < size);
		//return data[i];
		return data.find(i) != data.end();
	}

	void set(int i) {
		data.insert(i);
	}

	bool is_set(int i) const {
		return (data.find(i) != data.end());
	}

	unsigned int s_intersect(const BinarySet & other) const {
		unsigned int int_s = 0;
		auto it = data.cbegin();
		for (; it != data.cend(); ++it) {
			if (other.is_set(*it))
				int_s++;
		}
		return int_s;
	}

	unsigned int s_difference(const BinarySet & other) const {
		unsigned int int_s = 0;
		// A - B
		auto it = data.cbegin();
		for (; it != data.cend(); ++it) {
			if (!other.is_set(*it))
				int_s++;
		}
		return int_s;
	}

};


class RealVector {
private:
	//bool * data;
	std::unordered_map<int, double> data;
	//std::set<int> data;

public:
	size_t size;

	RealVector(size_t _size) {
		size = _size;
		//data = new bool[size];
		//memset(data, 0, size * sizeof(bool));
	}

	~RealVector() {
		//delete[] data;
	}

	void clear(){
		data.clear();
	}

	void set(int i, double v) {
		data[i] = v;
	}
    
    const double dot(const double * v) const{
        double result = 0.0;
        std::unordered_map<int, double>::const_iterator it = data.cbegin();
        for(; it != data.cend(); ++it){
            result += v[it->first] * it->second;
            //std::cout << v[it->first] <<"*"<< it->second << "("<< result<<") ";
        }            
        //std::cout << "\n";
        return result;
    }

	double operator [](int i) const {
		//assert(i < size);
		//return data[i];
		std::unordered_map<int, double>::const_iterator item = data.find(i);
		if( item != data.end())
			return item->second;
		else
			return 0.0;
	}
	
	void get_double_vector(double *p){
		for(size_t i=0; i<size; i++){
			p[i] = 0.0;
		}
		std::unordered_map<int, double>::const_iterator item = data.cbegin();
		for(; item!=data.cend(); ++item)
			p[item->first] = item->second;
	}
};

class RandomHyperplaneGenerator {
private:
	double ** rhs; // 1 random hyperplane vector per column. (rv_dim x n_hashes)
    // std::hash<std::string> band_hash; // Employed by the commented get_signature
	
public:
	int nbands, r; // nr of bands and integers per band
	int n_hashes; // nr of minwise estimates.
	size_t rv_dim;
	
	RandomHyperplaneGenerator(int _nbands, int _r, size_t _rv_dim) :
			nbands(_nbands), r(_r), rv_dim(_rv_dim) {
		n_hashes = nbands * r;
		// compute and store hyperplanes

        boost::variate_generator<boost::mt19937, boost::normal_distribution<> > generator(boost::mt19937(time(0)), boost::normal_distribution<>());

        
        /*
		unsigned vseed = std::chrono::system_clock::now().time_since_epoch().count();
		std::mt19937 rgen(vseed);
		std::normal_distribution<double> rndnorm (0.0,1.0);
		*/		
		rhs = new double * [n_hashes];
		for(int i=0; i<n_hashes; i++){
			rhs[i] = new double[rv_dim];
			for(size_t j=0; j<rv_dim; j++){
                double value = generator();
                if(value >= 0)
                    rhs[i][j] = 1;
                else
                    rhs[i][j] = -1;
            }
		}
        std::cout << "Random vectors initialized.\n";
	}
	
	~RandomHyperplaneGenerator(){
		for(int i=0; i<n_hashes; i++)
			delete[] rhs[i];
		delete[] rhs;
	}
	/**
	Computes the signature for the input vector and completely stores it
	into the signature pointer (nbands*r size).
	**/
	void get_complete_signature(const RealVector * in_vec, uint8_t *signature) {
		size_t in_size = in_vec->size;
		assert(in_size == rv_dim);
		
		for(int i=0; i<n_hashes; i++){
			double iprod = 0.0;
			for(size_t j=0; j<rv_dim; j++){
				iprod += rhs[i][j] * in_vec->operator[](j);
            }
            
			if(iprod >= 0)
				signature[i] = 1;
			else
				signature[i] = 0;
            //std::cout << "[("<< iprod<< ")"<< signature[i] << "]";
		}
        //std::cout << "\n";
		
		// signature must have n_hashes components !!		
	}
    
	/**
	Computes the signature for the input vector and stores it
	into the signature pointer (nbands size).
	**/
	void get_signature(const RealVector * in_vec, uint32_t *signature) {
		size_t in_size = in_vec->size;
		assert(in_size == rv_dim);
		
        char buffer[r+1];
        buffer[r] = '\0';

        int block = 0;
        int band_cnt = 0;
        for(int i=0; i<n_hashes; i++){
            double iprod = in_vec->dot(rhs[i]);
			if(iprod >= 0)
				buffer[block] = '1';
			else
				buffer[block] = '0';
            block++;
            
            if(block == r){
                block = 0;
                // Two options for generating int for band: string-hash and binary-to-integer.
                //uint32_t hband = band_hash(std::string(buffer)); // option A
                uint32_t hband = (uint32_t)strtol(&buffer[0], nullptr, 2);// option B
                
                //printf("\n\t_band %d : %s (%u)\n",band_cnt+1, buffer, hband );
                signature[band_cnt] = hband; // (Prime number 2^32 - 5)
                band_cnt++;      
            }            
        }
	}   
    
	/**
	Computes the signature for the input vector and stores it
	into the signature pointer (nbands size). The third parameter is
    employed as a dummy parameter in order to keep compliance with existing calls.
	**/
	void get_signature(const RealVector * in_vec, uint32_t *signature, uint32_t max_val) {
        get_signature(in_vec, signature);
	}   
    
    /*
	void get_signature(const RealVector * in_vec, uint32_t *signature, uint32_t max_hash) {
		size_t in_size = in_vec->size;
		assert(in_size == rv_dim);
		
        char buffer[r+1];
        buffer[r] = '\0';

        int block = 0;
        int band_cnt = 0;
        for(int i=0; i<n_hashes; i++){
            double iprod = in_vec->dot(rhs[i]);
			if(iprod >= 0)
				buffer[block] = '1';
			else
				buffer[block] = '0';
            block++;
            
            if(block == r){
                block = 0;
                uint32_t hband = band_hash(std::string(buffer));
                //printf("\n\t_band %d : %s (%u)\n",band_cnt+1, buffer, hband );
                signature[band_cnt] = hband % max_hash; // (Prime number 2^32 - 5)
                band_cnt++;      
            }            
        }
        assert( band_cnt == nbands);
        assert( block == 0);
		// signature must have nband components !!
	}    
    */
};

class MinwiseGenerator {
private:
	// Overall parameters
	//const int bigprime = 0x8000000b; //prime next to 2**31 ... When Murmur hashing was integrated this value became stale!
	// These two vectors shall be used jointly with the universal hashing function.
	// This will enable the generation of the permutations for each set representing a document.
//	int * a;
//	int * b;
	int * seeds;

	int compare_band(const UINT_T* sign1, const UINT_T* sign2) {
		/*
		 for (int i = 0; i < r; i++) {
		 std::cout << sign1[i] << "|" << sign2[i] << "  ";
		 }
		 std::cout << "\n";
		 */
		for (int i = 0; i < r; i++) {
			if (sign1[i] != sign2[i]) {
				return 0;
			}
		}
		return 1;
	}

public:
	int nbands, r; // nr of bands and integers per band
    unsigned setdim;
	int n_hashes; // nr of minwise estimates.
	UINT_T * buffer;

	MinwiseGenerator(int _nbands, int _r, unsigned _setdim) :
			nbands(_nbands), r(_r), setdim(_setdim) {
		n_hashes = nbands * r;

		std::cout << "Minwise Generator with parameters: " << "b=" << nbands
				<< " r=" << r << " num. permutations=" << n_hashes << "\n";

		//a = new int[n_hashes];
		//b = new int[n_hashes];
		seeds = new int[2*n_hashes];//stores a's and b's
		// filling the vectors.
		//unsigned vseed = std::chrono::system_clock::now().time_since_epoch().count();		
		//std::mt19937 rgen(vseed);
		//std::uniform_int_distribution<int> adist(0, bigprime - 1);
		//std::uniform_int_distribution<int> bdist(1, bigprime - 1);
        
        boost::variate_generator<boost::mt19937, boost::random::uniform_int_distribution<> > 
            generator( boost::mt19937(time(0)), boost::random::uniform_int_distribution<>(0, setdim - 1) );

        //int init_seed = rand();

		for (int i = 0; i < 2*n_hashes; i++) {
			//a[i] = adist(rgen);
			//b[i] = bdist(rgen);
            //seeds[i] = rgen();
            seeds[i] = generator();
		}

		buffer = new UINT_T[n_hashes];
	}

	

	/**
	This functions performs the same operation than the previous one,
	but it returns the complete 'n_hashes'-dimensional signature,
	instead of the 'nbands'-dimensional one.
	Note: Useful when the hash functions over each band is going to
	be performed outside this class.
	*/
	void getCompleteSignature(const BinarySet * in_set, UINT_T *signature) {
		size_t in_size = in_set->size;
		// signature must have n_hashes components !!

		for (int i = 0; i < n_hashes; i++)
			signature[i] = std::numeric_limits<int>::max();

		for (size_t i = 0; i < in_size; i++) { // for each feature in the universal set.
			if (in_set->operator[](i)) { // continues only if the feature is present in the current input set.
				for (int j = 0; j < n_hashes; j++) {
					//UINT_T h_j = hash(i, a[j], b[j], bigprime) % in_size; // recall that's a permutation--> [0, max_feature_index]
					UINT_T h_j = hash(i, seeds[j], seeds[j+n_hashes]) % in_size; // recall that's a permutation--> [0, max_feature_index]
                    // recall that if in_size is prime, then the permutation is perfect.
					signature[j] = std::min(signature[j], h_j);
				}
			}
		}
	}

    /**
     * Computes the minwise signature, given an input set and a previously allocated signature array.
     * in_set: binary array that holds the input set.
     * signature: array having nbands components.
     *
     * Note: This function resets every value of the signature array given.
     */
    void get_signature(const BinarySet * in_set, UINT_T *signature) {
    	size_t in_size = in_set->size;

    	for (int i = 0; i < n_hashes; i++)
    		buffer[i] = std::numeric_limits<int>::max();

    	for (size_t i = 0; i < in_size; i++) { // for each feature in the universal set.
    		if (in_set->operator[](i)) { // continues only if the feature is present in the current input set.
    			for (int j = 0; j < n_hashes; j++) {
    				//UINT_T h_j = hash(i, a[j], b[j], bigprime) % in_size; // recall that's a permutation--> [0, max_feature_index]
    				UINT_T h_j = hash(i, seeds[j], seeds[j+n_hashes]) % in_size; // recall that's a permutation--> [0, max_feature_index]
                       // recall that if in_size is prime, then the permutation is perfect.
    				buffer[j] = std::min(buffer[j], h_j);
    			}
    		}
    	}

    	//Process the buffer and fill in the signature.
    	for (int i = 0; i < nbands; i++) {
    		signature[i] = hash(&buffer[i * r], r);
    		//std::cout << "signature["<< i << "] "<< signature[i]<<"\n ";
    	}
    	//print_signature(signature);
    	//print_buffer();
    }


    /**
     * Computes the minwise signature, given an input set and a previously allocated signature array and
     * outputs into the signture array MOD maxval.
     * in_set: binary array that holds the input set.
     * signature: array having nbands components.
     * maxval: uint32_t number denoting the max hash value allowed foe every coordinate. 
     *
     * Note: This function resets every value of the signature array given.
     */
    void get_signature(const BinarySet * in_set, UINT_T *signature, uint32_t maxval) {
    	size_t in_size = in_set->size;

    	for (int i = 0; i < n_hashes; i++)
    		buffer[i] = std::numeric_limits<int>::max();

    	for (size_t i = 0; i < in_size; i++) { // for each feature in the universal set.
    		if (in_set->operator[](i)) { // continues only if the feature is present in the current input set.
    			for (int j = 0; j < n_hashes; j++) {
    				//UINT_T h_j = hash(i, a[j], b[j], bigprime) % in_size; // recall that's a permutation--> [0, max_feature_index]
    				UINT_T h_j = hash(i, seeds[j], seeds[j+n_hashes]) % in_size; // recall that's a permutation--> [0, max_feature_index]
                       // recall that if in_size is prime, then the permutation is perfect.
    				buffer[j] = std::min(buffer[j], h_j);
    			}
    		}
    	}

    	//Process the buffer and fill in the signature.
    	for (int i = 0; i < nbands; i++) {
    		signature[i] = hash(&buffer[i * r], r) % (maxval + 1);
    		//std::cout << "signature["<< i << "] "<< signature[i]<<"\n ";
    	}
    	//print_signature(signature);
    	//print_buffer();
    }


	
	/**
	 * compare bands and returns an approximate of the jaccard c.
	 */
	double compare_signatures(const UINT_T* sign1, const UINT_T* sign2) {
		// traverse each signature per band
		double matches = 0;

		/*
		 for (int bi = 0; bi < n_hashes; bi += r) {
		 matches += compare_band(&sign1[bi], &sign2[bi]);
		 }
		 */
		for (int bi = 0; bi < nbands; bi++) {
			matches += (sign1[bi] == sign2[bi]) ? 1 : 0;
		}

		/*
		 std::cout << "Signature...\n";
		 for (int i = 0; i < n_hashes; i++)
		 std::cout << " " << sign1[i];
		 std::cout << "\n";
		 for (int i = 0; i < n_hashes; i++)
		 std::cout << " " << sign2[i];
		 std::cout << "\n";

		 std::cout << "#matches: " << matches <<"/" << nbands << "("<< (matches/nbands)<<")"<< "\n";
		 */
		return matches / nbands;
	}

	void print_signature(const UINT_T* sign) {
		for (int bi = 0; bi < nbands; bi++) {
			std::cout << "|" << sign[bi] << " ";
		}
		std::cout << "|\n";
	}

	void print_buffer() {
		for (int bi = 0; bi < n_hashes; bi++) {
			if (bi % r == 0)
				std::cout << "| ";
			std::cout << buffer[bi] << " ";
		}
		std::cout << "|\n";
	}

	~MinwiseGenerator() {
		//delete[] a;
		//delete[] b;
		delete[] seeds;
		delete[] buffer;
	}
};

/**
 * Penalized Random Hyperplanes
 * */
class PRandomHyperplaneGenerator {
private:
	double ** rhs; // 1 random hyperplane vector per column. (rv_dim x n_hashes)
    // std::hash<std::string> band_hash; // Employed by the commented get_signature
	
public:
	int nbands; // nr of bands and integers per band
	int n_hashes; // nr of minwise estimates.
	size_t rv_dim;
	
	PRandomHyperplaneGenerator(int _nbands, size_t _rv_dim) :
			nbands(_nbands), rv_dim(_rv_dim) {
		n_hashes = nbands;
		// compute and store hyperplanes

        boost::variate_generator<boost::mt19937, boost::normal_distribution<> > generator(boost::mt19937(time(0)), boost::normal_distribution<>());

        
        /*
		unsigned vseed = std::chrono::system_clock::now().time_since_epoch().count();
		std::mt19937 rgen(vseed);
		std::normal_distribution<double> rndnorm (0.0,1.0);
		*/		
		rhs = new double * [n_hashes];
		for(int i=0; i<n_hashes; i++){
			rhs[i] = new double[rv_dim];
			for(size_t j=0; j<rv_dim; j++){
                double value = generator();
                rhs[i][j] = SGN(value); // +1 or -1
                //std::cout << rhs[i][j] << " ";
            }
            //std::cout<< "\n";
		}
        //std::cout << "Random vectors initialized.\n";
	}
	
	~PRandomHyperplaneGenerator(){
		for(int i=0; i<n_hashes; i++)
			delete[] rhs[i];
		delete[] rhs;
	}

	/**
	Computes the signature for the input vector and stores it
	into the signature pointer (nbands size).
	**/
	void get_signature(const RealVector * in_vec, double *signature) {
		size_t in_size = in_vec->size;
		assert(in_size == rv_dim);
		
        //std::cout << "Getting signature for input doc.\n";
        for(int i=0; i<n_hashes; i++){
            signature[i] = in_vec->dot(rhs[i]);
            //std::cout << signature[i] << " ";
        }
        //std::cout << "\n";
		// signature must have nband components !!
	}   
    
    /*
	void get_signature(const RealVector * in_vec, uint32_t *signature, uint32_t max_hash) {
		size_t in_size = in_vec->size;
		assert(in_size == rv_dim);
		
        char buffer[r+1];
        buffer[r] = '\0';

        int block = 0;
        int band_cnt = 0;
        for(int i=0; i<n_hashes; i++){
            double iprod = in_vec->dot(rhs[i]);
			if(iprod >= 0)
				buffer[block] = '1';
			else
				buffer[block] = '0';
            block++;
            
            if(block == r){
                block = 0;
                uint32_t hband = band_hash(std::string(buffer));
                //printf("\n\t_band %d : %s (%u)\n",band_cnt+1, buffer, hband );
                signature[band_cnt] = hband % max_hash; // (Prime number 2^32 - 5)
                band_cnt++;      
            }            
        }
        assert( band_cnt == nbands);
        assert( block == 0);
		// signature must have nband components !!
	}    
    */
};

double compute_jaccard(const int * s1, const int * s2, size_t ssize) {
	double len_inters = 0;
	double len_union = 0;
	//compute length of intersection
	for (size_t i = 0; i < ssize; i++) {
		if (s1[i] == 1 && s2[i] == 1) {
			len_inters++;
			len_union++;
		}
		if (s1[i] != s2[i] && (s1[i] == 1 || s2[i] == 1))
			len_union++;
	}

	return len_inters / len_union;
}

double compute_jaccard(const BinarySet & s1, const BinarySet & s2) {
	assert(s1.size == s2.size);
	double len_inters = s1.s_intersect(s2);
	double len_union = len_inters + s1.s_difference(s2) + s2.s_difference(s1);
	/*
	 size_t ssize = s1.size;

	 //compute length of intersection
	 for (int i = 0; i < ssize; i++) {
	 if (s1[i] && s2[i]) {
	 len_inters++;
	 len_union++;
	 }
	 if (s1[i] != s2[i] && (s1[i] || s2[i]))
	 len_union++;
	 }
	 */

	return len_inters / len_union;
}

/**
 * Wrapper class for a binary matrix (contained as attribute data)
 */
class Setdata {
	int nsets; // Nr. of sets grouped into the structure.
	int nelems; // Nr. of elements of the universal set.

public:
	//bool ** data;
	std::vector<BinarySet*> data; // matrix encoded as a vector of pointers to a Binary set.
	std::vector<std::string> labels; // matrix encoded as a vector of pointers to a Binary set.

	Setdata() {
		nsets = 0;
		nelems = 0;
	}

	int dims(){
		return nelems;
	}
	
	// keeped because of compatibility with older code.
	int sets(){
		return instances();
	}

	int instances() const{
		return nsets;
	}
	
	
	Setdata(int _nsets, int _nelems) {
		initialize(_nsets, _nelems);
	}

	void initialize(int _nsets, int _nelems) {
		nsets = _nsets;
		nelems = _nelems;

        labels.reserve(nsets);
		for (int i = 0; i < nsets; i++) {
			data.push_back(new BinarySet(nelems));
		}

		/*
		 if (nsets != 0 && nelems != 0) {

		 //data = new bool*[nsets];

		 for (int i = 0; i < nsets; i++) {
		 data[i] = new bool[nelems];
		 memset(data[i], 0, nelems * sizeof(bool));
		 }
		 }
		 */
	}

	const BinarySet & operator [](int i) const {
		return *(data[i]);
	}

	BinarySet & operator [](int i) {
		return *(data[i]);
	}

	//int& operator[](unsigned int i){ return pA[i];}
	//const int& operator[](unsigned int i)const{ return pA[i];}

	void print_sets() {
		std::cout << "N.SETS: " << nsets << " N.ELEMS: " << nelems << "\n";
		if (nsets != 0 && nelems != 0) {
			for (int i = 0; i < nsets; i++) {
				std::cout << "{ ";
				for (int j = 0; j < nelems; j++)
					std::cout << data[i]->operator[](j) << " ";
				std::cout << "}\n";
			}
		}
	}

	~Setdata() {
		for (int i = 0; i < nsets; i++)
			delete data[i];
		/*
		 if (nsets != 0 && nelems != 0) {
		 for (int i = 0; i < nsets; i++)
		 delete[] data[i];
		 delete[] data;
		 }
		 */
	}

};

/**
 * DistanceMatrix class for binary data representing sets.
 *
 */
class SetSimilarityMatrix {
private:
	double ** dmat;
	const Setdata * data;

	void compute() {
		for (size_t rowid = 0; rowid < msize - 1; rowid++) {
			for (size_t colid = rowid + 1; colid < msize; colid++) {
				double dval = compute_jaccard(data->operator[](rowid),
						data->operator[](colid));
				dmat[rowid][colid] = dval;
				dmat[colid][rowid] = dval;
			}
		}
	}

public:
	size_t msize; //square matrix

	SetSimilarityMatrix(const Setdata * data_) {
		data = data_;
		msize = data->instances();

		dmat = new double*[msize];
		for (size_t i = 0; i < msize; i++)
			dmat[i] = new double[msize];
		compute();
	}

	~SetSimilarityMatrix() {
		//deallocating
		for (size_t i = 0; i < msize; i++)
			delete[] dmat[i];
		delete[] dmat;
	}

	void write_for_cluto(std::string fname) {
		FILE *spMatFp;
		spMatFp = fopen(fname.c_str(), "w");
		//header
		fprintf(spMatFp, "%zu\n", msize);

		for (size_t lnum = 0; lnum < msize; lnum++) {
			for (size_t cnum = 0; cnum < msize; cnum++) {
				if (lnum == cnum)
					fprintf(spMatFp, "1.0 ");
				else {
					fprintf(spMatFp, "%f ", dmat[lnum][cnum]);
				}
			}
			fprintf(spMatFp, "\n");
		}
		fclose(spMatFp);
	}
};

/**
 * Function that reads a CSV file containing a set expressed in a compressed way.
 * The file header contains the nr. of rows and the nr. of elements in the universe set separated by space or comma.
 * E.g. A line containing '7,1, 2, 4, 8, 9, 11, 14' denotes a set having 7 elements:
 * these are: 1,2,4,8,9,11,14 .
 * Note. The first line is employed as a header containig the number of sets in the file and
 * 			the cardinality of the Universe set.
 *
 */
/*
 void read_set_data(const char * fname, Setdata & sdata) {
 std::ifstream inf(fname, std::iostream::in);
 char buf[256];

 //reading header
 int NSETS, NELEMS;
 //inf >> NSETS >> NELEMS;
 //std::cout << "N.SETS: " << NSETS << " N.ELEMS: " << NELEMS << "\n";

 inf.getline(buf, 256);
 char * pch;
 pch = strtok(buf, " ,");
 NSETS = atoi(pch);
 pch = strtok(NULL, " ,.-");
 NELEMS = atoi(pch);

 sdata.initialize(NSETS, NELEMS);

 int setid = 0;
 while (setid < NSETS) {
 inf.getline(buf, 256);
 pch = strtok(buf, " ,");
 pch = strtok(NULL, " ,");    //ignore the first field (set length)
 while (pch != NULL) {
 int n = atoi(pch);
 sdata.data[setid][n] = true;
 pch = strtok(NULL, " ,");
 }
 setid++;
 }
 inf.close();
 }
 */

/**
 * Wrapper class for a double matrix (contained as attribute data)
 */
class Realdata {
public:
	int ninstances; // Nr. of sets grouped into the structure.
	int nfeatures; // Nr. of elements of the universal set.
	//bool ** data;
	std::vector<RealVector*> data; // matrix encoded as a vector of pointers to a Binary set.
	std::vector<std::string> labels; // matrix encoded as a vector of pointers to a Binary set.

	Realdata() {
		ninstances = 0;
		nfeatures = 0;
	}

	int dims(){
		return nfeatures;
	}
	
	int instances(){
		return ninstances;
	}
	
	Realdata(int _ninstances, int _nfeatures) {
		initialize(_ninstances, _nfeatures);
	}

	void initialize(int _ninstances, int _nfeatures) {
		ninstances = _ninstances;
		nfeatures = _nfeatures;

        labels.reserve(ninstances);
		for (int i = 0; i < ninstances; i++) {
			data.push_back(new RealVector(nfeatures));
		}

		/*
		 if (nsets != 0 && nelems != 0) {

		 //data = new bool*[nsets];

		 for (int i = 0; i < nsets; i++) {
		 data[i] = new bool[nelems];
		 memset(data[i], 0, nelems * sizeof(bool));
		 }
		 }
		 */
	}

	const RealVector & operator [](int i) const {
		return *(data[i]);
	}

	RealVector & operator [](int i) {
		return *(data[i]);
	}

	//int& operator[](unsigned int i){ return pA[i];}
	//const int& operator[](unsigned int i)const{ return pA[i];}

	void print_vectors() {
		std::cout << "N.VECTORS: " << ninstances << " N.FEATS: " << nfeatures << "\n";
		if (ninstances != 0 && nfeatures != 0) {
			for (int i = 0; i < ninstances; i++) {
				std::cout << "{ ";
				for (int j = 0; j < nfeatures; j++)
					std::cout << data[i]->operator[](j) << " ";
				std::cout << "}\n";
			}
		}
	}

	~Realdata() {
		for (int i = 0; i < ninstances; i++)
			delete data[i];
		/*
		 if (nsets != 0 && nelems != 0) {
		 for (int i = 0; i < nsets; i++)
		 delete[] data[i];
		 delete[] data;
		 }
		 */
	}

};

/**
 * Function that reads a CSV file containing a set expressed as a sequence of 1's nad 0's.
 * The file header contains the nr. of rows and the nr. of elements in the universe set separated by space or comma.
 * E.g. A line containing '1, 0, 0, 0, 1, 1, 1' denotes a set having 4 elements(0,4,5,6) from
 * a universe set having size 7.
 *
 */
void read_csv_set_data(const char * fname, Setdata & sdata) {
	std::ifstream inf(fname, std::iostream::in);
	size_t MAXLINELEN = 2097152;// max line length (not including newline char)
	//MAXLINELEN = (size_t)ceil(log2(NELEMS));
	char buf[MAXLINELEN];

	//reading header
	int NSETS, NELEMS;
	//inf >> NSETS >> NELEMS;
	//std::cout << "N.SETS: " << NSETS << " N.ELEMS: " << NELEMS << "\n";

	inf.getline(buf, MAXLINELEN);
	char * pch;
	// counting the spaces.
	pch = strpbrk(buf, " ");
	int nspaces = 0;
	while (pch != NULL) {
		nspaces++;
		pch = strpbrk(pch + 1, " ");
	}
	//std::cout << "Nr. of spaces found:"<<nspaces<<"\n";
	// nspaces > 1  =>  sparse format.

	pch = strtok(buf, " ,");
	NSETS = atoi(pch);
	pch = strtok(NULL, " ,.-");
	NELEMS = atoi(pch);

	//std::cout << "N.SETS: " << NSETS << " N.ELEMS: " << NELEMS << "\n";

	sdata.initialize(NSETS, NELEMS);

	if (nspaces < 2) {
		int setid = 0;
		int items_counter;
		while (setid < NSETS) {
			inf.getline(buf, MAXLINELEN);
			items_counter = 0;
			assert((inf.rdstate() & std::ifstream::failbit) == 0); //limit is reached without finding the delimiting character, the failbit internal flag is set

			pch = strtok(buf, " ,");
			while (pch != NULL) {
				int n = atoi(pch);
				if (n == 1) {
					sdata[setid].set(items_counter);
					//sdata[setid][items_counter] = true;
				}
				items_counter++;
				pch = strtok(NULL, " ,");
			}
			setid++;
		}
	} else {
		//std::cout << "Sparse format reading not implemented yet!\n";
		//assert(nspaces < 2);

		int setid = 0;
		while (setid < NSETS) {
			inf.getline(buf, MAXLINELEN);
			assert((inf.rdstate() & std::ifstream::failbit) == 0); //MAXLINELEN is reached, increase it !

			pch = strtok(buf, " ,");
			while (pch != NULL) {
				int n = atoi(pch); //feature id (in clusto's format starts at 1.
				sdata[setid].set( n - 1 );
				//std::cout << "read feature: "<<n<<"\n";
				pch = strtok(NULL, " ,"); // feature value (not used unless other value different than 1 appears.
				//double fvalue = atof(pch);
				pch = strtok(NULL, " ,"); // next feature id.
			}
			setid++;
		}
	}
	inf.close();
}


/**
 * Function that reads a CSV file containing a set expressed as a sequence of 1's nad 0's or in sparse format.
 * The file header contains the nr. of rows and the nr. of elements in the universe set separated by space or comma.
 * E.g. A line containing '1, 0, 0, 0, 1, 1, 1' denotes a set having 4 elements(0,4,5,6) from
 * a universe set having size 7.
 * Also a label file is passed as a parameter containing a str label for each row in the input set file.
 */
void read_csv_set_data(const char * fname, const char * lbfname, Setdata & sdata) {
	std::ifstream lbinf(lbfname, std::iostream::in);
	std::ifstream inf(fname, std::iostream::in);
	size_t MAXLINELEN = 2097152;// max line length (not including newline char)
	//MAXLINELEN = (size_t)ceil(log2(NELEMS));
	char * buf = new char[MAXLINELEN];

	//reading header
	int NSETS, NELEMS;
	//inf >> NSETS >> NELEMS;
	//std::cout << "N.SETS: " << NSETS << " N.ELEMS: " << NELEMS << "\n";

	inf.getline(buf, MAXLINELEN);
	char * pch;
	// counting the spaces.
	pch = strpbrk(buf, " ");
	int nspaces = 0;
	while (pch != NULL) {
		nspaces++;
		pch = strpbrk(pch + 1, " ");
	}
	//std::cout << "Nr. of spaces found:"<<nspaces<<"\n";
	// nspaces > 1  =>  sparse format.

	pch = strtok(buf, " ,");
	NSETS = atoi(pch);
	pch = strtok(NULL, " ,.-");
	NELEMS = atoi(pch);

	//std::cout << "N.SETS: " << NSETS << " N.ELEMS: " << NELEMS << "\n";

	sdata.initialize(NSETS, NELEMS);

	if (nspaces < 2) {
		int setid = 0;
		int items_counter;
		while (setid < NSETS) {
			inf.getline(buf, MAXLINELEN);
			items_counter = 0;
			assert((inf.rdstate() & std::ifstream::failbit) == 0); //limit is reached without finding the delimiting character, the failbit internal flag is set

			pch = strtok(buf, " ,");
			while (pch != NULL) {
				int n = atoi(pch);
				if (n == 1) {
					sdata[setid].set(items_counter);
					//sdata[setid][items_counter] = true;
				}
				items_counter++;
				pch = strtok(NULL, " ,");
			}

            std::string lbl;
            lbinf >> lbl;
            sdata.labels.push_back(lbl);
            //std::cout << "LABEL("<< setid <<") "<< lbl <<"\n";
			setid++;
		}
	} else {
		//std::cout << "Sparse format reading not implemented yet!\n";
		//assert(nspaces < 2);

		int setid = 0;
		while (setid < NSETS) {
			inf.getline(buf, MAXLINELEN);
			assert((inf.rdstate() & std::ifstream::failbit) == 0); //MAXLINELEN is reached, increase it !

			pch = strtok(buf, " ,");
			while (pch != NULL) {
				int n = atoi(pch); //feature id (in clusto's format starts at 1.
				sdata[setid].set( n - 1 );
				//std::cout << "read feature: "<<n<<"\n";
				pch = strtok(NULL, " ,"); // feature value (not used unless other value different than 1 appears.
				//double fvalue = atof(pch);
				pch = strtok(NULL, " ,"); // next feature id.
			}
            std::string lbl;
            lbinf >> lbl;
            sdata.labels.push_back(lbl);
            //std::cout << "LABEL("<< setid <<") "<< lbl <<"\n";
			setid++;
		}
	}
	inf.close();
	lbinf.close();
	delete[] buf;
}


/**
 * Reads data in cluto's format (sparse or non sparse autodetection).
 * Stores the data and their labels into a Realdata object.
 * */
void read_csv_real_data(const char * fname, const char * lbfname, Realdata & sdata) {
	std::ifstream lbinf(lbfname, std::iostream::in);
	std::ifstream inf(fname, std::iostream::in);
	size_t MAXLINELEN = 2097152;// max line length (not including newline char)
	//MAXLINELEN = (size_t)ceil(log2(NELEMS));
	char * buf = new char[MAXLINELEN];

	//reading header
	int NINSTANCES, NFEATURES;
	//inf >> NFEATURES >> NFEATURES;
	//std::cout << "N.SETS: " << NINSTANCES << " N.ELEMS: " << NFEATURES << "\n";

	inf.getline(buf, MAXLINELEN);
	char * pch;
	// counting the spaces.
	pch = strpbrk(buf, " ");
	int nspaces = 0;
	while (pch != NULL) {
		nspaces++;
		pch = strpbrk(pch + 1, " ");
	}
	//std::cout << "Nr. of spaces found:"<<nspaces<<"\n";
	// nspaces > 1  =>  sparse format.

	pch = strtok(buf, " ,");
	NINSTANCES = atoi(pch);
	pch = strtok(NULL, " ,.-");
	NFEATURES = atoi(pch);

	//std::cout << "N.SETS: " << NINSTANCES << " N.ELEMS: " << NFEATURES << "\n";

	sdata.initialize(NINSTANCES, NFEATURES);

	if (nspaces < 2) { // i.e. Non-sparse file format found.
		int setid = 0;
		int items_counter;
		while (setid < NINSTANCES) {
			inf.getline(buf, MAXLINELEN);
			items_counter = 0;
			
			assert((inf.rdstate() & std::ifstream::failbit) == 0); //limit is reached without finding the delimiting character, the failbit internal flag is set

			pch = strtok(buf, " ,");
			while (pch != NULL) { // traversing the line content.
				double value = atof(pch);
				sdata[setid].set(items_counter, value);
				items_counter++;
				pch = strtok(NULL, " ,");
			}

            std::string lbl;
            lbinf >> lbl;
            sdata.labels.push_back(lbl);
            //std::cout << "LABEL("<< setid <<") "<< lbl <<"\n";
			setid++;
		}
	} else {
		//std::cout << "Sparse format reading not implemented yet!\n";
		//assert(nspaces < 2);
		int setid = 0;
		while (setid < NINSTANCES) {
			inf.getline(buf, MAXLINELEN);
			assert((inf.rdstate() & std::ifstream::failbit) == 0); //MAXLINELEN is reached, increase it !

			pch = strtok(buf, " ,");
			while (pch != NULL) {
				int n = atoi(pch); //feature id (in clusto's format starts at 1.
				//std::cout << "read feature: "<<n<<"\n";
				pch = strtok(NULL, " ,"); // feature value (not used unless other value different than 1 appears.
				double value = atof(pch);
				sdata[setid].set( n - 1, value );
				
				pch = strtok(NULL, " ,"); // next feature id.
			}
            std::string lbl;
            lbinf >> lbl;
            sdata.labels.push_back(lbl);
            //std::cout << "LABEL("<< setid <<") "<< lbl <<"\n";
			setid++;
		}
	}
	inf.close();
	lbinf.close();
	delete[] buf;
}

#endif /* MINWISE_LIB_HPP_ */
