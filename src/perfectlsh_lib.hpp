#ifndef PERFECTLSH_LIB_HPP_
#define PERFECTLSH_LIB_HPP_

#include <unordered_map>
#include <cinttypes>
#include <random>
#include <cstdlib>


template <typename IDTYPE, typename OBJECT>
class HashIndex
{
	// IDTYPE generally an integer
    typedef std::unordered_map<uint32_t, std::unordered_map<uint32_t, std::vector<IDTYPE>>> T;
	typedef std::unordered_map<uint32_t, std::vector<IDTYPE> > T2;
private:
    T table;
    uint32_t tableSize;
    size_t k;
    uint32_t seed;
    int32_t * rp; // parameters r' for t1 hash function
    int32_t * rpp;// parameters r'' for t2 hash function
    int64_t P = 4294967291; // 2^32 - 5 (recommended by Andoni and Indyk)


	/**
	1st level Hash function implementation suggested by Andoni.
	**/
    uint32_t t1b(const uint32_t * data){
    
        int64_t hi32 = 4294967295 << 32;
        int64_t lo32 = 4294967295;

        int64_t sum = 0;
        for(int i=0; i<k; i++){
            //std::cout << (int64_t)rp[i] <<"*"<< data[i] << " + ";
            //sum += ((int64_t)rp[i] * data[i]);
            int64_t mult = ((int64_t)rp[i] * data[i]);
            sum += ( (lo32 & mult) + 5*(hi32 & mult) ) % P;
        }
        //std::cout << " up & low " << (up32 & lo32) << "\n";

        //sum %= P;
        sum %= tableSize;
        return sum;
    }

	/**
	2nd level Hash function implementation suggested by Andoni.
	**/
    uint32_t t2b(const uint32_t * data){
    
        int64_t hi32 = 4294967295 << 32;
        int64_t lo32 = 4294967295;

        int64_t sum = 0;
        for(int i=0; i<k; i++){
            //std::cout << (int64_t)rp[i] <<"*"<< data[i] << " + ";
            //sum += ((int64_t)rp[i] * data[i]);
            int64_t mult = ((int64_t)rpp[i] * data[i]);
            sum += ( (lo32 & mult) + 5*(hi32 & mult) ) % P;
        }
        //std::cout << " up & low " << (up32 & lo32) << "\n";

        //sum %= P;
        //sum %= tableSize;
        return sum;
    }

    uint32_t t1b(const unsigned char * data){
		unsigned long hash = 0;
		int c;
		while ((c = *data++))
			hash = c + (hash << 6) + (hash << 16) - hash;
		
		return (uint32_t)(hash % tableSize);
	}

    uint32_t t2b(const unsigned char * data){
		unsigned long hash = 0;
		int c;
		while ((c = *data++))
			hash = c + (hash << 6) + (hash << 16) - hash;
		
		return (uint32_t)hash;
	}

	/**
	1st level hash function
	*/
    uint32_t t1(const int * data)
    {   
        //std::cout << "computing t1\n";
        int64_t sum = 0;
        for(int i=0; i<k; i++){
            //std::cout << (int64_t)rp[i] <<"*"<< data[i] << " + ";
            sum += ((int64_t)rp[i] * data[i]);
        }
        //std::cout << "\n";
        //std::cout << "sum prev. " << sum << " % "<<P <<" : " << (sum % P)<<" % " << tableSize << " : " << (sum % P)%tableSize<<"\n";
        sum %= P;
        sum %= tableSize;

        //if(sum < 0 )
        //  sum += tableSize;
        //std::cout << "t1 -->"<< sum << "\n";

        return sum;
    }

	/**
	2nd level hash function
	*/
    uint32_t t2(const int * data)
    {
        //std::cout << "computing t2\n";
        int64_t sum = 0;
        for(int i=0; i<k; i++) {
            //std::cout << rpp[i] <<"*"<< data[i] << " + ";
            sum += ((int64_t)rpp[i] * data[i]);
        }
        //std::cout << "\n";
        //std::cout << "sum prev. " << sum << " % "<<P <<" : " << (sum % P)<<"\n";
        sum %= P;
        //if(sum < 0 )
        //  sum += P;
        //std::cout << "t2 -->"<< sum << "\n";

        return sum;
    }


public:
	/**
	Constructor that enables the use of pointers to this class.
	*/
	HashIndex(){		
	}
	
	void initialize_structure(int _tableSize, size_t _k, int _seed){
		tableSize = _tableSize;
		k = _k;
		seed = _seed;
        // allocate memory for parameters.
        rp = new int32_t[k];
        rpp = new int32_t[k];
        // generate the random values
        //std::mt19937 unigen (seed);
        srand(seed);
        for(int i=0; i<k; i++) {
            rp[i] = rand();//(uint32_t)unigen();
            rpp[i] = rand();//(uint32_t)unigen();
        }		
	}
	
    HashIndex(int _tableSize, size_t _k, int _seed)
    {
		initialize_structure(_tableSize, _k, seed);
    }

	/*
	Compute both hash functions over the data attribute of
	object o.
	Note: This function does not store the o.data vector, it only
	computes the hash values and then store the id attribute of o.
	*/
    void index(const OBJECT & o)
    {
        // compute t1 and t2
        uint32_t t1v, t2v;
        t1v = t1b(o.data);
        t2v = t2b(o.data);
        table[t1v][t2v].push_back(o.id);
/*
        if(table.find(t1v) == table.end() )
            table[t1v][t2v].push_back(o.id);
        else if(table[t1v].find(t2v) == table[t1v].end() )
            table[t1v][t2v].push_back(o.id);
*/

    }

    std::vector<IDTYPE> * get_bucket(const OBJECT & o)
    {
        uint32_t t1v, t2v;
        t1v = t1b(o.data);
        t2v = t2b(o.data);

        if(table.find(t1v) == table.end() )
            return NULL;
        if(table[t1v].find(t2v) == table[t1v].end() )
            return NULL;
        return &table[t1v][t2v];
    }

	void traverse_table(std::vector<std::vector<IDTYPE>> & out){
		// for each key in the 1st level table
		for(typename T::iterator it=table.begin(); it != table.end(); ++it){
			// for each key in the 2nd level table (corresponding to key it->first)
			for(typename T2::iterator it2=it->second.begin(); it2 != it->second.end(); ++it2){
				//it2->second is a vector<IDTYPE>
				out.push_back(it2->second);
			}
		}
	}
	
    ~HashIndex()
    {
        // destroy memory allocated for parameters
        delete[] rp;
        delete[] rpp;
    }
};


class Point
{
public:
    int id;
    uint32_t * data;
    size_t dim;

    Point(int _id, int _dim):id(_id), dim(_dim)
    {
        data = new uint32_t[dim];
    }

    ~Point()
    {
        delete[] data;
        //std::cout << "Deallocating point\n";
    }
};

class BinarySignature
{
    typedef char T;
public:
    int id;
    T * data;
    size_t dim;

    BinarySignature(int _id, int _dim):id(_id), dim(_dim)
    {
        data = new T[dim + 1];
        data[dim] = '\0';
    }

    ~BinarySignature()
    {
        delete[] data;
        //std::cout << "Deallocating point\n";
    }
};

#endif /* PERFECTLSH_LIB_HPP_ */
