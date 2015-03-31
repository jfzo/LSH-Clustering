// perfectlshmwise_custom.cpp
//
#include <iostream>
#include <vector>
#include <set>
#include "perfectlsh_lib.hpp"
#include "minwise_lib.hpp"
#include <cmath>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <vector>

#include <functional>

#define ESIMCMP 1 // If set estimated similarity matrix will be computed.
#define SIMESTMOD 2 // Define the estimation mode (1 or 2)

void write_signatures(unsigned **signatures, int ninstances, int nbands, const char * outfname){
	FILE *spMatFp;	
	spMatFp = fopen(outfname, "w");
	//header
	fprintf(spMatFp, "%d %d\n", ninstances, nbands);	
	for(int i=0; i<ninstances; i++){
		for(int h=0; h<nbands; h++){
			//write signatures[i][h]
			fprintf(spMatFp, "%u ", signatures[i][h]);
		}
		fprintf(spMatFp, "\n");
	}
	fclose(spMatFp);
}


void write_simmatrix(double **data, int ninstances, const char * outfname){
	FILE *spMatFp;
	//std::cout << "Similarity matrix will be written to:"<< newname << "\n";
	spMatFp = fopen(outfname, "w");
	//header
	fprintf(spMatFp, "%d\n", ninstances);
	for (int lnum = 0; lnum <ninstances; lnum++) {
		for (int cnum = 0; cnum < ninstances; cnum++) {
			fprintf(spMatFp, "%.6f ", data[lnum][cnum] );
		}
		fprintf(spMatFp, "\n");
	}
	fclose(spMatFp);		
}


int main(int argc, char* argv[])
{
	if(argc < 7){
		std::cout << "Creates a [1-delta][1-s]-Neighbor structure from input data, i.e. for a query point,\n";
		std::cout << "returns all their neighbors at distance [1-s] with probability (1-delta).\n\n";
		std::cout << "Usage: "<<argv[0]<<" <input_dat> <input_labels> <fn-rate(delta)> <s-threshold([1-s]-neighbors)> <rows-per-band> <#bands>\n";
		return -1;
	}

    Setdata data;
    read_csv_set_data(argv[1], argv[2], data); // data file path, label file path and Setdata structure.
    std::cout << argv[1] << " labels: "<< argv[2] << "\n";
    std::cout << "NDataPoints:"<< data.instances() << " Dim:"<< data.dims() <<"\n";

#if ESIMCMP
    std::cout << "Approximate similarity matrix will be computed and stored.\n";
#endif

    double delta = atof(argv[3]); // error rate (rate of false positive)
    double s = atof(argv[4]); // similarity threshold / (1-s)-Neighbors will be returned with Pr. (1-delta)
    int r = atoi(argv[5]); // nr of rows per band
    int nbands = ceil(log10(delta)/log10(1-pow(s, r)));
    std::cout << "Recommended NBands : " << nbands<<"\n";
    
    nbands = atoi(argv[6]);
    
    std::cout << "NBands set to " << nbands<<" and band size manually set to " << r<<"\n";
    
    MinwiseGenerator gen(nbands, r, data.dims());
	//RandomHyperplaneGenerator gen(nbands, r, data.dims());

    /*
	HashIndex<int, BinarySignature> * tables = new HashIndex<int, BinarySignature>[nbands];
	for(int i=0; i<nbands; i++)
		tables[i].initialize_structure(data.instances(), r, 1);
	*/
    
	//unsigned char * buffer = new unsigned char[gen.nbands * r];
    //uint32_t * buffer = new uint32_t[gen.nbands];
	//BinarySignature inpt(-1, r); // new Point to hold every band of the signature
    
	// storing signatures
	uint32_t **sgns = new uint32_t * [data.instances()];
	for(int i=0; i<data.instances(); i++)
		sgns[i] = new uint32_t[gen.nbands];
    
    
	// apply minwise to each point in data.	
    for (int i = 0; i < data.instances(); i++) {
		//inpt.id = i;
		// create a Point object for the signature.
		
        gen.get_signature(data.data[i], sgns[i]); // unbounded hash value 
		//gen.get_signature(data.data[i], sgns[i], data.instances() * 10); // max hash value set to data.instances() * 2
        
        //std::cout << "Signature for instance "<<(i+1)<<"\n";
        
        
		//for(int j=0; j<nbands; j++){
        //    printf("BAND %d : %u\n", (j+1), sgns[i][j] % (data.instances()*2) );
            // write it
            /*
			// copy each of the gen.nbands r-dimensional signatures into a point.
			for(int h=0; h<r; h++){
				inpt.data[h] = buffer[j*r + h];
			}
            uint32_t hband = band_hash(std::string(inpt.data)) % 4294967291; // % 4294967291
            printf("BAND %d : %s (%d)\n", (j+1), inpt.data, hband);
            */
			/*
            // Index it.
			tables[j].index(inpt);
            */
		//}
	}

	char newname[250];
	
	// writting signatures
	strcpy(newname, argv[1]);
	strcat(newname, "_");
	strcat(newname, argv[5]);//r
	strcat(newname, "-");
	strcat(newname, argv[6]);//nbands
	strcat(newname, "_mw.signatures");	
	std::cout << "Signature matrix will be written to:"<< newname << "\n";
	write_signatures(sgns, data.instances(), nbands, newname);
    

#if ESIMCMP

#if SIMESTMOD == 1
    std::cout << "Estimating similarity in mode 1\n";
#endif
#if SIMESTMOD == 2
    std::cout << "Estimating similarity in mode 2\n";
#endif

	// (estimated) similarity matrix
	double ** simhat = new double* [data.instances()];
	for(int i=0; i<data.instances(); i++){
		simhat[i] = new double[data.instances()];
		for(int j=0; j<data.instances(); j++)
			simhat[i][j] = 0.0;
	}
	
	for(int i=0; i < data.instances()-1; i++){
		simhat[i][i] = 1.0;
		for(int j=i+1; j < data.instances(); j++){
			// comute nr. of matches between the two signatures.
			double matches = 0.0;
			for(int h=0; h<nbands; h++){
				if(sgns[i][h] == sgns[j][h])
					matches += 1.0;
			}
			//std::cout << "angle estim. " << 180-pow(matches / nbands, 1.0/r)*180 << "\n";
			double est_p1 = matches / nbands;

#if SIMESTMOD == 1
			simhat[i][j] = pow( 1 - pow( 1-est_p1 , 1.0/nbands) , 1.0/r); // --> MOD 1
#endif
#if SIMESTMOD == 2
			simhat[i][j] = pow(est_p1 , 1.0/(matches*r) ); // --> MOD 2
#endif
			simhat[j][i] = simhat[i][j];
		}
	}

	strcpy(newname, argv[1]);
	strcat(newname, "_");
	strcat(newname, argv[5]);//r
	strcat(newname, "-");
	strcat(newname, argv[6]);//nbands
	strcat(newname, "_mw.simestimate");
	std::cout << "Estimated similarity matrix will be written to:"<< newname << "\n";	
	write_simmatrix(simhat, data.instances(), newname);
		

	
	for(int i=0; i<data.instances(); i++)
		delete[] simhat[i];
	delete[] simhat;
#endif	

	for(int i=0; i<data.instances(); i++)
		delete[] sgns[i];
	delete[] sgns;
    /*
	delete[] tables;
    */

    
    	return 0;
}
