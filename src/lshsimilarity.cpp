// lshsimilarity.cpp
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
#include <ctime>

#include <functional>

#include <unordered_map>

#define ESIMCMP 0 // If set estimated similarity matrix will be computed.
#define SIMESTMOD 1 // Define the estimation mode (1 or 2)


typedef std::vector<int> bucket;
typedef std::unordered_map<uint32_t, bucket> table;

void write_simmatrix(double **data, int ninstances, const char * outfname)
{
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
    if(argc < 7) {
        std::cout << "Creates a [1-delta][1-s]-Neighbor structure from input data, i.e. for a query point,\n";
        std::cout << "returns all their neighbors at distance [1-s] with probability (1-delta).\n\n";
        std::cout << "Usage: "<<argv[0]<<" <input_dat> <input_labels> <fn-rate(delta)> <s-threshold([1-s]-neighbors)> <rows-per-band> <#bands>\n";
        return -1;
    }


    Setdata data;
    read_csv_set_data(argv[1], argv[2], data); // data file path, label file path and Setdata structure.
    std::cout << "Data file read.\n";
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
    
    std::vector<table> T(nbands);

    uint32_t *sgns = new uint32_t [gen.nbands];
    // apply minwise to each point in data.
    for (int dataid = 0; dataid < data.instances(); dataid++) {
        gen.get_signature(data.data[dataid], sgns, 10*data.instances()); // unbounded hash value
        /*
        for(int si=0; si<nbands; si++)
            std::cout << sgns[si] << " ";
        std::cout << "\n";
        */
        for(int j=0; j<nbands; j++) {
            // index point i into T[j] with key equal to sgns[j]
            T[j][sgns[j]].push_back( dataid );
        }
    }
    std::cout << "Indexing stage done.\n";


    // Similarity estimation stage.
    //

    double ** simhat = new double* [data.instances()];
    for(int i=0; i<data.instances(); i++) {
        simhat[i] = new double[data.instances()];
        for(int j=0; j<data.instances(); j++)
            simhat[i][j] = 0.0;
    }

    for (int i = 0; i < nbands; i++) {
        // traverse buckets in table T[i]
        auto it = T[i].begin();
        for(; it!=T[i].end(); ++it){
            // it->first (bucket-id) 
            // it->second (vector with documents)
            bucket& curr_bkt = it->second;

            size_t ndocs_inbkt = curr_bkt.size();

            for(auto dnr=0; dnr < ndocs_inbkt; dnr++){
                //assert(bucket[i] < data.instances() );
                //assert( curr_bkt[dnr] < data.instances() );
                for(auto dnraux=dnr+1; dnraux < ndocs_inbkt; dnraux++){
                    //assert( curr_bkt[dnraux] < data.instances() );
                    // increase similarity between documents curr_bkt[dnr] and curr_bkt[dnraux]
                    simhat[curr_bkt[dnr]][curr_bkt[dnraux]] += 1;
                    simhat[curr_bkt[dnraux]][curr_bkt[dnr]]= simhat[curr_bkt[dnr]][curr_bkt[dnraux]];
                }
            }
        }
    }


    // correcting the estimated similarity.
    for(int i=0; i<data.instances() - 1; i++) {
        simhat[i][i] = 1.0;
        for(int j=i+1; j<data.instances(); j++) {
            if(simhat[i][j] == 0 )
                continue;
            //std::cout << "Converting to "<< "pow("<<simhat[i][j]/nbands <<", 1.0/("<<simhat[i][j]*r<<") )" << "\n";
            simhat[i][j] = pow(simhat[i][j]/nbands , 1.0/(simhat[i][j]*r) );// MOD - 2

#if SIMESTMOD == 1
			simhat[i][j] = pow( 1 - pow( 1-simhat[i][j]/nbands , 1.0/nbands) , 1.0/r); // --> MOD 1
#endif
#if SIMESTMOD == 2
            simhat[i][j] = pow(simhat[i][j]/nbands , 1.0/(simhat[i][j]*r) );// --> MOD 2
#endif
            simhat[j][i] = simhat[i][j];
        }
    }
    std::cout << "Similarity estimation stage done.\n";

    char newname[250];
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

    delete[] sgns;
    return 0;
}
