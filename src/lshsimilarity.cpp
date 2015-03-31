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

#include <functional>

#define ESIMCMP 0 // If set estimated similarity matrix will be computed.
#define SIMESTMOD 1 // Define the estimation mode (1 or 2)


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

    HashIndex<int, IPoint> * tables = new HashIndex<int, IPoint>[nbands];
    for(int i=0; i<nbands; i++)
        tables[i].initialize_structure(data.instances(), r, 1); // table size, input dim and seed.

    IPoint inpt(-1);

    uint32_t *sgns = new uint32_t [gen.nbands];
    // apply minwise to each point in data.
    for (int i = 0; i < data.instances(); i++) {
        inpt.id = i; // Assign the id for the current point.
        // create a Point object for the signature.
        gen.get_signature(data.data[i], sgns); // unbounded hash value

        //std::cout << "Signature for instance "<<(i+1)<<"\n";
        for(int j=0; j<nbands; j++) {
            inpt.data = sgns[i];
            // Index it.
            tables[j].index(inpt);
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

    //std::cout<< "Displaying tables and buckets (ommiting buckets having less than 1 doc.)\n";
    for(int t=0; t<nbands; t++) {
        std::vector< std::vector<int>> out;
        tables[t].traverse_table(out);	// get all buckets from table[i]
        //traverse vector
        //std::cout << "Table "<<i<<"("<< out.size()<<" bcks.)"<<"\n";
        for(size_t bucketid=0; bucketid < out.size(); bucketid++) {
            // for each bucket retrieved
            std::vector<int>& bucket = out[bucketid];
            int bucket_size = bucket.size();

            if( bucket_size < 2 )
                continue;

            //std::cout << "Bucket nr "<< bucketid+1 <<" ("<< bucket_size<<" docs.)"<<"\n";
            for(int i=0; i<bucket_size; i++) {
                //std::cout << bucket[i] << " ";
                //assert(bucket[i] < data.instances() );
                for(int j=i+1; j<bucket_size; j++) {
                    //std::cout << "Before "<< simhat[bucket[i]][bucket[j]]<<" and after ";
                    simhat[bucket[i]][bucket[j]] += 1.0;
                    //std::cout << simhat[bucket[i]][bucket[j]] << "\n";
                    simhat[bucket[j]][bucket[i]] = simhat[bucket[i]][bucket[j]]; // symmetric condition.

                }
            }
            //std::cout <<"\n";
        }
    }

    /*
    // correcting the estimated similarity.
    for(int i=0; i<data.instances() - 1; i++) {
        simhat[i][i] = 1.0;
        for(int j=i+1; j<data.instances(); j++) {
            if(simhat[i][j] == 0 )
                continue;
            //std::cout << "Converting to "<< "pow("<<simhat[i][j]/nbands <<", 1.0/("<<simhat[i][j]*r<<") )" << "\n";
            simhat[i][j] = pow(simhat[i][j]/nbands , 1.0/(simhat[i][j]*r) );// MOD - 2
            simhat[j][i] = simhat[i][j];
        }
    }
    */
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
    return 0;
}
