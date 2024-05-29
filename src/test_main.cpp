#include <iostream>
#include <fstream>
#include <gmp.h>
#include <memory>
#include <stdexcept>
#include <nlohmann/json.hpp>
#include <chrono>


#include <alt_bn128.hpp>
#include "binfile_utils.hpp"
#include "zkey_utils.hpp"
#include "wtns_utils.hpp"
#include "groth16.hpp"

using json = nlohmann::json;


int main(int argc, char **argv)
{
    if (argc != 5) {
        std::cerr << "Invalid number of parameters:\n";
        std::cerr << "Usage: prover <circuit.zkey> <witness.wtns> <proof.json> <public.json>\n";
        return EXIT_FAILURE;
    }

    mpz_t altBbn128r;

    mpz_init(altBbn128r);
    mpz_set_str(altBbn128r, "21888242871839275222246405745257275088548364400416034343698204186575808495617", 10);

    try {
        std::string zkeyFilename = argv[1];
        std::string wtnsFilename = argv[2];
        std::string proofFilename = argv[3];
        std::string publicFilename = argv[4];

        auto zkey = BinFileUtils::openExisting(zkeyFilename, "zkey", 1);
        auto zkeyHeader = ZKeyUtils::loadHeader(zkey.get());

        std::string proofStr;
        if (mpz_cmp(zkeyHeader->rPrime, altBbn128r) != 0) {
            throw std::invalid_argument( "zkey curve not supported" );
        }

        auto wtns = BinFileUtils::openExisting(wtnsFilename, "wtns", 2);
        auto wtnsHeader = WtnsUtils::loadHeader(wtns.get());

        if (mpz_cmp(wtnsHeader->prime, altBbn128r) != 0) {
            throw std::invalid_argument( "different wtns curve" );
        }


        auto pointsA = zkey->getSectionData(5);
        AltBn128::FrElement *wtnsData = (AltBn128::FrElement *)wtns->getSectionData(2);

        AltBn128::Engine E;

        typename AltBn128::Engine::G1Point pi_a;

        auto start = std::chrono::high_resolution_clock::now();

        #pragma omp parallel for 
        for (int round = 0; round < 24; round++) {
          E.g1.multiMulByScalar(pi_a, (AltBn128::Engine::G1PointAffine *)pointsA, (uint8_t*)wtnsData, sizeof((AltBn128::Engine::FrElement *)wtnsData), 40000);
        }
        auto end = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double, std::milli> duration = end - start;

        // Display the duration in milliseconds
        std::cout << "Function execution time: " << duration.count() << " ms" << std::endl;



    } catch (std::exception* e) {
        mpz_clear(altBbn128r);
        std::cerr << e->what() << '\n';
        return EXIT_FAILURE;
    } catch (std::exception& e) {
        mpz_clear(altBbn128r);
        std::cerr << e.what() << '\n';
        return EXIT_FAILURE;
    }

    mpz_clear(altBbn128r);
    exit(EXIT_SUCCESS);
}
