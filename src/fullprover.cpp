#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>
#include <fstream>
#include <chrono>


#include "fullprover.hpp"
#include "fr.hpp"

#include "logging.hpp"
#include "wtns_utils.hpp"


std::string getfilename(std::string path)
{
    path = path.substr(path.find_last_of("/\\") + 1);
    size_t dot_i = path.find_last_of('.');
    return path.substr(0, dot_i);
}

FullProver::FullProver(const char *_zkeyFileName, const char *_witnessBinaryPath) : witnessBinaryPath(_witnessBinaryPath) {
    mpz_init(altBbn128r);
    mpz_set_str(altBbn128r, "21888242871839275222246405745257275088548364400416034343698204186575808495617", 10);

    circuit = getfilename(_zkeyFileName);
    zKey = BinFileUtils::openExisting(_zkeyFileName, "zkey", 1);
    zkHeader = ZKeyUtils::loadHeader(zKey.get());

    std::string proofStr;
    if (mpz_cmp(zkHeader->rPrime, altBbn128r) != 0) {
        throw std::invalid_argument( "zkey curve not supported" );
    }
    
    std::ostringstream ss1;
    ss1 << "circuit: " << circuit;
    LOG_DEBUG(ss1);

    prover = Groth16::makeProver<AltBn128::Engine>(
        zkHeader->nVars,
        zkHeader->nPublic,
        zkHeader->domainSize,
        zkHeader->nCoefs,
        zkHeader->vk_alpha1,
        zkHeader->vk_beta1,
        zkHeader->vk_beta2,
        zkHeader->vk_delta1,
        zkHeader->vk_delta2,
        zKey->getSectionData(4),    // Coefs
        zKey->getSectionData(5),    // pointsA
        zKey->getSectionData(6),    // pointsB1
        zKey->getSectionData(7),    // pointsB2
        zKey->getSectionData(8),    // pointsC
        zKey->getSectionData(9)     // pointsH1
    );
}

FullProver::~FullProver() {
    mpz_clear(altBbn128r);
}

ProverResponse::ProverResponse(ProverError _error) :
  type(ProverResponseType::ERROR), raw_json(""), error(_error), metrics(ProverResponseMetrics()) {}

ProverResponse::ProverResponse(const char *_raw_json, ProverResponseMetrics _metrics) :
  type(ProverResponseType::SUCCESS), raw_json(_raw_json), error(ProverError::NONE), metrics(_metrics) {}

ProverResponse FullProver::prove(const char *input) {
    LOG_TRACE("FullProver::prove begin");
    LOG_DEBUG(input);
    std::lock_guard<std::mutex> guard(mtx);
    
    // Generate witness
    json j = json::parse(input);

    std::string inputFile("/tmp/rapidsnark_input.json");
    std::string witnessFile("/tmp/rapidsnark_witness.wtns");
    
    std::ofstream file(inputFile);
    file << j;
    file.close();

    std::string command(witnessBinaryPath + "/" + circuit + " " + inputFile + " " + witnessFile);
    LOG_TRACE(command);
    std::array<char, 128> buffer;
    std::string result;

    auto start = std::chrono::high_resolution_clock::now();
    // std::cout << "Opening reading pipe" << std::endl;
    FILE* pipe = popen(command.c_str(), "r");
    if (!pipe)
    {
        return ProverResponse(ProverError::WITNESS_GENERATION_BINARY_PROBLEM);
    }
    while (fgets(buffer.data(), 128, pipe) != NULL) {
        // std::cout << "Reading..." << std::endl;
        result += buffer.data();
    }
    auto returnCode = pclose(pipe);
    auto end = std::chrono::high_resolution_clock::now();
    auto witness_generation_duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    if (returnCode != 0) {
        return ProverResponse(ProverError::WITNESS_GENERATION_BINARY_PROBLEM);
    }

    LOG_DEBUG(result);
    {
      std::stringstream ss;
      ss << "return code: " << returnCode;
      LOG_DEBUG(ss.str().data());
    }
    {
      std::stringstream ss;
      ss << "Time taken for witness generation: " << witness_generation_duration.count() << " milliseconds";
      std::cout << "Time taken for witness generation: " << witness_generation_duration.count() << " milliseconds" << std::endl;
      LOG_INFO(ss.str().data());
    }
    
    // Load witness
    auto wtns = BinFileUtils::openExisting(witnessFile, "wtns", 2);
    auto wtnsHeader = WtnsUtils::loadHeader(wtns.get());
            
    if (mpz_cmp(wtnsHeader->prime, altBbn128r) != 0) {
        LOG_ERROR("The generated witness file uses a different curve than bn128, which is currently the only supported curve.");
        return ProverResponse(ProverError::WITNESS_GENERATION_INVALID_CURVE);
    }

    AltBn128::FrElement *wtnsData = (AltBn128::FrElement *)wtns->getSectionData(2);

    start = std::chrono::high_resolution_clock::now();
    json proof = prover->prove(wtnsData)->toJson();
    end = std::chrono::high_resolution_clock::now();
    auto prover_duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    {
      std::stringstream ss;
      ss << "Time taken for Groth16 prover: " << prover_duration.count() << " milliseconds";
      std::cout << "Time taken for Groth16 prover: " << prover_duration.count() << " milliseconds" << std::endl;
      LOG_INFO(ss.str().data());
    }

    LOG_TRACE("FullProver::prove end");

    ProverResponseMetrics metrics;
    metrics.prover_time = prover_duration.count();
    metrics.witness_generation_time = witness_generation_duration.count();

    return ProverResponse(proof.dump().c_str(), metrics);
}

