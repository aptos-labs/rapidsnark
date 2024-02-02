#include <stdexcept>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>
#include <fstream>
#include <chrono>


#include "fullprover.hpp"
#include "fr.hpp"

#include "logging.hpp"
#include "nlohmann/json.hpp"
#include "wtns_utils.hpp"

#include <mutex>
#include "alt_bn128.hpp"
#include "groth16.hpp"
#include "binfile_utils.hpp"
#include "zkey_utils.hpp"


class FullProverImpl {
    bool unsupported_zkey_curve;

    std::string circuit;
    std::string witnessBinaryPath;

    std::unique_ptr<Groth16::Prover<AltBn128::Engine>> prover;
    std::unique_ptr<ZKeyUtils::Header> zkHeader;
    std::unique_ptr<BinFileUtils::BinFile> zKey;

    mpz_t altBbn128r;

  public:
    FullProverImpl(const char *_zkeyFileName, const char *_witnessBinaryPath);
    ~FullProverImpl();
    ProverResponse prove(const char *input);
};



FullProver::FullProver(const char *_zkeyFileName, const char *_witnessBinaryPath) {
  std::cout << "in FullProver constructor" << std::endl;
  try {
    std::cout << "try" << std::endl;
    impl = new FullProverImpl(_zkeyFileName, _witnessBinaryPath);
    state = FullProverState::OK;
  } catch (std::invalid_argument e) {
    std::cout << "caught" << std::endl;
    state = FullProverState::UNSUPPORTED_ZKEY_CURVE;
    impl = 0;
  } catch (std::system_error e) {
    std::cout << "caught 2" << std::endl;
    state = FullProverState::ZKEY_FILE_LOAD_ERROR;
    impl = 0;
  }

}

FullProver::~FullProver() {
  delete impl;
}

ProverResponse FullProver::prove(const char *input) {
  std::cout << "in FullProver::prove" << std::endl;
  if (state != FullProverState::OK) {
    return ProverResponse(ProverError::PROVER_NOT_READY);
  } else {
    return impl->prove(input);
  }
}









// FULLPROVERIMPL


std::string getfilename(std::string path)
{
    path = path.substr(path.find_last_of("/\\") + 1);
    size_t dot_i = path.find_last_of('.');
    return path.substr(0, dot_i);
}

FullProverImpl::FullProverImpl(const char *_zkeyFileName, const char *_witnessBinaryPath) : witnessBinaryPath(_witnessBinaryPath) {
  std::cout << "in FullProverImpl constructor" << std::endl;
    mpz_init(altBbn128r);
  std::cout << "in FullProverImpl constructor" << std::endl;
    mpz_set_str(altBbn128r, "21888242871839275222246405745257275088548364400416034343698204186575808495617", 10);
  std::cout << "in FullProverImpl constructor" << std::endl;

  std::cout << "in FullProverImpl constructor" << std::endl;
    circuit = getfilename(_zkeyFileName);
  std::cout << "in FullProverImpl constructor" << std::endl;
    zKey = BinFileUtils::openExisting(_zkeyFileName, "zkey", 1);
  std::cout << "in FullProverImpl constructor" << std::endl;
    zkHeader = ZKeyUtils::loadHeader(zKey.get());
  std::cout << "in FullProverImpl constructor" << std::endl;

  std::cout << "in FullProverImpl constructor" << std::endl;
    std::string proofStr;
    if (mpz_cmp(zkHeader->rPrime, altBbn128r) != 0) {
      unsupported_zkey_curve = true;
        throw std::invalid_argument( "zkey curve not supported" );
    }
    
  std::cout << "in FullProverImpl constructor" << std::endl;
    std::ostringstream ss1;
    ss1 << "circuit: " << circuit;
    LOG_DEBUG(ss1);

  std::cout << "in FullProverImpl constructor" << std::endl;
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
  std::cout << "in FullProverImpl constructor" << std::endl;
}

FullProverImpl::~FullProverImpl() {
    mpz_clear(altBbn128r);
}

ProverResponse::ProverResponse(ProverError _error) :
  type(ProverResponseType::ERROR), raw_json(""), error(_error), metrics(ProverResponseMetrics()) {}

ProverResponse::ProverResponse(const char *_raw_json, ProverResponseMetrics _metrics) :
  type(ProverResponseType::SUCCESS), raw_json(_raw_json), error(ProverError::NONE), metrics(_metrics) {}

ProverResponse FullProverImpl::prove(const char *input) {
  std::cout << "starting prove" << std::endl;
    LOG_TRACE("FullProverImpl::prove begin");
    LOG_DEBUG(input);
    
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

    LOG_TRACE("FullProverImpl::prove end");

    ProverResponseMetrics metrics;
    metrics.prover_time = prover_duration.count();
    metrics.witness_generation_time = witness_generation_duration.count();
    
    const char *proof_raw = strdup(proof.dump().c_str());

    return ProverResponse(proof_raw, metrics);
}

