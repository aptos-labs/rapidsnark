#include <stdexcept>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <sstream>


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


std::string getFormattedTimestamp() {
    // Get the current time point
    auto now = std::chrono::system_clock::now();
    
    // Convert the time point to a time_t object
    std::time_t now_time_t = std::chrono::system_clock::to_time_t(now);

    // Convert to std::tm for formatting
    std::tm* now_tm = std::gmtime(&now_time_t);

    // Extract milliseconds from the time point
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch()) % 1000;

    // Create a stringstream to format the timestamp
    std::stringstream ss;
    ss << std::put_time(now_tm, "%Y-%m-%dT%H:%M:%S") << '.' << std::setfill('0') << std::setw(3) << ms.count();
    ss << "Z";

    // Return the formatted timestamp as a string
    return ss.str();
}

void log(string level,string msg) {

  std::cout << "{\"timestamp\":\""
            << getFormattedTimestamp()
            << "\",\"level\":\""
            << level
            << "\",\"message\":\""
            << msg 
            << "\",\"target\":\"prover_service::rapidsnark\"}"
            << std::endl;
}

void log_info(string msg) { log("INFO", msg); }
void log_debug(string msg) { log("DEBUG", msg); }
void log_error(string msg) { log("ERROR", msg); }


class FullProverImpl {
    bool unsupported_zkey_curve;

    std::string circuit;

    std::unique_ptr<Groth16::Prover<AltBn128::Engine>> prover;
    std::unique_ptr<ZKeyUtils::Header> zkHeader;
    std::unique_ptr<BinFileUtils::BinFile> zKey;

    mpz_t altBbn128r;

  public:
    FullProverImpl(const char *_zkeyFileName);
    ~FullProverImpl();
    ProverResponse prove(const char *input) const;
};




FullProver::FullProver(const char *_zkeyFileName) {
  std::cout << "in FullProver constructor" << std::endl;
  try {
    std::cout << "try" << std::endl;
    impl = new FullProverImpl(_zkeyFileName);
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
  std::cout << "in FullProver destructor" << std::endl;
  delete impl;
}

ProverResponse FullProver::prove(const char *input) const {
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

FullProverImpl::FullProverImpl(const char *_zkeyFileName) {
  std::cout << "in FullProverImpl constructor" << std::endl;
    mpz_init(altBbn128r);
    mpz_set_str(altBbn128r, "21888242871839275222246405745257275088548364400416034343698204186575808495617", 10);

    circuit = getfilename(_zkeyFileName);
    zKey = BinFileUtils::openExisting(_zkeyFileName, "zkey", 1);
    zkHeader = ZKeyUtils::loadHeader(zKey.get());

    std::string proofStr;
    if (mpz_cmp(zkHeader->rPrime, altBbn128r) != 0) {
      unsupported_zkey_curve = true;
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

FullProverImpl::~FullProverImpl() {
    mpz_clear(altBbn128r);
}

ProverResponse::ProverResponse(ProverError _error) :
  type(ProverResponseType::ERROR), raw_json(""), error(_error), metrics(ProverResponseMetrics()) {}

ProverResponse::ProverResponse(const char *_raw_json, ProverResponseMetrics _metrics) :
  type(ProverResponseType::SUCCESS), raw_json(_raw_json), error(ProverError::NONE), metrics(_metrics) {}

ProverResponse FullProverImpl::prove(const char *witness_file_path) const {
  log_info("FullProverImpl::prove begin");
  log_debug(std::string(witness_file_path));



    std::string witnessFile(witness_file_path);
    
    // Load witness
    auto wtns = BinFileUtils::openExisting(witnessFile, "wtns", 2);
    auto wtnsHeader = WtnsUtils::loadHeader(wtns.get());
    log_info("Loaded witness file");
            
    if (mpz_cmp(wtnsHeader->prime, altBbn128r) != 0) {
        log_error("The generated witness file uses a different curve than bn128, which is currently the only supported curve.");
        return ProverResponse(ProverError::WITNESS_GENERATION_INVALID_CURVE);
    }

    AltBn128::FrElement *wtnsData = (AltBn128::FrElement *)wtns->getSectionData(2);

    auto start = std::chrono::high_resolution_clock::now();
    json proof = prover->prove(wtnsData)->toJson();
    auto end = std::chrono::high_resolution_clock::now();
    auto prover_duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    log_info("finished proof computation");

    {
      std::stringstream ss;
      ss << "Time taken for Groth16 prover: " << prover_duration.count() << " milliseconds";
      std::cout << "Time taken for Groth16 prover: " << prover_duration.count() << " milliseconds" << std::endl;
      log_info(ss.str().data());
    }


    log_info("constructing metrics struct");
    ProverResponseMetrics metrics;
    metrics.prover_time = prover_duration.count();

    const char *proof_raw = strdup(proof.dump().c_str());
    if (proof_raw == nullptr) {
      throw std::runtime_error("Failed to allocate memory for proof_raw");
    }

    log_info("FullProverImpl::prove end");
    return ProverResponse(proof_raw, metrics);
}

