#include <gmp.h>
#include <memory>
#include <string>
#include <cstring>
#include <stdexcept>
#include <alt_bn128.hpp>
#include <nlohmann/json.hpp>

#include "prover.h"
#include "zkey_utils.hpp"
#include "wtns_utils.hpp"
#include "groth16.hpp"
#include "binfile_utils.hpp"

using json = nlohmann::json;

static size_t ProofBufferMinSize()
{
    return 810;
}

static size_t PublicBufferMinSize(size_t count)
{
    return count * 82 + 4;
}

/**
 * Verifies the given prime numbers.
 *
 * @param zkey_prime A pointer to the first prime number.
 * @param wtns_prime A pointer to the second prime number.
 */
static void VerifyPrimes(mpz_srcptr zkey_prime, mpz_srcptr wtns_prime)
{
    mpz_t altBbn128r;

    mpz_init(altBbn128r);
    mpz_set_str(altBbn128r, "21888242871839275222246405745257275088548364400416034343698204186575808495617", 10);

    if (mpz_cmp(zkey_prime, altBbn128r) != 0) {
        throw std::invalid_argument( "zkey curve not supported" );
    }

    if (mpz_cmp(wtns_prime, altBbn128r) != 0) {
        throw std::invalid_argument( "different wtns curve" );
    }

    mpz_clear(altBbn128r);
}

/**
 * Builds a JSON string representation of the public data.
 * 
 * @param wtnsData Pointer to an array of AltBn128::FrElement representing the public data.
 * @param nPublic The number of elements in the public data array.
 * @return A JSON string representation of the public data.
 */
std::string BuildPublicString(AltBn128::FrElement *wtnsData, size_t nPublic)
{
    json jsonPublic;
    AltBn128::FrElement aux;
    for (u_int32_t i=1; i<= nPublic; i++) {
        AltBn128::Fr.toMontgomery(aux, wtnsData[i]);
        jsonPublic.push_back(AltBn128::Fr.toString(aux));
    }

    return jsonPublic.dump();
}

/**
 * Calculates the size of the public buffer required for the given zkey buffer.
 *
 * @param zkey_buffer A pointer to the zkey buffer.
 * @param zkey_size The size of the zkey buffer.
 * @return The size of the public buffer required.
 */
unsigned long CalcPublicBufferSize(const void *zkey_buffer, unsigned long zkey_size) {
    try {
        BinFileUtils::BinFile zkey(zkey_buffer, zkey_size, "zkey", 1);
        auto zkeyHeader = ZKeyUtils::loadHeader(&zkey);
        return PublicBufferMinSize(zkeyHeader->nPublic);
    } catch (...) {
    }

    return 0;
}

/**
 * Proves a Groth16 proof given the necessary inputs.
 *
 * @param zkey_buffer       Pointer to the buffer containing the zkey data.
 * @param zkey_size         Size of the zkey buffer.
 * @param wtns_buffer       Pointer to the buffer containing the wtns data.
 * @param wtns_size         Size of the wtns buffer.
 * @param proof_buffer      Pointer to the buffer to store the generated proof.
 * @param proof_size        Pointer to the size of the proof buffer. On input, it should contain the available size of the buffer. On output, it will be updated with the actual size of the generated proof.
 * @param public_buffer     Pointer to the buffer to store the generated public data.
 * @param public_size       Pointer to the size of the public buffer. On input, it should contain the available size of the buffer. On output, it will be updated with the actual size of the generated public data.
 * @param error_msg         Pointer to the buffer to store any error message that occurs during the proving process.
 * @param error_msg_maxsize Maximum size of the error message buffer.
 *
 * @return                  Returns PROVER_OK if the proving process is successful. If there is an error, it returns an appropriate error code.
 */
int
groth16_prover(const void *zkey_buffer,   unsigned long  zkey_size,
               const void *wtns_buffer,   unsigned long  wtns_size,
               char       *proof_buffer,  unsigned long *proof_size,
               char       *public_buffer, unsigned long *public_size,
               char       *error_msg,     unsigned long  error_msg_maxsize)
{
    try {
        BinFileUtils::BinFile zkey(zkey_buffer, zkey_size, "zkey", 1);
        auto zkeyHeader = ZKeyUtils::loadHeader(&zkey);

        BinFileUtils::BinFile wtns(wtns_buffer, wtns_size, "wtns", 2);
        auto wtnsHeader = WtnsUtils::loadHeader(&wtns);

        if (zkeyHeader->nVars != wtnsHeader->nVars) {
            snprintf(error_msg, error_msg_maxsize,
                     "Invalid witness length. Circuit: %u, witness: %u",
                     zkeyHeader->nVars, wtnsHeader->nVars);
            return PROVER_INVALID_WITNESS_LENGTH;
        }

        size_t proofMinSize  = ProofBufferMinSize();
        size_t publicMinSize = PublicBufferMinSize(zkeyHeader->nPublic);

        if (*proof_size < proofMinSize || *public_size < publicMinSize) {

            *proof_size  = proofMinSize;
            *public_size = publicMinSize;

            return PROVER_ERROR_SHORT_BUFFER;
        }

        VerifyPrimes(zkeyHeader->rPrime, wtnsHeader->prime);

        auto prover = Groth16::makeProver<AltBn128::Engine>(
            zkeyHeader->nVars,
            zkeyHeader->nPublic,
            zkeyHeader->domainSize,
            zkeyHeader->nCoefs,
            zkeyHeader->vk_alpha1,
            zkeyHeader->vk_beta1,
            zkeyHeader->vk_beta2,
            zkeyHeader->vk_delta1,
            zkeyHeader->vk_delta2,
            zkey.getSectionData(4),    // Coefs
            zkey.getSectionData(5),    // pointsA
            zkey.getSectionData(6),    // pointsB1
            zkey.getSectionData(7),    // pointsB2
            zkey.getSectionData(8),    // pointsC
            zkey.getSectionData(9)     // pointsH1
        );
        AltBn128::FrElement *wtnsData = (AltBn128::FrElement *)wtns.getSectionData(2);
        auto proof = prover->prove(wtnsData);

        std::string stringProof = proof->toJson().dump();
        std::string stringPublic = BuildPublicString(wtnsData, zkeyHeader->nPublic);

        size_t stringProofSize  = stringProof.length();
        size_t stringPublicSize = stringPublic.length();

        if (*proof_size < stringProofSize || *public_size < stringPublicSize) {

            *proof_size  = stringProofSize;
            *public_size = stringPublicSize;

            return PROVER_ERROR_SHORT_BUFFER;
        }

        std::strncpy(proof_buffer, stringProof.data(), *proof_size);
        std::strncpy(public_buffer, stringPublic.data(), *public_size);

    } catch (std::exception& e) {

        if (error_msg) {
            strlcpy(error_msg, e.what(), error_msg_maxsize);
        }
        return PROVER_ERROR;

    } catch (std::exception *e) {

        if (error_msg) {
            strlcpy(error_msg, e->what(), error_msg_maxsize);
        }
        delete e;
        return PROVER_ERROR;

    } catch (...) {
        if (error_msg) {
            strlcpy(error_msg, "unknown error", error_msg_maxsize);
        }
        return PROVER_ERROR;
    }

    return PROVER_OK;
}
