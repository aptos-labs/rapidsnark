#ifndef FULLPROVER_H
#define FULLPROVER_H



class FullProverImpl;

enum ProverResponseType {
  SUCCESS,
  ERROR
};

enum ProverError {
  NONE,
  INVALID_INPUT,
  WITNESS_GENERATION_BINARY_PROBLEM,
  WITNESS_GENERATION_INVALID_CURVE
};

struct ProverResponseMetrics {
  int prover_time;
  int witness_generation_time;

};

struct ProverResponse {
  ProverResponseType type;
  const char *raw_json;
  ProverError error;
  ProverResponseMetrics metrics;

  public:
    ProverResponse(ProverError _error);
    ProverResponse(const char *_raw_json, ProverResponseMetrics _metrics);

};



class FullProver {

  FullProverImpl *impl;


public: 
    FullProver(const char *_zkeyFileName, const char *_witnessBinaryPath);
    ~FullProver();
    ProverResponse prove(const char *input);


};

#endif // FULLPROVER_H
