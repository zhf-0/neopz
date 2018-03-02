//
// Created by Francisco Teixeira Orlandini on 2/8/18.
//

#ifndef PZ_SPZMODALANALYSISDATAREADER_H
#define PZ_SPZMODALANALYSISDATAREADER_H

#include <slepceps.h>
#include <pzreal.h>
#include <pzvec.h>
#include "parameter_handler.h"

class SPZModalAnalysisData;

struct SPZModalAnalysisDataReader{
private:
    ParameterHandler &prm;
    std::string path;
public:
    const std::string &GetPath() const;

private:
    SPZModalAnalysisDataReader();
    inline bool FileExists(const std::string &) const;
    void DeclareParameters();
    void PrintUsageMessage();
    void ParseCommandLine(const int argc, char *const *argv);
    void ReadComplexVector(const int &nEntries,const std::string &rName, const std::string &iName, const std::string &condName,TPZVec<STATE> &dest);
public:
    SPZModalAnalysisDataReader(ParameterHandler &prm, const int  argc, char *const *argv);
  void ReadParameters(SPZModalAnalysisData &data);
};

#endif //PZ_SPZMODALANALYSISDATAREADER_H
