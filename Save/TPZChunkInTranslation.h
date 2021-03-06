/* 
 * File:   TPZChunkInTranslator.h
 * Author: quinelato
 *
 * Created on September 18, 2017, 3:54 PM
 */

#ifndef TPZCHUNKINTRANSLATION_H
#define TPZCHUNKINTRANSLATION_H

#include "TPZContBufferedStream.h"


class TPZChunkInTranslation {
public:
    TPZChunkInTranslation(const int64_t &objId, const int &classId, TPZStream &stream, const size_t &chunkSize, std::map<std::string, uint64_t> &versionInfo);
    TPZChunkInTranslation(const TPZChunkInTranslation& orig);
    virtual ~TPZChunkInTranslation();
    int64_t GetObjId() const;
    int GetClassId() const;
private:
    void ReadFromStream(TPZStream &stream, const size_t nBytes);
private:
    TPZContBufferedStream mOldStream;
    TPZContBufferedStream mNewStream;
    
    std::map<std::string, uint64_t> mOldVersion;
    std::map<std::string, uint64_t> mNewVersion;
    
    int64_t mObjId;
    int mClassId;
    
    TPZManVector<int64_t, 2> mNewObjIds;
    
    friend TPZPersistenceManager;
};

#endif /* TPZCHUNKINTRANSLATION_H */

