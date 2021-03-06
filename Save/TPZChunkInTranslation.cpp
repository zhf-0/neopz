/* 
 * File:   TPZChunkInTranslator.cpp
 * Author: quinelato
 * 
 * Created on September 18, 2017, 3:54 PM
 */

#include "TPZChunkInTranslation.h"

TPZChunkInTranslation::TPZChunkInTranslation(const int64_t &objId, const int &classId, TPZStream &stream, const size_t &chunkSize, std::map<std::string, uint64_t> &versionInfo) :
mObjId(objId),
mClassId(classId),
mNewVersion(versionInfo) {
    this->ReadFromStream(stream, chunkSize);
}

TPZChunkInTranslation::TPZChunkInTranslation(const TPZChunkInTranslation& orig) {
}

TPZChunkInTranslation::~TPZChunkInTranslation() {
}

void TPZChunkInTranslation::ReadFromStream(TPZStream &stream, const size_t nBytes) {
    char *temp = new char[nBytes];
    stream.Read(temp, nBytes);
    mNewStream.Write(temp, nBytes);
	delete[] temp;
}

int64_t TPZChunkInTranslation::GetObjId() const {
    return mObjId;
}

int TPZChunkInTranslation::GetClassId() const{
    return mClassId;
}

