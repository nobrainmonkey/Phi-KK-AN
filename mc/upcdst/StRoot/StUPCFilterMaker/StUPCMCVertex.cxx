#include "StUPCMCVertex.h"
#include "StUPCMCEvent.h"
#include "StUPCEvent.h"
#include <algorithm>

ClassImp(StUPCMCVertex);

StUPCMCVertex::StUPCMCVertex(): TObject() { Clear(); }

void StUPCMCVertex::Clear(Option_t *) {
  mId = mNDaughters = mIntermediateFlag = mParentTrackId = 0;
  mPosX = mPosY = mPosZ = mTime = 0;
  mParentTrackIndex = -1;
  mEvent = 0;
  mRecoVertexIndices.clear();
}

StUPCMCTrack *StUPCMCVertex::getParentTrack() const {
  return mEvent ? mEvent->getTrack(mParentTrackIndex) : 0;
}

void StUPCMCVertex::addRecoVertexIndex(Int_t index) {
  if(index >= 0 && std::find(mRecoVertexIndices.begin(),mRecoVertexIndices.end(),index)==mRecoVertexIndices.end())
    mRecoVertexIndices.push_back(index);
}

Int_t StUPCMCVertex::getRecoVertexIndex(Int_t match) const {
  return match >= 0 && match < getNumberOfRecoVertices() ? mRecoVertexIndices[match] : -1;
}

StUPCVertex *StUPCMCVertex::getRecoVertex(Int_t match) const {
  StUPCEvent *reco = mEvent ? mEvent->getRecoEvent() : 0;
  return reco ? reco->getVertex(getRecoVertexIndex(match)) : 0;
}
