
//_____________________________________________________________________________
//    Class for UPC data
//    Author: Jaroslav Adam
//
//    contains parameters of primary vertex relevant for UPC analysis
//_____________________________________________________________________________

//local headers
#include "StUPCVertex.h"
#include "StUPCEvent.h"

ClassImp(StUPCVertex);

//_____________________________________________________________________________
StUPCVertex::StUPCVertex(): TObject(),
  mPosX(0), mPosY(0), mPosZ(0), mErrX(0), mErrY(0), mErrZ(0),
  mNPrimaryTracks(0), mNTracksUsed(0), mId(0),
  mIdTruth(0), mQaTruth(0), mIdParent(0), mMCVertexIndex(-1), mEvent(0)
{
  //default constructor

}//StUPCVertex

void StUPCVertex::Clear(Option_t *) {
  mPosX = mPosY = mPosZ = mErrX = mErrY = mErrZ = 0;
  mNPrimaryTracks = mNTracksUsed = mId = 0;
  mIdTruth = mQaTruth = mIdParent = 0;
  mMCVertexIndex = -1;
  mEvent = 0;
}

StUPCMCVertex *StUPCVertex::getMCVertex() const {
  return mEvent ? mEvent->getMCVertex(mMCVertexIndex) : 0;
}


