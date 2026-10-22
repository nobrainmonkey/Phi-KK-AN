#include "StUPCMCEvent.h"
#include "StUPCEvent.h"
#include "StUPCMCTrack.h"
#include "StUPCMCVertex.h"
#include "TClonesArray.h"

ClassImp(StUPCMCEvent);

StUPCMCEvent::StUPCMCEvent(): TObject(),
  mRunNumber(0), mEventNumber(0), mSimulationType(kUnknown),
  mHasDSMData(kFALSE), mLastDSM0(0), mLastDSM1(0), mLastDSM3(0),
  mTracks(new TClonesArray("StUPCMCTrack")),
  mVertices(new TClonesArray("StUPCMCVertex")),
  mNTracks(0), mNVertices(0), mRecoEvent(0) {
  clearEventInfo();
  mTracks->SetOwner(kTRUE);
  mVertices->SetOwner(kTRUE);
}

StUPCMCEvent::~StUPCMCEvent() {
  delete mTracks;
  delete mVertices;
}

void StUPCMCEvent::Clear(Option_t *) {
  mTracks->Clear("C");
  mVertices->Clear("C");
  mNTracks = mNVertices = 0;
  mSimulationType = kUnknown;
  clearEventInfo();
}

StUPCMCTrack *StUPCMCEvent::addTrack() {
  StUPCMCTrack *track = dynamic_cast<StUPCMCTrack*>(mTracks->ConstructedAt(mNTracks++));
  if(track) track->setEvent(this);
  return track;
}

StUPCMCVertex *StUPCMCEvent::addVertex() {
  StUPCMCVertex *vertex = dynamic_cast<StUPCMCVertex*>(mVertices->ConstructedAt(mNVertices++));
  if(vertex) vertex->setEvent(this);
  return vertex;
}

Int_t StUPCMCEvent::getNumberOfTracks() const { return mTracks->GetEntriesFast(); }
Int_t StUPCMCEvent::getNumberOfVertices() const { return mVertices->GetEntriesFast(); }

StUPCMCTrack *StUPCMCEvent::getTrack(Int_t index) const {
  if(index < 0 || index >= getNumberOfTracks()) return 0;
  StUPCMCTrack *track = dynamic_cast<StUPCMCTrack*>(mTracks->At(index));
  if(track) track->setEvent(const_cast<StUPCMCEvent*>(this));
  return track;
}

StUPCMCVertex *StUPCMCEvent::getVertex(Int_t index) const {
  if(index < 0 || index >= getNumberOfVertices()) return 0;
  StUPCMCVertex *vertex = dynamic_cast<StUPCMCVertex*>(mVertices->At(index));
  if(vertex) vertex->setEvent(const_cast<StUPCMCEvent*>(this));
  return vertex;
}

void StUPCMCEvent::clearEventInfo() {
  mHasEventInfo = mHasDSMData = kFALSE;
  mNRecoTracks = mNRecoVertices = mNClusters = mNHits = 0;
  mTrgIDs.Set(0);
  mRunNumber = 0;
  mEventNumber = 0;
  mFillNum = 0;
  mbCrossId = 0;
  mbCrossId7bit = 0;
  mMagField = 0;
  mZdcEastRate = 0;
  mZdcWestRate = 0;
  mZdcCoincRate = 0;
  mLastDSM0 = 0;
  mLastDSM1 = 0;
  mLastDSM3 = 0;
  mZdcEastUA = 0;
  mZdcWestUA = 0;
  mZdcEastTDC = 0;
  mZdcWestTDC = 0;
  mZdcTimeDiff = 0;
  mZdcVertexZ = 0;
  mBBCSmallEast = 0;
  mBBCSmallWest = 0;
  mBBCLargeEast = 0;
  mBBCLargeWest = 0;
  mVPDSumEast = 0;
  mVPDSumWest = 0;
  mVPDTimeDiff = 0;
  mTofMult = 0;
  mBemcMult = 0;
  mNGlobTracks = 0;
  mNPrimTracks = 0;
  mNPrimVertices = 0;
  for(Int_t i=0; i<mNZdcPmt; ++i) mZdcEastADC[i] = 0;
  for(Int_t i=0; i<mNZdcPmt; ++i) mZdcWestADC[i] = 0;
  for(Int_t i=0; i<mNZdcSmd; ++i) mZdcSmdEast[i] = 0;
  for(Int_t i=0; i<mNZdcSmd; ++i) mZdcSmdWest[i] = 0;
}

void StUPCMCEvent::copyEventInfo(const StUPCEvent &event, Bool_t hasTriggerData) {
  mTrgIDs = event.mTrgIDs;
  mRunNumber = event.mRunNum;
  mEventNumber = event.mEvtNum;
  mFillNum = event.mFillNum;
  mbCrossId = event.mbCrossId;
  mbCrossId7bit = event.mbCrossId7bit;
  mMagField = event.mMagField;
  mZdcEastRate = event.mZdcEastRate;
  mZdcWestRate = event.mZdcWestRate;
  mZdcCoincRate = event.mZdcCoincRate;
  mLastDSM0 = event.mLastDSM0;
  mLastDSM1 = event.mLastDSM1;
  mLastDSM3 = event.mLastDSM3;
  mZdcEastUA = event.mZdcEastUA;
  mZdcWestUA = event.mZdcWestUA;
  mZdcEastTDC = event.mZdcEastTDC;
  mZdcWestTDC = event.mZdcWestTDC;
  mZdcTimeDiff = event.mZdcTimeDiff;
  mZdcVertexZ = event.mZdcVertexZ;
  mBBCSmallEast = event.mBBCSmallEast;
  mBBCSmallWest = event.mBBCSmallWest;
  mBBCLargeEast = event.mBBCLargeEast;
  mBBCLargeWest = event.mBBCLargeWest;
  mVPDSumEast = event.mVPDSumEast;
  mVPDSumWest = event.mVPDSumWest;
  mVPDTimeDiff = event.mVPDTimeDiff;
  mTofMult = event.mTofMult;
  mBemcMult = event.mBemcMult;
  mNGlobTracks = event.mNGlobTracks;
  mNPrimTracks = event.mNPrimTracks;
  mNPrimVertices = event.mNPrimVertices;
  for(Int_t i=0; i<mNZdcPmt; ++i) mZdcEastADC[i] = event.mZdcEastADC[i];
  for(Int_t i=0; i<mNZdcPmt; ++i) mZdcWestADC[i] = event.mZdcWestADC[i];
  for(Int_t i=0; i<mNZdcSmd; ++i) mZdcSmdEast[i] = event.mZdcSmdEast[i];
  for(Int_t i=0; i<mNZdcSmd; ++i) mZdcSmdWest[i] = event.mZdcSmdWest[i];
  mNRecoTracks = event.getNumberOfTracks();
  mNRecoVertices = event.getNumberOfVertices();
  mNClusters = event.getNumberOfClusters();
  mNHits = event.getNumberOfHits();
  mHasEventInfo = kTRUE;
  mHasDSMData = hasTriggerData;
  if(!hasTriggerData) mLastDSM0 = mLastDSM1 = mLastDSM3 = 0;
}

Bool_t StUPCMCEvent::isTrigger(Int_t id) const {
  for(Int_t i=0; i<mTrgIDs.GetSize(); ++i) if(mTrgIDs[i]==id) return kTRUE;
  return kFALSE;
}
