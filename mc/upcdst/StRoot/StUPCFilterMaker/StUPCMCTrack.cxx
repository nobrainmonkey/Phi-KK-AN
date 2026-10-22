#include "StUPCMCTrack.h"
#include "StUPCMCVertex.h"
#include "StUPCMCEvent.h"
#include "StUPCEvent.h"
#include <algorithm>
#include <cmath>

ClassImp(StUPCMCTrack);

StUPCMCTrack::StUPCMCTrack(): TParticle() { Clear(); }

void StUPCMCTrack::Clear(Option_t *) {
  TParticle::operator=(TParticle());
  mId = mGeantId = mCharge = 0;
  mIsShower = kFALSE;
  for(Int_t i=0; i<kNumberOfDetectors; ++i) mHits[i] = 0;
  mPt = mEta = mPtot = mRapidity = 0;
  mStartVertexId = mStopVertexId = mIntermediateVertexId = 0;
  mStartVertexIndex = mStopVertexIndex = mIntermediateVertexIndex = -1;
  mEvent = 0;
  mRecoTrackIndices.clear();
}

void StUPCMCTrack::setNhits(Int_t detector, UChar_t hits) {
  if(detector >= 0 && detector < kNumberOfDetectors) mHits[detector] = hits;
}

UChar_t StUPCMCTrack::getNhits(Int_t detector) const {
  return detector >= 0 && detector < kNumberOfDetectors ? mHits[detector] : 0;
}

Int_t StUPCMCTrack::getNhits() const {
  Int_t sum = 0;
  for(Int_t i=0; i<kNumberOfDetectors; ++i) sum += mHits[i];
  return sum;
}

void StUPCMCTrack::setSourceKinematics(Float_t pt, Float_t eta, Float_t p, Float_t rapidity) {
  mPt = pt; mEta = eta; mPtot = p; mRapidity = rapidity;
}

Double_t StUPCMCTrack::getPhi() const { return std::atan2(Py(), Px()); }

StUPCMCVertex *StUPCMCTrack::getStartVertex() const {
  return mEvent ? mEvent->getVertex(mStartVertexIndex) : 0;
}

StUPCMCVertex *StUPCMCTrack::getStopVertex() const {
  return mEvent ? mEvent->getVertex(mStopVertexIndex) : 0;
}

StUPCMCVertex *StUPCMCTrack::getIntermediateVertex() const {
  return mEvent ? mEvent->getVertex(mIntermediateVertexIndex) : 0;
}

StUPCMCTrack *StUPCMCTrack::getParentTrack() const {
  StUPCMCVertex *vertex = getStartVertex();
  return vertex ? vertex->getParentTrack() : 0;
}

void StUPCMCTrack::addRecoTrackIndex(Int_t index) {
  if(index >= 0 && std::find(mRecoTrackIndices.begin(),mRecoTrackIndices.end(),index)==mRecoTrackIndices.end())
    mRecoTrackIndices.push_back(index);
}

Int_t StUPCMCTrack::getRecoTrackIndex(Int_t match) const {
  return match >= 0 && match < getNumberOfRecoTracks() ? mRecoTrackIndices[match] : -1;
}

StUPCTrack *StUPCMCTrack::getRecoTrack(Int_t match) const {
  Int_t index = getRecoTrackIndex(match);
  StUPCEvent *reco = mEvent ? mEvent->getRecoEvent() : 0;
  return reco && index >= 0 && index < reco->getNumberOfTracks() ? reco->getTrack(index) : 0;
}
