#ifndef StUPCMCTrack_h
#define StUPCMCTrack_h

#include "TParticle.h"
#include <vector>

class StUPCMCEvent;
class StUPCMCVertex;
class StUPCTrack;

// MuDst MC truth, independent of whether a primary track was reconstructed.
// TParticle supplies momentum, energy, PDG code and production position.
class StUPCMCTrack: public TParticle {
public:
  // Same detector ordering as SL23e StMuMcTrack::EHIT.
  enum Detector {kTpc, kSvt, kSsd, kCtb, kEem, kEmc, kEsm, kFtp,
    kGem, kHpd, kIst, kIgt, kFst, kFgt, kFpd, kMwc, kPgc, kPmd,
    kSmd, kPix, kTof, kVpd, kNumberOfDetectors};

  StUPCMCTrack();
  virtual ~StUPCMCTrack() {}
  void Clear(Option_t *option="");

  void setId(Int_t id) { mId = id; }
  void setGeantId(Int_t id) { mGeantId = id; }
  void setCharge(Short_t charge) { mCharge = charge; }
  void setIsShower(Bool_t shower) { mIsShower = shower; }
  void setNhits(Int_t detector, UChar_t hits);
  void setSourceKinematics(Float_t pt, Float_t eta, Float_t p, Float_t rapidity);
  void setStartVertexId(Int_t id) { mStartVertexId = id; }
  void setStopVertexId(Int_t id) { mStopVertexId = id; }
  void setIntermediateVertexId(Int_t id) { mIntermediateVertexId = id; }
  void setStartVertexIndex(Int_t index) { mStartVertexIndex = index; }
  void setStopVertexIndex(Int_t index) { mStopVertexIndex = index; }
  void setIntermediateVertexIndex(Int_t index) { mIntermediateVertexIndex = index; }
  void setEvent(StUPCMCEvent *event) { mEvent = event; }
  void addRecoTrackIndex(Int_t index);

  Int_t getId() const { return mId; }
  Int_t getGeantId() const { return mGeantId; }
  Int_t getPdgCode() const { return GetPdgCode(); }
  Short_t getCharge() const { return mCharge; }
  Bool_t getIsShower() const { return mIsShower; }
  UChar_t getNhits(Int_t detector) const;
  Int_t getNhits() const;
  UChar_t getNhitsTpc() const { return getNhits(kTpc); }
  // Preserve the original MuDst scalars, including its low-pT sentinel.
  Float_t getPt() const { return mPt; }
  Float_t getEta() const { return mEta; }
  Float_t getPtot() const { return mPtot; }
  Float_t getRapidity() const { return mRapidity; }
  Double_t getPhi() const;
  Int_t getStartVertexId() const { return mStartVertexId; }
  Int_t getStopVertexId() const { return mStopVertexId; }
  Int_t getIntermediateVertexId() const { return mIntermediateVertexId; }
  Int_t getStartVertexIndex() const { return mStartVertexIndex; }
  Int_t getStopVertexIndex() const { return mStopVertexIndex; }
  Int_t getIntermediateVertexIndex() const { return mIntermediateVertexIndex; }
  StUPCMCVertex *getStartVertex() const;
  StUPCMCVertex *getStopVertex() const;
  StUPCMCVertex *getIntermediateVertex() const;
  StUPCMCTrack *getParentTrack() const;
  StUPCMCEvent *getEvent() const { return mEvent; }
  const std::vector<Int_t>& getRecoTrackIndices() const { return mRecoTrackIndices; }
  Int_t getNumberOfRecoTracks() const { return mRecoTrackIndices.size(); }
  Int_t getRecoTrackIndex(Int_t match) const;
  StUPCTrack *getRecoTrack(Int_t match=0) const;

private:
  Int_t mId; // original MC track ID, scoped to this event
  Int_t mGeantId; // raw Geant3 ID, retained even when PDG conversion is unknown
  Short_t mCharge;
  Bool_t mIsShower;
  UChar_t mHits[kNumberOfDetectors];
  Float_t mPt;
  Float_t mEta;
  Float_t mPtot;
  Float_t mRapidity;
  Int_t mStartVertexId;
  Int_t mStopVertexId;
  Int_t mIntermediateVertexId;
  Int_t mStartVertexIndex; // output MC vertex index, -1 if absent/unresolved
  Int_t mStopVertexIndex;
  Int_t mIntermediateVertexIndex;
  std::vector<Int_t> mRecoTrackIndices; // indices of ALL matched saved primary UPC tracks
  StUPCMCEvent *mEvent; //! local navigation, restored by StUPCMCEvent::getTrack

  ClassDef(StUPCMCTrack, 1)
};
#endif
