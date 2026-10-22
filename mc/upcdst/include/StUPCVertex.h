#ifndef StUPCVertex_h
#define StUPCVertex_h

//_____________________________________________________________________________
//    Class for UPC data
//    Author: Jaroslav Adam
//_____________________________________________________________________________

#include "TObject.h"

class StUPCEvent;
class StUPCMCVertex;

class StUPCVertex: public TObject
{

public:

  StUPCVertex();
  ~StUPCVertex() {}
  void Clear(Option_t *option="");

  //setters
  void setPosX(Float_t x) { mPosX = x; }
  void setPosY(Float_t y) { mPosY = y; }
  void setPosZ(Float_t z) { mPosZ = z; }

  void setErrX(Float_t x) { mErrX = x; }
  void setErrY(Float_t y) { mErrY = y; }
  void setErrZ(Float_t z) { mErrZ = z; }

  void setNPrimaryTracks(Int_t ntrk) { mNPrimaryTracks = ntrk; }
  void setNTracksUsed(UShort_t ntrk) { mNTracksUsed = ntrk; }

  void setId(UInt_t id) { mId = id; }
  void setIdTruth(Int_t id) { mIdTruth = id; }
  void setQaTruth(Int_t quality) { mQaTruth = quality; }
  void setIdParent(Int_t id) { mIdParent = id; }
  void setMCVertexIndex(Int_t index) { mMCVertexIndex = index; }
  void setEvent(StUPCEvent *event) { mEvent = event; }

  //getters
  Float_t getPosX() const { return mPosX; }
  Float_t getPosY() const { return mPosY; }
  Float_t getPosZ() const { return mPosZ; }

  Float_t getErrX() const { return mErrX; }
  Float_t getErrY() const { return mErrY; }
  Float_t getErrZ() const { return mErrZ; }

  Int_t getNPrimaryTracks() const { return mNPrimaryTracks; }
  Int_t getNTracksUsed() const { return mNTracksUsed; }

  UInt_t getId() const { return mId; }
  Int_t getIdTruth() const { return mIdTruth; }
  Int_t getQaTruth() const { return mQaTruth; }
  Int_t getIdParent() const { return mIdParent; }
  Int_t getMCVertexIndex() const { return mMCVertexIndex; }
  StUPCMCVertex *getMCVertex() const;

private:

  StUPCVertex(const StUPCVertex &o); // not implemented
  StUPCVertex &operator=(const StUPCVertex &o); // not implemented

  Float_t mPosX; // position x
  Float_t mPosY; // position y
  Float_t mPosZ; // position z

  Float_t mErrX; // error in position x
  Float_t mErrY; // error in position y
  Float_t mErrZ; // error in position z

  Int_t mNPrimaryTracks; // number of all primary tracks associated with the vertex
  UShort_t mNTracksUsed; // number of tracks used for the vertex

  UInt_t mId; // original vertex ID

  Int_t mIdTruth; // raw StMuPrimaryVertex MC vertex ID
  Int_t mQaTruth; // raw StMuPrimaryVertex truth quality
  Int_t mIdParent; // raw StMuPrimaryVertex MC parent-track ID
  Int_t mMCVertexIndex; // output MC vertex index; -1 if no validated association
  StUPCEvent *mEvent; //! restored by StUPCEvent::getVertex

  ClassDef(StUPCVertex, 2)

};

#endif
