#ifndef StUPCMCVertex_h
#define StUPCMCVertex_h

#include "TObject.h"
#include <vector>

class StUPCMCEvent;
class StUPCMCTrack;
class StUPCVertex;

class StUPCMCVertex: public TObject {
public:
  StUPCMCVertex();
  virtual ~StUPCMCVertex() {}
  void Clear(Option_t *option="");

  void setId(Int_t id) { mId = id; }
  void setPosition(Float_t x, Float_t y, Float_t z) { mPosX=x; mPosY=y; mPosZ=z; }
  void setTime(Float_t time) { mTime = time; }
  void setNDaughters(Int_t n) { mNDaughters = n; }
  void setIntermediateFlag(Int_t flag) { mIntermediateFlag = flag; }
  void setParentTrackId(Int_t id) { mParentTrackId = id; }
  void setParentTrackIndex(Int_t index) { mParentTrackIndex = index; }
  void setEvent(StUPCMCEvent *event) { mEvent = event; }
  void addRecoVertexIndex(Int_t index);

  Int_t getId() const { return mId; }
  Float_t getPosX() const { return mPosX; }
  Float_t getPosY() const { return mPosY; }
  Float_t getPosZ() const { return mPosZ; }
  Float_t getTime() const { return mTime; }
  Int_t getNDaughters() const { return mNDaughters; }
  Int_t getIntermediateFlag() const { return mIntermediateFlag; }
  Int_t getParentTrackId() const { return mParentTrackId; }
  Int_t getParentTrackIndex() const { return mParentTrackIndex; }
  StUPCMCTrack *getParentTrack() const;
  StUPCMCEvent *getEvent() const { return mEvent; }
  const std::vector<Int_t>& getRecoVertexIndices() const { return mRecoVertexIndices; }
  Int_t getNumberOfRecoVertices() const { return mRecoVertexIndices.size(); }
  Int_t getRecoVertexIndex(Int_t match) const;
  StUPCVertex *getRecoVertex(Int_t match=0) const;

private:
  Int_t mId; // original MC vertex ID, scoped to this event
  Float_t mPosX; // cm, copied from StMuMcVertex
  Float_t mPosY;
  Float_t mPosZ;
  Float_t mTime; // copied without conversion from StMuMcVertex::Time()
  Int_t mNDaughters;
  Int_t mIntermediateFlag;
  Int_t mParentTrackId; // original MC parent-track ID
  Int_t mParentTrackIndex; // output MC track index, -1 if absent/unresolved
  std::vector<Int_t> mRecoVertexIndices; // indices of ALL matched saved reco UPC vertices
  StUPCMCEvent *mEvent; //! local navigation, restored by StUPCMCEvent::getVertex

  ClassDef(StUPCMCVertex, 1)
};
#endif
