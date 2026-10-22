#ifndef StUPCMCEvent_h
#define StUPCMCEvent_h

#include "TObject.h"
#include "TArrayI.h"

class TClonesArray;
class StUPCEvent;
class StUPCMCTrack;
class StUPCMCVertex;

// Complete MuDst MC collections, owned by the corresponding reconstructed UPC event.
class StUPCMCEvent: public TObject {
public:
  enum SimulationType {kUnknown=0, kStarsim=1, kEmbedding=2};

  StUPCMCEvent();
  virtual ~StUPCMCEvent();
  void Clear(Option_t *option="");
  void setRunNumber(Int_t run) { mRunNumber = run; }
  void setEventNumber(Int_t event) { mEventNumber = event; }
  void setSimulationType(Int_t type) { mSimulationType = type; }
  // Copy the input event's trigger words; these are not generated MC truth.
  void setLastDSM(UShort_t dsm0, UShort_t dsm1, UShort_t dsm3) {
    mLastDSM0 = dsm0; mLastDSM1 = dsm1; mLastDSM3 = dsm3;
    mHasDSMData = kTRUE;
  }
  // Snapshot the completed input-event metadata; do not change MC truth or links.
  void copyEventInfo(const StUPCEvent &event, Bool_t hasTriggerData);
  Bool_t hasEventInfo() const { return mHasEventInfo; }
  Bool_t hasTriggerData() const { return mHasDSMData; }
  Bool_t isTrigger(Int_t id) const;
  const TArrayI &getTriggerIds() const { return mTrgIDs; }
  Int_t getFillNumber() const { return mFillNum; }
  UInt_t getBunchCrossId() const { return mbCrossId; }
  UInt_t getBunchCrossId7bit() const { return mbCrossId7bit; }
  Double_t getMagneticField() const { return mMagField; }

  Double_t getZDCEastRate() const { return mZdcEastRate; }
  Double_t getZDCWestRate() const { return mZdcWestRate; }
  Double_t getZDCCoincRate() const { return mZdcCoincRate; }


  UShort_t getZDCUnAttEast() const { return mZdcEastUA; }
  UShort_t getZDCUnAttWest() const { return mZdcWestUA; }
  UShort_t getZDCEastADC(Int_t pmt) const { return mZdcEastADC[pmt-1]; }
  UShort_t getZDCWestADC(Int_t pmt) const { return mZdcWestADC[pmt-1]; }
  UShort_t getZDCEastTDC() const { return mZdcEastTDC; }
  UShort_t getZDCWestTDC() const { return mZdcWestTDC; }
  UShort_t getZDCTimeDiff() const { return mZdcTimeDiff; }
  Float_t getZdcVertexZ() const { return mZdcVertexZ/10.; } // convert from mm to cm

  UShort_t getZdcSMDEast(Int_t verthori, Int_t strip) const { return mZdcSmdEast[strip-1 + verthori*8]; }
  UShort_t getZdcSMDWest(Int_t verthori, Int_t strip) const { return mZdcSmdWest[strip-1 + verthori*8]; }

  UInt_t getBBCSmallEast() const { return mBBCSmallEast; }
  UInt_t getBBCSmallWest() const { return mBBCSmallWest; }
  UInt_t getBBCLargeEast() const { return mBBCLargeEast; }
  UInt_t getBBCLargeWest() const { return mBBCLargeWest; }

  UShort_t getVPDSumEast() const { return mVPDSumEast; }
  UShort_t getVPDSumWest() const { return mVPDSumWest; }
  UShort_t getVPDTimeDiff() const { return mVPDTimeDiff; }

  UShort_t getTOFMultiplicity() const { return mTofMult; }
  UInt_t getBEMCMultiplicity() const { return mBemcMult; }

  Int_t getNGlobTracks() const { return mNGlobTracks; }
  Int_t getNPrimTracks() const { return mNPrimTracks; }
  UInt_t getNPrimVertices() const { return mNPrimVertices; }


  Int_t getNumberOfRecoTracks() const { return mNRecoTracks; }
  Int_t getNumberOfRecoVertices() const { return mNRecoVertices; }
  Int_t getNumberOfClusters() const { return mNClusters; }
  Int_t getNumberOfHits() const { return mNHits; }

  void setRecoEvent(StUPCEvent *event) { mRecoEvent = event; }
  StUPCMCTrack *addTrack();
  StUPCMCVertex *addVertex();

  Int_t getRunNumber() const { return mRunNumber; }
  Int_t getEventNumber() const { return mEventNumber; }
  Int_t getSimulationType() const { return mSimulationType; }
  Bool_t hasDSMData() const { return mHasDSMData; }
  UShort_t getLastDSM0() const { return mLastDSM0; }
  UShort_t getLastDSM1() const { return mLastDSM1; }
  UShort_t getLastDSM3() const { return mLastDSM3; }
  Int_t getNumberOfTracks() const;
  Int_t getNumberOfVertices() const;
  StUPCMCTrack *getTrack(Int_t index) const;
  StUPCMCVertex *getVertex(Int_t index) const;
  StUPCEvent *getRecoEvent() const { return mRecoEvent; }

private:
  StUPCMCEvent(const StUPCMCEvent&); // not implemented: owns its collections
  StUPCMCEvent& operator=(const StUPCMCEvent&); // not implemented

  void clearEventInfo();
  Bool_t mHasEventInfo; // full header snapshot available (schema >= 3)
  TArrayI mTrgIDs; // fired trigger IDs

  Int_t mFillNum; // beam fill number
  UInt_t mbCrossId; // bunch crossing ID
  UInt_t mbCrossId7bit; // bunch crossing ID 7bit
  Double32_t mMagField; // magnetic field

  Double32_t mZdcEastRate; // ZDC rate, east
  Double32_t mZdcWestRate; // ZDC rate, west
  Double32_t mZdcCoincRate; // ZDC rate, east&&west coincidence


  UShort_t mZdcEastUA; // ZDC unattenuated signal, east
  UShort_t mZdcWestUA; // ZDC unattenuated signal, west
  static const Int_t mNZdcPmt = 3;
  UShort_t mZdcEastADC[mNZdcPmt]; // ZDC 3 PMT ADCs, east
  UShort_t mZdcWestADC[mNZdcPmt]; // ZDC 3 PMT ADCs, west
  UShort_t mZdcEastTDC; // ZDC TDC, east
  UShort_t mZdcWestTDC; // ZDC TDC, west
  UShort_t mZdcTimeDiff; // ZDC time difference, east-west
  Float_t mZdcVertexZ; // ZDC vertex z position, mm
  static const Int_t mNZdcSmd = 16; // number of ZDC SMD channels, east and west
  UShort_t mZdcSmdEast[mNZdcSmd]; // ZDC SMD data east
  UShort_t mZdcSmdWest[mNZdcSmd]; // ZDC SMD data west

  UInt_t mBBCSmallEast; // BBC truncated sum, small tiles, east
  UInt_t mBBCSmallWest; // BBC truncated sum, small tiles, west
  UInt_t mBBCLargeEast; // BBC truncated sum, large tiles, east
  UInt_t mBBCLargeWest; // BBC truncated sum, large tiles, west

  UShort_t mVPDSumEast; // VPD ADC sum east
  UShort_t mVPDSumWest; // VPD ADC sum west
  UShort_t mVPDTimeDiff; // VPD time difference

  UShort_t mTofMult; // TOF multiplicity from StTriggerData
  UInt_t mBemcMult; // BEMC multiplicity, number of all BEMC clusters in event

  Int_t mNGlobTracks; // number of global tracks
  Int_t mNPrimTracks; // number of primary tracks
  UInt_t mNPrimVertices; // number of primary vertices

  Int_t mNRecoTracks; // saved primary tracks, distinct from MC count
  Int_t mNRecoVertices; // saved reconstructed vertices
  Int_t mNClusters; // saved BEMC clusters
  Int_t mNHits; // saved TOF hits

  Int_t mRunNumber;
  Int_t mEventNumber;
  Int_t mSimulationType;
  Bool_t mHasDSMData; // input trigger data was available; false for old files
  UShort_t mLastDSM0; // input StTriggerData::lastDSM(0)
  UShort_t mLastDSM1; // input StTriggerData::lastDSM(1)
  UShort_t mLastDSM3; // input StTriggerData::lastDSM(3)
  TClonesArray *mTracks; //-> complete StUPCMCTrack collection
  TClonesArray *mVertices; //-> complete StUPCMCVertex collection
  Int_t mNTracks; //! local append counter
  Int_t mNVertices; //! local append counter
  StUPCEvent *mRecoEvent; //! restored by StUPCEvent::getMCEvent

  ClassDef(StUPCMCEvent, 3)
};
#endif
