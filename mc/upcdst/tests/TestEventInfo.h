#ifndef TestEventInfo_h
#define TestEventInfo_h
#include <iostream>
#include "StUPCEvent.h"
#include "StUPCMCEvent.h"

// Compare every persisted event header field after ROOT serialization.
inline bool sameEventInfo(const StUPCMCEvent &mc, const StUPCEvent &rc) {
#define CHECK_EVENT_FIELD(method) if(mc.method()!=rc.method()) { std::cerr << #method << " mismatch" << std::endl; return false; }
  CHECK_EVENT_FIELD(getRunNumber)
  CHECK_EVENT_FIELD(getEventNumber)
  CHECK_EVENT_FIELD(getFillNumber)
  CHECK_EVENT_FIELD(getBunchCrossId)
  CHECK_EVENT_FIELD(getBunchCrossId7bit)
  CHECK_EVENT_FIELD(getMagneticField)
  CHECK_EVENT_FIELD(getZDCEastRate)
  CHECK_EVENT_FIELD(getZDCWestRate)
  CHECK_EVENT_FIELD(getZDCCoincRate)
  CHECK_EVENT_FIELD(getLastDSM0)
  CHECK_EVENT_FIELD(getLastDSM1)
  CHECK_EVENT_FIELD(getLastDSM3)
  CHECK_EVENT_FIELD(getZDCUnAttEast)
  CHECK_EVENT_FIELD(getZDCUnAttWest)
  CHECK_EVENT_FIELD(getZDCEastTDC)
  CHECK_EVENT_FIELD(getZDCWestTDC)
  CHECK_EVENT_FIELD(getZDCTimeDiff)
  CHECK_EVENT_FIELD(getZdcVertexZ)
  CHECK_EVENT_FIELD(getBBCSmallEast)
  CHECK_EVENT_FIELD(getBBCSmallWest)
  CHECK_EVENT_FIELD(getBBCLargeEast)
  CHECK_EVENT_FIELD(getBBCLargeWest)
  CHECK_EVENT_FIELD(getVPDSumEast)
  CHECK_EVENT_FIELD(getVPDSumWest)
  CHECK_EVENT_FIELD(getVPDTimeDiff)
  CHECK_EVENT_FIELD(getTOFMultiplicity)
  CHECK_EVENT_FIELD(getBEMCMultiplicity)
  CHECK_EVENT_FIELD(getNGlobTracks)
  CHECK_EVENT_FIELD(getNPrimTracks)
  CHECK_EVENT_FIELD(getNPrimVertices)
#undef CHECK_EVENT_FIELD
  if(mc.getTriggerIds().GetSize()!=rc.getTriggerIds().GetSize()) return false;
  for(int i=0;i<rc.getTriggerIds().GetSize();++i)
    if(mc.getTriggerIds()[i]!=rc.getTriggerIds()[i] || !mc.isTrigger(rc.getTriggerIds()[i])) return false;
  for(int i=1;i<=3;++i)
    if(mc.getZDCEastADC(i)!=rc.getZDCEastADC(i) || mc.getZDCWestADC(i)!=rc.getZDCWestADC(i)) return false;
  for(int h=0;h<2;++h) for(int i=1;i<=8;++i)
    if(mc.getZdcSMDEast(h,i)!=rc.getZdcSMDEast(h,i) || mc.getZdcSMDWest(h,i)!=rc.getZdcSMDWest(h,i)) return false;
  return mc.getNumberOfRecoTracks()==rc.getNumberOfTracks()
    && mc.getNumberOfRecoVertices()==rc.getNumberOfVertices()
    && mc.getNumberOfClusters()==rc.getNumberOfClusters()
    && mc.getNumberOfHits()==rc.getNumberOfHits();
}

inline void fillEventInfoFixture(StUPCEvent &event) {
  event.setRunNumber(19); event.setEventNumber(23);
  event.setFillNumber(456); event.setBunchCrossId(0xfedcba98U); event.setBunchCrossId7bit(123);
  event.setMagneticField(-4.75);
  event.addTriggerId(600001); event.addTriggerId(600009);
  event.setZDCEastRate(123.25); event.setZDCWestRate(456.5); event.setZDCCoincRate(78.75);
  event.setLastDSM0(0xffff); event.setLastDSM1(0x8000); event.setLastDSM3(0x1234);
  event.setZDCUnAttEast(100); event.setZDCUnAttWest(101);
  event.setZDCEastTDC(102); event.setZDCWestTDC(103); event.setZDCTimeDiff(104);
  event.setZdcVertexZ(-1234.567f);
  for(int i=1;i<=3;++i) { event.setZDCEastADC(200+i,i); event.setZDCWestADC(300+i,i); }
  for(int h=0;h<2;++h) for(int i=1;i<=8;++i) {
    event.setZdcSMDEast(h,i,400+8*h+i); event.setZdcSMDWest(h,i,500+8*h+i);
  }
  event.setBBCSmallEast(0xfffffffeU); event.setBBCSmallWest(600);
  event.setBBCLargeEast(601); event.setBBCLargeWest(602);
  event.setVPDSumEast(603); event.setVPDSumWest(604); event.setVPDTimeDiff(605);
  event.setTOFMultiplicity(606); event.setBEMCMultiplicity(607);
  event.setNGlobTracks(608); event.setNPrimTracks(609); event.setNPrimVertices(610);
}
#endif
