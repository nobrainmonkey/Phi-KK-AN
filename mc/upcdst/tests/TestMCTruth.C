#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TDatabasePDG.h"
#include "TBufferFile.h"
#include "StUPCEvent.h"
#include "StUPCTrack.h"
#include "StUPCVertex.h"
#include "StUPCMCTrack.h"
#include "StUPCMCVertex.h"
#include "StUPCMCEvent.h"
#include "StChain/StChain.h"
#include "StMuDSTMaker/COMMON/StMuDstMaker.h"
#include "StMuDSTMaker/COMMON/StMuDst.h"
#include "StMuDSTMaker/COMMON/StMuEvent.h"
#include "StMuDSTMaker/COMMON/StMuTrack.h"
#include "StMuDSTMaker/COMMON/StMuPrimaryVertex.h"
#include "StMuDSTMaker/COMMON/StMuMcTrack.h"
#include "StMuDSTMaker/COMMON/StMuMcVertex.h"
#include "StUPCFilterMaker.h"
#include "StEvent/StTriggerData.h"
#include "TestEventInfo.h"

namespace {
void require(bool condition, const char *message) {
  if(!condition) throw std::runtime_error(message);
}
void same(double actual, double expected, const char *message) {
  if(actual == expected || std::fabs(actual-expected) <= 1.e-12+1.e-6*std::fabs(expected)) return;
  std::cerr << std::setprecision(16) << message << ": saved=" << actual << " source=" << expected << std::endl;
  require(false, message);
}
double maxPidRounding = 0;
void samePID(float actual, float source, const char *message) {
  // StUPCTrack already stores nSigma as Float16_t. Test its exact ROOT encoding,
  // rather than accepting an arbitrary tolerance or requiring lossless floats.
  TBufferFile buffer(TBuffer::kWrite);
  buffer.WriteFloat16(&source);
  buffer.SetReadMode();
  buffer.SetBufferOffset(0);
  float expected=0;
  buffer.ReadFloat16(&expected);
  same(actual,expected,message);
  if(std::fabs(actual-source)>maxPidRounding) maxPidRounding=std::fabs(actual-source);
}
int indexOf(const std::map<int,int>& ids, int id) {
  std::map<int,int>::const_iterator it = ids.find(id);
  return id > 0 && it != ids.end() ? it->second : -1;
}

void testIO(const char *output, const char *legacy) {
  // Independent objects exercise ownership and empty/data defaults.
  StUPCEvent event, second;
  second.addTrack()->setPtEtaPhi(9,0,0);
  event.setIsMC();
  event.getMCEvent()->setRunNumber(19);
  event.getMCEvent()->setEventNumber(23);
  event.getMCEvent()->setSimulationType(StUPCMCEvent::kEmbedding);
  require(!event.getMCEvent()->hasDSMData(),"DSM unavailable by default");
  event.getMCEvent()->setLastDSM(0xffff,0x8000,0x1234);
  StUPCMCVertex *start = event.addMCVertex();
  start->setId(91); start->setPosition(1,2,3); start->setTime(0.25);
  start->setNDaughters(2); start->setIntermediateFlag(3);
  StUPCMCVertex *decay = event.addMCVertex();
  decay->setId(17); decay->setPosition(4,5,6); decay->setParentTrackId(42);
  decay->setParentTrackIndex(0);
  StUPCMCTrack *parent = event.addMCTrack();
  parent->setId(42); parent->setGeantId(11); parent->setCharge(1);
  parent->SetPdgCode(321); parent->SetMomentum(0.1,0.2,0.3,0.8);
  parent->SetProductionVertex(1,2,3,0.25);
  parent->setSourceKinematics(-999,0.9,0.4,0.5);
  parent->setIsShower(kTRUE);
  parent->setStartVertexId(91); parent->setStartVertexIndex(0);
  parent->setStopVertexId(17); parent->setStopVertexIndex(1);
  parent->setIntermediateVertexId(17); parent->setIntermediateVertexIndex(1);
  for(int d=0; d<StUPCMCTrack::kNumberOfDetectors; ++d) parent->setNhits(d,200+d);
  StUPCMCTrack *child = event.addMCTrack();
  child->setId(7); child->setStartVertexId(17); child->setStartVertexIndex(1);
  for(int i=0; i<2; ++i) {
    StUPCTrack *rc = event.addTrack();
    rc->setIdTruth(42); rc->setQaTruth(97-i); rc->setIdParentVx(91);
    rc->setMCTrackIndex(0); rc->setVertexId(5);
    parent->addRecoTrackIndex(i);
  }
  StUPCTrack *unmatched = event.addTrack();
  unmatched->setIdTruth(65535); unmatched->setQaTruth(80);
  StUPCVertex *rcv = event.addVertex();
  rcv->setId(5); rcv->setIdTruth(91); rcv->setQaTruth(99);
  rcv->setIdParent(42); rcv->setMCVertexIndex(0);
  start->addRecoVertexIndex(0);
  require(event.getRecoTrackIndices(0).size()==2, "multiple reco matches");
  require(child->getParentTrack()==parent, "MC ancestry navigation");
  require(event.getTrack(0)->getVertex()->getMCVertex()==start, "reco vertex navigation");
  require(!unmatched->getMCTrack(), "unmatched track must not resolve");
  require(!event.getMCTrack(-1) && !event.getMCVertex(99), "invalid indices must return null");

  fillEventInfoFixture(event);
  event.addCluster(); event.addHit();
  event.getMCEvent()->copyEventInfo(event,kTRUE);
  require(event.getMCEvent()->hasEventInfo() && sameEventInfo(*event.getMCEvent(),event),"copy all event metadata");
  require(!event.getMCEvent()->isTrigger(123),"absent trigger ID");

  TFile file(output,"RECREATE");
  require(!file.IsZombie(), "create I/O test file");
  TTree tree("mUPCTree","MC schema test");
  StUPCEvent *ptr = &event;
  tree.Branch("mUPCEvent", &ptr);
  tree.Fill();
  event.clearEvent();
  require(second.getNumberOfTracks()==1, "event arrays must not be shared");
  require(event.getNumberOfMCTracks()==0 && event.getNumberOfMCVertices()==0, "clear truth collections");
  require(!event.getMCEvent()->hasDSMData() && event.getMCEvent()->getLastDSM0()==0
          && event.getMCEvent()->getLastDSM1()==0 && event.getMCEvent()->getLastDSM3()==0,"clear DSM availability and words");
  require(!event.getMCEvent()->hasEventInfo(),"clear event metadata availability");
  StUPCEvent empty;
  require(sameEventInfo(*event.getMCEvent(),empty),"clear every copied metadata field");
  StUPCMCEvent emptyCopy;
  emptyCopy.copyEventInfo(event,kFALSE);
  require(emptyCopy.hasEventInfo() && !emptyCopy.hasTriggerData() && sameEventInfo(emptyCopy,empty),"clear every reco header field before missing-trigger event");
  StUPCTrack *fresh = event.addTrack();
  require(fresh->getIdTruth()==0 && fresh->getQaTruth()==0 && fresh->getMCTrackIndex()==-1, "reset reco truth");
  rcv = event.addVertex();
  require(rcv->getMCVertexIndex()==-1 && rcv->getIdTruth()==0, "reset vertex truth");
  event.getMCEvent()->copyEventInfo(event,kFALSE);
  tree.Fill();
  event.clearEvent();
  parent = event.addMCTrack();
  require(parent->getStartVertexIndex()==-1 && parent->getId()==0 && parent->getNhits()==0, "reset MC truth");
  parent->setId(1001);
  start = event.addMCVertex();
  require(start->getParentTrackIndex()==-1 && start->getTime()==0, "reset MC vertex");
  start->setId(2001);
  event.getMCEvent()->copyEventInfo(event,kTRUE); // available all-zero data differs from missing data
  tree.Fill();
  tree.Write();
  tree.SetDirectory(0);
  file.Close();

  TFile read(output);
  TTree *saved = dynamic_cast<TTree*>(read.Get("mUPCTree"));
  require(saved && saved->GetEntries()==3, "read I/O test entries");
  StUPCEvent *restored = 0;
  saved->SetBranchAddress("mUPCEvent", &restored);
  require(saved->GetEntry(0)>0, "read truth event");
  require(restored->getMCEvent() && restored->getMCEvent()->getRunNumber()==19
          && restored->getMCEvent()->getEventNumber()==23
          && restored->getMCEvent()->getSimulationType()==StUPCMCEvent::kEmbedding,"persist MC event metadata");
  require(restored->getMCEvent()->hasDSMData() && restored->getMCEvent()->getLastDSM0()==0xffff
          && restored->getMCEvent()->getLastDSM1()==0x8000 && restored->getMCEvent()->getLastDSM3()==0x1234,"persist full 16-bit DSM words");
  require(restored->getMCEvent()->hasEventInfo() && sameEventInfo(*restored->getMCEvent(),*restored),"persist all event metadata");
  require(restored->getMCEvent()->getFillNumber()==456 && restored->getMCEvent()->getBBCSmallEast()==0xfffffffeU
          && restored->getMCEvent()->getBunchCrossId()==0xfedcba98U,"persist full-width header fields");
  require(restored->getMCEvent()->getZdcVertexZ()==Float_t(-1234.567f/10.),"preserve ZDC mm-to-cm units");
  parent = restored->getMCTrack(0);
  require(parent && parent->getId()==42 && parent->getGeantId()==11, "persist MC identity");
  require(parent->getCharge()==1 && parent->getIsShower(), "persist MC flags");
  same(parent->Px(),0.1,"persist truth momentum");
  same(parent->Energy(),0.8,"persist truth energy");
  same(parent->getPt(),-999,"preserve source pT sentinel");
  same(parent->getEta(),0.9,"persist source eta");
  same(parent->getPtot(),0.4,"persist source momentum");
  same(parent->getRapidity(),0.5,"persist source rapidity");
  for(int d=0; d<StUPCMCTrack::kNumberOfDetectors; ++d) require(parent->getNhits(d)==200+d,"persist detector hit counts");
  require(parent->getStartVertex()->getId()==91 && parent->getStopVertex()->getId()==17, "persist MC vertex links");
  require(parent->getIntermediateVertex()->getId()==17,"persist intermediate vertex link");
  require(restored->getMCTrack(1)->getParentTrack()==parent,"persist ancestry");
  require(restored->getTrack(0)->getMCTrack()==parent,"persist reco track link");
  require(restored->getTrack(0)->getQaTruth()==97 && restored->getTrack(0)->getIdParentVx()==91,"persist reco truth values");
  require(restored->getRecoTrackIndices(0).size()==2,"persist multiple matches");
  require(parent->getNumberOfRecoTracks()==2 && parent->getRecoTrackIndex(1)==1
          && parent->getRecoTrack(1)==restored->getTrack(1),"persist MC-to-reco track links");
  require(!parent->getRecoTrack(-1) && !parent->getRecoTrack(2),"invalid reverse track index");
  require(restored->getRecoVertexIndices(0).size()==1,"persist reverse vertex link");
  require(parent->getStartVertex()->getRecoVertex()==restored->getVertex(0),"persist MC-to-reco vertex navigation");
  require(restored->getVertex(0)->getIdParent()==42 && restored->getVertex(0)->getQaTruth()==99,"persist reco vertex truth");
  require(restored->getTrack(2)->getIdTruth()==65535 && !restored->getTrack(2)->getMCTrack(),"persist unmatched marker");
  require(restored->getMCParticle(0)==parent,"legacy particle accessor");
  saved->GetEntry(1);
  require(restored->getMCEvent()->hasEventInfo() && sameEventInfo(*restored->getMCEvent(),*restored),"persist missing-trigger event metadata");
  require(!restored->getMCEvent()->hasDSMData() && restored->getMCEvent()->getLastDSM0()==0
          && restored->getMCEvent()->getLastDSM1()==0 && restored->getMCEvent()->getLastDSM3()==0,"no DSM leak between entries");
  require(restored->getNumberOfMCTracks()==0 && !restored->getTrack(0)->getMCTrack(),"no association leak between entries");
  saved->GetEntry(2);
  require(restored->getMCEvent()->hasEventInfo() && sameEventInfo(*restored->getMCEvent(),*restored),"persist zero-trigger metadata");
  require(restored->getMCEvent()->hasDSMData() && restored->getMCEvent()->getLastDSM0()==0
          && restored->getMCEvent()->getLastDSM1()==0 && restored->getMCEvent()->getLastDSM3()==0,"persist valid zero DSM words");
  require(restored->getNumberOfMCTracks()==1 && restored->getMCTrack(0)->getId()==1001,"truth collection reuse");
  require(!restored->getMCTrack(0)->getStartVertex(),"no MC vertex link leak");
  require(restored->getMCTrack(0)->getNumberOfRecoTracks()==0
          && restored->getMCVertex(0)->getNumberOfRecoVertices()==0,"no reverse link leak");
  saved->ResetBranchAddresses();
  delete restored;
  read.Close();

  if(legacy && legacy[0]) {
    TFile old(legacy);
    TTree *oldTree = dynamic_cast<TTree*>(old.Get("mUPCTree"));
    require(oldTree!=0,"open legacy test tree");
    restored = 0;
    oldTree->SetBranchAddress("mUPCEvent",&restored);
    oldTree->GetEntry(0);
    require(restored->getNumberOfMCParticles()==1 && restored->getMCParticle(0)->GetPdgCode()==321,"read legacy TParticle");
    require(restored->getNumberOfMCTracks()==0 && !restored->getMCTrack(0),"do not fabricate missing legacy truth metadata");
    require(restored->getNumberOfMCVertices()==0,"legacy has no MC vertices");
    require(!restored->getMCEvent(),"legacy has no full MC event");
    require(restored->getTrack(0)->getIdTruth()==0 && restored->getTrack(0)->getMCTrackIndex()==-1,"legacy track defaults");
    require(restored->getVertex(0)->getMCVertexIndex()==-1,"legacy vertex defaults");
    same(restored->getTrack(0)->getPt(),0.4,"legacy reco kinematics");
    oldTree->GetEntry(1);
    require(restored->getNumberOfMCParticles()==0,"read legacy empty truth event");
    oldTree->ResetBranchAddresses();
    delete restored;
  }
  std::cout << "MC_IO_PASS" << std::endl;
}

void convert(const char *input, const char *output) {
  StChain *chain = new StChain;
  StMuDstMaker *reader = new StMuDstMaker(0,0,"",input,"",100000);
  StUPCFilterMaker *filter = new StUPCFilterMaker(reader,output);
  filter->setIsMC(2);
  require(chain->Init()==0,"initialize conversion chain");
  Long64_t entries = reader->chain()->GetEntries();
  require(entries>0,"nonempty embedding input");
  for(Long64_t i=0; i<entries; ++i) {
    chain->Clear();
    require(chain->Make()==0,"convert embedding event");
  }
  chain->Finish();
  delete chain;
  std::cout << "MC_CONVERT_PASS events=" << entries << std::endl;
}

void writeMuEfficiency(StMuDst *mu, Long64_t entry, std::ostream& out) {
  // Build the reference directly from MuDst IDs, without UPC indices or getters.
  TClonesArray *tracks=mu->mcArray(1), *vertices=mu->mcArray(0);
  std::map<int,StMuMcVertex*> vertexIds;
  for(int i=0;i<vertices->GetEntriesFast();++i) {
    StMuMcVertex *v=(StMuMcVertex*)vertices->At(i);
    vertexIds[v->Id()]=v;
  }
  std::map<int,std::vector<StMuTrack*> > all, kept;
  std::map<const StMuTrack*,int> vertexTruth;
  for(unsigned int iv=0;iv<mu->numberOfPrimaryVertices();++iv) {
    StMuDst::setVertexIndex(iv);
    TObjArray *primary=mu->primaryTracks();
    int count=0;
    for(int i=0;i<primary->GetEntriesFast();++i) if(primary->At(i)) ++count;
    for(int i=0;i<primary->GetEntriesFast();++i) {
      StMuTrack *rc=(StMuTrack*)primary->At(i);
      if(!rc) continue;
      all[rc->idTruth()].push_back(rc);
      if(count<=4) kept[rc->idTruth()].push_back(rc);
      vertexTruth[rc]=mu->primaryVertex()->idTruth();
    }
  }
  for(int i=0;i<tracks->GetEntriesFast();++i) {
    StMuMcTrack *mc=(StMuMcTrack*)tracks->At(i);
    if(mc->GePid()!=11 && mc->GePid()!=12) continue;
    StMuMcVertex *v=vertexIds[mc->IdVx()];
    if(!v || v->IdParTrk()!=0) continue;
    const std::vector<StMuTrack*>& matches=kept[mc->Id()];
    StMuTrack *best=0;
    int quality=0, selected=0, correctVertex=0, tof=0;
    for(size_t j=0;j<matches.size();++j) {
      StMuTrack *rc=matches[j];
      if(!best || rc->qaTruth()>best->qaTruth()
         || (rc->qaTruth()==best->qaTruth() && rc->nHitsFit()>best->nHitsFit())) best=rc;
      if(rc->qaTruth()>=50) ++quality;
      if(rc->qaTruth()<50 || rc->nHitsFit()<15 || rc->nHitsPoss()<=0
         || double(rc->nHitsFit())/rc->nHitsPoss()<0.52 || std::fabs(rc->eta())>=1.0 || rc->pt()<0.1) continue;
      ++selected;
      if(vertexTruth[rc]==mc->IdVx()) ++correctVertex;
      if(rc->btofPidTraits().matchFlag()!=0) ++tof;
    }
    out << entry << '\t' << mu->event()->runNumber() << '\t' << mu->event()->eventNumber()
        << '\t' << mc->Id() << '\t' << (mc->GePid()==11 ? 321 : -321)
        << '\t' << mc->pT() << '\t' << mc->Eta() << '\t' << std::atan2(mc->Pxyz().y(),mc->Pxyz().x())
        << '\t' << v->XyzV().z() << '\t' << int(mc->NoHits(StMuMcTrack::ktpc))
        << '\t' << matches.size() << '\t' << quality << '\t' << selected << '\t' << correctVertex << '\t' << tof
        << '\t' << (best ? best->pt() : -999) << '\t' << (best ? best->eta() : -999)
        << '\t' << (best ? best->phi() : -999) << '\t' << (best ? best->qaTruth() : 0)
        << '\t' << (best ? int(best->nHitsFit()) : 0) << '\t' << (best ? int(best->nHitsPoss()) : 0)
        << '\t' << (best ? int(best->charge()) : 0) << '\t' << (best ? vertexTruth[best] : 0)
        << '\t' << all[mc->Id()].size() << '\n';
  }
}

void validate(const char *input, const char *output, const char *reference) {
  std::ofstream efficiency;
  if(reference && reference[0]) {
    efficiency.open(reference);
    require(efficiency.good(),"open MuDst efficiency reference");
    efficiency << std::setprecision(12)
      << "entry\trun\tevent\tmc_id\tpdg\tmc_pt\tmc_eta\tmc_phi\tmc_vz\tmc_tpc_hits"
      << "\tn_match\tn_quality\tn_selected\tn_correct_vertex\tn_tof\tbest_rc_pt\tbest_rc_eta\tbest_rc_phi"
      << "\tbest_qa\tbest_nhitsfit\tbest_nhitsmax\tbest_charge\tbest_vertex_truth\tn_before_skim\n";
  }
  StChain *chain = new StChain;
  StMuDstMaker *reader = new StMuDstMaker(0,0,"",input,"",100000);
  reader->SetStatus("*",0);
  reader->SetStatus("MuEvent",1);
  reader->SetStatus("PrimaryVertices",1);
  reader->SetStatus("PrimaryTracks",1);
  reader->SetStatus("GlobalTracks",1);
  reader->SetStatus("MCAll",1);
  require(chain->Init()==0,"initialize validation reader");
  TFile file(output);
  TTree *tree = dynamic_cast<TTree*>(file.Get("mUPCTree"));
  Long64_t entries=reader->chain()->GetEntries();
  require(tree && tree->GetEntries()==entries,"preserve input events including no-primary events");
  StUPCEvent *event=0;
  tree->SetBranchAddress("mUPCEvent",&event);
  Long64_t nMC=0,nMCV=0,nRC=0,nRCV=0,nMatches=0,nVertexMatches=0,nNoReco=0,nDSM=0,nNonzeroDSM=0;
  for(Long64_t entry=0; entry<entries; ++entry) {
    chain->Clear();
    require(chain->Make()==0 && tree->GetEntry(entry)>0,"read corresponding events");
    StMuDst *mu=reader->muDst();
    require(event->getEventNumber()==mu->event()->eventNumber() && event->getRunNumber()==mu->event()->runNumber(),"event alignment");
    require(event->getMCEvent() && event->getMCEvent()->getRunNumber()==event->getRunNumber()
            && event->getMCEvent()->getEventNumber()==event->getEventNumber()
            && event->getMCEvent()->getSimulationType()==2,"MC/reco event alignment");
    const StTriggerData *trigger=mu->event()->triggerData();
    StUPCMCEvent *mcEvent=event->getMCEvent();
    require(mcEvent->hasEventInfo() && sameEventInfo(*mcEvent,*event),"all MC event metadata matches reco event");
    if(mcEvent->getLastDSM0() || mcEvent->getLastDSM1() || mcEvent->getLastDSM3()) ++nNonzeroDSM;
    require(mcEvent->hasDSMData()==(trigger!=0),"MC DSM availability matches input trigger data");
    require(mcEvent->getLastDSM0()==(trigger ? trigger->lastDSM(0) : 0)
            && mcEvent->getLastDSM1()==(trigger ? trigger->lastDSM(1) : 0)
            && mcEvent->getLastDSM3()==(trigger ? trigger->lastDSM(3) : 0),"MC DSM words match MuDst");
    if(trigger) {
      ++nDSM;
      require(mcEvent->getLastDSM0()==event->getLastDSM0() && mcEvent->getLastDSM1()==event->getLastDSM1()
              && mcEvent->getLastDSM3()==event->getLastDSM3(),"MC and reco DSM words agree");
    }
    TClonesArray *tracks=mu->mcArray(1), *vertices=mu->mcArray(0);
    require(event->getNumberOfMCTracks()==tracks->GetEntriesFast(),"preserve every MC track");
    require(event->getNumberOfMCVertices()==vertices->GetEntriesFast(),"preserve every MC vertex");
    std::map<int,int> trackIds,vertexIds;
    for(int i=0;i<tracks->GetEntriesFast();++i) trackIds[((StMuMcTrack*)tracks->At(i))->Id()]=i;
    for(int i=0;i<vertices->GetEntriesFast();++i) vertexIds[((StMuMcVertex*)vertices->At(i))->Id()]=i;
    for(int i=0;i<vertices->GetEntriesFast();++i) {
      StMuMcVertex *src=(StMuMcVertex*)vertices->At(i);
      StUPCMCVertex *dst=event->getMCVertex(i);
      require(dst && dst->getId()==src->Id(),"MC vertex ID");
      require(dst->getEvent()==event->getMCEvent(),"MC vertex-to-event navigation");
      same(dst->getPosX(),src->XyzV().x(),"MC vertex x");
      same(dst->getPosY(),src->XyzV().y(),"MC vertex y");
      same(dst->getPosZ(),src->XyzV().z(),"MC vertex z");
      same(dst->getTime(),src->Time(),"MC vertex time");
      require(dst->getNDaughters()==src->NoDaughters() && dst->getIntermediateFlag()==src->IsIntermedate(),"MC vertex flags");
      require(dst->getParentTrackId()==src->IdParTrk(),"MC parent ID");
      require(dst->getParentTrackIndex()==indexOf(trackIds,src->IdParTrk()),"MC parent index");
      if(dst->getParentTrackIndex()>=0) require(dst->getParentTrack()->getId()==src->IdParTrk(),"MC parent navigation");
      ++nMCV;
    }
    for(int i=0;i<tracks->GetEntriesFast();++i) {
      StMuMcTrack *src=(StMuMcTrack*)tracks->At(i);
      StUPCMCTrack *dst=event->getMCTrack(i);
      require(dst && dst->getId()==src->Id(),"MC track ID");
      require(dst->getEvent()==event->getMCEvent() && dst->getEvent()->getRecoEvent()==event,"MC track-to-event navigation");
      require(dst->getGeantId()==src->GePid() && dst->getPdgCode()==TDatabasePDG::Instance()->ConvertGeant3ToPdg(src->GePid()),"MC particle identity");
      require(dst->getCharge()==src->Charge() && dst->getIsShower()==src->IsShower(),"MC charge/shower flag");
      same(dst->Px(),src->Pxyz().x(),"MC px"); same(dst->Py(),src->Pxyz().y(),"MC py"); same(dst->Pz(),src->Pxyz().z(),"MC pz");
      same(dst->Energy(),src->E(),"MC energy"); same(dst->getPt(),src->pT(),"MC pT"); same(dst->getEta(),src->Eta(),"MC eta");
      same(dst->getPtot(),src->Ptot(),"MC p magnitude"); same(dst->getRapidity(),src->Rapidity(),"MC rapidity");
      for(int d=0;d<StUPCMCTrack::kNumberOfDetectors;++d) require(dst->getNhits(d)==src->NoHits(d),"MC detector hit count");
      require(dst->getStartVertexId()==src->IdVx() && dst->getStopVertexId()==src->IdVxEnd() && dst->getIntermediateVertexId()==src->ItrmdVertex(),"MC vertex IDs");
      require(dst->getStartVertexIndex()==indexOf(vertexIds,src->IdVx()),"MC start index");
      require(dst->getStopVertexIndex()==indexOf(vertexIds,src->IdVxEnd()),"MC stop index");
      require(dst->getIntermediateVertexIndex()==indexOf(vertexIds,src->ItrmdVertex()),"MC intermediate index");
      if(dst->getStartVertex()) {
        same(dst->Vx(),dst->getStartVertex()->getPosX(),"production position");
        same(dst->Vy(),dst->getStartVertex()->getPosY(),"production y");
        same(dst->Vz(),dst->getStartVertex()->getPosZ(),"production z");
        same(dst->T(),dst->getStartVertex()->getTime(),"production time");
      }
      require(event->getMCParticle(i)==dst,"compatibility particle accessor");
      ++nMC;
    }
    int irc=0,ivrc=0;
    std::map<int,int> matchedTracks,matchedVertices;
    for(unsigned int iv=0;iv<mu->numberOfPrimaryVertices();++iv) {
      StMuDst::setVertexIndex(iv);
      TObjArray *primary=mu->primaryTracks();
      int n=0;
      for(int it=0;it<primary->GetEntriesFast();++it) if(primary->At(it)) ++n;
      if(n>4) continue; // unchanged production skim selection
      StMuPrimaryVertex *srcv=mu->primaryVertex();
      StUPCVertex *dstv=event->getVertex(ivrc++);
      require(dstv && dstv->getId()==iv,"saved primary vertex order");
      same(dstv->getPosX(),srcv->position().x(),"reco vertex x");
      same(dstv->getPosY(),srcv->position().y(),"reco vertex y");
      same(dstv->getPosZ(),srcv->position().z(),"reco vertex z");
      same(dstv->getErrX(),srcv->posError().x(),"reco vertex x error");
      same(dstv->getErrY(),srcv->posError().y(),"reco vertex y error");
      same(dstv->getErrZ(),srcv->posError().z(),"reco vertex z error");
      require(dstv->getNTracksUsed()==srcv->nTracksUsed() && dstv->getNPrimaryTracks()==n,"reco vertex track counts");
      require(dstv->getIdTruth()==srcv->idTruth() && dstv->getQaTruth()==srcv->qaTruth() && dstv->getIdParent()==srcv->idParent(),"reco vertex truth values");
      int vm=indexOf(vertexIds,srcv->idTruth());
      require(dstv->getMCVertexIndex()==vm,"reco vertex truth link");
      if(vm>=0) { require(dstv->getMCVertex()==event->getMCVertex(vm),"reco vertex navigation"); ++matchedVertices[vm]; ++nVertexMatches; }
      for(int it=0;it<primary->GetEntriesFast();++it) {
        StMuTrack *src=(StMuTrack*)primary->At(it);
        if(!src) continue;
        StUPCTrack *dst=event->getTrack(irc++);
        require(dst!=0,"saved primary track order");
        same(dst->getPt(),src->pt(),"reco pT"); same(dst->getEta(),src->eta(),"reco eta"); same(dst->getPhi(),src->phi(),"reco phi");
        require(dst->getCharge()==src->charge() && dst->getNhits()==src->nHits()
                && dst->getNhitsFit()==src->nHitsFit() && dst->getNhitsMax()==src->nHitsPoss()
                && dst->getNhitsDEdx()==src->nHitsDedx(),"reco charge and hit counts");
        require(dst->getFlag(StUPCTrack::kPrimary),"primary track flag");
        same(dst->getChi2(),src->chi2(),"reco chi2");
        same(dst->getDcaXY(),src->dcaGlobal().perp(),"reco DCA transverse");
        same(dst->getDcaZ(),src->dcaGlobal().z(),"reco DCA z");
        same(dst->getDEdxSignal(),src->dEdx(),"reco dEdx");
        samePID(dst->getNSigmasTPCElectron(),src->nSigmaElectron(),"reco electron PID");
        samePID(dst->getNSigmasTPCPion(),src->nSigmaPion(),"reco pion PID");
        samePID(dst->getNSigmasTPCKaon(),src->nSigmaKaon(),"reco kaon PID");
        samePID(dst->getNSigmasTPCProton(),src->nSigmaProton(),"reco proton PID");
        require(dst->getIdTruth()==src->idTruth() && dst->getQaTruth()==src->qaTruth() && dst->getIdParentVx()==src->idParentVx(),"primary truth fields");
        int tm=indexOf(trackIds,src->idTruth());
        require(dst->getMCTrackIndex()==tm,"primary truth index incl. background markers");
        require(dst->getVertex()==dstv,"primary-to-reco-vertex navigation");
        if(tm>=0) { require(dst->getMCTrack()==event->getMCTrack(tm),"primary MC navigation"); ++matchedTracks[tm]; ++nMatches; }
        else require(!dst->getMCTrack(),"unmatched primary navigation");
      }
    }
    require(irc==event->getNumberOfTracks() && ivrc==event->getNumberOfVertices(),"preserve reco skim selection");
    for(int i=0;i<event->getNumberOfMCTracks();++i) {
      StUPCMCTrack *mc=event->getMCEvent()->getTrack(i);
      require(mc->getNumberOfRecoTracks()==matchedTracks[i],"stored reverse track associations");
      for(int k=0;k<mc->getNumberOfRecoTracks();++k)
        require(mc->getRecoTrack(k) && mc->getRecoTrack(k)->getMCTrackIndex()==i,"reciprocal track association");
    }
    for(int i=0;i<event->getNumberOfMCVertices();++i) {
      StUPCMCVertex *mc=event->getMCEvent()->getVertex(i);
      require(mc->getNumberOfRecoVertices()==matchedVertices[i],"stored reverse vertex associations");
      for(int k=0;k<mc->getNumberOfRecoVertices();++k)
        require(mc->getRecoVertex(k) && mc->getRecoVertex(k)->getMCVertexIndex()==i,"reciprocal vertex association");
    }
    nRC+=irc; nRCV+=ivrc; if(!mu->numberOfPrimaryVertices()) ++nNoReco;
    if(efficiency.is_open()) writeMuEfficiency(mu,entry,efficiency);
  }
  tree->ResetBranchAddresses();
  delete event;
  file.Close();
  chain->Finish();
  delete chain;
  std::cout << "MC_VALIDATE_PASS events=" << entries << " mc_tracks=" << nMC << " mc_vertices=" << nMCV
            << " reco_tracks=" << nRC << " reco_vertices=" << nRCV << " track_matches=" << nMatches
            << " vertex_matches=" << nVertexMatches << " events_without_reco_vertex=" << nNoReco
            << " max_pid_float16_rounding=" << maxPidRounding << " events_with_dsm=" << nDSM << " events_with_nonzero_dsm=" << nNonzeroDSM << std::endl;
}

void testLegacyMC(const char *input, bool withDSM=false) {
  TFile file(input);
  TTree *tree=dynamic_cast<TTree*>(file.Get("mUPCTree"));
  require(tree && tree->GetEntries()>0,"open pre-DSM MC event fixture");
  StUPCEvent *event=0;
  tree->SetBranchAddress("mUPCEvent",&event);
  for(Long64_t i=0;i<tree->GetEntries();++i) {
    require(tree->GetEntry(i)>0,"read pre-DSM entry");
    StUPCMCEvent *mc=event->getMCEvent();
    require(mc && !mc->hasEventInfo(),"old MC schema has no full metadata snapshot");
    require(mc->getFillNumber()==0 && mc->getTriggerIds().GetSize()==0 && mc->getNumberOfRecoTracks()==0,"new metadata defaults in old MC files");
    require(mc->hasDSMData()==withDSM,"old DSM availability preserved");
    require(mc->getLastDSM0()==(withDSM ? event->getLastDSM0() : 0)
            && mc->getLastDSM1()==(withDSM ? event->getLastDSM1() : 0)
            && mc->getLastDSM3()==(withDSM ? event->getLastDSM3() : 0),"old DSM words preserved");
  }
  tree->ResetBranchAddresses();
  delete event;
  std::cout << "MC_LEGACY_DSM_PASS events=" << tree->GetEntries() << std::endl;
}
}

int TestMCTruth(const char *mode, const char *input, const char *output, const char *legacy) {
  try {
    if(std::string(mode)=="io") testIO(output,legacy);
    else if(std::string(mode)=="convert") convert(input,output);
    else if(std::string(mode)=="validate") validate(input,output,legacy);
    else if(std::string(mode)=="legacy_metadata") testLegacyMC(input,true);
    else if(std::string(mode)=="legacy_mc") testLegacyMC(input);
    else throw std::runtime_error("unknown test mode");
    return 0;
  } catch(const std::exception& error) {
    std::cerr << "MC_TEST_FAIL: " << error.what() << std::endl;
    return 1;
  }
}
