#include "MLAnalyzer/RecHitAnalyzer/interface/RecHitAnalyzer.h"

const unsigned nJets = 2;
TH1D *h_jet_pT;
TH1D *h_jet_E;
TH1D *h_jet_eta;
TH1D *h_jet_m0;
TH1D *h_jet_nJet;
TH1D *h_nGG;
TH1D *h_nQQ;
float vJet_m0_[nJets];
float vJet_pt_[nJets];
float vJetIsQuark_[nJets];
float vJetIds_[nJets];
std::vector<float> vSubJetE_[nJets];
std::vector<float> vSubJetPx_[nJets];
std::vector<float> vSubJetPy_[nJets];
std::vector<float> vSubJetPz_[nJets];
std::vector<int> vJetIds;

// Initialize branches _____________________________________________________//
void RecHitAnalyzer::branchesEvtSel_dijet_gg_qq ( TTree* tree, edm::Service<TFileService> &fs ) {

  h_jet_pT    = fs->make<TH1D>("h_jet_pT"  , "p_{T};p_{T};Particles", 100,  0., 500.);
  h_jet_E     = fs->make<TH1D>("h_jet_E"   , "E;E;Particles"        , 100,  0., 800.);
  h_jet_eta   = fs->make<TH1D>("h_jet_eta" , "#eta;#eta;Particles"  , 100, -5., 5.);
  h_jet_nJet  = fs->make<TH1D>("h_jet_nJet", "nJet;nJet;Events"     ,  10,  0., 10.);
  h_jet_m0    = fs->make<TH1D>("h_jet_m0"  , "m0;m0;Events"         ,  60, 50., 110.);
  h_nGG       = fs->make<TH1D>("h_nGG"     , "nGG;nGG;Events"       ,   3,  0.,   3.);
  h_nQQ       = fs->make<TH1D>("h_nQQ"     , "nQQ;nQQ;Events"       ,   3,  0.,   3.);

  char hname[50];
  for ( unsigned iJ = 0; iJ != nJets; iJ++ ) {
    sprintf(hname, "jetM%d",  iJ);
    tree->Branch(hname,            &vJet_m0_[iJ]);
    sprintf(hname, "jetPt%d", iJ);
    tree->Branch(hname,            &vJet_pt_[iJ]);
    sprintf(hname, "jetIsQuark%d", iJ);
    tree->Branch(hname,            &vJetIsQuark_[iJ]);
    sprintf(hname, "jetIds%d", iJ);
    tree->Branch(hname,            &vJetIds_[iJ]);

    sprintf(hname, "subJet%d_E", iJ);
    tree->Branch(hname,            &vSubJetE_[iJ]);
    sprintf(hname, "subJet%d_Px", iJ);
    tree->Branch(hname,            &vSubJetPx_[iJ]);
    sprintf(hname, "subJet%d_Py", iJ);
    tree->Branch(hname,            &vSubJetPy_[iJ]);
    sprintf(hname, "subJet%d_Pz", iJ);
    tree->Branch(hname,            &vSubJetPz_[iJ]);
  }

} // branchesEvtSel_jet()

bool RecHitAnalyzer::runEvtSel_dijet_gg_qq( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{

  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken( genParticleCollectionT_, genParticles );

  edm::Handle<reco::PFJetCollection> jets;
  iEvent.getByToken(jetCollectionT_, jets);
  float dR;

  vGoodJetIdxs.clear();
  vJetIds.clear();
  /*
  edm::Handle<reco::GenJetCollection> genJets;
  iEvent.getByToken(genJetCollectionT_, genJets);
  std::vector<float> vJetFakePhoIdxs;
  */

  unsigned nGG = 0;
  unsigned nQQ = 0;
  //int nQcdIn = 0;
  //int nQcdOut = 0;

  if ( debug ) std::cout << " >>>>>>>>>>>>>>>>>>>> evt:" << std::endl;
  for (reco::GenParticleCollection::const_iterator iGen = genParticles->begin();
      iGen != genParticles->end();
      ++iGen) {

    if ( iGen->numberOfMothers() != 2 ) continue;
    //if ( iGen->status() != 3 ) continue;
    if ( iGen->status() != 23 ) continue;
    if ( debug ) std::cout << " >> id:" << iGen->pdgId() << " status:" << iGen->status() << " nDaught:" << iGen->numberOfDaughters() << " pt:"<< iGen->pt() << " eta:" <<iGen->eta() << " phi:" <<iGen->phi() << " nMoms:" <<iGen->numberOfMothers()<< std::endl;

    // Loop over jets
    for ( unsigned iJ(0); iJ != jets->size(); ++iJ ) {
      //if ( debug ) std::cout << " >>>>>> jet[" << iJ << "]" << std::endl;
      reco::PFJetRef iJet( jets, iJ );
      if ( std::abs(iJet->pt())  < 70. ) continue;
      if ( std::abs(iJet->eta()) > 2. ) continue;
      dR = reco::deltaR( iJet->eta(),iJet->phi(), iGen->eta(),iGen->phi() );
      if ( debug ) {
        std::cout << " >>>>>> jet[" << iJ << "] Pt:" << iJet->pt() << " jetEta:" << iJet->eta() << " jetPhi:" << iJet->phi()
        << " dR:" << dR << std::endl;
      }
      if ( dR > 0.4 ) continue;
      vGoodJetIdxs.push_back( iJ );
      vJetIds.push_back( std::abs(iGen->pdgId()) );
      if ( debug ) {
        std::cout << " >>>>>> DR matched: jet[" << iJ << "] pdgId:" << std::abs(iGen->pdgId()) << std::endl;
      }
      break;
    } // reco jets

  } // gen particles
  h_jet_nJet->Fill( vGoodJetIdxs.size() );
  if ( vGoodJetIdxs.size() != nJets ) return false;

  for ( unsigned iJ(0); iJ != nJets; ++iJ ) {
    if ( vJetIds[iJ] < 4 ) nQQ++;
    else if ( vJetIds[iJ] == 21 ) nGG++;
  }

  if ( vGoodJetIdxs[0] == vGoodJetIdxs[1] ) return false; // protect against double counting: only valid for nJets==2
  if ( nQQ+nGG != nJets ) return false; // require dijet
  if ( nQQ != nJets && nGG != nJets ) return false; // require gg or qq final state
  //if ( nQQ != 1 && nGG != 1 ) return false; // require qg final state

  for ( unsigned iJ(0); iJ != nJets; ++iJ ) {
    reco::PFJetRef iJet( jets, vGoodJetIdxs[iJ] );
    // Plot kinematics
    h_jet_pT->Fill( std::abs(iJet->pt()) );
    h_jet_eta->Fill( iJet->eta() );
    h_jet_E->Fill( iJet->energy() );
    h_jet_m0->Fill( iJet->mass() );

    // Write out event ID
    vJet_pt_[iJ] = iJet->pt();
    vJet_m0_[iJ] = iJet->mass();
    vJetIsQuark_[iJ] = vJetIds[iJ] < 4 ? true : false;
    vJetIds_[iJ] = vJetIds[iJ];

    vSubJetE_[iJ].clear();
    vSubJetPx_[iJ].clear();
    vSubJetPy_[iJ].clear();
    vSubJetPz_[iJ].clear();
    //std::vector<reco::PFCandidatePtr> jetConstituents = iJet->getPFConstituents();
    unsigned int nConstituents = iJet->getPFConstituents().size();
    for ( unsigned int j = 0; j < nConstituents; j++ ) {
      const reco::PFCandidatePtr subJet = iJet->getPFConstituent( j );
      std::cout << " >> " << j << ": E:" << subJet->energy() << " px:" << subJet->px() << " py:" << subJet->py() << " pz:" << subJet->pz() << std::endl;
      vSubJetE_[iJ].push_back( subJet->energy() );
      vSubJetPx_[iJ].push_back( subJet->px() );
      vSubJetPy_[iJ].push_back( subJet->py() );
      vSubJetPz_[iJ].push_back( subJet->pz() );
    }
  }

  h_nGG->Fill( nGG );
  h_nQQ->Fill( nQQ );
    /*
    for ( unsigned int k = 0; k < genJets->size(); k++ ) {
      reco::GenJetRef iJet( genJets, k );
      if ( std::abs(iJet->pt()) < 50. ) continue;
      dR = reco::deltaR( iJet->eta(),iJet->phi(), iGen->eta(),iGen->phi() );
      if ( debug ) std::cout << " >>>>>> genJet[" << k << "]Pt:" << iJet->pt() << " jetEta:" << iJet->eta() << " jetPhi:" << iJet->phi() << " dR:" << dR << std::endl;
    } // genJets
    */
    /*
    // Get incoming QCD partons
    for ( unsigned int i = 0; i < iGen->numberOfDaughters(); i++ ) {
      const reco::Candidate* ppParton = iGen->daughter(i);
      if ( ppParton->status() != 3 ) continue;
      if ( debug ) std::cout << " >>>> id:" << ppParton->pdgId() << " status:" << ppParton->status() << " nDaught:" << ppParton->numberOfDaughters() << " pt:" << ppParton->pt() << " eta:" <<ppParton->eta() << " phi:" <<ppParton->phi() <<std::endl;
      // Get outgoing QCD partons
      for ( unsigned int j = 0; j < ppParton->numberOfDaughters(); j++ ) {
        const reco::Candidate* qcdIn = ppParton->daughter(j);
        if ( qcdIn->status() != 3 ) continue;
        if ( debug ) std::cout << " >>>>>> id:" << qcdIn->pdgId() << " status:" << qcdIn->status() << " nDaught:" << qcdIn->numberOfDaughters() << " pt:" << qcdIn->pt() << " eta:" <<qcdIn->eta() << " phi:" <<qcdIn->phi() << std::endl;
        nQcdIn++;
        for ( unsigned int k = 0; k < qcdIn->numberOfDaughters(); k++ ) {
          const reco::Candidate* qcdOut = qcdIn->daughter(k);
          if ( qcdOut->status() != 3 ) continue;
          if ( debug ) std::cout << " >>>>>>>> id:" << qcdOut->pdgId() << " status:" << qcdOut->status() << " nDaught:" << qcdOut->numberOfDaughters() << " pt:" << qcdOut->pt() << " eta:" <<qcdOut->eta() << " phi:" <<qcdOut->phi() << std::endl;
          //nQcdOut++;
          if ( std::abs(qcdOut->pdgId()) != 21 ) continue;
          nGG++;

        // Loop over jets
        for ( unsigned iJ(0); iJ != jets->size(); ++iJ ) {

          //if ( debug ) std::cout << " >>>>>> jet[" << iJ << "]" << std::endl;
          reco::PFJetRef iJet( jets, iJ );
          if ( std::abs(iJet->pt()) < 50. ) continue;
          //if ( std::abs(iJet->eta()) > 2. ) continue;
          dR = reco::deltaR( iJet->eta(),iJet->phi(), qcdOut->eta(),qcdOut->phi() );
          if ( debug ) {
            std::cout << " >>>>>> jet[" << iJ << "]Pt:" << iJet->pt() << " jetEta:" << iJet->eta() << " jetPhi:" << iJet->phi()
            << " dR:" << dR << std::endl;
          }
        } // PFjets
        for ( unsigned int k = 0; k < genJets->size(); k++ ) {
          reco::GenJetRef iJet( genJets, k );
          if ( std::abs(iJet->pt()) < 50. ) continue;
          dR = reco::deltaR( iJet->eta(),iJet->phi(), qcdOut->eta(),qcdOut->phi() );
          if ( debug ) std::cout << " >>>>>> genJet[" << k << "]Pt:" << iJet->pt() << " jetEta:" << iJet->eta() << " jetPhi:" << iJet->phi() << " dR:" << dR << std::endl;
        } // genJets

        } // final
      } // outgoing partons
    } // incoming partons

  } // gen particles
    */

  /*
  if ( debug ) std::cout << " >> genJets.size: " << genJets->size() << std::endl;
  for ( unsigned int i = 0; i < genJets->size(); i++ ) {
    reco::GenJetRef iJet( genJets, i );
    if ( debug ) std::cout << " >> genJet[" << i << "]Pt:" << iJet->pt() << " jetEta:" << iJet->eta() << " jetPhi:" << iJet->phi() << std::endl;
    unsigned int nConstituents = iJet->getGenConstituents().size();
    for ( unsigned int j = 0; j < nConstituents; j++ ) {
      const reco::GenParticle* subJet = iJet->getGenConstituent( j );
      //if ( std::abs(subJet->pdgId()) != 21 ) continue;
      dR = reco::deltaR( iJet->eta(),iJet->phi(), subJet->eta(),subJet->phi() );
      if ( debug ) {
        std::cout << " >> pdgId:" << subJet->pdgId() << " status:" << subJet->status()
        << " pT:" << subJet->pt() << " eta:" << subJet->eta() << " phi: " << subJet->phi() << " E:" << subJet->energy();
        std::cout << " moth:" << subJet->mother()->pdgId();
        std::cout << " dR: " << dR << std::endl;
      }
    }
  }
  if ( debug ) std::cout << " >> nQcdIn:" << nQcdIn << " nQcdOut:" << nQcdOut << " nGG:" << nGG << std::endl;
  if ( nGG != 2 ) return false;
  if ( nQcdIn !=2 || nQcdOut !=2 ) return false;
  */

  if ( debug ) std::cout << " >> has_dijet: passed" << std::endl;
  return true;

} // has_dijet
