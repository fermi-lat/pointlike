/** @file PhotonList.h
@brief implementation of PhotonList

*/
    

#include "PhotonList.h"
#include "TTree.h"
#include "TFile.h"
#include "TLeaf.h"
#include "TIterator.h"
#include "TString.h"
#include "TKey.h"
#include "TSystem.h"

#include <stdexcept>



namespace {
#ifdef WIN32
// Prevent problems with Root's dynamic loader in 4.0x.yy by instantiating a tree.
// This may not be necessary after Root 4.02.00.
TTree g_PhotonList_windows_dynamic_loader_bug_4_02_00;
#endif

void resetSigHandlers() {
    gSystem->ResetSignal(kSigBus);
    gSystem->ResetSignal(kSigSegmentationViolation);
    gSystem->ResetSignal(kSigSystem);
    gSystem->ResetSignal(kSigPipe);
    gSystem->ResetSignal(kSigIllegalInstruction);
    gSystem->ResetSignal(kSigQuit);
    gSystem->ResetSignal(kSigInterrupt);
    gSystem->ResetSignal(kSigWindowChanged);
    gSystem->ResetSignal(kSigAlarm);
    gSystem->ResetSignal(kSigChild);
    gSystem->ResetSignal(kSigUrgent);
    gSystem->ResetSignal(kSigFloatingException);
    gSystem->ResetSignal(kSigTermination);
    gSystem->ResetSignal(kSigUser1);
    gSystem->ResetSignal(kSigUser2);
  }

}

PhotonList::PhotonList(const std::string infile, 
                       const std::string treename)
{

    resetSigHandlers();

    m_file = new TFile(infile.c_str(), "readonly");
    if( !m_file->IsOpen()) throw std::invalid_argument("Could not open file "+infile);


    if( ! treename.empty() ){
        m_tree = dynamic_cast<TTree*>( m_file->Get(treename.c_str()));
    }else{
        TIter nextTopLevelKey(m_file->GetListOfKeys());
        TKey *key;

        // loop on keys, get the first top-level TTree 
        while  ( (key=(TKey*)nextTopLevelKey()) ) {
            TString className(key->GetClassName());
            if( className.CompareTo("TTree")==0 )  {
                // Found It
                m_tree = dynamic_cast<TTree*>(m_file->Get(key->GetName()));
            }
        }
    }
    if( m_tree==0) {
        throw std::invalid_argument(std::string("Could not find tree ")+treename);
    }

    int n=0;
    static std::string names[]={"ra", "dec", "energy", "event_class" ,"layer", "source"};
    m_leaf_ra = getLeaf(names[n++]);
    m_leaf_dec = getLeaf(names[n++]);
    m_leaf_energy = getLeaf(names[n++]);

    m_leaf_event_class = getLeaf(names[n++], -1);// optional, default -1

    m_leaf_layer = getLeaf(names[n++], -1);
    m_leaf_source = getLeaf(names[n++], -1);
}
RootLeaf* PhotonList::getLeaf(std::string name)
{
    TLeaf * leaf = m_tree->GetLeaf(name.c_str());
    if( leaf==0 ){
        throw std::invalid_argument("PhotonList::getLeaf: Did not find leaf "+name);
    }
    return new RootLeaf(leaf);
}
RootLeaf* PhotonList::getLeaf(std::string name, double defaul)
{
    try{ 
        return getLeaf(name);
    }catch(...){ 
        return new RootLeaf(defaul); 
    }
}

RootLeaf::operator double()const
{
    if( m_leaf==0) return m_default;
    return m_leaf->GetValue();
}

//~~~~~~~~~~~~~~~~~~~~~~~ 

size_t PhotonList::size()const{ return m_tree->GetEntries();}

skymaps::Photon PhotonList::operator[](int index)const
{
    m_tree->GetEntry(index);
    return skymaps::Photon(astro::SkyDir(*m_leaf_ra, *m_leaf_dec),
        *m_leaf_energy,*m_leaf_event_class);
}

skymaps::Photon PhotonListIterator::operator*()const
{
    // note: possibility of selection here: test, then increment index  
    return m_list[m_index];
}


PhotonListIterator PhotonList::begin()
{ 
    return PhotonListIterator(*this, 0);
}

PhotonListIterator PhotonList::end()
{
    return PhotonListIterator(*this, size());
}

