/** @file PhotonList.h
@brief definition of PhotonList

*/

#ifndef tools_PhotonList_h
#define tools_PhotonList_h

#include "astro/Photon.h"
#include <string>

// ROOT forward declarations
class TTree;
class TFile;
class TLeaf;



    // forward declarations 
    class PhotonListIterator;
    class RootLeaf;
/**
    @class PhotonList
    @container that adapts a root file containing photon data, returned as astro::Photon objects 

    Usage:
*/
class PhotonList
{
public:

    /** @brief ctor sets up container
        @param infile input file name
        @param table_name  if null, take first found
    */
    PhotonList(const std::string infile, 
        const std::string table_name="");

    astro::Photon operator[](int index)const;

    size_t size() const;

    PhotonListIterator begin();
    PhotonListIterator end();

private:
    // ROOT stuff
    RootLeaf * getLeaf(std::string name);
    RootLeaf * getLeaf(std::string name, double def);
    TFile* m_file;
    TTree* m_tree;
    RootLeaf* m_leaf_ra, * m_leaf_dec, * m_leaf_energy;
    RootLeaf* m_leaf_event_class;
    RootLeaf* m_leaf_layer;
    RootLeaf* m_leaf_source;
};


// make it a container by implementing a forward iterator
class PhotonListIterator {
public:
    PhotonListIterator(const PhotonList& list, long index=0)
        :m_list(list), m_index(index){}
        astro::Photon operator*()const;             ///< dereference
        operator int()const{return m_index;} ///< allows comparison
        long operator++(){return ++m_index;} ///< increment operator
private:
    const PhotonList& m_list;
    long m_index;
};


// wrap a ROOT TLeaf, allow simple default if it does not exist
class RootLeaf {
public:
    RootLeaf(TLeaf* leaf):m_leaf(leaf){};
    RootLeaf(double defaul):m_leaf(0), m_default(defaul){};
    operator double()const; ///< return the value, or default
private:
    TLeaf* m_leaf;
    double m_default;
};


#endif
