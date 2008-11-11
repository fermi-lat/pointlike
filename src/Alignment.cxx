/** @file Alignment.cxx
@brief implementation of Alignment

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/src/Alignment.cxx,v 1.3 2008/09/10 21:49:49 burnett Exp $

*/


#include "pointlike/Alignment.h"
#include <stdexcept>

using namespace pointlike;

namespace{
// wired-in from Marshall's data for now.
    struct { double t; double x; double y; double z;}
    data[]=
    {
        {0,             0,     0,     0},
        {236511638,  -183,  -169,  -497}, // first light stuff 
                                             // afterward
        {242053458,    19,    50,    8}, // correction now made in data
        {999999999,     0,     0,     0}  // end of time
    };

}

Alignment::Alignment()
: m_active(false)
{}

Alignment::Alignment(const std::string& calibration)
: m_active(true)
, m_start(0)
, m_end(0)
{
    if(calibration != "default"){
        throw std::invalid_argument("Alignment data file not implemented yet");
        /// @TODO: implement this to 
    }
}
  
Alignment::Alignment(double x, double y, double z)
: m_start(0), m_end(999999999)
{
    set_rot(x,y,z);
}

CLHEP::Hep3Vector Alignment::rotate(const CLHEP::Hep3Vector& in, double time)const
{
    if( !active()) return in;
    return rotation(time)*in;

}
void Alignment::set_rot(double arcsecx,double arcsecy,double arcsecz) const
{
    //N.B. the signs are reversed.
    static double c(M_PI/648000);
    m_rot = CLHEP::HepRotationX(-arcsecx*c)
           *CLHEP::HepRotationY(-arcsecy*c)
           *CLHEP::HepRotationZ(-arcsecz*c);

}
const CLHEP::HepRotation& Alignment::rotation(double time)const
{
    if( !in_range(time, m_start, m_end)) {
        // update
        int i =0; 
        while( !in_range(time, data[i].t, data[i+1].t) ){++i;}
        m_start=data[i].t;
        m_end = data[i+1].t;
        set_rot(data[i].x, data[i].y, data[i].z);
        
    }
    return m_rot; // note that the constants are misalignment
}
