/** @file Alignment.h 
    @brief declaration of the Alignment wrapper class

    $Header: /nfs/slac/g/glast/ground/cvs/pointlike/pointlike/Alignment.h,v 1.29 2008/07/22 03:58:10 burnett Exp $
*/


#ifndef pointlike_Alignment_h
#define pointlike_Alignment_h

#include "CLHEP/Vector/Rotation.h"

#include <map>
#include <vector>

namespace pointlike{

    class Alignment {

    public:
        /// @brief default ctor does nothing
        Alignment();

        /// @brief ctor loads the data. default uses wired-in
        Alignment(const std::string& calibration);

        /// @brief ctor with fixed values, to be always valid
        Alignment(double x, double y, double z);
        
        /// @brief the appropriate rotation matrix for the given time
        const CLHEP::HepRotation& rotation(double time)const;

        /// @brief apply the saved rotation
        CLHEP::Hep3Vector rotate(const CLHEP::Hep3Vector& in, double time)const;

        /// @brief if the alignment has been set
        bool active()const{return m_active;}

    private:
        bool m_active;
        mutable double m_start, m_end; ///< current start, end times
        mutable CLHEP::HepRotation m_rot; ///< current rotation
        std::map<double, std::vector<double> > m_calibdata;

        void set_rot(double arcsecx,double arcsecy,double arcsecz)const ;
        bool in_range(double t, double a, double b)const {return t>=a && t<b;}
    };


} // pointlike


#endif
