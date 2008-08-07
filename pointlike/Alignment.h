/** @file Alignment.h 
    @brief declaration of the Alignment calibration management class

    $Header: /nfs/slac/g/glast/ground/cvs/pointlike/pointlike/Alignment.h,v 1.2 2008/07/29 19:13:59 burnett Exp $
*/


#ifndef pointlike_Alignment_h
#define pointlike_Alignment_h

#include "CLHEP/Vector/Rotation.h"


namespace pointlike{

    /** @class Alignment
        @brief manage LAT/SC alignment constants


    */
    class Alignment {

    public:
        /// @brief default ctor does nothing, defines default identity transform
        Alignment();

        /// @brief ctor loads the data. set name to "default" to uses wired-in set
        /// 
        /// The format of the file containing the data is four columns: MET,x,y,z,
        /// @param MET start time for valididty. end time is the next entry.
        /// @param x,y,z Misalignment values in arcsec.

        Alignment(const std::string& calibration);

        /// @brief ctor with fixed values, to be always valid
        /// @param x,y,z Values in arcsec for successive rotations about x,y,z axes. 
        ///  These are the observed misalignment
        Alignment(double x, double y, double z);
        
        /// @brief the appropriate rotation matrix for the given time
        const CLHEP::HepRotation& rotation(double time)const;

        /// @brief apply the saved rotation
        /// @param in  direction unit vector in instrument coords
        /// @param out MET time
        CLHEP::Hep3Vector rotate(const CLHEP::Hep3Vector& in, double time)const;

        /// @brief if the alignment has been set
        bool active()const{return m_active;}

    private:
        bool m_active;
        mutable double m_start, m_end; ///< current start, end times
        mutable CLHEP::HepRotation m_rot; ///< current rotation

        void set_rot(double arcsecx,double arcsecy,double arcsecz)const ;
        bool in_range(double t, double a, double b)const {return t>=a && t<b;}
    };


} // pointlike


#endif
