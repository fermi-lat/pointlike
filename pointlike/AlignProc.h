/** @file AlignProc.h
@brief definition of AlignProc

*/

#ifndef pointlike_AlignProc_h
#define pointlike_AlignProc_h

#include "astro/Photon.h"
#include "pointlike/RotationInfo.h"
#include "pointlike/PointSourceLikelihood.h"
#include <vector>


namespace pointlike{

    //! class AlignProc
    //! Tries to characterize misalignment of GLAST coordinates by taking pointing information along with
    //! known bright source positions. The alignment is done with a likelihood analysis over a grid of rotations about each axis,
    //! followed by a least squares quadratic fit for the likelihood surface
    class AlignProc{
    public:

    //! class Photon
    //! Wrapper class for astro::Photon with added pointing information
        class Photona : public astro::Photon {
        public:
            Photona(const astro::SkyDir& dir, double energy, 
                double time, int event_class, 
                const astro::SkyDir& scz, 
                const astro::SkyDir& scx)
                : astro::Photon(dir, energy, time, event_class, 0)
            {
                CLHEP::Hep3Vector scy (scz().cross(scx()));
                CLHEP::Hep3Vector sz = scz();
                CLHEP::Hep3Vector sx = scx();
                //unused double dot = scz().dot(scx());
                m_rot = CLHEP::HepRotation(scx(), scy, scz());
                m_gdir = dir();
            }
            CLHEP::Hep3Vector GDir() {return m_gdir;}
            CLHEP::HepRotation Rot() {return m_rot;}

        private:
            CLHEP::HepRotation m_rot;
            CLHEP::Hep3Vector m_gdir;
        };


        //! constructor takes an array of source positions and creates an grid of rotations
        //! @param sources an array of source positions
        //! @param files an array of either root merit or FT1 fits files 
        //! @param arcmins determines the resolution of the grid about the identity rotation (0,0,0)
        AlignProc(std::vector<astro::SkyDir>& sources, std::vector<std::string>& files, double arcmins, double offx, double offy, double offz, int start=-1, int stop=-1); 

        //! takes event information and returns a photon with pointing information
        pointlike::AlignProc::Photona event(std::vector<float>& row);

        //! accumulates likelihood information from each alignment photon (contains pointing information)
        int add(pointlike::AlignProc::Photona& p);

        //! returns the loglikelihood of one of the 27 grid points
        //! @param x +1,0,-1
        //! @param y +1,0,-1
        //! @param z +1,0,-1
        long double loglikelihood(int x, int y, int z) {return m_roti.likelihood(x,y,z);}

        //! returns the fit parameters for the likelihood surface
        //! surface is presumed to be a 3-D paraboloid
        //! Loglike(x,y,z) = a[0]*x**2+a[1]*x+a[2]*y**2+a[3]*y+a[4]*z**2+a[5]*z+a[6]+a[7]*x*y+a[8]*x*z+a[9]*y*z
        //! where 'a' is the returned array
        //! curvaturex = (-2*a[0])**(-0.5) etc
        std::vector<double> fitparam();

        //! Rotation for performing iterative correction to misalignment
        static HepRotation s_hr;

        //! Corrects previous rotation by adding in new rotation
        //! @param x rotation about x-axis in radians
        //! @param y rotation about y-axis in radians
        //! @param z rotation about z-axis in radians
        static void addRot(double x,double y, double z) { HepRotationX mx(x); HepRotationY my(y); HepRotationZ mz(z); s_hr = mx*my*mz*s_hr;}

        //maximum scaled deviation
        static int s_umax;
        static int setumax(int umax) { int t = s_umax; s_umax = umax; return t;}
        static double s_scale;
        static double setscale(double scale) {double t = s_scale; s_scale = scale; return t;}

    private:
        void loadroot(const std::string& file);
        void loadfits(const std::string& file);
        int m_photons; //number of photons
        int m_start; //start time 0... (-1 for all)
        int m_stop; //end time ...stop
        double m_arcmin; //resolution of the grid
        RotationInfo m_roti; //contains likelihood data for each rotation grid point
        const std::vector<astro::SkyDir> m_skydir; //source positions
    };
}
#endif
