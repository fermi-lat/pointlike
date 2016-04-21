/** @file AlignProc.h
@brief definition of AlignProc

*/

#ifndef pointlike_AlignProc_h
#define pointlike_AlignProc_h

#include "pointlike/RotationInfo.h"
#include "pointlike/PointSourceLikelihood.h"
#include "skymaps/PhotonBinner.h"
#include "TMatrixD.h"
#include <vector>

namespace astro {
    class PointingHistory;
}

namespace pointlike{

    //! class AlignProc
    //! Tries to characterize misalignment of GLAST coordinates by taking pointing information along with
    //! known bright source positions. The alignment is done with a likelihood analysis over a grid of rotations about each axis,
    //! followed by a least squares quadratic fit for the likelihood surface
    class AlignProc{
    public:

        //! class LikeSurface
        //! Functor to the likelihood surface L(x,y,z) where x,y,z are the rotations in arcseconds of RotX(x)*RotY(y)*RotZ(Z)
        
        class LikeSurface {
        public:
            LikeSurface(std::vector<double>& fit): m_params(fit){}
            
            //returns the ith element of the parameter array where the parameter array is described below
            //Loglike(x,y,z) = a[0]*x**2+a[1]*x+a[2]*y**2+a[3]*y+a[4]*z**2+a[5]*z+a[6]+a[7]*x*y+a[8]*x*z+a[9]*y*z
            double operator()(int i) {return m_params[i];}

            //returns the log-likelihood of the fit surface with the rotations in arcseconds x,y,z
            double operator()(double x,double y, double z) {return m_params[0]*x*x+m_params[1]*x+m_params[2]*y*y+m_params[3]*y
                +m_params[4]*z*z+m_params[5]*z+m_params[6]+m_params[7]*x*y+m_params[8]*x*z+m_params[9]*y*z;}

            //returns the curvature of the likelihood surface to determine the deviation
            //sigma = (-2*curv)**-0.5
            std::vector<double> curvature() {
                TMatrixD fij(3,3);
                Double_t *det=0;
                std::vector<double> curv;
                fij[0][0]=2*m_params[0];
                fij[1][1]=2*m_params[2];
                fij[2][2]=2*m_params[4];
                fij[0][1]=m_params[7];
                fij[0][2]=m_params[8];
                fij[1][2]=m_params[9];
                fij[1][0]=fij[0][1];
                fij[2][0]=fij[0][2];
                fij[2][1]=fij[1][2];
                fij=fij.Invert(det);
                for(int i(0);i<3;++i) {
                    double temp = fij[i][i];
                    curv.push_back(sqrt(fabs(temp)));
                }

                return curv;
            }

        private:
            std::vector<double> m_params;
        };

        //! class Photon
        //! Wrapper class for astro::Photon with added pointing information
		// Inherit from skymaps::Photon for compatibility with new PhotonBinner - EEW
        class Photona : public skymaps::Photon {
        public:
            Photona(const astro::SkyDir& dir, double energy, 
                double time, int event_class, 
                const astro::SkyDir& scz, 
                const astro::SkyDir& scx,
                int id=0)
                : skymaps::Photon(dir, energy, time, event_class, id)
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

            /// make transformation in GLAST frame
            /// @param corr matrix that corrects direction in the GLAST frame
            skymaps::Photon transform(const CLHEP::HepRotation& corr)const
            {
                CLHEP::Hep3Vector 
                    local( m_rot.inverse()*dir()),
                    transformed( m_rot * corr * local );

                // account for worse energy resolution of back events by simply scaling up the energy.
                int evtclass(eventClass());
                //assert( evtclass>=0 && evtclass<3); // just a check: should be 0 or 1 if not DC2
                return skymaps::Photon(SkyDir(transformed), energy(),time(),evtclass, source());
            }

        private:
            CLHEP::HepRotation m_rot;
            CLHEP::Hep3Vector m_gdir;
        };


        //! constructor takes an array of source positions and creates an grid of rotations
        //! @param sources an array of source positions
        //! @param files an array of either root merit file
        //! @param arcmins determines the resolution of the grid about the identity rotation (0,0,0)
        AlignProc(std::vector<astro::SkyDir>& sources, std::vector<std::map<std::pair<int,int>,double> >& alphas ,std::vector<std::string>& rootfiles, 
            const skymaps::PhotonBinner& pb, double arcsecs, double start=-1.0, double stop=-1.0,int event_class=-1); 

        AlignProc(std::vector<astro::SkyDir>& sources, std::vector<std::map<std::pair<int,int>,double> >& alphas ,std::vector<std::string>& fitsfiles, const std::string& ft2file, 
            const skymaps::PhotonBinner& pb, double arcsecs, double start=-1.0, double stop=-1.0,int event_class=-1);
        //! takes event information and returns a photon with pointing information
        pointlike::AlignProc::Photona events(std::vector<float>& row);

        //! accumulates likelihood information from each alignment photon (contains pointing information)
        int add(pointlike::AlignProc::Photona& p);

        //! returns the loglikelihood of one of the (2*RotInfo::s_points+1)**3 grid points
        //! @param x +RotInfo::s_points...,0,....-RotInfo::s_points
        //! @param y +RotInfo::s_points...,0,....-RotInfo::s_points
        //! @param z +RotInfo::s_points...,0,....-RotInfo::s_points
        double loglikelihood(int x, int y, int z) {return m_roti.likelihood(x,y,z);}

        //! returns the fit parameters for the likelihood surface
        //! surface is presumed to be a 3-D second order equation
        LikeSurface fitsurface();

        //! returns misalignment in form {RotX(x[0])*RotY(x[1])*RotZ(x[2])}
        std::vector<double> alignment();

        //! Corrects alignment by applying rotation
        //! @param x rotation about x-axis in arcseconds
        //! @param y rotation about y-axis in arcseconds
        //! @param z rotation about z-axis in arcseconds
        static void addRot(double x,double y, double z); 

        //maximum scaled deviation
        
        static int setumax(int umax) { int t = s_umax; s_umax = umax; return t;}
        static int umax() {return s_umax;}
        static void set_rot(const CLHEP::HepRotation& R);//{s_hr=R;}

        static int class_level() {return s_classlevel;}
        static int set_class_level(int level){int t = s_classlevel; s_classlevel = level; return t;}

    private:

        bool addgti(const std::string& s);
        static int s_umax;
        static HepRotation s_hr;
        static int s_classlevel;
        int m_event_class;
        void loadroot(const std::string& file);
        void loadfits(const std::string& file);
        int m_photons; //number of photons
        astro::PointingHistory* m_ph;
        double m_start; //start time 0... (-1 for all)
        double m_stop; //end time ...stop
        double m_arcsec; //resolution of the grid
        RotationInfo m_roti; //contains likelihood data for each rotation grid point
        const std::vector<astro::SkyDir> m_skydir; //source positions
        const std::vector<std::map<std::pair<int,int>,double> > m_alphas; // source signal fraction
        const skymaps::PhotonBinner m_binner;
    };
}
#endif
