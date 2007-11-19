/** @file Convolution.cxx
@brief Convolves healpix maps with a sky function

@author M. Roth 

$Header: /nfs/slac/g/glast/ground/cvs/healpix/healpix/Map.h,v 1.4 2007/06/04 22:14:25 mar0 Exp $
*/
#include "healpix/Convolution.h"
#include "healpix/Map.h"
#include "healpix/AlmOp.h"
#include "healpix/HealPixel.h"
#include "base/alm_filter_tools.h"
#include "pointlike/PointSourceLikelihood.h"
#include <iostream>

using namespace healpix;

namespace{
    void ShowPercent(int sofar, int total, int found)
    {
        static int toskip(50), skipped(0);
        if(++skipped<toskip) return; skipped=0;
        static int lastpercent(-1);
        int percent( static_cast<int>(100 * sofar / total +0.5) );
        if( percent==lastpercent) return;
        lastpercent=percent;
        char   s[50];
        sprintf(s, "%d%%, %d found.", percent, found);
        std::cout << s;
        if (sofar < total)
        {
            for (size_t j = 0; j < strlen(s); ++j)
                std::cout << "\b";
        }
        else
            std::cout << std::endl;
    }

    //klugy scale factor to visualize pmap
    double klug = 1e10;
}

namespace healpix {
    Convolution::Convolution (const astro::SkyFunction &sf, const astro::SkyFunction &ker, int level) {
        Map<double> sfm(level);
        Map<double> kerm(level);
        AlmOp<xcomplex<double> > sfh(2*sfm.map()->Nside(),2*sfm.map()->Nside());
        AlmOp<xcomplex<double> > kerh(2*sfm.map()->Nside(),2*sfm.map()->Nside());
        for(int i(0);i<sfm.map()->Npix();++i) {
            HealPixel hp(i,level);
            (*sfm.map())[sfm.map()->nest2ring(i)] = sf(hp);
            (*kerm.map())[kerm.map()->nest2ring(i)]= ker(hp);
        }
        map2alm_iter(*sfm.map(),*sfh.Alms(),0);
        map2alm_iter(*kerm.map(),*kerh.Alms(),0);
        AlmOp<xcomplex<double> > conv = sfh*kerh;
        alm2map(*conv.Alms(),*sfm.map());
        m_map = sfm.pmap();
    }

    Convolution::Convolution(const astro::SkyFunction &sf, int level){
        std::cout << "*********************************************************" << std::endl;
        std::cout << "                 Convolution: level " << level << std::endl << std::endl;
        pointlike::PsfFunction psf(2.25);
        std::cout << "          Allocating map and harmonic storage...";
        Map<double> sfm(level);
        Map<double> psfm(level);
        //point spread function parameter for energy band
        double sigma = pointlike::PointSourceLikelihood::sigma_level[level]/pow(2.,1.*(level-6))/2.5*M_PI/180;
        //Nyquist frequency for spherical harmonics is half the number of discrete latitudes
        AlmOp<xcomplex<double> > sfh(2*sfm.map()->Nside(),2*sfm.map()->Nside());
        AlmOp<xcomplex<double> > psfh(2*psfm.map()->Nside(),2*psfm.map()->Nside());
        std::cout << "done!" << std::endl;
        std::cout << "          Populating sky and psf maps...";
        //put psf at north pole (azimuthal symmetry)
        //iterate over each latitudnal ring
        for(int i = 1;i<4*sfm.map()->Nside();i++) {
            int startpix,ringpix;
            double theta;
            bool shift;
            sfm.map()->get_ring_info2(i,startpix,ringpix,theta,shift);
            double u = 0.5*(theta*theta)/(sigma*sigma);
            u = psf(u);
            for(int j=0;j<ringpix;++j) {
                ShowPercent(j+startpix,sfm.map()->Npix(),j+startpix);
                HealPixel hp(sfm.map()->ring2nest(j+startpix),level+1);
                (*psfm.map())[startpix+j] = u;
                (*sfm.map())[startpix+j] = sf(hp)*klug;
            }
        }
        std::cout << "done!" << std::endl;
        std::cout << "          Calculating sky harmonics...";
        map2alm_iter(*sfm.map(),*sfh.Alms(),0);
        std::cout << "done!" << std::endl;
        std::cout << "          Calculating psf harmonics...";
        map2alm_iter(*psfm.map(),*psfh.Alms(),0);
        std::cout << "done!" << std::endl;
        mf_constantnoise(*sfh.Alms(),*psfh.Alms());
        std::cout << "          Performing convolution...";
        alm2map(*sfh.Alms(),*sfm.map());
        std::cout << "done!" << std::endl;
        m_map = sfm.pmap();
    }
}