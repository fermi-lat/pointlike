/** @file Convolution.cxx
@brief Convolves healpix maps with a sky function

@author M. Roth 

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/src/Convolution.cxx,v 1.3 2007/11/22 02:38:49 burnett Exp $
*/
#include "pointlike/Convolution.h"
#include "healpix/AlmOp.h"
#include "healpix/HealPixel.h"
#include "src/base/alm_filter_tools.h"
#include "pointlike/PointSourceLikelihood.h"
#include <map>

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
    double klug = 1;
    int minE = 555;
}

namespace pointlike {
    Convolution::Convolution (const pointlike::SkySpectrum &sf, const pointlike::SkySpectrum &ker, double energy, int level) : m_level(level) {
        createConv(sf,ker,energy);

    }

    Convolution::Convolution(const pointlike::SkySpectrum &sf, double energy, int level) : m_level(level) {
        createConv(sf,energy);
    }

    void Convolution::createConv(const pointlike::SkySpectrum& sf, double energy) {
        iterator it = this->find((int)energy);
        if(energy<minE&&it==end()) {
            std::cout << "*********************************************************" << std::endl;
            std::cout << "     Convolution: level " << m_level << " and energy " << energy << std::endl << std::endl;
            pointlike::PsfFunction psf(2.25);
            std::cout << "          Allocating map and harmonic storage...";
            Map<double> sfm(m_level);
            Map<double> psfm(m_level);
            double p0 = 0.058;
            double p1 = 0.000377;
            double sigmasq = 0.404*((p0*p0*pow(energy/100,-1.6))+p1*p1)*M_PI/180;
            //point spread function parameter for energy band
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
                double u = 0.5*(theta*theta)/(sigmasq);
                u = psf(u);
                for(int j=0;j<ringpix;++j) {
                    ShowPercent(j+startpix,sfm.map()->Npix(),j+startpix);
                    HealPixel hp(sfm.map()->ring2nest(j+startpix),m_level);
                    (*psfm.map())[startpix+j] = u;
                    (*sfm.map())[startpix+j] = sf.value(hp,energy)*klug;
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
            this->insert(std::map<int,healpix::Map<double> >::value_type((int)energy,sfm));
        }
    }

    void Convolution::createConv(const pointlike::SkySpectrum& sf, const pointlike::SkySpectrum& ker, double energy) {
        iterator it = this->find((int)energy);
        if(it==end()) {
            std::cout << "*********************************************************" << std::endl;
            std::cout << "     Convolution: level " << m_level << " and energy " << energy << std::endl << std::endl;
            pointlike::PsfFunction psf(2.25);
            std::cout << "          Allocating map and harmonic storage...";
            Map<double> sfm(m_level);
            Map<double> kerm(m_level);
            AlmOp<xcomplex<double> > sfh(2*sfm.map()->Nside(),2*sfm.map()->Nside());
            AlmOp<xcomplex<double> > kerh(2*sfm.map()->Nside(),2*sfm.map()->Nside());
            std::cout << "done!" << std::endl;
            std::cout << "          Populating sky and psf maps...";
            for(int i(0);i<sfm.map()->Npix();++i) {
                ShowPercent(i,sfm.map()->Npix(),i);
                HealPixel hp(i,m_level);
                (*sfm.map())[sfm.map()->nest2ring(i)] = sf(hp,energy);
                (*kerm.map())[kerm.map()->nest2ring(i)]= ker(hp,energy);
            }
            map2alm_iter(*sfm.map(),*sfh.Alms(),0);
            map2alm_iter(*kerm.map(),*kerh.Alms(),0);
            std::cout << "done!" << std::endl;
            std::cout << "          Calculating sky harmonics...";
            map2alm_iter(*sfm.map(),*sfh.Alms(),0);
            std::cout << "done!" << std::endl;
            std::cout << "          Calculating psf harmonics...";
            map2alm_iter(*kerm.map(),*kerh.Alms(),0);
            std::cout << "done!" << std::endl;
            mf_constantnoise(*sfh.Alms(),*kerh.Alms());
            std::cout << "          Performing convolution...";
            alm2map(*sfh.Alms(),*sfm.map());
            std::cout << "done!" << std::endl;
            this->insert(std::map<int,healpix::Map<double> >::value_type((int)energy,sfm));
        }
    }

    double Convolution::value(const astro::SkyDir& dir,double e) const{
        double e1(layer(e,true)), e2(layer(e,false));
        int pixel_index = find(e1)->second.cmap().nest2ring(healpix::HealPixel(dir,m_level).index());
        double f1((find(e1)->second.cmap())[pixel_index]),f2((find(e2)->second.cmap())[pixel_index]);
        double alpha ( log(f1/f2)/log(e2/e1) );
        return f1*pow( e1/e, alpha);
    }

    double Convolution::integral(const astro::SkyDir& dir, double a, double b) const{
        double fa(value(dir, a)), fb(value(dir,b));
        double q ( 1. - log(fa/fb)/log(b/a) );
        return fa* a * (pow(b/a, q)-1)/q;
    }

    int Convolution::layer(double e, bool front) const{
        const_iterator it = front?begin():end();
        const_iterator stop = front?begin():end();
        double laste = it->first;
        for(; it!=stop; ++it){
            if( front?(it->first > e):(it->first <e) ) break;
            laste = it->first;
        }
        return (int)laste;
    }
}
