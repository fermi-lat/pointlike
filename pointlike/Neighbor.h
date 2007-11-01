/** @file Neighbor.h
    @brief declare class Neighbor

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/pointlike/DiffuseFunction.h,v 1.8 2007/10/29 16:41:34 burnett Exp $

*/
#ifndef pointlike_Neighbor_h
#define pointlike_Neighbor_h

#include "pointlike/SkySpectrum.h"


namespace pointlike {
class PointSourceLikelihood;

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/** @class DiffuseFunction
    @brief a SkyFunction that adapts a diffuse map. also includes extragal diffuse

*/

class Neighbor: public SkySpectrum{
public:
    Neighbor(PointSourceLikelihood& near):m_near(near){}

    ///@brief return differential value 
    ///@param e energy in MeV
    virtual double value(const astro::SkyDir& dir, double e)const;
   
    ///@brief integral for the energy limits, in the given direction
    virtual double integral(const astro::SkyDir& dir, double a, double b)const;


private:
    pointlike::PointSourceLikelihood& m_near;
    double energy;

};

}
#endif
