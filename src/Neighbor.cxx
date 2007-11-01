/** @file Neighbor.cxx
    @brief implement class Neighbor

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/pointlike/DiffuseFunction.h,v 1.8 2007/10/29 16:41:34 burnett Exp $

*/

#include "pointlike/Neighbor.h"
#include "pointlike/PointSourceLikelihood.h"

using namespace pointlike;


double Neighbor::value(const astro::SkyDir& dir, double e) const
{
    return 0; // placeholder
}
    ///@brief integral for the energy limits, in the given direction
double Neighbor::integral(const astro::SkyDir& dir, double a, double b)const
{
    return 0; // placeholder
   

}


