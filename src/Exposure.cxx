/** @file Exposure.cxx

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/src/Exposure.cxx,v 1.1 2007/11/21 07:00:39 burnett Exp $
*/

#include "pointlike/Exposure.h"

#include <stdexcept>

using namespace pointlike;


Exposure::Exposure()
{
}

Exposure::~Exposure()
{}


double Exposure::value(const astro::SkyDir& dir, double)const
{
    return 0;
}

///@brief integral for the energy limits, in the given direction
double Exposure::integral(const astro::SkyDir& dir, double , double )const
{
    return 0;
}

std::string Exposure::name()const
{
    return "Exposure";
}
