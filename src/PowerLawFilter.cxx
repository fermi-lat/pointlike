/** @file PowerLawFilter.cxx
@brief implementation of PowerLawFilter

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/src/PowerLawFilter.cxx,v 1.2 2007/08/27 23:24:00 mar0 Exp $
*/

#include "pointlike/PowerLawFilter.h"

#include <vector>
#include <cmath>

using namespace pointlike;


PowerLawFilter::PowerLawFilter(const std::vector<std::pair<double, double> > & values)
: m_slope(0.0)
, m_constant(0.0)
, m_metric(0.0)
, m_chi_sq(0.0)
{
    if (values.size() < 2) return;  // No pattern if less than two points.

    std::vector<std::pair<double, double> > log_values;


    log_values.clear();
    for (std::vector<std::pair<double, double> >::const_iterator it = values.begin();
        it != values.end(); ++it)
    {
        double ln_x = (it->first > 0)? log(it->first): -300.0,
            ln_y = (it->second > 0)? log(it->second): -300.0;
        log_values.push_back(std::make_pair(ln_x, ln_y));
    }

    // Now remove entries from the back end of the vector with log values that are very small, coresponding to "near zero" input values.
    while (log_values.size() > 1 && log_values.back().second < -299.0)
    {
	    log_values.pop_back();
    }

    if (log_values.size() < 2) return;  // No pattern if less than two points.

    double sum_x_sqrd(0), sum_y_sqrd(0), sum_xy(0), x_ave(0), y_ave(0);

    // reference http://mathworld.wolfram.com/LeastSquaresFitting.html

    for (std::vector<std::pair<double, double> >::const_iterator it = log_values.begin();
        it != log_values.end(); ++it)
    {
        sum_x_sqrd += it->first * it->first;
        sum_y_sqrd += it->second * it->second;
        sum_xy += it->first * it->second;
        x_ave += it->first;
        y_ave += it->second;
    }

    int n = log_values.size();
    x_ave /= n;
    y_ave /= n;

    double ss_xy = sum_xy - (n * x_ave * y_ave),
        ss_xx = sum_x_sqrd - (n * x_ave * x_ave),
        ss_yy = sum_y_sqrd - (n * y_ave * y_ave);

    m_slope = (ss_xx < 1e-10)? 0: ss_xy / ss_xx;
    m_constant = y_ave - (m_slope * x_ave);
    m_metric = ((ss_xx * ss_yy) < 1e-10)? 0: (ss_xy * ss_xy) / (ss_xx * ss_yy);

    // Calculate chi squared
    for (std::vector<std::pair<double, double> >::const_iterator it = log_values.begin();
        it != log_values.end(); ++it)
    {
        double chi = it->second - m_constant - (m_slope * it->first);
        m_chi_sq += chi * chi;
    }

    n = n; // for debug
}



