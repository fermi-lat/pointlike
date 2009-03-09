/** @file ConfidenceLevel.h 
    @brief declaration of the ConfidenceLevel goodness-of-fit evaluation class

    $Header$
*/

#ifndef pointlike_ConfidenceLevel_h
#define pointlike_ConfidenceLevel_h

#include "astro/SkyDir.h"
#include <vector>



namespace pointlike{

    class SimpleLikelihood;

    /** @class ConfidenceLevel
        @brief Evaluate goodness-of-fit confidence level for a SimpleLikelihood 
        analysis



    */
    class ConfidenceLevel {

    public:
        /// @brief ctor
        /// @param slike a simple likelihood object to check
        /// @param the direction to use
        /// @param realizations[100] number of realizations to generate
        ConfidenceLevel(const SimpleLikelihood& slike, const astro::SkyDir& dir, 
            int realizations=100);

        /// @brief return the confidence level value
        double operator()() const;

        /// @brief the negative log-likelihood, as reevaluated
        /// Note that there is an offset from the SimpleLikelihood value
        double likelihood(double alpha, bool save=false) const;

        /// @brief realize an experiment, assuming likelihood called with save=true
        /// @param n[0] number of experiments: if 0, use default
        double realize()const;

    private:
       
        const SimpleLikelihood& m_slike;  ///< object to manage.
        int m_events;
        std::vector<double> m_signal, m_back;
        mutable std::vector<double> m_predict;
        std::vector< std::pair<astro::SkyDir,int> > m_pixels;

        void setup(const astro::SkyDir& dir);

        int m_realizations; ///< number of realizations
    };


} // pointlike


#endif
