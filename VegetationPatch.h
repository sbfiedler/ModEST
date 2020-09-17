//!  VegetationPatch.
/*!
  This class handles all vegetation processes and parameters relevant on the patch scale.
*/

#ifndef VEGETATIONPATCH_H
#define VEGETATIONPATCH_H

#include <vector>

using namespace std;

//Forward declaration
class Input;
class Patch;
class Plant;

class VegetationPatch
{
    public:
        /*!
         * \brief A constructor.
         */
        VegetationPatch(int xCor, int yCor, const Input& input, const Patch& _patch);

        //! A destructor.
        virtual ~VegetationPatch();

        /*!
         * \brief Delete the list of plants.
         */
        void deletePlantIDList();

        //Member variables
        const Patch& patch;     /*!< Reference to Patch object */
        double xCor;     /*!< x-coordinate of the patch */
        double yCor;     /*!< y-coordinate of the patch */

        vector< pair<int, double> > presentPlants; /*!< List of all present plants in the patch */
        vector<double> relativeRoots;   /*!< Relative roots in the patch per soil layer [-] */
        double absVegCover;   /*!< Total absolute vegetation cover in the patch [cm2] */
        double relVegCover;   /*!< Total relative vegetation cover in the patch [-] */
        double relFPC; /*!< the relative area of ground covered by foliage above it [-] */
        double abovegroundCMass;   /*!< Aboveground alive carbon mass in the patch [kg] */
        double abovegroundNMass;   /*!< Aboveground alive nitrogen mass in the patch [kg] */
        vector<double> belowgroundCMass; /*!< List of alive carbon mass per layer of the patch [kg] */
        vector<double> belowgroundNMass; /*!< List of alive nitrogen mass per layer of the patch [kg] */
        double totalAliveCMass; /*!< Total alive carbon mass (incl. heartwood mass) in the patch [kg] */
        double standingDeadCMass; /*!< Total standing dead carbon mass (sapwood and heartwood) in the patch from dead plants [kg] */
        double standingDeadNMass; /*!< Total standing dead nitrogen mass (sapwood and heartwood) in the patch from dead plants [kg] */
        double abovegroundLitterCMass;  /*!< Dead aboveground carbon mass in the patch [kg] */
        double deadAboveLitterCMass;  /*!< Dead aboveground carbon mass in the patch from dead plants [kg] */
        double abovegroundLitterNMass;  /*!< Dead aboveground nitrogen mass in the patch [kg] */
        double deadAboveLitterNMass;  /*!< Dead aboveground nitrogen mass in the patch from dead plants [kg] */
        vector<double> belowgroundLitterCMass; /*!< List of all dead carbon mass per layer of the patch [kg] */
        vector<double> deadBelowLitterCMass; /*!< List of all dead carbon mass per layer of the patch from dead plants [kg] */
        vector<double> belowgroundLitterNMass; /*!< List of all dead nitrogen mass per layer of the patch [kg] */
        vector<double> deadBelowLitterNMass; /*!< List of all dead nitrogen mass per layer of the patch from dead plants [kg] */
        vector<double> fracRelWaterWP;   /*!< fractional relative water content at wilting point for each soil layer [-] */
};

#endif // VEGETATIONPATCH_H
