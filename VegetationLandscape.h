//!  VegetationLandscape.
/*!
  This class handles all water processes and parameters relevant on the landscape scale.
*/

#ifndef VEGETATIONLANDSCAPE_H
#define VEGETATIONLANDSCAPE_H

#include <vector>
#include <string>
#include <chrono>

using namespace std;

//Forward declaration
class Input;
class Patch;
class Plant;
class Landscape;

class VegetationLandscape
{
    public:
        /*!
         * \brief A constructor.
         * \param input pointer to the Input object
         */
        VegetationLandscape(const Input& _input, vector< vector<Patch*> >& _grid, Landscape* _land);

        //! A destructor.
        virtual ~VegetationLandscape();

        /*!
         * \brief Initialise plants within the landscape.
         */
        void distributePlants();

        /*!
         * \brief Calculate which plants are in which patch. This is important to calculate the intersection area between patches and crown or root areas
         */
        void plantsPerPatch();

        /*!
         * \brief patchesPerPlant
         */
        void patchesPerPlant();

        /*!
         * \brief calculatePatchVariables
         */
        void calculatePatchVariables();

        /*!
         * \brief
         */
        void coverPerPlantAndPatch();

        /*!
         * \brief
         */
        void rootsPerPatch();

        /*!
         * \brief
         */
        void variablesPerPatch();

        /*!
         * \brief calculateProcesses
         */
        void calculateProcesses(int day, double dayLength);
        void vegMortality(int day, int plant);
        void deleteDeadPlants();
        void dispersal(int day, int plant);
        /*!
         * \brief This process is newly included and
         * calculates first the intersected area of overgrowing plants
         * and afterwards the APAR that will be used by these larger plants.
         */
        void shadingAPAR(double extraRadiation);
        double calculateIntersectionArea(double myRadius, double neighborRadius, double distance);
        void calculateAPAR(double extraRadiation, double shadedArea, double totalShadedArea, int plant);

        //Public member variables
        vector<Plant*> plantList; /*!< List of all plants in the landscpae */

    private:
        //Private member variables
        const Input& input; /*!< Reference to the input object */
        const int xCells;   /*!< Number of x-cells */
        const int yCells;   /*!< Number of y-cells */
        double xLength; /*!< maximum length of the landscape in x direction [m] */
        double yLength; /*!< maximum length of the landscape in y direction [m] */
        int plantIDs;   /*!< maximum number of assigned plant IDs [-] */
        int newIndividuals; /*!< number new germinated individuals [-] */

        vector<int> deadPlants; /*!< List dead plants in this time step */
        vector< vector<string> > iniVegetation; /*!< Plants within the landscape (coordinate and initial properties of the plant) */
        vector< vector<Patch*> >& grid; /*!< Reference to two-dimensional landscape (grid) of the type "Patch" for each element (grid cell) */
        Landscape* land;    /*!< Pointer to Landscape object */
};

#endif // VEGETATIONLANDSCAPE_H
