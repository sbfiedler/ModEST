//!  Landscape.
/*!
  This class handles all processes and parameters relevant on the landscape scale.
*/

#ifndef LANDSCAPE_H
#define LANDSCAPE_H

#include <vector>

using namespace std;

//Forward declaration
class Input;
class Patch;
class WaterLandscape;
class NutrientLandscape;
class NutrientPatch;
class VegetationLandscape;
class Output;

class Landscape
{
    public:
        /*!
         * \brief A constructor.
         * \param _input pointer to the Input object
         */
        Landscape(const Input& _input);

        //! A destructor.
        virtual ~Landscape();

        //Member functions
        /*!
         * \brief Calls all processes within the landscape in a logical order.
         * \param day current day
         */
        void calculateProcesses(int day);

        /*!
         * \brief
         * \param day
         */
        double extraterrestrialRadiation(int day);

        /*!
         * \brief
         * \param day
         */
        double dayLength(int day);

        //Member variables
        vector< vector<Patch*> > grid;  /*!< Two-dimensional landscape (grid) of the type "Patch" for each element (grid cell) */
        WaterLandscape* water;   /*!< Pointer to WaterLandscape object */
        NutrientLandscape* nutrient;   /*!< Pointer to NutrientLandscape object */
        NutrientPatch* nutrientPatch;   /*!< Pointer to NutrientPatch object */
        VegetationLandscape* vegetation;   /*!< Pointer to VegetationLandscape object */
        Output* output;  /*!< Pointer to Output object */
        const Input& input; /*!< Reference to Input object */
};

#endif // LANDSCAPE_H
