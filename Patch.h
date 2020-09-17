//!  Patch.
/*!
  This class calls patch processes.
*/

#ifndef PATCH_H
#define PATCH_H

#include "WaterPatch.h"
#include "NutrientPatch.h"
#include "VegetationPatch.h"

//Forward declaration
class Input;

class Patch
{
    public:
        /*!
         * \brief A constructor.
         * \param index position of soil type in the soil type vector
         * \param xCor x-coordinate of the patch
         * \param yCor y-coordinate of the patch
         * \param input pointer to the Input object
         */
        Patch(int index, int xCor, int yCor, const Input& input);

        //! A destructor.
        virtual ~Patch();

        //Member variables
        WaterPatch* waterPatch;   /*!< Define new WaterPatch object */
        NutrientPatch* nutrientPatch;   /*!< Define new WaterPatch object */
        VegetationPatch* vegetationPatch;   /*!< Define new WaterPatch object */
};

#endif // PATCH_H
