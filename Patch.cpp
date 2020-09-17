#include "Patch.h"

/////////////////////////////////////////////////////////
Patch::Patch(int index, int xCor, int yCor, const Input& input)
{
    waterPatch = new WaterPatch(index, xCor, yCor, input, *this);
    nutrientPatch = new NutrientPatch(index, xCor, yCor, input, *this);
    vegetationPatch = new VegetationPatch(xCor, yCor, input, *this);
}

/////////////////////////////////////////////////////////
Patch::~Patch()
{
    //Delete patches
    delete nutrientPatch;
    nutrientPatch = 0;
    delete vegetationPatch;
    vegetationPatch = 0;
    delete waterPatch;
    waterPatch = 0;
}
