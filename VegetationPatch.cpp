#include "VegetationPatch.h"
#include "Input.h"
#include "Patch.h"
#include "Plant.h"

/////////////////////////////////////////////////////////
VegetationPatch::VegetationPatch(int xCor, int yCor, const Input& _input, const Patch& _patch):
    xCor(((xCor+1.0)*_input.cellSize)-(_input.cellSize/2.0)), yCor(((yCor+1.0)*_input.cellSize)-(_input.cellSize/2.0)), patch(_patch)
{
    relativeRoots.resize(patch.waterPatch->waterLayer.size(), 0.0);
    abovegroundLitterCMass = 0.0;
    abovegroundLitterNMass = 0.0;
    standingDeadCMass = 0.0;
    standingDeadNMass = 0.0;
	deadAboveLitterCMass = 0.0;
    deadAboveLitterNMass = 0.0;
    belowgroundLitterCMass.resize(patch.waterPatch->waterLayer.size(), 0.0);
    belowgroundLitterNMass.resize(patch.waterPatch->waterLayer.size(), 0.0);
    belowgroundCMass.resize(patch.waterPatch->waterLayer.size(), 0.0);
    deadBelowLitterCMass.resize(patch.waterPatch->waterLayer.size(), 0.0);
    belowgroundNMass.resize(patch.waterPatch->waterLayer.size(), 0.0);
    deadBelowLitterNMass.resize(patch.waterPatch->waterLayer.size(), 0.0);
    fracRelWaterWP.resize(patch.waterPatch->waterLayer.size(), 0.0);
}

/////////////////////////////////////////////////////////
void VegetationPatch::deletePlantIDList()
{
    presentPlants.erase(presentPlants.begin(),presentPlants.end());
}

/////////////////////////////////////////////////////////
VegetationPatch::~VegetationPatch()
{

}
