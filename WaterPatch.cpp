#include "WaterPatch.h"
#include "Patch.h"
#include "Input.h"
#include <cmath>    //to use mathematical functions

/////////////////////////////////////////////////////////
void WaterPatch::fastInfiltration(int layer)
{
    /// \todo check and adapt for more than two layers, maybe ask Gregor (Gregor schaut wohin der größte Gradient ist)
    /// \todo compared to EcoHyD fast infiltratrion can also occur in the first layer now

    ///Saturation function g1 of vegetation roots in the respective layer
    //Originally: Britta used just shrub cover instead of roots (eqn. 3)
    //Here: According to Tong's version, where relative roots of shrubs and perennials will be used (relative roots per PFT in respective layer weighed by the relative cover they have)
    double saturationFunction; //function to describe the functional response of infiltration to vegetation based on relative root proportion in soil layer [0, 1]
    saturationFunction = ( infRateBare + ( patch.vegetationPatch->relativeRoots[layer] ) ) / (1 + ( patch.vegetationPatch->relativeRoots[layer] ));

    ///Fast, preferential infiltration into lower layers via macropores due to cracks and roots (eqn. 2)
    waterLayer[layer].infWater = min( surfaceWater * soilProp.infRate * saturationFunction , soilProp.maxTotalInf );

    //Update surface and layer water after infiltration
    surfaceWater -= waterLayer[layer].infWater;
    waterLayer[layer].absWater += waterLayer[layer].infWater;
    waterLayer[layer].relWater = waterLayer[layer].absWater / waterLayer[layer].depthLayer;
}

/////////////////////////////////////////////////////////
void WaterPatch::slowInfiltration(int layer)
{
    //fillable pores
    double fillablePores;    // space in pores, which can be filled (n_a) [vol%]
    fillablePores = soilProp.fieldCap - waterLayer[layer].relWater; //// \todo Britta originally uses mean relative water of the first layer in the ladscape! Why? Also: In paper there is totally different calculation for it?

    //Holling type II function: the more vegetation the bigger saturationFunction
    //Originally: Britta used just vegetation cover instead of roots (eqn. 5)
    //until saturation occurs, the strange value is chosen in a way that for a cover of cGrass+cShrub=1 the function equals 0.95 and without vegetation the function equals 0.2
    //// \todo make the following constants generic?
    //0.013 for saturationFunction(0.0, 0.0) = 0.2; ... at least some infiltration, if there is no vegetation
    //0.067 for saturationFunction(0.5, 0.5) = 0.95; ... almost full infiltration at a total cover of 1

    //Here: According to Tong's version, where relative roots of shrubs and perennials will be used (relative roots per PFT in respective layer weighed by the relative cover they have)
    double saturationFunction; //function to describe the functional response of infiltration to vegetation based on relative root proportion in soil layer [0, 1]
    saturationFunction = ( 0.013 + ( patch.vegetationPatch->relativeRoots[layer] ) ) / ( 0.067 + ( patch.vegetationPatch->relativeRoots[layer] ) );

    //no ponding if hydraulic conductivity is bigger than available water (divided by time step), available surface water infiltrates in the upper soil layer (eqn. 6)
    if (waterLayer[layer].absWater <= geomHydrCon) {
        waterLayer[layer].infWaterSat = waterLayer[layer].absWater;
        waterLayer[layer].timeSat = -1.0;
        waterLayer[layer].infWater = waterLayer[layer].infWaterSat;
        waterLayer[layer].wettingFrontFlag = false;
    }

    //else ponding can theoretically occur
    else {
        //there was no ponding in the previous time step (eqn. 4)
        if (waterLayer[layer].wettingFrontFlag == false) {
            //// \todo check eqn 4 with Rawls et al. 1992
            waterLayer[layer].timeSat = ( ( (fillablePores * soilProp.suction) / (surfaceWater / geomHydrCon - 1) / surfaceWater) ) * saturationFunction;
            waterLayer[layer].infWaterSat = min(1.0, waterLayer[layer].timeSat) * surfaceWater;
        }
        //there was already ponding in the previous time step
        else {
            waterLayer[layer].timeSat = 0.0;
            waterLayer[layer].infWaterSat = 0.0;
        }
    }

    //if ponding occurs in this time step (eqn. 7)
    if (waterLayer[layer].timeSat < 1.0) { // TODO: Katja found problem. This might be solution: if ((ts >= 0) && (ts < 1)) { //previous ts < 1
        //this routine is taken from WaSiM-ETH
        //after Peschke (in Dyck, Peschke 1989)
        double A = geomHydrCon * (1 - waterLayer[layer].timeSat);
        double B = max(waterLayer[layer].infWaterSat + 2 * fillablePores * soilProp.suction, 0.0);
        double C = A * A / 4.0 + A * B + waterLayer[layer].infWaterSat * waterLayer[layer].infWaterSat;
        waterLayer[layer].infWater = min( ( A / 2.0 + (A * A + 2 * A * B) / (4.0 * sqrt(C)) ) * saturationFunction + waterLayer[layer].infWaterSat, surfaceWater); //// \todo differs from eqn. 7
        waterLayer[layer].wettingFrontFlag = true;
    }

    //if ponding does not occur
    else {
        waterLayer[layer].infWater = waterLayer[layer].infWaterSat;
        waterLayer[layer].wettingFrontFlag = false;
    }

    surfaceWater -= waterLayer[layer].infWater;
    waterLayer[layer].absWater += waterLayer[layer].infWater;
    waterLayer[layer].relWater = waterLayer[layer].absWater / waterLayer[layer].depthLayer;
}

/////////////////////////////////////////////////////////
void WaterPatch::drainage(int layer)
{
    //if the field capacity of the layer is exceeded, there is water flow to the layer below
    if (layer < int(waterLayer.size() - 1) && waterLayer[layer].relWater > soilProp.fieldCap) {
        waterLayer[layer].drainedWater = (waterLayer[layer].relWater - soilProp.fieldCap) * waterLayer[layer].depthLayer;
        waterLayer[layer+1].absWater += waterLayer[layer].drainedWater;
        waterLayer[layer].relWater = soilProp.fieldCap;
        waterLayer[layer].absWater = waterLayer[layer].relWater * waterLayer[layer].depthLayer;
        waterLayer[layer+1].relWater = waterLayer[layer+1].absWater / waterLayer[layer+1].depthLayer;
    }

    //if the field capacity of the last layer is exceeded, there is deep drainage
    if (layer == int(waterLayer.size() - 1) && waterLayer[layer].relWater > soilProp.fieldCap) {
        waterLayer[layer].drainedWater = (waterLayer[layer].relWater - soilProp.fieldCap) * waterLayer[layer].depthLayer;
        deepDrainedWater += waterLayer[layer].drainedWater;
        waterLayer[layer].relWater = soilProp.fieldCap;
        waterLayer[layer].absWater = waterLayer[layer].relWater * waterLayer[layer].depthLayer;
    }
}

/////////////////////////////////////////////////////////
void WaterPatch::diffusion(int layer)
{
    //water transport between layers (darcy's law)
    if(layer < int(waterLayer.size() - 1)) {
        if (waterLayer[layer].relWater > waterLayer[layer+1].relWater) {
            //eqn 9
            waterLayer[layer].diffusedWater = soilProp.diffSpeed * geomHydrCon * ((waterLayer[layer].relWater - waterLayer[layer+1].relWater) * waterLayer[layer].depthLayer) / ((waterLayer[layer].depthLayer + waterLayer[layer+1].depthLayer) / 2.0);	//Darcy's Law, Maidment 5.17
            waterLayer[layer+1].diffusedWater = 0.0;
            waterLayer[layer].absWater -= waterLayer[layer].diffusedWater;
            waterLayer[layer+1].absWater += waterLayer[layer].diffusedWater;
            waterLayer[layer].relWater = waterLayer[layer].absWater / waterLayer[layer].depthLayer;
            waterLayer[layer+1].relWater = waterLayer[layer+1].absWater / waterLayer[layer+1].depthLayer;
        }
        else if (waterLayer[layer+1].relWater > waterLayer[layer].relWater) {
            //eqn 9
            waterLayer[layer].diffusedWater =  soilProp.diffSpeed * geomHydrCon * ((waterLayer[layer+1].relWater - waterLayer[layer].relWater) * waterLayer[layer+1].depthLayer) / ((waterLayer[layer].depthLayer + waterLayer[layer+1].depthLayer) / 2.0);
            waterLayer[layer+1].diffusedWater = 0.0;
            waterLayer[layer+1].absWater -= waterLayer[layer].diffusedWater;
            waterLayer[layer].absWater += waterLayer[layer].diffusedWater;
            waterLayer[layer+1].relWater = waterLayer[layer+1].absWater / waterLayer[layer+1].depthLayer;
            waterLayer[layer].relWater = waterLayer[layer].absWater / waterLayer[layer].depthLayer;
        }
    }
}

/////////////////////////////////////////////////////////
void WaterPatch::surfaceEvaporation(double potEP)
{
    //potential evaporation [mm/h] dependent on temperature, extraterrestrial radiation, slope and aspect factor, eqn. 10
    potEP = potEP * aspectFactor * cos(inclination);
    //acual evaporation depends on vegetation cover (shading effect), eqn. 12
    potEP = potEP * (1.0 - patch.vegetationPatch->relVegCover / 2.0);

    if (surfaceWater > 0) {
        //if there is more surface water than there could potentially evaporate some surface water is left
        if (surfaceWater - potEP > 0) {
            surfaceWater -= potEP;
        }
        //otherwise all surface water evaporates
        else {
            surfaceWater = 0;
        }
    }
}

/////////////////////////////////////////////////////////
void WaterPatch::soilEvaporation(double potEP)
{
    //decreasing function with increasing cover
    double vegfunctionL1 = 1.2 - 0.2 * patch.vegetationPatch->relVegCover; //linear decreasing function to estimate transpiration vs evaporation for the upper soil layer
    double actEP = 0; //actual evapotranspiration [mm]

    potEP = potEP * aspectFactor;

    //this routine is oriented according to
    //Porporato's models but also to the HBV model
    //it is flattened between 0 and wp and constant above the wp

    //evapotranpiration from the upper layer
    if(waterLayer[0].relWater > patch.vegetationPatch->fracRelWaterWP[0]){
        //constant value evaporates above wilting point, but a residual water content must always be there
        actEP = min(potEP * cos(inclination) * vegfunctionL1 * soilProp.evapoFactor,(waterLayer[0].relWater - soilProp.residualWater) * waterLayer[0].depthLayer);
    }
    else
        actEP = max(0.0, min(potEP * pow(waterLayer[0].relWater / patch.vegetationPatch->fracRelWaterWP[0], 2) * cos(inclination) * vegfunctionL1 * soilProp.evapoFactor,(waterLayer[0].relWater - soilProp.residualWater) * waterLayer[0].depthLayer));

    //NEW: extract the transpired water from the transpiration routine to only get the evaporated water
    actEP = max(actEP - waterLayer[0].transpiredAbsWater, 0.0);
    waterLayer[0].transpiredAbsWater = 0.0;

    //absolute moisture
    waterLayer[0].absWater = max(waterLayer[0].absWater - actEP, 0.0);
    //relative moisture
    waterLayer[0].relWater = waterLayer[0].absWater / waterLayer[0].depthLayer;
}

/////////////////////////////////////////////////////////
void WaterPatch::layerEvapotranspiration(int layer, double potEP)
{
    //// \todo this function should be renewed completely (transpiration should go into the veg model, new routine for evaporation from the layers)
    //decreasing function with increasing cover
    double vegfunctionL1 = 1.2 - 0.2 * patch.vegetationPatch->relVegCover; //linear decreasing function to estimate transpiration vs evaporation for the upper soil layer
    double WP = 0.041; ///// \todo Wilting point should come from the vegetation model, even better: transpiration should be calculated in the veg model
    double actEP = 0; //actual evapotranspiration [mm]

    potEP = potEP * aspectFactor;

    //this routine is oriented according to
    //Porporato's models but also to the HBV model
    //it is flattened between 0 and wp and constant above the wp

    //evapotranpiration from the upper layer
    //// \todo roots should be included here as well for transpiration, not only cover for shade
    if(layer == 0) {
        if(waterLayer[layer].relWater > WP){
            //constant value evaporates above wilting point, but a residual water content must always be there
            actEP = min(potEP * cos(inclination) * vegfunctionL1 * soilProp.evapoFactor,(waterLayer[layer].relWater - soilProp.residualWater) * waterLayer[layer].depthLayer);
        }
        else
            actEP = max(0.0, min(potEP * pow(waterLayer[layer].relWater / WP,2) * cos(inclination) * vegfunctionL1 * soilProp.evapoFactor,(waterLayer[layer].relWater - soilProp.residualWater) * waterLayer[layer].depthLayer));
    }

    //transpiration from the lower layers
    else {
        if (waterLayer[layer].relWater <= WP)
            actEP = max(0.0,min(potEP * pow(waterLayer[layer].relWater / WP,2) * patch.vegetationPatch->relativeRoots[layer] * cos(inclination),(waterLayer[layer].relWater - soilProp.residualWater) * waterLayer[layer].depthLayer));

        else
            actEP = min(potEP * cos(inclination) * patch.vegetationPatch->relativeRoots[layer],(waterLayer[layer].relWater - soilProp.residualWater) * waterLayer[layer].depthLayer);
    }

    //absolute moisture
    waterLayer[layer].absWater = max(waterLayer[layer].absWater - actEP, 0.0);
    //relative moisture
    waterLayer[layer].relWater = waterLayer[layer].absWater / waterLayer[layer].depthLayer;
}

/////////////////////////////////////////////////////////
WaterPatch::WaterPatch(int index, int _xCor, int _yCor, const Input& input, const Patch& _patch):
    surfaceWater(input.iniSurfaceWater), elevation(input.elevation[_xCor][_yCor]), xCor(_xCor), yCor(_yCor), patch(_patch)
{
    soilProp =  {
                 input.soilProperties[index][0], atof(input.soilProperties[index][1].c_str()), atof(input.soilProperties[index][2].c_str())*24.0, atof(input.soilProperties[index][3].c_str()),
                 atof(input.soilProperties[index][4].c_str()), atof(input.soilProperties[index][5].c_str()), atof(input.soilProperties[index][6].c_str()), atof(input.soilProperties[index][7].c_str())*24.0,
                 atof(input.soilProperties[index][8].c_str()), atof(input.soilProperties[index][9].c_str()), atof(input.soilProperties[index][10].c_str()),  atof(input.soilProperties[index][11].c_str()),
                 atof(input.soilProperties[index][12].c_str()), atof(input.soilProperties[index][13].c_str())
                };

    /// Calculation of unsaturated hydraulic conductivity
    // based on saturated and unsaturated conductivity
    // see equation in: Kemp et al. 1997, p.81, following Gardner 1958
    // K = Ks / (1 + Psi / Psi*)^p
    unsatHydrCon = soilProp.satHydrCon / (pow ( 1.0 +( -30.0/-26.0 ), 3.0 ) );

    /// Calculation of geometric mean hydraulic conductivity
    // based on saturated and unsaturated conductivity
    geomHydrCon = sqrt(unsatHydrCon * soilProp.satHydrCon);

    ///Infiltration without vegetation is determined by hydraulic conductivity
    /// \todo Why these levels? Should be adapted to daily values, since geomHyrCon is now in another range. Values given here were set when there was still an hourly calculation.
    if (geomHydrCon < 1) infRateBare = 0.1;
    else if (geomHydrCon < 10) infRateBare = 0.5;
    else infRateBare = 0.8;

    //Create water layers per patch
    waterLayer.resize(input.nSoilLayers);
    WaterLayer helpLayer;
    for(unsigned int i = 0; i < waterLayer.size(); i++) {
        helpLayer = {input.depthLayers[i], input.iniWaterLayer[i] * input.depthLayers[i], input.iniWaterLayer[i], 0, 0, 0, 0, 0, -1, false};
        waterLayer[i] = helpLayer;
    }

    deepDrainedWater = 0;
}

/////////////////////////////////////////////////////////
WaterPatch::~WaterPatch()
{

}
