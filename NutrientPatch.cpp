#include "NutrientPatch.h"
#include "VegetationPatch.h"
#include "Input.h"
#include "Patch.h"
#include <cmath>
#include <vector>
#include <algorithm>
#include <cstdlib>


/////////////////////////////////////////////////////////
void NutrientPatch::calcSoilSurfTemperature(int day)
{
    /*!< \todo using albedo as 0.154 for sandy-loam soils, possibly add albedo to soil properties */
    double albedo = 0.154;
    double MJm2d2Radiation;
    MJm2d2Radiation = input.solRadiation[day] * 24 * 60 * 60 / 1000000;
    double radiationTerm;
    radiationTerm = (MJm2d2Radiation * (1 - albedo ) - 14)/20;
    TsoilSurfBare = input.tempMean[day] + radiationTerm * (input.tempMax[day] + input.tempMin[day])/2;
    double totalAliveCMassKgHec = patch.vegetationPatch->totalAliveCMass * 10000 / pow(input.cellSize, 2);
    TsoilCoverF = totalAliveCMassKgHec / (totalAliveCMassKgHec + exp(7.563 - 1.297 * pow(10, -4) * totalAliveCMassKgHec));

    TsoilSurf = (TsoilCoverF * nutrientLayer[0].soilTemp) + ((1 - TsoilCoverF) * TsoilSurfBare);
}

/////////////////////////////////////////////////////////
void NutrientPatch::calcSoilTemperature(int layer, int day)
{
    nutrientLayer[layer].soilWaterScalingFactor = patch.waterPatch->waterLayer[layer].absWater / ((0.356 - 0.144 * soilProp.bulkDensity) * nutrientLayer[layer].maxDepthLayer);

    nutrientLayer[layer].dampingDepth = dampingDepthMaximum * exp(log(500 / dampingDepthMaximum) * pow((1.0 - nutrientLayer[layer].soilWaterScalingFactor) / (1.0 + nutrientLayer[layer].soilWaterScalingFactor), 2));
    nutrientLayer[layer].dampingDepthRatio = nutrientLayer[layer].meanDepthLayer / nutrientLayer[layer].dampingDepth;
    nutrientLayer[layer].soilDepthTF = nutrientLayer[layer].dampingDepthRatio / (nutrientLayer[layer].dampingDepthRatio + exp(-0.867 - (2.078 * nutrientLayer[layer].dampingDepthRatio))); // df in SWAT
    nutrientLayer[layer].soilTemp = (TsoilLag * nutrientLayer[layer].soilTemp) + (1 - TsoilLag) * ((nutrientLayer[layer].soilDepthTF * (input.MAT[input.yearOfDay[day]] - TsoilSurf)) + TsoilSurf);
}

/////////////////////////////////////////////////////////
void NutrientPatch::incomeRes(int layer)
{
    double belowgroundDeadCMass_hectar = 0;
    double abovegroundDeadCMass_hectar = 0;
    double belowgroundDeadNMass_hectar = 0;
    double abovegroundDeadNMass_hectar = 0;
    belowgroundDeadCMass_hectar = (patch.vegetationPatch->belowgroundLitterCMass[layer]) * 10000 / pow(input.cellSize, 2); //to kg ha-1;
    abovegroundDeadCMass_hectar = (patch.vegetationPatch->abovegroundLitterCMass) * 10000 / pow(input.cellSize, 2);  //kg ha-1;
    belowgroundDeadNMass_hectar = (patch.vegetationPatch->belowgroundLitterNMass[layer]) * 10000 / pow(input.cellSize, 2);  //kg ha-1;
    abovegroundDeadNMass_hectar = (patch.vegetationPatch->abovegroundLitterNMass) * 10000 / pow(input.cellSize, 2);  //kg ha-1;
    if(layer < 1)
    {
        nutrientLayer[layer].CresNew = abovegroundDeadCMass_hectar + belowgroundDeadCMass_hectar;
        nutrientLayer[layer].NresNew = abovegroundDeadNMass_hectar + belowgroundDeadNMass_hectar;
    }
    else
    {
        nutrientLayer[layer].CresNew = belowgroundDeadCMass_hectar;
        nutrientLayer[layer].NresNew = belowgroundDeadNMass_hectar;
    }

    patch.vegetationPatch->belowgroundLitterCMass[layer] = 0.0;
    patch.vegetationPatch->abovegroundLitterCMass = 0.0;
    patch.vegetationPatch->belowgroundLitterNMass[layer] = 0.0;
    patch.vegetationPatch->abovegroundLitterNMass = 0.0;
}

/////////////////////////////////////////////////////////
void NutrientPatch::setFactorsNcycle(int layer)
{
    //nutrientLayer[layer].NresComposF = min(1.0, exp(-0.693 * ((nutrientLayer[layer].CNres-25)/25)));
    //nutrientLayer[layer].NcyclTF = (0.9 * nutrientLayer[layer].soilTemp / (nutrientLayer[layer].soilTemp + exp(9.93 - (0.312 * nutrientLayer[layer].soilTemp)))) + 0.1;
    //nutrientLayer[layer].NcyclWF = patch.waterPatch->waterLayer[layer].absWater / nutrientLayer[layer].fieldCapmm;
    //mineralizationRateCoef = 0.05; // as default in SWAT, varies from 0.02-0.1
    //nutrientLayer[layer].resDecayFactor = mineralizationRateCoef * nutrientLayer[layer].NresComposF * pow(nutrientLayer[layer].NcyclTF * nutrientLayer[layer].NcyclWF,0.5);
    double xx2 = 0;
    double fc = nutrientLayer[layer].fc + nutrientLayer[layer].wp;
    double wc = patch.waterPatch->waterLayer[layer].absWater + nutrientLayer[layer].wp;
    if (wc <= nutrientLayer[layer].wp) xx2 = 0.4 * wc / nutrientLayer[layer].wp;
    else if(wc <= fc) xx2 = 0.4 + 0.6 * (wc - nutrientLayer[layer].wp) / (fc - nutrientLayer[layer].wp);
    else xx2 = 1;
    nutrientLayer[layer].NcyclWF =  max( min((1 + ( 1 - xx2)/(0.25)) * pow(xx2, 4), 1.0 ), 0.0 );

    //nutrientLayer[layer].voidSoil = nutrientLayer[layer].soilPorosity * (1. - nutrientLayer[layer].awc / nutrientLayer[layer].fc); // in SWAT the last fc is actually sat, but here fc = sat
    nutrientLayer[layer].voidSoil = nutrientLayer[layer].soilPorosity * (1. - patch.waterPatch->waterLayer[layer].absWater / nutrientLayer[layer].fc); // modified from SWAT, using absolute water instead of

    double xx3 = 0;
    if (nutrientLayer[layer].voidSoil >= 0.1) xx3 = 0.2 + (nutrientLayer[layer].voidSoil - 0.1) / (nutrientLayer[layer].soilPorosity - 0.1);
    else xx3 = 0.2 * nutrientLayer[layer].voidSoil / 0.1;
    nutrientLayer[layer].NcyclOF = 0.5 + 0.5 * xx3 / (xx3 + exp(-20 * xx3));
    double tn = -5;
    double top = 35;
    double tx = 50;
    double qq = (tn - top)/(top - tx);
    nutrientLayer[layer].NcyclTF = pow((nutrientLayer[layer].soilTemp - tn), qq) * (tx - nutrientLayer[layer].soilTemp) / (pow((top - tn), qq) * (tx - top));
    nutrientLayer[layer].NcyclTF = max( min(nutrientLayer[layer].NcyclTF, 1.0 ), 0.0 );
    nutrientLayer[layer].NcyclEF = pow(nutrientLayer[layer].NcyclTF * nutrientLayer[layer].NcyclOF * nutrientLayer[layer].NcyclWF, 0.67);
    nutrientLayer[layer].NcyclEF = max( min(nutrientLayer[layer].NcyclEF, 1.0 ), 0.0 );
    // combined temperature and humidity factor 0-1
    nutrientLayer[layer].combNcyclTFWF = nutrientLayer[layer].NcyclWF * nutrientLayer[layer].NcyclTF;
}

/////////////////////////////////////////////////////////
void NutrientPatch::denitrification(int layer)
{
    double vof = 1 / (1 + pow(nutrientLayer[layer].voidSoil/0.04, 5));
    nutrientTest[layer][9] = vof;
    nutrientTest[layer][10] = nutrientLayer[layer].voidSoil;
    nutrientTest[layer][11] = -1;
    nutrientLayer[layer].denit = 0;
    if (nutrientLayer[layer].NO3 > 0) {
        nutrientLayer[layer].denit = nutrientLayer[layer].NO3 * min((1 - exp(-1.4 * nutrientLayer[layer].NcyclTF * vof * nutrientLayer[layer].CsomPerc)), 1.0);
        nutrientTest[layer][11] = nutrientLayer[layer].CsomPerc;
    }
}

/////////////////////////////////////////////////////////
void NutrientPatch::SOCdecomp(int layer) // fsol_cdec - ks Sc in Kemanian 2011
{
    double decF = min(pow(nutrientLayer[layer].CsomPerc / nutrientLayer[layer].CsomSat, 0.5), 1.0);
    double kS = turnoverRate * nutrientLayer[layer].combNcyclTFWF * decF;
    nutrientLayer[layer].decompSOC = kS * nutrientLayer[layer].Csom;
    //nutrientLayer[layer].decompSON = turnoverRate * decF * nutrientLayer[layer].combNcyclTFWF * nutrientLayer[layer].Nsom;
}

/////////////////////////////////////////////////////////
void NutrientPatch::ResDecomp(int layer) // rdc
{
    nutrientLayer[layer].decompRes = 0.05 * nutrientLayer[layer].combNcyclTFWF * nutrientLayer[layer].Cres; // kR = 0.05; optimum decomposition factor

}

/////////////////////////////////////////////////////////
void NutrientPatch::newCNsom(int layer)
{
    nutrientLayer[layer].CNres = nutrientLayer[layer].Cres / (nutrientLayer[layer].Nres + nutrientLayer[layer].NO3); // updating CNres: as in SWAT, NH4 is not taken into account*/

    //nutrientLayer[layer].CNres = max(min(nutrientLayer[layer].Cres / (nutrientLayer[layer].Nres + nutrientLayer[layer].NO3), 80.0), 5.0) ;
    nutrientLayer[layer].CNres = max(min(nutrientLayer[layer].Cres / nutrientLayer[layer].Nres, 80.0), 5.0); // limits set to 5-80
    nutrientLayer[layer].CNsom = max(min(nutrientLayer[layer].Csom / nutrientLayer[layer].Nsom, 80.0), 5.0);

    double mineralFraction = (nutrientLayer[layer].NO3 + nutrientLayer[layer].NH4) / (nutrientLayer[layer].hectarWeight * 1000000);
    nutrientLayer[layer].CNsomNew = 8.5 + 2.7 * (1 - (1 / (1 + pow(nutrientLayer[layer].CNres/110, 3)))) * (1 + (1 / (1 + pow((mineralFraction)/8, 3))));

    nutrientLayer[layer].CNsomNew = max( min( nutrientLayer[layer].CNsomNew, 14.0), 8.5); // must be within 8.5-14 Kemanian 2011
}

/////////////////////////////////////////////////////////
void NutrientPatch::ResHumific(int layer)
{
    double humificRateMax; // hx in paper
    humificRateMax = 0.09 * (2 - exp(-5.5 * (nutrientLayer[layer].pClay)));
    double humificF;
    if (nutrientLayer[layer].CsomPerc > nutrientLayer[layer].CsomSat) humificF = 0;
    else humificF = 1 - pow((nutrientLayer[layer].CsomPerc / nutrientLayer[layer].CsomSat), 6);
    nutrientLayer[layer].humificRate = 0;
    nutrientLayer[layer].humificRate = humificF * humificRateMax; // rhc

    nutrientLayer[layer].humificRes = nutrientLayer[layer].humificRate * nutrientLayer[layer].decompRes;
    nutrientLayer[layer].respRes = (1-nutrientLayer[layer].humificRate) * nutrientLayer[layer].decompRes;
    nutrientLayer[layer].humificResN = nutrientLayer[layer].humificRes / nutrientLayer[layer].CNsomNew;
}

/////////////////////////////////////////////////////////
// net mineralization
void NutrientPatch::netMineralization(int layer) // also calculates immobilization if the value is negative
{
    nutrientLayer[layer].minNetRes = 0;

    nutrientLayer[layer].minNetRes = nutrientLayer[layer].decompRes  * max( min( ((1/nutrientLayer[layer].CNres)-(nutrientLayer[layer].humificRate/nutrientLayer[layer].CNsomNew)) , 1.0), -1.0);
    nutrientLayer[layer].minSOM = (nutrientLayer[layer].decompSOC/nutrientLayer[layer].Csom) * nutrientLayer[layer].Nsom;
    // adding ks SN tp netMin

    if (nutrientLayer[layer].minNetRes < 0.0) nutrientLayer[layer].minSOM = 0.0;
}

/////////////////////////////////////////////////////////
void NutrientPatch::checkNavail(int layer)
{
    double adjustFactor = 1.0;
    double useNH4 = 0.0;
    double useNO3 = 0.0;

    if(nutrientLayer[layer].minNetRes > 0 &&
       nutrientLayer[layer].humificResN + abs(nutrientLayer[layer].minNetRes)  > nutrientLayer[layer].Nres) // mineralization, use only N res as source
    {
        adjustFactor = max( min( nutrientLayer[layer].Nres / (nutrientLayer[layer].humificResN + abs(nutrientLayer[layer].minNetRes)), 1.0), 0.0);
    }

    if(nutrientLayer[layer].minNetRes < 0 &&
         nutrientLayer[layer].humificResN + abs(nutrientLayer[layer].minNetRes)  > nutrientLayer[layer].Nres &&
         nutrientLayer[layer].humificResN + abs(nutrientLayer[layer].minNetRes) - (nutrientLayer[layer].NO3 + nutrientLayer[layer].NH4)  > nutrientLayer[layer].Nres) // immobilization larger than available N
    {
        adjustFactor = max( min( (nutrientLayer[layer].Nres + nutrientLayer[layer].NO3 + nutrientLayer[layer].NH4) / (nutrientLayer[layer].humificResN + abs(nutrientLayer[layer].minNetRes)), 1.0), 0.0);
        useNH4 = nutrientLayer[layer].NH4;
        useNO3 = nutrientLayer[layer].NO3;
    }

    if(nutrientLayer[layer].minNetRes < 0 &&
       nutrientLayer[layer].humificResN + abs(nutrientLayer[layer].minNetRes)  > nutrientLayer[layer].Nres &&
       nutrientLayer[layer].humificResN + abs(nutrientLayer[layer].minNetRes) - (nutrientLayer[layer].NO3 + nutrientLayer[layer].NH4)  <= nutrientLayer[layer].Nres) // immobilization smaller than available N, but needing mineral N.
    {
        useNH4 = min(nutrientLayer[layer].humificResN + abs(nutrientLayer[layer].minNetRes) - nutrientLayer[layer].Nres, nutrientLayer[layer].NH4);
        useNO3 = min(max(nutrientLayer[layer].humificResN + abs(nutrientLayer[layer].minNetRes) - nutrientLayer[layer].Nres - nutrientLayer[layer].NH4, 0.0), nutrientLayer[layer].NO3);
    }


   if(nutrientLayer[layer].Nres == 0)
    {
        adjustFactor = 0;
    }

    nutrientLayer[layer].humificResN *= adjustFactor;
    nutrientLayer[layer].humificRes *= adjustFactor;
    nutrientLayer[layer].respRes *= adjustFactor;
    nutrientLayer[layer].decompRes *= adjustFactor;
    nutrientLayer[layer].minNetRes *= adjustFactor;

    nutrientLayer[layer].NH4 -= useNH4;
    nutrientLayer[layer].NO3 -= useNO3;
    nutrientLayer[layer].Nsom += useNH4 + useNO3;

    nutrientTest[layer][0] = useNO3;
    nutrientTest[layer][1] = nutrientLayer[layer].minNetRes;
    nutrientTest[layer][7] = adjustFactor;
    nutrientTest[layer][8] = nutrientLayer[layer].humificResN;

    nutrientLayer[layer].debitNetMin = 0.0;
    nutrientLayer[layer].debitNetMin = useNH4 + useNO3;
}

/////////////////////////////////////////////////////////
void NutrientPatch::nitrificationVolatilization(int layer)
{
    nutrientLayer[layer].nitrified = 0;
    nutrientLayer[layer].volatilized = 0;
    double nitriTfactor = 0.41 * ((nutrientLayer[layer].soilTemp - 5) / 10);
    //nutrientLayer[layer].NH4 = 10;
    if (nutrientLayer[layer].NH4 > 0 && nitriTfactor >= 0.001 && nutrientLayer[layer].soilTemp > 5)
    {
        double nitriWfactor = 0;

        if (patch.waterPatch->waterLayer[layer].absWater < 0.25 * nutrientLayer[layer].humificResN - 0.75 * nutrientLayer[layer].wp)
        {
            nitriWfactor = max((patch.waterPatch->waterLayer[layer].absWater - nutrientLayer[layer].wp) / (0.25 * (nutrientLayer[layer].fc - nutrientLayer[layer].wp)), 0.0);
        }
        else nitriWfactor = 1;

        double volatDfactor = 1 - (nutrientLayer[layer].meanDepthLayer / (nutrientLayer[layer].meanDepthLayer + exp(4.706 - 0.0305 * nutrientLayer[layer].meanDepthLayer)));
        double nitrReg = nitriTfactor * nitriWfactor;
        double volatReg = nitriTfactor * volatDfactor * 0.15; // cation exchange factor = 0.15
        double NnitrVolat = nutrientLayer[layer].NH4 * (1 - exp(-nitrReg - volatReg));
        double fracNitr = 1 - exp(-nitrReg);
        double fracVolat = 1 - exp(-volatReg);
        if(fracNitr > 0.0 || fracVolat > 0.0) {
            nutrientLayer[layer].nitrified = (fracNitr/(fracNitr + fracVolat)) * NnitrVolat;
            nutrientLayer[layer].volatilized = (fracVolat/(fracNitr + fracVolat)) * NnitrVolat;
        }
        else {
            nutrientLayer[layer].nitrified = 0.0;
            nutrientLayer[layer].volatilized = 0.0;
        }
    }
}

/////////////////////////////////////////////////////////
// dissolved NO3 available for leaching
void NutrientPatch::leachingNO3(int layer)
{
    double saturationH20 = soilProp.fieldCap * nutrientLayer[layer].thicknessLayer;
    if (layer < 1)
    {
        nutrientLayer[layer].mobileH2O = patch.waterPatch->totalRunOff + max(0.0, patch.waterPatch->waterLayer[layer].diffusedWater);
        //nutrientLayer[layer].mobileH2O = max(0.0, patch.waterPatch->waterLayer[layer].diffusedWater);
    }
    else nutrientLayer[layer].mobileH2O = max(0.0, patch.waterPatch->waterLayer[layer].drainedWater);

    if (nutrientLayer[layer].mobileH2O > 0) nutrientLayer[layer].concAvailNO3 = ((nutrientLayer[layer].NO3) * (1 - exp(-nutrientLayer[layer].mobileH2O / ((1-0.5) * saturationH20) ))) / nutrientLayer[layer].mobileH2O;
    else nutrientLayer[layer].concAvailNO3 = 0;

    if (layer < 1) {
            nutrientLayer[layer].leaching_out = nutrientLayer[layer].mobileH2O * nutrientLayer[layer].concAvailNO3;
            if(patch.waterPatch->waterLayer[layer].diffusedWater < 0) {
                nutrientLayer[layer].leaching_in = patch.waterPatch->waterLayer[layer].diffusedWater * nutrientLayer[layer + 1].concAvailNO3;
            }
    }
    else {
        nutrientLayer[layer].leaching_out = min(patch.waterPatch->waterLayer[layer-1].diffusedWater * nutrientLayer[layer].concAvailNO3, 0.0) +
        patch.waterPatch->waterLayer[layer].drainedWater * nutrientLayer[layer].concAvailNO3;

        nutrientLayer[layer].leaching_in = min(patch.waterPatch->waterLayer[layer-1].diffusedWater * nutrientLayer[layer-1].concAvailNO3, 0.0);
    }
}

/////////////////////////////////////////////////////////
void NutrientPatch::processesToPools(int layer)
{
    nutrientTest[layer][2] = nutrientLayer[layer].leaching_in;
    nutrientTest[layer][3] = nutrientLayer[layer].leaching_out;
    nutrientTest[layer][4] = nutrientLayer[layer].denit;
    nutrientTest[layer][5] = nutrientLayer[layer].nitrified;
    nutrientTest[0][6] = dailyNO3deposition;

    nutrientLayer[layer].Cres += nutrientLayer[layer].CresNew - nutrientLayer[layer].respRes - nutrientLayer[layer].humificRes;
    nutrientLayer[layer].Csom +=  nutrientLayer[layer].humificRes - nutrientLayer[layer].decompSOC;

    nutrientLayer[layer].Nres += nutrientLayer[layer].NresNew - nutrientLayer[layer].humificResN - abs(nutrientLayer[layer].minNetRes);

    if (nutrientLayer[layer].minNetRes < 0) // immobilization
    {
        nutrientLayer[layer].Nsom += abs(nutrientLayer[layer].minNetRes) + nutrientLayer[layer].humificResN;
    }

      if (nutrientLayer[layer].minNetRes > 0) // mineralization
    {
        nutrientLayer[layer].Nsom += nutrientLayer[layer].humificResN - nutrientLayer[layer].minSOM;
        nutrientLayer[layer].NH4 += nutrientLayer[layer].minSOM + abs(nutrientLayer[layer].minNetRes);

    }

    // denitrification
    nutrientLayer[layer].NO3 = max(nutrientLayer[layer].NO3 - nutrientLayer[layer].denit, 0.0);
    if (nutrientLayer[layer].denit > nutrientLayer[layer].NO3) nutrientLayer[layer].denit = nutrientLayer[layer].NO3;

    // nitrification and volatilization
    nutrientLayer[layer].NH4 -= nutrientLayer[layer].nitrified - nutrientLayer[layer].volatilized;
    nutrientLayer[layer].NO3 += nutrientLayer[layer].nitrified;

    // leaching
    nutrientLayer[layer].NO3 += max(nutrientLayer[layer].leaching_in - nutrientLayer[layer].leaching_out, 0.0);

    //daily deposition
    nutrientLayer[0].NO3 += dailyNO3deposition;
    nutrientLayer[0].NH4 += dailyNH4deposition;

    //nutrientLayer[layer].CNres = max(min(nutrientLayer[layer].Cres / (nutrientLayer[layer].Nres + nutrientLayer[layer].NO3), 80.0), 5.0) ;
    nutrientLayer[layer].CNres = max(min(nutrientLayer[layer].Cres / nutrientLayer[layer].Nres, 100.0), 5.0) ;
    nutrientLayer[layer].CNsom = max(min(nutrientLayer[layer].Csom / nutrientLayer[layer].Nsom, 100.0), 5.0);

    nutrientLayer[layer].CsomPerc = nutrientLayer[layer].Csom * 100 / nutrientLayer[layer].hectarWeight;
    nutrientLayer[layer].NsomPerc = nutrientLayer[layer].Nsom * 100 / nutrientLayer[layer].hectarWeight;

    nutrientLayer[layer].NO3_weight = nutrientLayer[layer].NO3 * 0.0001 * pow(input.cellSize, 2);
    nutrientLayer[layer].NH4_weight = nutrientLayer[layer].NH4 * 0.0001 * pow(input.cellSize, 2);

    nutrientLayer[layer].minNetRes += nutrientLayer[layer].debitNetMin;

}

/////////////////////////////////////////////////////////
NutrientPatch::NutrientPatch(int index, int _xCor, int _yCor, const Input& _input, const Patch& _patch) : xCor(_xCor), yCor(_yCor), patch(_patch), input(_input)
{
    soilProp =  {
                 input.soilProperties[index][0], atof(input.soilProperties[index][1].c_str()), atof(input.soilProperties[index][2].c_str())*24.0, atof(input.soilProperties[index][3].c_str()),
                 atof(input.soilProperties[index][4].c_str()), atof(input.soilProperties[index][5].c_str()), atof(input.soilProperties[index][6].c_str()), atof(input.soilProperties[index][7].c_str())*24.0,
                 atof(input.soilProperties[index][8].c_str()), atof(input.soilProperties[index][9].c_str()), atof(input.soilProperties[index][10].c_str()),  atof(input.soilProperties[index][11].c_str()),
                 atof(input.soilProperties[index][12].c_str()), atof(input.soilProperties[index][13].c_str())
                };

    dampingDepthMaximum = 1000 + ((2500 * soilProp.bulkDensity) / (soilProp.bulkDensity + 686 * exp(-5.63 * soilProp.bulkDensity)));
    double TD = 0;
    for (int layer = 0; layer < input.nSoilLayers; layer++) {
            TD += input.depthLayers[layer];
    }
    totalDepth = TD;
    TsoilSurfBare = input.MAT[0];
    TsoilSurf = input.MAT[0];
    TsoilLag = 0.80;
    turnoverRate = 0.000123; // daily turnover rate, approx 4.5% annual turnover
    //residue = 0;
    dailyNO3deposition = 0.5 * input.dailyNdeposition;
    dailyNH4deposition = 0.5 * input.dailyNdeposition;

    //Create nutrient layers per patch
    nutrientLayer.resize(input.nSoilLayers);
    NutrientLayer helpNLayer;
    for(unsigned int i = 0; i < nutrientLayer.size(); i++) {
        helpNLayer = {(int)i,
                      input.depthLayers[i]*1.0,
                      0,
                      0};

        for (int layer = 0; layer <= helpNLayer.myLayer; layer++)
        {
            helpNLayer.maxDepthLayer += input.depthLayers[layer];
        }

        helpNLayer.meanDepthLayer = helpNLayer.maxDepthLayer - (helpNLayer.thicknessLayer / 2);
        helpNLayer.soilWaterScalingFactor = 0;
        helpNLayer.soilDepthTF = 0;
        helpNLayer.soilTemp = input.MAT[0];
        helpNLayer.dampingDepth = 0;
        helpNLayer.dampingDepthRatio = 0;
        helpNLayer.weight = pow(input.cellSize, 2) * (helpNLayer.thicknessLayer/1000) * (soilProp.bulkDensity*1000); // soil layer weight in kg
        helpNLayer.hectarWeight = helpNLayer.weight * (pow(100,2) / pow(input.cellSize, 2));
        helpNLayer.Cres = 1*60.0;
        helpNLayer.Nres = 1;
        helpNLayer.CNres = 60.0;
        helpNLayer.CNsomNew = 0;
        helpNLayer.Csom = (3.22/100) * helpNLayer.hectarWeight; // results from soil sampling (percentage), upscalled to kg ha-1
        helpNLayer.Nsom = (0.26/100) * helpNLayer.hectarWeight; // results from soil sampling (percentage), upscalled to kg ha-1
        helpNLayer.CsomPerc = 0;
        helpNLayer.NsomPerc = 0;
        helpNLayer.CNsom = helpNLayer.Csom/helpNLayer.Nsom;
        helpNLayer.fc = soilProp.fieldCap * helpNLayer.thicknessLayer;
        helpNLayer.soilPorosity = 1 - soilProp.bulkDensity/2.65;
        helpNLayer.pClay = soilProp.clayContent;
        helpNLayer.WPcalc = 0.4 * helpNLayer.pClay * soilProp.bulkDensity;
        if(helpNLayer.WPcalc <= 0) helpNLayer.WPcalc = 0.005;
        helpNLayer.AWC = soilProp.fieldCap - helpNLayer.WPcalc;
        helpNLayer.awc = helpNLayer.AWC * helpNLayer.thicknessLayer;
        helpNLayer.wp = helpNLayer.WPcalc * helpNLayer.thicknessLayer;
        helpNLayer.CsomSat = helpNLayer.hectarWeight * (2.11 + 0.0375 * helpNLayer.pClay) / 100; //
        nutrientLayer[i] = helpNLayer;
    }

    // correcting C and N values according to the FAO dataset
    nutrientLayer[0].Csom = (0.9/100) * helpNLayer.hectarWeight;
    nutrientLayer[0].Nsom = nutrientLayer[0].Csom/14;
    nutrientLayer[1].Csom = (0.3/100) * helpNLayer.hectarWeight;
    nutrientLayer[1].Nsom = nutrientLayer[1].Csom/14;
    nutrientLayer[0].CsomPerc = 0.9; // as percentage, like in SWAT 5.2.2018
    nutrientLayer[0].NsomPerc = nutrientLayer[0].CsomPerc/14;
    nutrientLayer[1].CsomPerc = 0.3;
    nutrientLayer[1].NsomPerc = nutrientLayer[1].CsomPerc/14;
    nutrientLayer[0].NH4 = input.iniNH4 * nutrientLayer[0].hectarWeight / 1000000; // from mg/kg to kg/ha
    nutrientLayer[1].NH4 = input.iniNH4 * nutrientLayer[0].hectarWeight / 1000000; // from mg/kg to kg/ha
    nutrientLayer[0].NO3 = input.iniNO3 * nutrientLayer[0].hectarWeight / 1000000; // from mg/kg to kg/ha
    nutrientLayer[1].NO3 = input.iniNO3 * nutrientLayer[1].hectarWeight / 1000000; // from mg/kg to kg/ha
    nutrientLayer[0].NO3_weight = nutrientLayer[0].NO3 * 0.0001 * pow(input.cellSize, 2); // from kg/kg to kg/cell
    nutrientLayer[0].NH4_weight = nutrientLayer[1].NH4 * 0.0001 * pow(input.cellSize, 2); // from kg/kg to kg/cell
    nutrientLayer[1].NO3_weight = nutrientLayer[0].NO3 * 0.0001 * pow(input.cellSize, 2); // from kg/kg to kg/cell
    nutrientLayer[1].NH4_weight = nutrientLayer[1].NH4 * 0.0001 * pow(input.cellSize, 2); // from kg/kg to kg/cell
    nutrientTest.resize(2, vector<double>(12));
}

/////////////////////////////////////////////////////////
NutrientPatch::~NutrientPatch()
{

};


