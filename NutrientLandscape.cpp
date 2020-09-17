#include "NutrientLandscape.h"
#include "NutrientPatch.h"
#include "Patch.h"
#include "Input.h"
#include "cmath"
#include <iostream>

/////////////////////////////////////////////////////////
/*void NutrientLandscape::potSolarRadiation(int day)
{
    //from SWAT
    double relativeDistanceFromSun = 1.0 + 0.033 * cos( 2.0 * M_PI  * (day+1) / 365);
    double radLatitude = input.latitude * M_PI / 180;
    double solarDeclination = asin(0.4 * sin(((day+1) - 82) * 2 * M_PI /365));

    //0.4093 * sin( 2.0 * M_PI / 365 * (day+1) - 1.405)

    double avEr = 0.2618; // angular velocity of Earth rotation
    double sunRiseTime = acos(-tan(solarDeclination) * tan(radLatitude)) / avEr;

    potSolRadiation = 30 * relativeDistanceFromSun * (avEr * sunRiseTime * sin(solarDeclination) * sin (radLatitude) + cos(solarDeclination) * cos(radLatitude) * sin(avEr * sunRiseTime));

    //return potSolRadiation;
}

void NutrientLandscape::actSolarRadiation(int year, int day, double potSolRadiation)
{
    // equations from Liu and Scott 2001 Agricultural and Forest Meteorology https://doi.org/10.1016/S0168-1923(00)00173-8
    double TempEffect = 0.0;

    if (day < 365) TempEffect = input.tempMax[year][day] - ( (input.tempMin[year][day] + input.tempMin[year][day+1]) / 2 );
    else TempEffect = input.tempMax[year][day] - input.tempMin[year][day];

    double a = 0.717;
    double b = 0.028;
    double c= 1.67;

    actSolRad = potSolRadiation * a * (1-exp(-b * pow(TempEffect, c)));

    //return actSolRad;
}*/

/////////////////////////////////////////////////////////
void NutrientLandscape::temperatureSurf(int day)
{
   for(int x = 0; x < xCells; x++) {
        for(int y = 0; y < yCells; y++) {
                grid[x][y]->nutrientPatch->calcSoilSurfTemperature(day);
        }
    }
}

/////////////////////////////////////////////////////////
void NutrientLandscape::temperatureGrid(int day)
{
   for(int x = 0; x < xCells; x++) {
        for(int y = 0; y < yCells; y++) {
            for(int layer = 0; layer < input.nSoilLayers; layer++) {
                grid[x][y]->nutrientPatch->calcSoilTemperature(layer, day);
            }
        }
    }
}

/////////////////////////////////////////////////////////
void NutrientLandscape::residualIncomeGrid()
{
   for(int x = 0; x < xCells; x++) {
        for(int y = 0; y < yCells; y++) {
            for(int layer = 0; layer < input.nSoilLayers; layer++) {
                grid[x][y]->nutrientPatch->incomeRes(layer);
            }
        }
    }
}

/////////////////////////////////////////////////////////
void NutrientLandscape::setFactorsNcycleGrid()
{
   for(int x = 0; x < xCells; x++) {
        for(int y = 0; y < yCells; y++) {
            for(int layer = 0; layer < input.nSoilLayers; layer++) {
                grid[x][y]->nutrientPatch->setFactorsNcycle(layer);
            }
        }
    }
}

/////////////////////////////////////////////////////////
void NutrientLandscape::denitrificationGrid()
{
   for(int x = 0; x < xCells; x++) {
        for(int y = 0; y < yCells; y++) {
            for(int layer = 0; layer < input.nSoilLayers; layer++) {
                grid[x][y]->nutrientPatch->denitrification(layer);
            }
        }
    }
}

/////////////////////////////////////////////////////////
void NutrientLandscape::decompositionGrid()
{
   for(int x = 0; x < xCells; x++) {
        for(int y = 0; y < yCells; y++) {
            for(int layer = 0; layer < input.nSoilLayers; layer++) {
                grid[x][y]->nutrientPatch->SOCdecomp(layer);
                grid[x][y]->nutrientPatch->ResDecomp(layer);
            }
        }
    }
}

/////////////////////////////////////////////////////////
void NutrientLandscape::newCNsomGrid()
{
   for(int x = 0; x < xCells; x++) {
        for(int y = 0; y < yCells; y++) {
            for(int layer = 0; layer < input.nSoilLayers; layer++) {
                grid[x][y]->nutrientPatch->newCNsom(layer);
            }
        }
    }
}

/////////////////////////////////////////////////////////
void NutrientLandscape::ResHumificGrid()
{
   for(int x = 0; x < xCells; x++) {
        for(int y = 0; y < yCells; y++) {
            for(int layer = 0; layer < input.nSoilLayers; layer++) {
                grid[x][y]->nutrientPatch->ResHumific(layer);
            }
        }
    }
}

/////////////////////////////////////////////////////////
void NutrientLandscape::netMineralizationGrid()
{
   for(int x = 0; x < xCells; x++) {
        for(int y = 0; y < yCells; y++) {
            for(int layer = 0; layer < input.nSoilLayers; layer++) {
                grid[x][y]->nutrientPatch->netMineralization(layer);
            }
        }
    }
}

/////////////////////////////////////////////////////////
void NutrientLandscape::checkNavailGrid()
{
   for(int x = 0; x < xCells; x++) {
        for(int y = 0; y < yCells; y++) {
            for(int layer = 0; layer < input.nSoilLayers; layer++) {
                grid[x][y]->nutrientPatch->checkNavail(layer);
            }
        }
    }
}

/////////////////////////////////////////////////////////
void NutrientLandscape::nitrificationVolatilizationGrid()
{
   for(int x = 0; x < xCells; x++) {
        for(int y = 0; y < yCells; y++) {
            for(int layer = 0; layer < input.nSoilLayers; layer++) {
                grid[x][y]->nutrientPatch->nitrificationVolatilization(layer);
            }
        }
    }
}

/////////////////////////////////////////////////////////
void NutrientLandscape::leachingNO3Grid()
{
   for(int x = 0; x < xCells; x++) {
        for(int y = 0; y < yCells; y++) {
            for(int layer = 0; layer < input.nSoilLayers; layer++) {
                grid[x][y]->nutrientPatch->leachingNO3(layer);
            }
        }
    }
}

/////////////////////////////////////////////////////////
void NutrientLandscape::processesToPoolsGrid()
{
   for(int x = 0; x < xCells; x++) {
        for(int y = 0; y < yCells; y++) {
            for(int layer = 0; layer < input.nSoilLayers; layer++) {
                grid[x][y]->nutrientPatch->processesToPools(layer);
            }
        }
    }
}

//////////////////////////////////////////////////
/// \brief NutrientLandscape::NutrientLandscape
/// \param input
/// \param _grid
NutrientLandscape::NutrientLandscape(const Input& _input, vector< vector<Patch*> >& _grid):
    xCells(_input.xCells), yCells(_input.yCells), grid(_grid), input(_input)
{

}

/////////////////////////////////////////////////////////
NutrientLandscape::~NutrientLandscape()
{

}
