#include "Landscape.h"
#include "Input.h"
#include "Patch.h"
#include "Input.h"
#include "WaterLandscape.h"
#include "NutrientLandscape.h"
#include "NutrientPatch.h"
#include "VegetationLandscape.h"
#include "Output.h"
#include "cmath"

/////////////////////////////////////////////////////////
void Landscape::calculateProcesses(int day)
{
    vegetation->shadingAPAR(extraterrestrialRadiation(day));
    vegetation->calculatePatchVariables();
    water->precipitation(day);
    water->infiltrationAndDrainage();
    water->runoff();
    water->diffusion();
    water->evaporation(day);
    nutrient->temperatureSurf(day);
    nutrient->temperatureGrid(day);
    nutrient->residualIncomeGrid();
    nutrient->setFactorsNcycleGrid();
    nutrient->denitrificationGrid();
    nutrient->decompositionGrid();
    nutrient->newCNsomGrid();
    nutrient->ResHumificGrid();
    nutrient->netMineralizationGrid();
    nutrient->checkNavailGrid();
    nutrient->nitrificationVolatilizationGrid();
    nutrient->leachingNO3Grid();
    nutrient->processesToPoolsGrid();
    vegetation->calculateProcesses(day, dayLength(day));
    water->meanLandscape();
    vegetation->calculatePatchVariables();
    output->writeFiles(day, *vegetation, *nutrientPatch, *water); // vegetation patch variables not updated to this point (beginning of next day instead)
    vegetation->deleteDeadPlants();
}

/////////////////////////////////////////////////////////
double Landscape::extraterrestrialRadiation(int day)
{
    double extraRadiation;
    double conversionFactor = 2450000.0/86400.0;

    if(input.measuredSolRadiation)
        extraRadiation = input.solRadiation[day] / conversionFactor; // used measured solar radiation and transform from W*m-2 (J*m-2*s-1) into mm*day-1

    else {
        //calculation of evaporation equivalent of extraterrestrial evaporation (after Tietjen et al. 2009)
        //routine found in Maidment "Handbook of Hydrology", page 4.31
        double distEarthSun = 1.0 + 0.033 * cos( 2.0 * M_PI / 365 * input.dayOfYear[day-1]);
        double solDeclination = 0.4093 * sin( 2.0 * M_PI / 365 * input.dayOfYear[day-1] - 1.405);
        double sunSetHour = acos( -1.0 * tan( input.latitude * M_PI/180.0 ) * tan( solDeclination ) );
        extraRadiation = 15.392 * distEarthSun * ( sunSetHour * sin( input.latitude * M_PI/180.0 ) * sin( solDeclination ) + cos(input.latitude * M_PI/180.0 ) * cos( solDeclination ) * sin( sunSetHour ) );
    }

    return extraRadiation;
}

/////////////////////////////////////////////////////////
double Landscape::dayLength(int day)
{
    //Sitch et al. 2000, Appendix eqn. 8-10
    double solDeclination = 0.4093 * sin( 2.0 * M_PI / 365 * input.dayOfYear[day-1] - 1.405); // (after Tietjen et al. 2009)
    double u = sin(input.latitude * M_PI/180.0) * sin(solDeclination);
    double v = cos(input.latitude * M_PI/180.0) * cos(solDeclination);
    double dayLen = 24/M_PI * acos(-u/v); // day length in [hours]

    return dayLen;
}

/////////////////////////////////////////////////////////
Landscape::Landscape(const Input& _input):
     input(_input)
{
    //Create grid
    grid.resize(input.xCells, vector<Patch*>(input.yCells));
    for(unsigned int x = 0; x < grid.size(); x++) {
        for(unsigned int y = 0; y < grid.at(x).size(); y++) {
            for(unsigned int i = 0; i < input.soilProperties.at(i).size(); i++) {
                if(input.soilProperties[i][0] == input.soilType[x][y]) {
                   grid[x][y] = new Patch(i, x, y, input);
                   break;
                }
            }
        }
    }

    //Create object
    water = new WaterLandscape(input, grid, this);
    nutrient = new NutrientLandscape(input, grid);
    vegetation = new VegetationLandscape(input, grid, this);
    output = new Output(input, grid);
}

/////////////////////////////////////////////////////////
Landscape::~Landscape()
{
    delete water;
    water = 0;
    delete nutrient;
    nutrient = 0;
    delete vegetation;
    vegetation = 0;
    delete output;
    output = 0;

    //Delete objects
    for(unsigned int x = 0; x < grid.size(); x++) {
        for(unsigned int y = 0; y < grid.at(x).size(); y++) {
            delete grid[x][y];
            grid[x][y] = 0;
        }
    }
}
