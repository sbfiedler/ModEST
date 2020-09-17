#include "Output.h"
#include "Input.h"
#include "WaterLandscape.h"
#include "NutrientLandscape.h"
#include "VegetationLandscape.h"
#include "Plant.h"
#include "NutrientPatch.h"
#include "Patch.h"
#include <string>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <iomanip>

/////////////////////////////////////////////////////////
void Output::initialiseWaterFile()
{
    string dailyWaterFileName;

    //create and open result file
    if(exe)
        dailyWaterFileName = "Output" + delimiter + "DailyWater.txt";
    else
        dailyWaterFileName = ".." + delimiter + "ModEST" + delimiter + "Output" + delimiter + "DailyWater.txt";

    cout << "Write daily water into: " << dailyWaterFileName << endl;
        cout << endl;

    dailyWaterFile.open(dailyWaterFileName.c_str(), ios::out | ios::trunc);

    //if something goes wrong...
    if(!dailyWaterFile){
        cerr << dailyWaterFileName
             << " could not be opened!\n";
        exit(-1);
    }

    dailyWaterFile 	<< "Date" << "\t" << "xCor" << "\t" << "yCor" << "\t" << "SurfaceWater[mm]";

    for(int i = 0; i < input.nSoilLayers; i++) {
       dailyWaterFile << "\t" << "RelWaterL" + to_string(i+1) + "[vol%]" << "\t" << "AbsWaterL" + to_string(i+1) + "[mm]";
    }

    dailyWaterFile << "\n";
}

/////////////////////////////////////////////////////////
void Output::writeWaterFile(int day, const WaterLandscape& water)
{

    //dailyWaterFile 	<< input.date[day] << "\t" << water.meanAbsSurfaceWater;

    for(int x = 0; x < input.xCells; x++) {
        for(int y = 0; y < input.yCells; y++) {
            dailyWaterFile << input.date[day] << "\t" << x << "\t" << y << "\t" << grid[x][y]->waterPatch->surfaceWater;
            for(int i = 0; i < input.nSoilLayers; i++) {
                //dailyWaterFile << "\t" << water.meanRelWater[i] << "\t" << water.meanAbsWater[i];
                dailyWaterFile << "\t" << grid[x][y]->waterPatch->waterLayer[i].relWater << "\t" << grid[x][y]->waterPatch->waterLayer[i].absWater;
            }
            dailyWaterFile << "\n";
        }
    }

    //dailyWaterFile << "\t" << water.meanDeepDrainedWater << "\n";

    if(day == input.simDays) dailyWaterFile.close();
}

/////////////////////////////////////////////////////////
void Output::initialiseNutrientFile()
{
    string dailyNutrientFileName;

    //create and open result file
    if(exe)
        dailyNutrientFileName = "Output" + delimiter + "DailyNutrient.txt";
    else
        dailyNutrientFileName = ".." + delimiter + "ModEST" + delimiter + "Output" + delimiter + "DailyNutrient.txt";

    cout << "Write daily nutrients into: " << dailyNutrientFileName << endl;
        cout << endl;

    dailyNutrientFile.open(dailyNutrientFileName.c_str(), ios::out | ios::trunc);

    //if something goes wrong...
    if(!dailyNutrientFile){
        cerr << dailyNutrientFileName
             << " could not be opened!\n";
        exit(-1);
    }

    dailyNutrientFile 	<< "Date" << "\t" << "CoordX" << "\t" << "CoordY" << "\t" << "SurfTempBareSoil" << "\t" << "SurfTempSoil";

//    dailyNutrientFile << "\t" << "Cresidual[kg/ha]" << "\t" << "Nresidual[kg/ha]" << "\t" << "Denitrified[kg/ha]"  << "\t" << "DecomposedResidual[kg/ha]" << "\t" << "DecomposedSOC[kg/ha]" << "\t" << "HumifiedResidual[kg/ha]" << "\t" << "RespirationResidual[kg/ha]" << "\t" << "NetMineralization[kg/ha]";

    for(int i = 0; i < input.nSoilLayers; i++) {
        dailyNutrientFile << "\t" <<
        "SoilTemp_L" + to_string(i+1) << "\t" <<
        "Csom_L" + to_string(i+1) << "\t" <<
        "Nsom_L" + to_string(i+1) << "\t" <<
        "NO3_L" + to_string(i+1) << "\t" <<
        "NH4_L" + to_string(i+1) << "\t" <<
        "Cresidual_L" + to_string(i+1) << "\t" <<
        "Nresidual_L" + to_string(i+1) << "\t" <<
        "Denitrified_L" + to_string(i+1) << "\t" <<
        "DecomposedResidual_L" + to_string(i+1) << "\t" <<
        "DecomposedSOC_L" + to_string(i+1) << "\t" <<
        "HumifiedResidual_L" + to_string(i+1) << "\t" <<
        "HumifiedResidualN_L" + to_string(i+1) << "\t" <<
        "RespirationResidual_L" + to_string(i+1) << "\t" <<
        "NetMineralization_L" + to_string(i+1) << "\t" <<
        "CNsomNew_L" + to_string(i+1) << "\t" <<
        "CNsom_L" + to_string(i+1) << "\t" <<
        "CNres_L" + to_string(i+1);

    }

    dailyNutrientFile 	<< "\n";
}

/////////////////////////////////////////////////////////
void Output::writeNutrientFile(int day, const NutrientPatch& nutrientPatch)
{
    double oldCsom;
    for(int x = 0; x < input.xCells; x++) {
        for(int y = 0; y < input.yCells; y++) {
            dailyNutrientFile 	<< input.date[day] << "\t" << x << "\t" << y << "\t" << grid[x][y]->nutrientPatch->TsoilSurfBare << "\t" << grid[x][y]->nutrientPatch->TsoilSurf;

            for(int i = 0; i < input.nSoilLayers; i++) {
                if(day == 0) oldCsom = grid[x][y]->nutrientPatch->nutrientLayer[i].Csom;
                double dailyCsomInc = oldCsom - grid[x][y]->nutrientPatch->nutrientLayer[i].Csom;
                oldCsom = grid[x][y]->nutrientPatch->nutrientLayer[i].Csom;

                dailyNutrientFile << "\t" <<
                grid[x][y]->nutrientPatch->nutrientLayer[i].soilTemp  << "\t" <<
                dailyCsomInc << "\t" <<
                grid[x][y]->nutrientPatch->nutrientLayer[i].Nsom << "\t" <<
                grid[x][y]->nutrientPatch->nutrientLayer[i].NO3_weight << "\t" <<
                grid[x][y]->nutrientPatch->nutrientLayer[i].NH4_weight << "\t" <<
                grid[x][y]->nutrientPatch->nutrientLayer[i].Cres << "\t" <<
                grid[x][y]->nutrientPatch->nutrientLayer[i].Nres << "\t" <<
                grid[x][y]->nutrientPatch->nutrientLayer[i].denit << "\t" <<
                grid[x][y]->nutrientPatch->nutrientLayer[i].decompRes << "\t" <<
                grid[x][y]->nutrientPatch->nutrientLayer[i].decompSOC << "\t" <<
                grid[x][y]->nutrientPatch->nutrientLayer[i].humificRes << "\t" <<
                grid[x][y]->nutrientPatch->nutrientLayer[i].humificResN << "\t" <<
                grid[x][y]->nutrientPatch->nutrientLayer[i].respRes << "\t" <<
                grid[x][y]->nutrientPatch->nutrientLayer[i].minNetRes << "\t" <<
                grid[x][y]->nutrientPatch->nutrientLayer[i].CNsomNew << "\t" <<
                grid[x][y]->nutrientPatch->nutrientLayer[i].CNsom << "\t" <<
                grid[x][y]->nutrientPatch->nutrientLayer[i].CNres;
            }
            dailyNutrientFile << "\n";
        }
    }
}

/////////////////////////////////////////////////////////
void Output::initialisePlantFile()
{
    string dailyPlantFileName;

    //create and open result file   
    if(exe)
        dailyPlantFileName = "Output" + delimiter + "DailyPlant.txt";
    else
       dailyPlantFileName = ".." + delimiter + "ModEST" + delimiter + "Output" + delimiter + "DailyPlant.txt";

    cout << "Write daily plants into: " << dailyPlantFileName << endl;
        cout << endl;

    dailyPlantFile.open(dailyPlantFileName.c_str(), ios::out | ios::trunc);

    //if something goes wrong...
    if(!dailyPlantFile){
        cerr << dailyPlantFileName
             << " could not be opened!\n";
        exit(-1);
    }

    dailyPlantFile 	<< "Date" << "\t" << "ID\t" << "PFT\t" << "xCor[m]\t" << "yCor[m]\t" << "GPP[kgC*day-1]\t" << "NPP[kgC*day-1]\t" << "height[m]\t" << "leafCMass[kg]\t" << "rootCMass[kg]\t" << "sapwoodCMass[kg]\t" << "heartwoodCMass[kg]\t" << "leafNMass[kg]\t" << "rootNMass[kg]\t" << "sapwoodNMass[kg]\t"  << "litterCAboveground[kg]\t" << "litterCBelowground[kg]\t" << "litterNAboveground[kg]\t" << "litterNBelowground[kg]\t" << "dailyCIncrement[kg]" << "\n";
}

/////////////////////////////////////////////////////////
void Output::writePlantFile(int day, const VegetationLandscape& plants)
{
    for(int i = 0; i < plants.plantList.size(); i++) {
        dailyPlantFile << input.date[day] << "\t" << plants.plantList[i]->myID << "\t" << plants.plantList[i]->PFTname << "\t" << plants.plantList[i]->xCor << "\t" << plants.plantList[i]->yCor;
        dailyPlantFile << "\t" << plants.plantList[i]->GPP << "\t" << plants.plantList[i]->NPP << "\t" << plants.plantList[i]->height << "\t" << plants.plantList[i]->leafCMass << "\t" << plants.plantList[i]->rootCMass << "\t" << plants.plantList[i]->sapwoodCMass << "\t" << plants.plantList[i]->heartwoodCMass << "\t" << plants.plantList[i]->leafNMass << "\t" << plants.plantList[i]->rootNMass << "\t" << plants.plantList[i]->sapwoodNMass << "\t" << plants.plantList[i]->abovegroundCLitter << "\t" << plants.plantList[i]->belowgroundCLitter << "\t" << plants.plantList[i]->abovegroundNLitter << "\t" << plants.plantList[i]->belowgroundNLitter << "\t" << plants.plantList[i]->dailyCInc << endl;
    }

    if(day == input.simDays) dailyPlantFile.close();
}

/////////////////////////////////////////////////////////
void Output::initialiseVegetationFile()
{
    string dailyVegetationFileName;

    //create and open result file
    if(exe)
        dailyVegetationFileName = "Output" + delimiter + "DailyVegetation.txt";
    else
        dailyVegetationFileName = ".." + delimiter + "ModEST" + delimiter + "Output" + delimiter + "DailyVegetation.txt";

    cout << "Write daily vegetation into: " << dailyVegetationFileName << endl;
        cout << endl;

    dailyVegetationFile.open(dailyVegetationFileName.c_str(), ios::out | ios::trunc);

    //if something goes wrong...
    if(!dailyVegetationFile){
        cerr << dailyVegetationFileName
             << " could not be opened!\n";
        exit(-1);
    }

    dailyVegetationFile << "Date" << "\t" << "xCor" << "\t" << "yCor" << "\t" << "vegetationCover[-]";
    dailyVegetationFile  << "\t" << "abovegroundCMass[kg]" << "\t";
    for(int i = 0; i < input.nSoilLayers; i++) {
       dailyVegetationFile  << "\t" << "belowgroundCMassL" + to_string(i+1) + "[kg]";
    }
    dailyVegetationFile  << "\t" << "abovegroundNMass[kg]";
    for(int i = 0; i < input.nSoilLayers; i++) {
       dailyVegetationFile  << "\t" << "belowgroundNMassL" + to_string(i+1) + "[kg]";
    }
    dailyVegetationFile << "\n";
}

/////////////////////////////////////////////////////////
void Output::writeVegetationFile(int day, const VegetationLandscape& vegetation)
{

    for(int x = 0; x < input.xCells; x++) {
        for(int y = 0; y < input.yCells; y++) {
            dailyVegetationFile << input.date[day] << "\t" << x << "\t" << y;
            dailyVegetationFile << "\t" << grid[x][y]->vegetationPatch->relVegCover;
            dailyVegetationFile << "\t" << grid[x][y]->vegetationPatch->abovegroundCMass;
            for(int i = 0; i < input.nSoilLayers; i++) {
                 dailyVegetationFile << "\t" << grid[x][y]->vegetationPatch->belowgroundCMass[i];
            }
            dailyVegetationFile  << "\t" << grid[x][y]->vegetationPatch->abovegroundNMass;
            for(int i = 0; i < input.nSoilLayers; i++) {
                dailyVegetationFile << "\t" << grid[x][y]->vegetationPatch->belowgroundNMass[i];
            }
            dailyVegetationFile << "\n";
        }
    }

    if(day == input.simDays) dailyVegetationFile.close();
}

/////////////////////////////////////////////////////////
void Output::initialiseDebugFile()
{
    string dailyDebugFileName;

    //create and open result file
    dailyDebugFileName = ".." + delimiter + "ModEST" + delimiter + "Output" + delimiter + "DailyDebug.txt";

    dailyDebugFile.open(dailyDebugFileName.c_str(), ios::out | ios::trunc);

    //if something goes wrong...
    if(!dailyDebugFile){
        cerr << dailyDebugFileName
             << " could not be opened!\n";
        exit(-1);
    }

    dailyDebugFile 	<< "Date" << "\t" << "CoordX" << "\t" << "CoordY" ;

    for(int i = 0; i < input.nSoilLayers; i++) {
        dailyDebugFile << "\t" <<

        "voidSoil_" + to_string(i+1) << "\t" <<
        "NcyclWF_" + to_string(i+1)<< "\t" <<
        "NcyclTF_" + to_string(i+1) << "\t" <<
        "CsomPerc_" + to_string(i+1) << "\t" <<
        "NO3_" + to_string(i+1) << "\t" <<
        "NH4_" + to_string(i+1) << "\t" <<
        "nitrified_" + to_string(i+1) << "\t" <<
        "volatilized_" + to_string(i+1) << "\t" <<
        "minNetRes_" + to_string(i+1) << "\t" <<
        "minSom_" + to_string(i+1) << "\t" <<
        "humifN_" + to_string(i+1) << "\t" <<
        "newCres_" + to_string(i+1) << "\t" <<
        "newNres_" + to_string(i+1) << "\t" <<
        "denit_" + to_string(i+1) << "\t" <<
        "leaching_in_" + to_string(i+1) << "\t" <<
        "leaching_out_" + to_string(i+1);

    }

    dailyDebugFile 	<< "\n";

}

/////////////////////////////////////////////////////////
void Output::writeDebugFile(int day, const NutrientPatch& nutrientPatch)
{
    for(int x = 0; x < input.xCells; x++) {
        for(int y = 0; y < input.yCells; y++) {
            dailyDebugFile 	<< input.date[day] << "\t" << x << "\t" << y ;

            for(int i = 0; i < input.nSoilLayers; i++) {
                dailyDebugFile << "\t" <<
                grid[x][y]->nutrientPatch->nutrientLayer[i].voidSoil  << "\t" <<
                grid[x][y]->nutrientPatch->nutrientLayer[i].NcyclWF << "\t" <<
                grid[x][y]->nutrientPatch->nutrientLayer[i].NcyclTF << "\t" <<
                grid[x][y]->nutrientPatch->nutrientLayer[i].CsomPerc << "\t" <<
                grid[x][y]->nutrientPatch->nutrientLayer[i].NO3 << "\t" <<
                grid[x][y]->nutrientPatch->nutrientLayer[i].NH4 << "\t" <<
                grid[x][y]->nutrientPatch->nutrientLayer[i].nitrified << "\t" <<
                grid[x][y]->nutrientPatch->nutrientLayer[i].volatilized << "\t" <<
                grid[x][y]->nutrientPatch->nutrientLayer[i].minNetRes << "\t" <<
                grid[x][y]->nutrientPatch->nutrientLayer[i].minSOM << "\t" <<
                grid[x][y]->nutrientPatch->nutrientLayer[i].humificResN << "\t" <<
                grid[x][y]->nutrientPatch->nutrientLayer[i].CresNew << "\t" <<
                grid[x][y]->nutrientPatch->nutrientLayer[i].NresNew << "\t" <<
                grid[x][y]->nutrientPatch->nutrientLayer[i].denit << "\t" <<
                grid[x][y]->nutrientPatch->nutrientLayer[i].leaching_in << "\t" <<
                grid[x][y]->nutrientPatch->nutrientLayer[i].leaching_out;
            }
            dailyDebugFile << "\n";
        }
    }
}

/////////////////////////////////////////////////////////
void Output::initialiseTestFile()
{
    string dailyTestFileName;

    //create and open result file
    if(exe)
        dailyTestFileName = "Output" + delimiter + "DailyTest.txt";
    else
        dailyTestFileName = ".." + delimiter + "ModEST" + delimiter + "Output" + delimiter + "DailyTest.txt";

    dailyTestFile.open(dailyTestFileName.c_str(), ios::out | ios::trunc);

    //if something goes wrong...
    if(!dailyTestFile){
        cerr << dailyTestFileName
             << " could not be opened!\n";
        exit(-1);
    }

    dailyTestFile 	<< "Date" << "\t" << "Plant" << "\t" << "leafCMass[kg]" << "\t" << "rootCMass[kg]" << "\t" << "sapwoodCMass[kg]" << "\t" << "heartwoodCMass[kg]" << "\t";
    dailyTestFile 	<< "height[m]" << "\t" << "lmToRmScal" << "\t" << "wscal" << "\t" << "nscal" << "\t" << "leafNDemand" << "\t" << "rootNDemand" << "\t";
    dailyTestFile 	<< "sapwoodNDemand" << "\t" << "storageNDemand" << "\t" << "availableNfunc" << "\t" << "soilTempfunc" << "\t" << "NC" << "\t" << "CNfunc" << "\t";
    dailyTestFile 	<< "avNitrogen" << "\t" << "totNdemand" << "\t" << "nitrogenAmount" << "\t" << "testAmount";
    dailyTestFile   << "\n";
}

/////////////////////////////////////////////////////////
void Output::writeTestFile(int day, const VegetationLandscape& plants)
{
    for(int i = 0; i < plants.plantList.size(); i++) {
        dailyTestFile << input.date[day] << "\t" << plants.plantList[i]->myID << "\t" << plants.plantList[i]->leafCMass << "\t" << plants.plantList[i]->rootCMass << "\t" << plants.plantList[i]->sapwoodCMass << "\t" << plants.plantList[i]->heartwoodCMass << "\t" ;
        dailyTestFile << plants.plantList[i]->height << "\t" << plants.plantList[i]->allocationTest[0] << "\t" << plants.plantList[i]->allocationTest[1] << "\t" << plants.plantList[i]->allocationTest[2] << "\t" << plants.plantList[i]->allocationTest[3] << "\t" << plants.plantList[i]->allocationTest[4] << "\t" ;
        dailyTestFile << plants.plantList[i]->allocationTest[5] << "\t" << plants.plantList[i]->allocationTest[6] << "\t" << plants.plantList[i]->allocationTest[7] << "\t" << plants.plantList[i]->allocationTest[8] << "\t" << plants.plantList[i]->allocationTest[9] << "\t" << plants.plantList[i]->allocationTest[10] << "\t" ;
        dailyTestFile << plants.plantList[i]->allocationTest[11] << "\t";
        dailyTestFile << plants.plantList[i]->allocationTest[12] << "\t";
        dailyTestFile << plants.plantList[i]->allocationTest[13] << "\t";
        dailyTestFile << plants.plantList[i]->allocationTest[14];
        dailyTestFile << "\n";
    }
}

/////////////////////////////////////////////////////////
void Output::initialiseTestFileB()
{
    string dailyTestFileBName;

    //create and open result file
    if(exe)
        dailyTestFileBName = "Output" + delimiter + "DailyTestB.txt";
    else
        dailyTestFileBName = ".." + delimiter + "ModEST" + delimiter + "Output" + delimiter + "DailyTestB.txt";

    dailyTestFileB.open(dailyTestFileBName.c_str(), ios::out | ios::trunc);

    //if something goes wrong...
    if(!dailyTestFileB){
        cerr << dailyTestFileBName
             << " could not be opened!\n";
        exit(-1);
    }

    dailyTestFileB 	<< "Date" << "\t" << "meauredRad[J*m-2*s-1]" << "\t" << "BrittaRad[mm*day-1]" << "\t" << "JosePotRad[MJ*m-2]" << "\t" << "JoseActRad[MJ*m-2]" << "\n";
}


/////////////////////////////////////////////////////////
void Output::writeTestFileB(int day)
{

    //BRITTA
    double distEarthSun = 1.0 + 0.033 * cos( 2.0 * M_PI / 365 * input.dayOfYear[day-1]);
    double solDeclination = 0.4093 * sin( 2.0 * M_PI / 365 * input.dayOfYear[day-1] - 1.405);
    double sunSetHour = acos( -1.0 * tan( input.latitude * M_PI/180.0 ) * tan( solDeclination ) );
    double extraRadiationBritta = 15.392 * distEarthSun * ( sunSetHour * sin( input.latitude * M_PI/180.0 ) * sin( solDeclination ) + cos(input.latitude * M_PI/180.0 ) * cos( solDeclination ) * sin( sunSetHour ) );

    //JOSE
    double relativeDistanceFromSun = 1.0 + 0.033 * cos( 2.0 * M_PI  * (day+1) / 365);
    double radLatitude = input.latitude * M_PI / 180;
    double solarDeclination = asin(0.4 * sin(((day+1) - 82) * 2 * M_PI /365));
    double avEr = 0.2618; // angular velocity of Earth rotation
    double sunRiseTime = acos(-tan(solarDeclination) * tan(radLatitude)) / avEr;
    double potSolRadiationJose = 40 * relativeDistanceFromSun * (avEr * sunRiseTime * sin(solarDeclination) * sin (radLatitude) + cos(solarDeclination) * cos(radLatitude) * sin(avEr * sunRiseTime));
    double TempEffect = input.tempMax[day] - input.tempMin[day];
    //if (day < 365) TempEffect = input.tempMax[day] - ( (input.tempMin[day] + input.tempMin[day+1]) / 2 );
    double a = 0.717;
    double b = 0.028;
    double c= 1.67;
    double actSolRadJose = potSolRadiationJose * a * (1-exp(-b * pow(TempEffect, c)));

    dailyTestFileB << input.date[day] << "\t";
    dailyTestFileB << input.solRadiation[day] << "\t";
    dailyTestFileB << extraRadiationBritta << "\t";
    dailyTestFileB << potSolRadiationJose << "\t";
    dailyTestFileB << actSolRadJose << "\n";
}

/////////////////////////////////////////////////////////
void Output::initialiseWaterOutputPaper2()
{
    string dailyWaterFileName;

    //create and open result file
    if(exe)
        dailyWaterFileName = "Output" + delimiter + "DailyWaterP2.txt";
    else
        dailyWaterFileName = ".." + delimiter + "ModEST" + delimiter + "Output" + delimiter + "DailyWaterP2.txt";

    cout << "Write daily water into: " << dailyWaterFileName << endl;
        cout << endl;

    dailyWaterFileP2.open(dailyWaterFileName.c_str(), ios::out | ios::trunc);

    //if something goes wrong...
    if(!dailyWaterFileP2){
        cerr << dailyWaterFileName
             << " could not be opened!\n";
        exit(-1);
    }

    dailyWaterFileP2 	<< "Date\t" << "xCor\t" << "yCor\t" << "DeepDrainage[mm]";

    dailyWaterFileP2 << "\n";
}

/////////////////////////////////////////////////////////
void Output::writeWaterOutputPaper2(int day, const WaterLandscape& water)
{
    for(int x = 0; x < input.xCells; x++) {
        for(int y = 0; y < input.yCells; y++) {
            dailyWaterFileP2 << input.date[day] << "\t" << x << "\t" << y << "\t" << grid[x][y]->waterPatch->waterLayer[1].drainedWater;
            dailyWaterFileP2 << "\n";
        }
    }

    if(day == input.simDays) dailyWaterFileP2.close();
}

/////////////////////////////////////////////////////////
void Output::initialisePlantOutputPaper2()
{
    string dailyPlantFileNameP2;

    //create and open result file
    if(exe)
        dailyPlantFileNameP2 = "Output" + delimiter + "DailyPlantP2.txt";
    else
       dailyPlantFileNameP2 = ".." + delimiter + "ModEST" + delimiter + "Output" + delimiter + "DailyPlantP2.txt";

    cout << "Write daily plants into: " << dailyPlantFileNameP2 << endl;
        cout << endl;

    dailyPlantFileP2.open(dailyPlantFileNameP2.c_str(), ios::out | ios::trunc);

    //if something goes wrong...
    if(!dailyPlantFileP2){
        cerr << dailyPlantFileNameP2
             << " could not be opened!\n";
        exit(-1);
    }

    dailyPlantFileP2 << "Date" << "\t" << "ID\t" << "PFT\t" << "xCor[m]\t" << "yCor[m]\t" << "NPP[kg/day]\t" << "NLitter[kg/day]\t" << "CLitter[kg/day]\t" << "totalCBiomass[kg]\t" << "totalBiomass[kg]\t" << "plantCover[m2]" << "\n";
}

/////////////////////////////////////////////////////////
void Output::writePlantOutputPaper2(int day, const VegetationLandscape& plants)
{
    for(int i = 0; i < plants.plantList.size(); i++) {
        int test = plants.plantList[i]->myID;
        dailyPlantFileP2 << input.date[day] << "\t" << plants.plantList[i]->myID << "\t" << plants.plantList[i]->PFTname << "\t" << plants.plantList[i]->xCor << "\t" << plants.plantList[i]->yCor;
        dailyPlantFileP2 << "\t" << plants.plantList[i]->NPP << "\t" << plants.plantList[i]->nitrogenLitter << "\t" << plants.plantList[i]->carbonLitter;
        dailyPlantFileP2 << "\t" << plants.plantList[i]->leafCMass + plants.plantList[i]->sapwoodCMass + plants.plantList[i]->rootCMass + plants.plantList[i]->heartwoodCMass + plants.plantList[i]->reproductiveMass;
        dailyPlantFileP2 << "\t" << plants.plantList[i]->leafCMass + plants.plantList[i]->sapwoodCMass + plants.plantList[i]->rootCMass + plants.plantList[i]->heartwoodCMass + plants.plantList[i]->leafNMass + plants.plantList[i]->sapwoodNMass + plants.plantList[i]->rootNMass;
        dailyPlantFileP2 << "\t" << plants.plantList[i]->crownArea << endl;
    }

    if(day == input.simDays) dailyPlantFileP2.close();
}

/////////////////////////////////////////////////////////
void Output::initialiseNutrientOutputPaper2()
{
    string dailyNutrientFileName;

    //create and open result file
    if(exe)
        dailyNutrientFileName = "Output" + delimiter + "DailyNutrientP2.txt";
    else
        dailyNutrientFileName = ".." + delimiter + "ModEST" + delimiter + "Output" + delimiter + "DailyNutrientP2.txt";

    cout << "Write daily nutrient into: " << dailyNutrientFileName << endl;
        cout << endl;

    dailyNutrientFileP2.open(dailyNutrientFileName.c_str(), ios::out | ios::trunc);

    //if something goes wrong...
    if(!dailyWaterFileP2){
        cerr << dailyNutrientFileName
             << " could not be opened!\n";
        exit(-1);
    }

    dailyNutrientFileP2 << "Date\t" << "xCor\t" << "yCor\t";
    dailyNutrientFileP2 << "soilCarbon[kg/ha]\t";
    dailyNutrientFileP2 << "soilAvailNPool[kg/ha]";
    dailyNutrientFileP2 << "\n";
}

/////////////////////////////////////////////////////////
void Output::writeNutrientOutputPaper2(int day, const NutrientPatch& nutrientPatch)
{
    for(int x = 0; x < input.xCells; x++) {
        for(int y = 0; y < input.yCells; y++) {
            dailyNutrientFileP2 << input.date[day] << "\t" << x << "\t" << y;
            dailyNutrientFileP2 << "\t" << grid[x][y]->nutrientPatch->nutrientLayer[0].Cres + grid[x][y]->nutrientPatch->nutrientLayer[1].Cres + grid[x][y]->nutrientPatch->nutrientLayer[0].Csom + grid[x][y]->nutrientPatch->nutrientLayer[1].Csom;
            dailyNutrientFileP2 << "\t" << grid[x][y]->nutrientPatch->nutrientLayer[0].NH4 + grid[x][y]->nutrientPatch->nutrientLayer[1].NH4 + grid[x][y]->nutrientPatch->nutrientLayer[0].NO3 + grid[x][y]->nutrientPatch->nutrientLayer[1].NO3;
            dailyNutrientFileP2 << "\n";
        }
    }

    if(day == input.simDays) dailyNutrientFileP2.close();
}

/////////////////////////////////////////////////////////
void Output::initialiseVegetationOutputPaper2()
{
    string dailyVegetationFileName;

    //create and open result file
    if(exe)
        dailyVegetationFileName = "Output" + delimiter + "DailyVegetationP2.txt";
    else
        dailyVegetationFileName = ".." + delimiter + "ModEST" + delimiter + "Output" + delimiter + "DailyVegetationP2.txt";

    cout << "Write daily vegetation into: " << dailyVegetationFileName << endl;
        cout << endl;

    dailyVegetationFileP2.open(dailyVegetationFileName.c_str(), ios::out | ios::trunc);

    //if something goes wrong...
    if(!dailyVegetationFileP2){
        cerr << dailyVegetationFileName
             << " could not be opened!\n";
        exit(-1);
    }

    dailyVegetationFileP2 << "Date\t" << "xCor\t" << "yCor\t" << "standingDeadCMass[kg]";
    dailyVegetationFileP2 << "\n";
}

/////////////////////////////////////////////////////////
void Output::writeVegetationOutputPaper2(int day, const VegetationLandscape& vegetation)
{

    for(int x = 0; x < input.xCells; x++) {
        for(int y = 0; y < input.yCells; y++) {
            dailyVegetationFileP2 << input.date[day] << "\t" << x << "\t" << y;
            dailyVegetationFileP2 << "\t" << grid[x][y]->vegetationPatch->standingDeadCMass << "\n";
        }
    }

    if(day == input.simDays) dailyVegetationFileP2.close();
}


/////////////////////////////////////////////////////////
void Output::writeFiles(int day, const VegetationLandscape& vegetation, const NutrientPatch& nutrientPatch, const WaterLandscape& water)
{
    if(paper) {
        writeWaterOutputPaper2(day, water);
        writeNutrientOutputPaper2(day, nutrientPatch);
        writePlantOutputPaper2(day, vegetation);
        //writeVegetationOutputPaper2(day, vegetation);
    }
    else {
        writeWaterFile(day, water);
        //writeNutrientFile(day, nutrientPatch);
        writePlantFile(day, vegetation);
        //writeVegetationFile(day, vegetation);
        //writeTestFile(day, vegetation);
        //writeTestFileB(day);
        //output->writeDebugFile(day, *nutrientPatch);
    }
}

/////////////////////////////////////////////////////////
Output::Output(const Input& _input, vector< vector<Patch*> >& _grid):
    input(_input), grid(_grid)
{
    if(paper) {
        initialiseWaterOutputPaper2();
        initialiseNutrientOutputPaper2();
        initialisePlantOutputPaper2();
        //initialiseVegetationOutputPaper2();
    }
    else {
        initialiseWaterFile();
        //initialiseNutrientFile();
        initialisePlantFile();
        //initialiseVegetationFile();
        //initialiseDebugFile();
        //initialiseTestFile();
        //initialiseTestFileB();
    }
}

/////////////////////////////////////////////////////////
Output::~Output()
{
    //dtor
}
