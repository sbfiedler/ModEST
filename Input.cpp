#include "Input.h"
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <random>

/////////////////////////////////////////////////////////
void Input::readMaster()
{
    string fileName;
    if(exe)
        fileName = "Input" + delimiter + "master.txt";
    else
        fileName = ".." + delimiter + "ModEST" + delimiter + "Input" + delimiter + "master.txt";

    ifstream masterFile;

    cout << "Read general model parameters from: " << fileName << endl;
    cout << endl;

    masterFile.open(fileName.c_str(), ios::binary|ios::in);

    if(!masterFile)
    {
        cout << "Master file could not be opened!";
    }

    int fsize = 10000;
	masterFile.ignore(fsize, ':');               //parses until end of file is reached (INT_MAX) OR until next ":" is found
	masterFile >> xCells;
	masterFile.ignore(fsize, ':');
	masterFile >> yCells;
	masterFile.ignore(fsize, ':');
	masterFile >> cellSize;
    masterFile.ignore(fsize, ':');
    masterFile >> nSoilLayers;

    depthLayers.resize(nSoilLayers);
    for(int i = 0; i < nSoilLayers; i++) {
        masterFile.ignore(fsize, ':');
        masterFile >> depthLayers[i];
    }

    masterFile.ignore(fsize, ':');
    masterFile >> latitude;

    masterFile.ignore(fsize, ':');
    masterFile >> simStartDate;
    masterFile.ignore(fsize, ':');
    masterFile >> simEndDate;
    masterFile.ignore(fsize, ':');
    masterFile >> measuredSolRadiation;

    masterFile.ignore(fsize, ':');
    masterFile >> iniSurfaceWater;

    iniWaterLayer.resize(nSoilLayers);
    for(int i = 0; i < nSoilLayers; i++) {
        masterFile.ignore(fsize, ':');
        masterFile >> iniWaterLayer[i];
    }

    masterFile.ignore(fsize, ':');
    masterFile >> iniNH4;
    masterFile.ignore(fsize, ':');
    masterFile >> iniNO3;
    masterFile.ignore(fsize, ':');
    masterFile >> dailyNdeposition;
    dailyNdeposition /= 1000.0;

    masterFile.ignore(fsize, ':');
    masterFile >> PFTs;
    masterFile.ignore(fsize, ':');
    masterFile >> spatialVegInput;
    if(!spatialVegInput) {
        masterFile.ignore(fsize, ':');
        masterFile >> nPFTs;
        pftNames.resize(nPFTs);
        nIndividuals.resize(nPFTs);
        meanHeight.resize(nPFTs);
        sdHeight.resize(nPFTs);
        totalIndividuals = 0;
        for(int pft = 0; pft < nPFTs; pft++) {
            //masterFile.ignore(fsize, '*');
            masterFile.ignore(fsize, ':');
            masterFile >> pftNames[pft];
            masterFile.ignore(fsize, ':');
            masterFile >> nIndividuals[pft];
            masterFile.ignore(fsize, ':');
            masterFile >> meanHeight[pft];
            masterFile.ignore(fsize, ':');
            masterFile >> sdHeight[pft];
            totalIndividuals += nIndividuals[pft];
        }
    }

    calculateSimDays();

    cout << "Latitude:\t\t" << latitude << " deg" << endl;
    cout << "Simulation start date:\t" << simStartDate << endl;
    cout << "Simulation end date:\t" << simEndDate << endl;
    cout << "Simulation time:\t" << simDays << " days" << endl;
    cout << "Grid scale:\t\t" << xCells << " x " << yCells << " grid cells" << endl;
    cout << "Cell size:\t\t" << cellSize << " m" << endl;
    cout << "Soil layers:\t\t" << nSoilLayers << " layers" << endl;
    cout << "Measured radiation?:\t" << measuredSolRadiation << endl;
    cout << endl;

    masterFile.close();
}

/////////////////////////////////////////////////////////
void Input::calculateSimDays()
{
    /////Calculate total simulation days "simDays"
    int simStartDay = stoi(simStartDate.substr(0,2));
    int simStartMonth = stoi(simStartDate.substr(3,2));
    int simStartYear= stoi(simStartDate.substr(6,4));

    int simEndDay = stoi(simEndDate.substr(0,2));
    int simEndMonth = stoi(simEndDate.substr(3,2));
    int simEndYear= stoi(simEndDate.substr(6,4));

    simDays = 0;
    for(int year = simStartYear; year <= simEndYear; year++) {
        int startMonth = 1;
        int endMonth = 12;
        if(year == simStartYear) startMonth = simStartMonth;
        if(year == simEndYear) endMonth = simEndMonth;
        for(int month = startMonth; month <= endMonth; month++) {

            int startDay = 1;
            int endDay;
            if(year == simStartYear && month == startMonth) startDay = simStartDay;

            if(month == 1) endDay = 31;
            if(month == 2) {
                endDay = 28;
                if(year % 4 == 0)
                    if(year % 100 != 0)
                        endDay = 29;
                if(year % 100 == 0 && year % 400 == 0)
                        endDay = 29;
            }
            if(month == 3) endDay = 31;
            if(month == 4) endDay = 30;
            if(month == 5) endDay = 31;
            if(month == 6) endDay = 30;
            if(month == 7) endDay = 31;
            if(month == 8) endDay = 31;
            if(month == 9) endDay = 30;
            if(month == 10) endDay = 31;
            if(month == 11) endDay = 30;
            if(month == 12) endDay = 31;
            if(year == simEndYear && month == endMonth) endDay = simEndDay;

            for(int day = startDay; day <= endDay; day++) {
                simDays++;
            }
        }
    }

    /////Calculate day of the year "dayOfYear" and year for each day
    dayOfYear.resize(simDays);
    yearOfDay.resize(simDays);
    date.resize(simDays);
    int count = 0;
    int thisYear = -1;
    int endDay;
    for(int year = simStartYear; year <= simEndYear; year++) {
        int startMonth = 1;
        int endMonth = 12;
        int doy = 0;
        thisYear++;
        if(year == simStartYear) {
            for(int month = startMonth; month <= simStartMonth; month++) {
                endDay = 0;
                if(month == 1) endDay = 31;
                if(month == 2) {
                    endDay = 28;
                    if(year % 4 == 0)
                        if(year % 100 != 0)
                            endDay = 29;
                    if(year % 100 == 0 && year % 400 == 0)
                            endDay = 29;
                }
                if(month == 3) endDay = 31;
                if(month == 4) endDay = 30;
                if(month == 5) endDay = 31;
                if(month == 6) endDay = 30;
                if(month == 7) endDay = 31;
                if(month == 8) endDay = 31;
                if(month == 9) endDay = 30;
                if(month == 10) endDay = 31;
                if(month == 11) endDay = 30;
                if(month == 12) endDay = 31;
                if(month == simStartMonth) endDay = simStartDay-1;

                doy += endDay;
            }
        }

        if(year == simStartYear) startMonth = simStartMonth;
        if(year == simEndYear) endMonth = simEndMonth;
        for(int month = startMonth; month <= endMonth; month++) {
            if(month == 1) endDay = 31;
            if(month == 2) {
                endDay = 28;
                if(year % 4 == 0)
                    if(year % 100 != 0)
                        endDay = 29;
                if(year % 100 == 0 && year % 400 == 0)
                        endDay = 29;
            }
            if(month == 3) endDay = 31;
            if(month == 4) endDay = 30;
            if(month == 5) endDay = 31;
            if(month == 6) endDay = 30;
            if(month == 7) endDay = 31;
            if(month == 8) endDay = 31;
            if(month == 9) endDay = 30;
            if(month == 10) endDay = 31;
            if(month == 11) endDay = 30;
            if(month == 12) endDay = 31;
            if(year == simEndYear && month == simEndMonth) endDay = simEndDay;

            int startDay = 1;
            if(year == simStartYear && month == startMonth) startDay = simStartDay;
            for(int day = startDay; day <= endDay; day++) {
                doy++;
                dayOfYear[count] = doy;
                yearOfDay[count] = thisYear;
                if(day < 10 && month < 10) date[count] = "0" + to_string(day) + ".0" + to_string(month) + "." + to_string(year);
                if(day < 10 && month >= 10) date[count] = "0" + to_string(day) + "." + to_string(month) + "." + to_string(year);
                if(day >= 10 && month < 10) date[count] = to_string(day) + ".0" + to_string(month) + "." + to_string(year);
                if(day >= 10 && month >= 10) date[count] = to_string(day) + "." + to_string(month) + "." + to_string(year);
                count++;
            }
        }
    }
}

/////////////////////////////////////////////////////////
void Input::readWeather()
{
    string fileName;
    if(exe)
        fileName = "Input" + delimiter + "weather.txt";
    else
        fileName = ".." + delimiter + "ModEST" + delimiter + "Input" + delimiter + "weather.txt";

    ifstream weatherFile;
    string dummy;

    cout << "Read weather from: " << fileName << endl;
        cout << endl;

    weatherFile.open(fileName.c_str(), ios::binary|ios::in);

    if(!weatherFile)
    {
        cout << "Weather file could not be opened!";
    }

    prec.resize(simDays);
    tempMean.resize(simDays);
    tempMin.resize(simDays);
    tempMax.resize(simDays);
    ambientCO2.resize(simDays);
    solRadiation.resize(simDays);

    if(measuredSolRadiation)
        for(int i = 0; i < 7; i++) {
             weatherFile >> dummy;
        }
    else
        for(int i = 0; i < 6; i++) {
             weatherFile >> dummy;
        }
    for(int day = 0; day < simDays; day++) {
        if(measuredSolRadiation)
            weatherFile >> dummy >> prec[day] >> tempMean[day] >> tempMin[day] >> tempMax[day] >> ambientCO2[day] >> solRadiation[day];
        else
           weatherFile >> dummy >> prec[day] >> tempMean[day] >> tempMin[day] >> tempMax[day] >> ambientCO2[day];
    }

    if(measuredSolRadiation) {
        cout << "Date\t\t" <<  "Prec[mm]\t" << "MeanTemp[degC]\t" << "MinTemp[degC]\t" << "MaxTemp[degC]\t" << "AmbientCO2[ppmv]" << "SolRad[Wm-2]"<< endl;
        cout << date[0] << "\t" <<  prec[0] << "\t\t" << tempMean[0] << "\t\t" << tempMin[0] << "\t\t" << tempMax[0] << "\t\t" << ambientCO2[0] << "\t\t" << solRadiation[0] << endl;
        cout << "...\t\t" <<  "...\t\t" <<  "...\t\t" << "...\t\t" << "...\t\t" << "...\t\t" << "..." << endl;
        cout << date[simDays-1] << "\t" <<  prec[simDays-1] <<  "\t\t" << tempMean[simDays-1] <<  "\t\t" << tempMin[simDays-1] <<  "\t\t" << tempMax[simDays-1] <<  "\t\t" << ambientCO2[simDays-1] << "\t\t" << solRadiation[simDays-1] << endl << endl;
    }
    else {
        cout << "Date\t\t" <<  "Prec[mm]\t" << "MeanTemp[degC]\t" << "MinTemp[degC]\t" << "MaxTemp[degC]\t" << endl;
        cout << date[0] << "\t" <<  prec[0] << "\t\t" << tempMean[0] << "\t\t" << tempMin[0] << "\t\t" << tempMax[0] << "\t\t" << ambientCO2[0] << endl;
        cout << "...\t\t" << "...\t\t" << "...\t\t" << "...\t\t" << "..." << endl;
        cout << date[simDays-1] << "\t" <<  prec[simDays-1] << "\t\t" << tempMean[simDays-1] << "\t\t" << tempMin[simDays-1] << "\t\t" << tempMax[simDays-1] << "\t\t" << ambientCO2[simDays-1] << endl << endl;
    }
    weatherFile.close();

    calculateMAT();
}

/////////////////////////////////////////////////////////
void Input::calculateMAT(){

    double sumMAT = 0.0;
    int begMissingDays = 0;
    int endMissingDays = 0;
    int begLastDay = 0;
    int endLastDay = 0;
    double begSumMat = 0.0;
    int year = -1;

    for(int day = 0; day < simDays; day++) {

        if(day == 0 && dayOfYear[day] > 1) {
            sumMAT += tempMean[day];
            MAT.push_back(0.0);
            begMissingDays = dayOfYear[day] - 1;
            year++;
        }

        else if(dayOfYear[day] == 1) {
            sumMAT = 0.0;
            sumMAT += tempMean[day];
            MAT.push_back(0.0);
            year++;
        }

        else {
            sumMAT += tempMean[day];
        }

        if(day == simDays-1 && dayOfYear[day] < 365) {
            endMissingDays = 365 - dayOfYear[day]; //does not account for leap years here!
            endLastDay = dayOfYear[day];
        }

        else if(dayOfYear[day] >= 365 && year == 0 && begMissingDays > 0) {
            begSumMat = sumMAT;
            begLastDay = day;
        }

        else if(dayOfYear[day] >= 365) MAT[year] = sumMAT/dayOfYear[day];
    }

    if(begMissingDays > 0) {
        int test = begMissingDays;
        for(int day = begLastDay+1; day < begLastDay+1+begMissingDays; day++) {
            begSumMat += tempMean[day];
            test--;
        }
        MAT[0] = begSumMat/(begLastDay+1+begMissingDays);
    }

    if(endMissingDays > 0) {
        for(int day = (simDays - endLastDay - endMissingDays); day < (simDays - endLastDay); day++) {
            sumMAT += tempMean[day];
        }
        MAT[year] = sumMAT/365.0; //does not account for leap years here!
    }
}

/////////////////////////////////////////////////////////
void Input::readTopography()
{
    string fileName;
    if(exe)
        fileName = "Input" + delimiter + "elevationMap.txt";
    else
        fileName = ".." + delimiter + "ModEST" + delimiter + "Input" + delimiter + "elevationMap.txt";

    string param;
    ifstream topoFile;

    cout << "Read topography from: " << fileName << endl;
        cout << endl;

    topoFile.open(fileName.c_str(), ios::binary|ios::in);

    if(!topoFile)
    {
        cout << "Topography file could not be opened!";
    }

    elevation.resize(xCells, vector<double>(yCells));

    // (0,0) in the bottom left corner, and (xCells-1,yCells-1) in the top right corner
    for(int j = (yCells -1); j >= 0; j--) {
        for(int i = 0; i < xCells; i++) {
            topoFile >> elevation[i][j];
		}
	}

    topoFile.close();
}

/////////////////////////////////////////////////////////
void Input::readSoilMap()
{
    string fileName;
    if(exe)
        fileName = "Input" + delimiter + "soilpropertyMap.txt";
    else
        fileName = ".." + delimiter + "ModEST" + delimiter + "Input" + delimiter + "soilpropertyMap.txt";

    string param;
    ifstream soiMapFile;

    cout << "Read soil map from: " << fileName << endl;
        cout << endl;

    soiMapFile.open(fileName.c_str(), ios::binary|ios::in);

    if(!soiMapFile)
    {
        cout << "Soilmap file could not be opened!";
    }

    soilType.resize(xCells, vector<string>(yCells));

    // (0,0) in the bottom left corner, and (xCells-1,yCells-1) in the top right corner
    for(int j = (yCells -1); j >= 0; j--) {
        for(int i = 0; i < xCells; i++) {
            soiMapFile >> soilType[i][j];
        }
    }

    soiMapFile.close();
}

/////////////////////////////////////////////////////////
void Input::readSoilProperties()
{
    string fileName;
    if(exe)
        fileName = "Input" + delimiter + "soilProperties.txt";
    else
        fileName = ".." + delimiter + "ModEST" + delimiter + "Input" + delimiter + "soilProperties.txt";

    int soilPropertiesFileLines;
    soilPropertiesFileLines = 0;
    ifstream soilPropertiesFile;
    string line;

    cout << "Read soil properties from: " << fileName << endl;
        cout << endl;

    soilPropertiesFile.open(fileName.c_str(), ios::binary|ios::in);

    if(!soilPropertiesFile)
    {
        cout << "Soilparameter file could not be opened!";
    }

    while(getline(soilPropertiesFile, line))
    {
        soilPropertiesFileLines++;
    }
    soilPropertiesFile.close();

    soilProperties.resize((soilPropertiesFileLines - 1), vector<string>(14));

    soilPropertiesFile.open(fileName.c_str(), ios::binary|ios::in);
    int fsize = 10000;
    soilPropertiesFile.ignore(fsize, '\n'); // Skipping the first line with headers

    for(int soiltype = 0; soiltype < (soilPropertiesFileLines-1); soiltype++) {
        for(int properties = 0; properties < 14; properties++) {
            soilPropertiesFile >> soilProperties[soiltype][properties];
        }
    }

    soilPropertiesFile.close();
}

/////////////////////////////////////////////////////////
void Input::readIniVegetation()
{
    if(spatialVegInput) {
        string fileName;
        if(exe)
            fileName = "Input" + delimiter + "vegetationMap.txt";
        else
            fileName = ".." + delimiter + "ModEST" + delimiter + "Input" + delimiter + "vegetationMap.txt";

        int fileLines;
        int dummy;
        fileLines = 0;
        ifstream iniVegetationFile;
        string line;

        cout << "Read initial vegetation map from: " << fileName << endl;
            cout << endl;

        iniVegetationFile.open(fileName.c_str(), ios::binary|ios::in);

        if(!iniVegetationFile)
        {
            cout << "Initial vegetation map could not be opened!";
        }

        while(getline(iniVegetationFile, line))
        {
            fileLines++;
        }
        iniVegetationFile.close();

        iniVegetation.resize((fileLines - 1), vector<string>(5));

        iniVegetationFile.open(fileName.c_str(), ios::binary|ios::in);
        int fsize = 10000;
        iniVegetationFile.ignore(fsize, '\n'); // Skipping the first line with headers

        for(int plantID = 0; plantID < (fileLines-1); plantID++) {
            iniVegetationFile >> iniVegetation[plantID][0]; // ID of the individual
            iniVegetationFile >> iniVegetation[plantID][1]; // PFT
            iniVegetationFile >> dummy; // RidgefieldRipline
            iniVegetationFile >> dummy; // RidgefieldPosition
            iniVegetationFile >> iniVegetation[plantID][2]; //xCor
            iniVegetationFile >> iniVegetation[plantID][3]; //yCor
            iniVegetationFile >> iniVegetation[plantID][4]; //height
        }

        iniVegetationFile.close();
    }
    else {
        iniVegetation.resize(totalIndividuals, vector<string>(5));
        vector<double> xCors;
        vector<double> yCors;
        int ID = 1;
        //set (random) seed and random generator depending on current time
        myclock::duration d = myclock::now() - begin;
        unsigned randomSeed = d.count();
        default_random_engine generator(randomSeed);
        for(int pft = 0; pft < nPFTs; pft++) {
            for(int ind = 0; ind < nIndividuals[pft]; ind++) {
                iniVegetation[ID-1][0] = to_string(ID); // ID of the individual
                iniVegetation[ID-1][1] = pftNames[pft]; // PFT

                double xCor, yCor;
                bool foundMyPlace = false;
                while(!foundMyPlace) {
                    uniform_real_distribution<double> xCorDistribution(0.0, xCells*cellSize);
                    xCor = xCorDistribution(generator);
                    uniform_real_distribution<double> yCorDistribution(0.0, yCells*cellSize);
                    yCor = yCorDistribution(generator);
                    if(xCors.size() == 0) foundMyPlace = true;
                    else {
                        foundMyPlace = true;
                        for(int c = 0; c < xCors.size(); c++) {
                            double distance = 1.5;
                            if(xCor > xCors[c] - distance && xCor < xCors[c] + distance && yCor > yCors[c] - distance && yCor < yCors[c] + distance) {
                                foundMyPlace = false;
                                break;
                            }
                        }
                    }
                }
                iniVegetation[ID-1][2] = to_string(xCor); //xCor
                xCors.push_back(xCor);
                iniVegetation[ID-1][3] = to_string(yCor); //yCor
                yCors.push_back(yCor);


                normal_distribution<double> heightDistribution(meanHeight[pft], sdHeight[pft]);
                double height = -1.0;
                while(height <= 0.0) {
                    height = heightDistribution(generator);
                }
                iniVegetation[ID-1][4] = to_string(height); //height
                ID++;
            }
        }
    }
}

/////////////////////////////////////////////////////////
void Input::readPlantTraits()
{
    string fileName;
    if(exe)
        fileName = "Input" + delimiter + "plantTraits.txt";
    else
        fileName = ".." + delimiter + "ModEST" + delimiter + "Input" + delimiter + "plantTraits.txt";
    int fileLines;
    string dummy;
    fileLines = 0;
    ifstream plantTraitFile;
    string line;

    cout << "Read plant traits from: " << fileName << endl;
        cout << endl;

    plantTraitFile.open(fileName.c_str(), ios::binary|ios::in);

    if(!plantTraitFile)
    {
        cout << "Plant trait file could not be opened!";
    }

    while(getline(plantTraitFile, line))
    {
        fileLines++;
    }
    plantTraitFile.close();

    plantTraits.resize(PFTs, vector<string>(fileLines));

    plantTraitFile.open(fileName.c_str(), ios::binary|ios::in);

    for(int properties = 0; properties < fileLines; properties++) {
        plantTraitFile >> dummy;
        for(int PFT = 0; PFT < PFTs; PFT++) {
            plantTraitFile >> plantTraits[PFT][properties];
        }
    }

    plantTraitFile.close();
}

/////////////////////////////////////////////////////////
Input::Input()
{
    begin = myclock::now();
    readMaster();
    readTopography();
    readWeather();
    readSoilProperties();
    readSoilMap();
    readIniVegetation();
    readPlantTraits();
}

/////////////////////////////////////////////////////////
Input::~Input()
{
    //dtor
}
