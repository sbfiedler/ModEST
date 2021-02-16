#include "Plant.h"
#include "Input.h"
#include "Landscape.h"
#include "Patch.h"
#include <cmath>
#include <iostream>
#include <random>

/////////////////////////////////////////////////////////
Plant::Plant(int index, const Input& _input, Landscape* _land) : myID(atof(_input.iniVegetation[index][0].c_str())), PFTname(_input.iniVegetation[index][1]), xCor(atof(_input.iniVegetation[index][2].c_str())),
  yCor(atof(_input.iniVegetation[index][3].c_str())), height(atof(_input.iniVegetation[index][4].c_str())), input(_input), land(_land)
{
    relRootsPerLayer.resize(input.nSoilLayers, 0);
    allocationTest.resize(15,0);

    for(int i = 0; i < input.plantTraits.size(); i++) {
        if(input.plantTraits[i][0].c_str() == PFTname) {
            phenologyType = input.plantTraits[i][1];
            GDD5 = atof(input.plantTraits[i][2].c_str());
            SLA = atof(input.plantTraits[i][3].c_str());
            laToSa = atof(input.plantTraits[i][4].c_str());
            lmToRm = atof(input.plantTraits[i][5].c_str());
            allom2 = atof(input.plantTraits[i][6].c_str());
            allom3 = atof(input.plantTraits[i][7].c_str());
            allom1 = atof(input.plantTraits[i][8].c_str());
            reinickeRp = atof(input.plantTraits[i][9].c_str());
            maxCrownArea = atof(input.plantTraits[i][10].c_str());
            photosyntheticPathway = input.plantTraits[i][11];
            CO2UptakeEfficiency = atof(input.plantTraits[i][12].c_str());
            temp1 = atof(input.plantTraits[i][13].c_str());
            temp2 = atof(input.plantTraits[i][14].c_str());
            temp3 = atof(input.plantTraits[i][15].c_str());
            temp4 = atof(input.plantTraits[i][16].c_str());
            optLambda = atof(input.plantTraits[i][17].c_str());
            gMin = atof(input.plantTraits[i][18].c_str());
            WP = atof(input.plantTraits[i][19].c_str());
            eMax = atof(input.plantTraits[i][20].c_str());
            b = atof(input.plantTraits[i][21].c_str());
            respirationCoeff = atof(input.plantTraits[i][22].c_str());
            fracGResp = atof(input.plantTraits[i][23].c_str());
            carbToNitroLeaves = atof(input.plantTraits[i][24].c_str());
            carbToNitroSapwood = atof(input.plantTraits[i][25].c_str());
            carbToNitroRoots = atof(input.plantTraits[i][26].c_str());
            minLeafCN = atof(input.plantTraits[i][27].c_str());
            maxLeafCN = atof(input.plantTraits[i][28].c_str());
            k = atof(input.plantTraits[i][29].c_str());
            kmax = atof(input.plantTraits[i][30].c_str())/1000.0;
            nUptakePerRootC = atof(input.plantTraits[i][31].c_str())/1000.0;
            leafTurnoverTime = atof(input.plantTraits[i][32].c_str());
            rootTurnoverTime = atof(input.plantTraits[i][33].c_str());
            sapwoodTurnoverTime = atof(input.plantTraits[i][34].c_str());
            seedMass = atof(input.plantTraits[i][35].c_str());
            meanDispersalDistance = atof(input.plantTraits[i][36].c_str())*100.0;
            sdDispersalDistance = atof(input.plantTraits[i][37].c_str())*100.0;
            germinationProb = atof(input.plantTraits[i][38].c_str());
            LAI = atof(input.plantTraits[i][39].c_str());
            woodDensity = atof(input.plantTraits[i][40].c_str());
            mortThreshold = atof(input.plantTraits[i][41].c_str());
            for(int j = 0; j < input.nSoilLayers; j++) relRootsPerLayer[j] = atof(input.plantTraits[i][42 + j].c_str());
            break;
       }
    }

    xPatch = min(input.xCells-1, int(xCor / input.cellSize));
    yPatch = min(input.yCells-1, int(yCor / input.cellSize));
    currentGDD5 = 0;
    abovegroundCLitter = 0.0;
    belowgroundCLitter = 0.0;
    abovegroundNLitter = 0.0;
    belowgroundNLitter = 0.0;
    waterDemand = 0.0;
    waterSupply = 0.0;
    dead = false;
    NPP = 0.0;
    nitrogenTurnover = 0.0;
    nScal = 1.0;
    wScal = 1.0;
    reproductiveMass = 0.0;
    shadedAPAR = 0.0;
    age = 0;
    transpiredWater = 0.0;
    uptakenNitrogen = 0.0;
    nitrogenLitter = 0.0;
    carbonLitter = 0.0;
    totalRespiration = 0.0;
    dailyCInc = 0.0;
    heartwoodCMass = 0.0;
    phenologyStatus = 1.0;

    beginning = myclock::now();

    calculateIniPlantProperties();
}

/////////////////////////////////////////////////////////
Plant::Plant(int _myID, string _PFTname, double _xCor, double _yCor, const Input& _input, Landscape* _land) : myID(_myID), PFTname(_PFTname), xCor(_xCor),
  yCor(_yCor), input(_input), land(_land)
{
    relRootsPerLayer.resize(input.nSoilLayers, 0);
    allocationTest.resize(15,0);

    for(int i = 0; i < input.plantTraits.size(); i++) {
        if(input.plantTraits[i][0].c_str() == PFTname) {
            phenologyType = input.plantTraits[i][1];
            GDD5 = atof(input.plantTraits[i][2].c_str());
            SLA = atof(input.plantTraits[i][3].c_str());
            laToSa = atof(input.plantTraits[i][4].c_str());
            lmToRm = atof(input.plantTraits[i][5].c_str());
            allom2 = atof(input.plantTraits[i][6].c_str());
            allom3 = atof(input.plantTraits[i][7].c_str());
            allom1 = atof(input.plantTraits[i][8].c_str());
            reinickeRp = atof(input.plantTraits[i][9].c_str());
            maxCrownArea = atof(input.plantTraits[i][10].c_str());
            photosyntheticPathway = input.plantTraits[i][11];
            CO2UptakeEfficiency = atof(input.plantTraits[i][12].c_str());
            temp1 = atof(input.plantTraits[i][13].c_str());
            temp2 = atof(input.plantTraits[i][14].c_str());
            temp3 = atof(input.plantTraits[i][15].c_str());
            temp4 = atof(input.plantTraits[i][16].c_str());
            optLambda = atof(input.plantTraits[i][17].c_str());
            gMin = atof(input.plantTraits[i][18].c_str());
            WP = atof(input.plantTraits[i][19].c_str());
            eMax = atof(input.plantTraits[i][20].c_str());
            b = atof(input.plantTraits[i][21].c_str());
            respirationCoeff = atof(input.plantTraits[i][22].c_str());
            fracGResp = atof(input.plantTraits[i][23].c_str());
            carbToNitroLeaves = atof(input.plantTraits[i][24].c_str());
            carbToNitroSapwood = atof(input.plantTraits[i][25].c_str());
            carbToNitroRoots = atof(input.plantTraits[i][26].c_str());
            minLeafCN = atof(input.plantTraits[i][27].c_str());
            maxLeafCN = atof(input.plantTraits[i][28].c_str());
            k = atof(input.plantTraits[i][29].c_str());
            kmax = atof(input.plantTraits[i][30].c_str());
            nUptakePerRootC = atof(input.plantTraits[i][31].c_str());
            leafTurnoverTime = atof(input.plantTraits[i][32].c_str());
            rootTurnoverTime = atof(input.plantTraits[i][33].c_str());
            sapwoodTurnoverTime = atof(input.plantTraits[i][34].c_str());
            seedMass = atof(input.plantTraits[i][35].c_str());
            meanDispersalDistance = atof(input.plantTraits[i][36].c_str())*100.0;
            sdDispersalDistance = atof(input.plantTraits[i][37].c_str())*100.0;
            germinationProb = atof(input.plantTraits[i][38].c_str());
            LAI = 1.5; /*!< \todo in plantTraits the values are for planted individuals in Ridgefield but should be for saplings if model runs with saplings from the beginning. 1.5 from LPJ  */
            woodDensity = atof(input.plantTraits[i][40].c_str());
            mortThreshold = atof(input.plantTraits[i][41].c_str());
            for(int j = 0; j < input.nSoilLayers; j++) relRootsPerLayer[j] = atof(input.plantTraits[i][42 + j].c_str());
            break;
       }
    }

    height = 0.1; /*!< \todo What initial height value? */
    xPatch = int(xCor / input.cellSize);
    yPatch = int(yCor / input.cellSize);
    currentGDD5 = 0;
    abovegroundCLitter = 0.0;
    belowgroundCLitter = 0.0;
    abovegroundNLitter = 0.0;
    belowgroundNLitter = 0.0;
    waterDemand = 0.0;
    waterSupply = 0.0;
    dead = false;
    NPP = 0.0;
    nitrogenTurnover = 0.0;
    nScal = 1.0;
    wScal = 1.0;
    reproductiveMass = 0.0;
    shadedAPAR = 0.0;
    age = 0;
    transpiredWater = 0.0;
    uptakenNitrogen = 0.0;
    nitrogenLitter = 0.0;
    carbonLitter = 0.0;
    totalRespiration = 0.0;
    dailyCInc = 0.0;
    heartwoodCMass = 0.0;
    phenologyStatus = 1.0;

    beginning = myclock::now();

    calculateIniPlantProperties();
}

/////////////////////////////////////////////////////////
void Plant::calculateIniPlantProperties()
{
    //Height-diameter relationship (Sitch et al. 2000, eqn. 3)
    stemDiameter = pow(height / allom2, 1.0 / allom3); /*!< \todo Use Jonas' equations */

    //Self-thinning constraint (Sitch et al. 2000, eqn. 4-8)
    crownArea = min(allom1 * pow(stemDiameter, reinickeRp), maxCrownArea);
    //SLA = 0.0002 * exp(6.15 - 0.46*log(leafTurnoverTime * 12)); // SLA is an input variable now
    leafCMass = (LAI * crownArea) / SLA; /*!< \todo Calculate from Maries allometric equation */
    FPC = (1 - exp(-0.5 * LAI)) * crownArea; // changed for an IBM approach

    //Crown radius (NEW)
    crownRadius = pow(crownArea / M_PI, 0.5);

    //Functional balance (Sitch et al. 2000, eqn. 2)
    rootCMass = leafCMass / lmToRm;

    //Initial total wood mass (equals sapwood mass so far) calculated via volume function for a cylinder and wood density
    double sapwoodVolume = ((pow(stemDiameter, 2) * M_PI) / 4) * height; // in m3
    sapwoodCMass = sapwoodVolume * woodDensity;
    sapwoodArea = sapwoodVolume / height;
    //The pipe model (Sitch et al. 2000, eqn. 1)
    leafArea = laToSa * sapwoodArea;

    //Inital values for nitrogen mass for each compartment (given that there is an optimal C:N)
    leafNMass = leafCMass / carbToNitroLeaves;
    sapwoodNMass = sapwoodCMass / carbToNitroSapwood;
    rootNMass = rootCMass / carbToNitroRoots;
    storageNMass = maxNitrogenStorage();

    debtMass = 0.0;
    leafCMassPreviousYear = leafCMass;
    rootCMassPreviousYear = rootCMass;
    sapwoodCMassPreviousYear = sapwoodCMass;
    totalCMassPreviousYear = leafCMass + rootCMass + sapwoodCMass;
    totalCMass = leafCMass + rootCMass + sapwoodCMass + heartwoodCMass;
    ANPP = 0.0;
}

/////////////////////////////////////////////////////////
void Plant::updatePlantProperties()
{
    double woodVolume = (sapwoodCMass+heartwoodCMass) / woodDensity; // in m³
    double sapwoodVolume = sapwoodCMass / woodDensity; // in m³

    //Height-diameter relationship (Sitch et al. 2000, eqn. 3) and volume function for a cylinder
    stemDiameter = pow((4*woodVolume) / (M_PI * allom2), 1/(allom3 + 2));
    height = allom2 * pow(stemDiameter, allom3);

    if(height > 0) sapwoodArea = sapwoodVolume / height;
    else sapwoodArea = 0.0;
    //The pipe model (Sitch et al. 2000, eqn. 1)
    leafArea = laToSa * sapwoodArea;

    //Self-thinning constraint (Sitch et al. 2000, eqn. 4-8)
    crownArea = min(allom1 * pow(stemDiameter, reinickeRp), maxCrownArea);
    if(crownArea > 0) LAI = (leafCMass * SLA) / (crownArea);
    else LAI = 0.0;
    FPC = (1 - exp(-0.5 * LAI)) * crownArea; // changed for an IBM approach

    //NEW
    crownRadius = pow(crownArea / M_PI, 0.5);

    dailyCInc = leafCMass + rootCMass + sapwoodCMass + heartwoodCMass - totalCMass;
    totalCMass = leafCMass + rootCMass + sapwoodCMass + heartwoodCMass;

    age++;

    if(sapwoodCMass <= 0.0) {
        inCaseOfDeath();
        dead = true;
    }
}

/////////////////////////////////////////////////////////
void Plant::tissueTurnover()
{
    //Sitch et al. 2000, eqn. 41 (adapted for individuals and daily turnover)
    abovegroundCLitter = abovegroundCLitter + leafCMass / (365*leafTurnoverTime);
    belowgroundCLitter = belowgroundCLitter + rootCMass / (365*rootTurnoverTime);

    //Smith et al. 2014, between eqn. C13 and C14
    nitrogenTurnover = (leafCMass / (leafTurnoverTime*365)) * 1.0/carbToNitroLeaves + (rootCMass / (rootTurnoverTime*365)) * 1.0/carbToNitroRoots + (sapwoodCMass / (sapwoodTurnoverTime*365)) * 1.0/carbToNitroSapwood;

    if(storageNMass <= maxNitrogenStorage()) {
        if((storageNMass + nitrogenTurnover) <= maxNitrogenStorage()) {
            storageNMass += nitrogenTurnover;
            nitrogenTurnover = 0.0;
        }
        else {
            double diff = maxNitrogenStorage() - storageNMass;
            storageNMass += diff;
            nitrogenTurnover -= diff;
            belowgroundNLitter += 0.5 * nitrogenTurnover;
            abovegroundNLitter += 0.5 * nitrogenTurnover;
            nitrogenTurnover = 0.0;
        }
    }
    else {
        //nitrogenTurnover will be allocated to the litter pools (50:50)
        belowgroundNLitter += 0.5 * nitrogenTurnover;
        abovegroundNLitter += 0.5 * nitrogenTurnover;
        nitrogenTurnover = 0.0;
    }

    //Sitch et al. 2000, eqn. 42 (adapted for individuals and daily turnover)
    leafCMass = leafCMass * (1.0 - (1.0 / (leafTurnoverTime*365)));
    heartwoodCMass = heartwoodCMass + (sapwoodCMass / (sapwoodTurnoverTime*365));
    sapwoodCMass = sapwoodCMass * (1.0 - 1.0 / (sapwoodTurnoverTime*365));
    rootCMass = rootCMass * (1.0 - 1.0 / (rootTurnoverTime*365));

    updatePlantProperties();
}

/////////////////////////////////////////////////////////
void Plant::photosynthesis(int day, double dayLength, double lambda)
{
    //changed after Schapoff et al., eqn 1, 26 (now: included shading by substracting the APAR that will be used by overgrown individuals = shadedAPAR)
    APAR = max(0.0, (0.5 * (land->extraterrestrialRadiation(day) * 2450000.0) * FPC * phenologyStatus) - shadedAPAR); //2450000.0 to transfer solar radiation mm*day-1 to J*m-2*day-1
    APAR *= 0.0000046; //convert from J*day-1 to mol*day-1 for solar radiation at 550 nm
    APAR *= cMass; //convert from mol*day-1 to gC*day-1

    //Temperature stress function after Sibyll Schapoff (Email from 15th March 2018) (cite: LPJml4) instead of Sitch et al. 2000, eqn. 17
    double k2 = (temp1+temp2)*0.5;
    double k1 = 2*log(1.0/0.99-1)/(temp1-temp2);
    double k3 = log(0.99/0.01)/(temp4-temp3);
    double low = 1.0/(1+exp(k1*(k2-input.tempMean[day])));
    double high = 1-0.01*exp(k3*(input.tempMean[day]-temp3));
    double tempStressFunction = low * high; /*!< PFT-specific temperature inhibition function limiting photosynthesis at low and high temperature */

    double c1, c2;

    if(photosyntheticPathway == "C3") {

        //Schapoff et al. 2000, eqn. 31a
        double tau = 2600.0 * pow(0.57, (input.tempMean[day] - 25.0) * 0.1);

        //Sitch et al. 2000, eqn. 19
        double CO2CompPoint  = ambientO2Pa / (2 * tau);    /*!< CO2 compensation point */

        //From LPJGuess code: Eqn 7, Haxeltine & Prentice 1996a
        double interCO2Pa = lambda * input.ambientCO2[day] * atmPa * CO2ConvFactor; //Intercellular partial pressure of CO2 given stomatal opening (Pa)

        //Schapoff et al. 2017, eqn. 29
        //c1 = CO2UptakeEfficiency * tempStressFunction * ((ambientPartCO2Pressure * lambda - CO2CompPoint) / (ambientPartCO2Pressure * lambda - CO2CompPoint));
        //Sitch et al. 2000, eqn. 15 (eqn adapted to equation in LPJGuess code)
        c1 = CO2UptakeEfficiency * tempStressFunction * ((interCO2Pa - CO2CompPoint) / (interCO2Pa + 2 * CO2CompPoint));

        //Sitch et al. 2000, eqn. 21
        c2 = ( interCO2Pa - CO2CompPoint ) / ( interCO2Pa + 30.0 * (1.0 + ambientO2Pa / 30000.0) );
    }

    if(photosyntheticPathway == "C4") {
        //Sitch et al. 2000, eqn. 16
        c1 = CO2UptakeEfficiency * tempStressFunction * lambda / 0.4;

        //Sitch et al. 2000, eqn. 22
        c2 = 1.0;
    }

    vmax(c1, c2, dayLength, day);

    //Sitch et al. 2000, eqn. 23
    double leafRespiration = b * Vmax;   /*!< daily leaf respiration in [gC * m-2 * day-1] */

    //Sitch et al. 2000, eqn. 14
    //double lightLimitedPhotosynthsisRate = c1 * 0.0000046 * 12.0107 * 1/dayLength * APAR;
    //Schapoff et al. 2017, eqn. 28
    double lightLimitedPhotosynthsisRate = c1 * 1/dayLength * APAR; // gC/hour

    //Sitch et al. 2000, eqn. 20
    double rubiscoLimitedPhotosynthsisRate = (double)1/24 * c2 * Vmax; // gC/hour

    //Sitch et al. 2000, eqn. 24
    GPP = (((lightLimitedPhotosynthsisRate + rubiscoLimitedPhotosynthsisRate - pow(pow(lightLimitedPhotosynthsisRate + rubiscoLimitedPhotosynthsisRate, 2) - 4.0 * theta * lightLimitedPhotosynthsisRate * rubiscoLimitedPhotosynthsisRate , 0.5)) / 2 * theta) * dayLength); /*!< daily gross photosynthesis in [gC * m-2 * day-1] */
    /*!< \todo Here can be negative. Is this valid? */
    if(GPP < 0) GPP = 0.0;

    //Sitch et al. 2000, eqn. 26
    NPP = GPP - dayLength / 24 * leafRespiration; /*!< net daytime photosynthesis in [gC/m2/day] */
    actLeafRespiration = (dayLength / 24 * leafRespiration) / 1000.0; // kg

    //Convert netDaytimePhotosynthesis from gC/m2/day to mm/m2/day using ideal gas equation
    NPPinMM = (NPP / cMass * 8.314 * (input.tempMean[day]+273.15)/atmPa*1000.0); //from LPJml code

    //Convert photosynthesis into kgC * m-2 * day-1
    NPP /= 1000.0;
    GPP /= 1000.0;

    //Sitch et al. 2000, eqn. 27
    if(NPPinMM <= 0)
        gPot = 0.0;
    else
        gPot = 1.6 / CO2ConvFactor / 3600 * NPPinMM / input.ambientCO2[day] / (1 - lambda) / dayLength; // in mm/s
}

/////////////////////////////////////////////////////////
void Plant::vmax(double c1, double c2, double dayLength, int day) {

    //Sitch et al. 2000, eqn. 25
    double s = (24.0 / dayLength) * b;
    double sigma = sqrt(max(0.0, 1.0 - ((c2 - s) / (c2 - theta * s))));
    Vmax = 1/b * c1/c2 * ((2.0 * theta - 1.0) * s - (2 * theta * s - c2) * sigma) * APAR * nScal; //nScal differently calculated in Smith et al. 2014
    if(Vmax < 0) Vmax = 0.0; //included since there were negative values that led to negative leaf respiration and thus to increasing NPP in photosynthesis function

//    //Smith et al. 2014, eqn. C10 and C11
//    /*!< \todo Too high values. Will always be corrected by if statements */
//    optLeafNMass = 2083.0 * Vmax * ( exp(-0.693*(input.tempMean[day]-25)) / (dayLength*3600) ) * exp(0.12 * LAI) + 7.15 * 0.001 * leafCMass; //0.001 from LPJGuess code
//    // Can not have higher nitrogen concentration than minimum leaf C:N ratio
//    if (leafCMass / optLeafNMass < minLeafCN) {
//        optLeafNMass = leafCMass / minLeafCN;
//    }
//    // Can not have lower nitrogen concentration than maximum leaf C:N ratio
//    else if (leafCMass / optLeafNMass > maxLeafCN) {
//        optLeafNMass = leafCMass / maxLeafCN;
//    }

//    //The following from LPJGuess code
//    // Calculate nitrogen-limited Vmax for current leaf nitrogen
//    // Haxeltine & Prentice 1996b Eqn 28
//    const double M = 25.0; // corresponds to parameter p in Eqn 28, Haxeltine & Prentice 1996b

//    // Conversion factor in calculation of leaf nitrogen: includes conversion of:
//    //		- Vm from gC/m2/day to umolC/m2/sec
//    //      - nitrogen from mg/m2 to kg/m2
//    double CN = 1.0 / (3600 * dayLength * CMASS);
//    double tfac = exp(-0.0693 * (input.tempMean[day] - 25.0));
//    double maxVmax = nactive / (M * CN * tfac);

//    // Calculate optimal leaf nitrogen based on [potential] Vmax (Eqn 28 Haxeltine & Prentice 1996b)
//    nactive_opt = M * vm * CN * tfac;

//    if (Vmax > maxVmax && ifnlimvmax) {
//        vmaxnlim = maxVmax / Vmax;	// Save vmax nitrogen limitation
//        Vmax = maxVmax;
//    }
//    else {
//        vmaxnlim = 1.0;
//    }
}

/////////////////////////////////////////////////////////
void Plant::calculatewScal(int day, double phen, double potLeafCMass, double potFPC)
{
    gAct = gPot + gMin;

    //// WATER DEMAND
    //Tietjen et al. 2009, eqn. 10
    /*!< \todo  aspect factor and inclination missing here */
    potEP = 0.0023 * (input.tempMean[day] + 17.8) * pow((input.tempMax[day] - input.tempMin[day]), 0.5) * land->extraterrestrialRadiation(day); //[mm/day]

    //waterDemand = potEP * FPC * alpha * (1 - exp(-gAct * phenologyStatus/scalingFactor)); //Sitch et al. 2000, eqn. 34
    //changed since potEP is area-based and there should be the fractional plant cover in the equation
    /*!< \todo potEP in mm/day, gAct and gMin in mm/s. Is this valid? */
    waterDemand = potEP * potFPC * alpha / (1 + scalingFactor/gAct); //Schaphoff et al. 2017, eqn. 114; looks different to the equation in source code (water_stressed.c, line 135)

    //// WATER SUPPLY
    //Schapoff et al. 2017, eqn. 109
    //changed since eMax is area-based and there should be the fractional plant cover in the equation (also added in the LPJml code). Also included, rootCMass/leafCMass factor have a feedback of higher root mass ratio on water yield
    waterSupply = 0.0;
    if(potLeafCMass > 0) waterSupply = eMax * potFPC * avRelSoilMoistureFC * phen * (rootCMass*lmToRm)/potLeafCMass; //[mm/day]

    //For rare cases if waterSupply is higher than actual absolute available water which would lead to negative soil moisture then
    if(waterSupply > avAbsSoilMoisture) waterSupply = avAbsSoilMoisture;

    //According to the description in Smith et al. 2001 for use in allocation and phenology routine (However not sure if this calculated as intended)
    wScal = min(1.0, waterSupply/waterDemand);
}

/////////////////////////////////////////////////////////
void Plant::transpiration(int day, double dayLength)
{
    transpiredWater = 0.0;

    //Calculate gAct, potEP, waterSupply, waterDemand and wScal
    calculatewScal(day, phenologyStatus, leafCMass, FPC);

    //Sitch et al. 2000, eqn. 32
    double ET = min(waterSupply, waterDemand);

    //New routine for extracting the water from the soil
    if(ET > 0.0 && avAbsSoilMoisture > 0.0) {
        for(int i = 0; i < intersectedPatches.size(); i++) {
            for(int j = 0; j < input.nSoilLayers; j++) {
                transpiredWater += ((intersectedPatches[i].second * relRootsPerLayer[j] * land->grid[intersectedPatches[i].first.first][intersectedPatches[i].first.second]->waterPatch->waterLayer[j].absWater) / avAbsSoilMoisture) * ET;
                land->grid[intersectedPatches[i].first.first][intersectedPatches[i].first.second]->waterPatch->waterLayer[j].transpiredAbsWater += ((intersectedPatches[i].second * relRootsPerLayer[j] * land->grid[intersectedPatches[i].first.first][intersectedPatches[i].first.second]->waterPatch->waterLayer[j].absWater) / avAbsSoilMoisture) * ET;
                land->grid[intersectedPatches[i].first.first][intersectedPatches[i].first.second]->waterPatch->waterLayer[j].absWater -= ((intersectedPatches[i].second * relRootsPerLayer[j] * land->grid[intersectedPatches[i].first.first][intersectedPatches[i].first.second]->waterPatch->waterLayer[j].absWater) / avAbsSoilMoisture) * ET;
                land->grid[intersectedPatches[i].first.first][intersectedPatches[i].first.second]->waterPatch->waterLayer[j].relWater = land->grid[intersectedPatches[i].first.first][intersectedPatches[i].first.second]->waterPatch->waterLayer[j].absWater / land->grid[intersectedPatches[i].first.first][intersectedPatches[i].first.second]->waterPatch->waterLayer[j].depthLayer;
            }
        }
    }

    //// WATER-STRESSED PHOTOSYNTHESIS
    //The following lines from LPJml code, water_stressed.c
    if(waterDemand > waterSupply) {
        gAct = 0;
        if(potEP > 0) {
            //gAct = scalingFactor * log(1-waterSupply / (potEP * alpha));    //Sitch et al. 2000, eqn. 35
            gAct = scalingFactor / (((potEP * FPC * alpha) / waterSupply) - 1.0); //Schaphoff et al. 2017, eqn. 114
            if(gAct < 0) gAct = 0;
        }
        double lambda = 0.0;
        double lamdaLow = 0.02;
        double lambdaHigh = optLambda + 0.05;
        double lambdaAcc = 0.0;
        double ymid = 0.0;
        double yacc = 0.1;
        double fac = (dayLength*3600 * (gAct - gMin * phenologyStatus)) / 1.6 * input.ambientCO2[day] * 1e-6;
        photosynthesis(day, dayLength, lamdaLow);
        double ylow = fac * (1 - lamdaLow) - NPPinMM;

        for(int i = 0; i < 10; i++) { /*!< \todo how many maximum iterations? 10 from LPJml4 code */
            lambda = (lamdaLow + lambdaHigh) * 0.5;
            if(lambdaHigh - lamdaLow < lambdaAcc) break;
            photosynthesis(day, dayLength, lambda);
            ymid = fac * (1 - lambda) - NPPinMM;
            if(abs(ymid) < yacc) break;
            if(ylow * ymid <= 0.0) lambdaHigh = lambda;
            else {
                lamdaLow = lambda;
                ylow = ymid;
            }
        }
        photosynthesis(day, dayLength, lambda);
    }
}

/////////////////////////////////////////////////////////
void Plant::phenology(int day)
{
    // Phenology status is expressed as fraction between 0 and 1, representing the fraction of full leaf coverage attained by the PFT (after Sitch et al. 2000)
    /*! \todo should be tested for summergreen and raingreen */

    phenologyStatus = 0.0;

    if(phenologyType == "evergreen") {
        phenologyStatus = 1.0; // status = 1.0 means full leaf coverage the entire year
    }

    if(phenologyType == "summergreen") {

        if(input.latitude >= 0) {
            if(input.dayOfYear[day] == 0) currentGDD5 = 0; //reset GDD 1st of January (because of northern hemisphere), differs in leap years
        }
        else
            if(input.dayOfYear[day] == 182) currentGDD5 = 0; //reset GDD around 1st of July (because of southern hemisphere)

        if(input.tempMean[day] > 5.0) currentGDD5++;
        phenologyStatus = double(currentGDD5)/double(GDD5);
        if(phenologyStatus > 1.0) phenologyStatus = 1.0;

        //Senescence
        if((input.latitude >= 0 && input.dayOfYear[day] >= 182) || (input.latitude < 0 && input.dayOfYear[day] >= 0)) {
            if(input.tempMean[day] < 5.0) { /*! \todo TEST THIS IF IT SENESCENCE IS MODELLED CORRECTLY. TEMP THRESHOLD SHOULD BE PFT-SPECIFIC */
                phenologyStatus = 0.0;
                double Ntransfer = 0.5 * leafCMass/carbToNitroLeaves; //Smith et al. 2014: The N store is replenished up to its maximum capacity with up to 50% of the N mass of shed leaves
                Ntransfer = min(leafNMass, Ntransfer);
                if((storageNMass + Ntransfer) > maxNitrogenStorage()) {
                    Ntransfer = maxNitrogenStorage() - storageNMass;
                }
                storageNMass += Ntransfer;
                abovegroundNLitter += leafNMass - Ntransfer;
                abovegroundCLitter += leafCMass;
                leafCMass = 0.0;
                leafNMass = 0.0;
            }
        }
    }

    if(phenologyType == "raingreen") {
        //First calculate wScal as if plant would be fully covered in leaves
        potLeafCMass = leafCMass;
        if(potLeafCMass == 0.0)
            potLeafCMass = rootCMass * lmToRm;
        potentialLAI = 0.0;
        if(crownArea > 0) potentialLAI = (potLeafCMass * SLA) / (crownArea);
        potentialFPC = (1 - exp(-0.5 * potentialLAI)) * crownArea; // changed for an IBM approach
        calculatewScal(day, 1.0, potLeafCMass, potentialFPC);
        wScal_phen = wScal;

        if(wScal > 0.35) { /*! \todo 0.35 should probably be a PFT-specific parameter */
            phenologyStatus = 1.0;
            /*! \todo How to better resprout leaves? */
            leafCMass += 0.0001; //add very little leaf C mass to beginn with
        } else { //Senescence
            phenologyStatus = 0.0;
            if(leafCMass > 0.0) {
                double Ntransfer = 0.5 * leafCMass/carbToNitroLeaves; //Smith et al. 2014: The N store is replenished up to its maximum capacity with up to 50% of the N mass of shed leaves
                Ntransfer = min(leafNMass, Ntransfer);
                if((storageNMass + Ntransfer) > maxNitrogenStorage()) {
                    Ntransfer = maxNitrogenStorage() - storageNMass;
                }
                storageNMass += Ntransfer;
                abovegroundNLitter += leafNMass - Ntransfer;
                abovegroundCLitter += leafCMass;
                leafCMass = 0.0;
                leafNMass = 0.0;
            }
        }
    }
}

/////////////////////////////////////////////////////////
void Plant::allocation()
{
    double leafMassInc = 0.0;
    double rootMassInc = 0.0;
    double sapwoodMassInc = 0.0;
    double NSEG = 20.0; //number of segments (parameter in numerical methods)
    ANPP += NPP;

    //Calculation the most limitating stress factor (either nitrogen or water stress) after Smith et al. 2014, eqn. C20
    if(NPP >= 0.0)
        lmToRmScal = lmToRm * min(wScal, nScal);
    else
        lmToRmScal = lmToRm;

    allocationTest[0] = lmToRmScal;
    allocationTest[1] = wScal;
    allocationTest[2] = nScal;

    //Calculation of leaf mass increment (leafMassInc) satisfying Sitch et al. 2000, eqn. 213 using bisection method (Press et al. 1986)
    double x1, x2;
    if(NPP > 0.0)
        x2 = NPP / (1.0 + 1.0/lmToRmScal); //new
    if(NPP < 0.0) {
        x2 = -NPP / (1.0 + 1.0/lmToRmScal); //new
        if(x2 > leafCMass) x2 = 0.0;
    }
    if(NPP == 0.0)
        x2 = 0.0;
    x1 = 0.0;

    //Bisection loop (after LPJml code)
    if(x1 == 0.0 && x2 == 0.0)
        leafMassInc = 0.0;
    else {
        if(x2 < x1) {
            double swap = x1;
            x1 = x2;
            x2 = swap;
        }

        double dx = (x2 - x1) / NSEG;
        double xmid;
        double ymid;
        if(allocationFunction(x1) < 0.0)
            for(xmid = x1 + dx; allocationFunction(xmid) < 0.0 && xmid <= x2 - dx; xmid += dx);
        else
            for(xmid = x1 + dx; allocationFunction(xmid) > 0.0 && xmid <= x2 - dx; xmid += dx);

        double xlow = xmid - dx;
        double xhigh = xmid;
        double xacc = 0.001; /*!< \todo what accuracy? 0.00001 from LPJml4 code */
        double yacc = 1.0e-10; /*!< \todo what accuracy? 0.0000000001 from LPJml4 code */
        double ylow = allocationFunction(xlow);

        for(int i = 0; i < 40; i++) { /*!< \todo how many maximum iterations? 40 from LPJml4 code */
            xmid = (xlow + xhigh) * 0.5;
            if(xhigh - xlow < xacc) break;
            ymid = allocationFunction(xmid);
            if(abs(ymid) < yacc) break;
            if(ylow * ymid <= 0.0) xhigh = xmid;
            else {
                xlow = xmid;
                ylow = ymid;
            }
        }
        leafMassInc = xmid;
    }
    if(NPP > 0.0) {
        if(leafMassInc != 0.0) rootMassInc = leafMassInc / lmToRmScal; //new
        else rootMassInc = 0.0;
    }
    if(NPP < 0.0) {
        leafMassInc *= -1.0;
        if(leafMassInc != 0.0) rootMassInc = leafMassInc / lmToRmScal; //new
        else rootMassInc = 0.0;
    }
    //Sitch et al. 2000, eqn. 1
    sapwoodMassInc = NPP - leafMassInc - rootMassInc;

    if(leafMassInc < 0.0) {
        if(leafMassInc > leafCMass) { //does this make sense? shouldnt it be (-leafMassInc > leafCMass)
            leafMassInc = leafCMass;
        }
        carbonLitter += -leafMassInc;
        abovegroundCLitter += -leafMassInc;
        //Smith et al. 2014, between eqn. C13 and C14
        nitrogenTurnover += -leafMassInc * 1.0/carbToNitroLeaves;
    }
    leafCMass += leafMassInc;

    if(rootMassInc < 0.0) {
        if(rootMassInc > rootCMass) {
            rootMassInc = rootCMass;
        }
        carbonLitter += -rootMassInc;
        belowgroundCLitter += -rootMassInc;
        //Smith et al. 2014, between eqn. C13 and C14
        nitrogenTurnover += -rootMassInc * 1.0/carbToNitroRoots;
    }
    rootCMass += rootMassInc;

    if(sapwoodMassInc < 0.0) {
        if(sapwoodMassInc > sapwoodCMass) {
            sapwoodMassInc = sapwoodCMass;
        }
        heartwoodCMass += -sapwoodMassInc;
        nitrogenTurnover += -sapwoodMassInc * 1.0/carbToNitroSapwood;
    }
    sapwoodCMass += sapwoodMassInc;

    if(storageNMass <= maxNitrogenStorage()) {
        if((storageNMass + nitrogenTurnover) <= maxNitrogenStorage()) {
            storageNMass += nitrogenTurnover;
            nitrogenTurnover = 0.0;
        }
        else {
            double diff = maxNitrogenStorage() - storageNMass;
            storageNMass += diff;
            nitrogenTurnover -= diff;
            nitrogenLitter += nitrogenTurnover;
            belowgroundNLitter += 0.5 * nitrogenTurnover;
            abovegroundNLitter += 0.5 * nitrogenTurnover;
            nitrogenTurnover = 0.0;
        }
    }
    else {
        //nitrogenTurnover will be allocated to the litter pools (50:50)
        nitrogenLitter += nitrogenTurnover;
        belowgroundNLitter += 0.5 * nitrogenTurnover;
        abovegroundNLitter += 0.5 * nitrogenTurnover;
        nitrogenTurnover = 0.0;
    }
    updatePlantProperties();
}

/////////////////////////////////////////////////////////
double Plant::allocationFunction(double leafMassInc)
{
    //Sitch et al. 2000, eqn. 43-45
    return pow(allom2, 2.0/allom3) * 4.0/M_PI * (sapwoodCMass + NPP - leafMassInc - (leafCMass + leafMassInc) / lmToRmScal + rootCMass + heartwoodCMass) / woodDensity -
            pow((sapwoodCMass + NPP - leafMassInc - (leafCMass + leafMassInc) / lmToRmScal + rootCMass) / (leafCMass + leafMassInc) / woodDensity / SLA * laToSa, 1 + 2/allom3);
}

/////////////////////////////////////////////////////////
void Plant::respiration(int day)
{
    totalRespiration = 0.0;

    //Sitch et al. 2000, eqn. 37 & 38
    const double k = 0.095218;   /*!< \todo constant only found in lpj code > canexch.cpp > respiration */
    double arrheniusFunction; /*!< Arrhenius temperature-respiration relationship */
    arrheniusFunction = exp( 308.56 * ( (1/56.02) - ( 1/(input.tempMean[day]+46.02) ) ) ); //Sitch et al. 2000, eqn. 38 with mistake airTemp + 46.02 and not airTemp - 46.02 (from LPJml4 description)

    //Sapwood autotrophic respiration
    double sapwoodRespiration = respirationCoeff * k * (sapwoodCMass / carbToNitroSapwood) * arrheniusFunction;

    //Root autotrophic respiration
    arrheniusFunction = exp( 308.56 * ( (1/56.02) - ( 1/(avSoilTemperature+46.02) ) ) );
    double rootRespiration = respirationCoeff * k * (rootCMass / carbToNitroRoots) * phenologyStatus * arrheniusFunction;

    //Total autotrophic respiration (in LPJ this is calculated yearly and for the whole grid), Sitch et al. 2000, eqn. 39
    double autotrophicRespiration = sapwoodRespiration + rootRespiration;

    //Update NPP (by extracting 25 percent for growh respiration), Sitch et al. 2000, eqn. 40
    NPPyesterday = NPP;
    if((NPP - autotrophicRespiration) > 0) {
        NPP = (1.0 - fracGResp) * (NPP - autotrophicRespiration); //new: only if there is enough NPP left
        totalRespiration += autotrophicRespiration + fracGResp * (NPP - autotrophicRespiration);
    }
    else {
        NPP -= autotrophicRespiration; //new
        totalRespiration += autotrophicRespiration;
    }

    totalRespiration += actLeafRespiration;
}

/////////////////////////////////////////////////////////
void Plant::mortality()
{
    double deltaCMass = (leafCMass + sapwoodCMass + rootCMass) / totalCMassPreviousYear;

    if(deltaCMass <= mortThreshold) {
        inCaseOfDeath();
        dead = true;
    }

    leafCMassPreviousYear = leafCMass;
    rootCMassPreviousYear = rootCMass;
    sapwoodCMassPreviousYear = sapwoodCMass;
    totalCMassPreviousYear = leafCMass + sapwoodCMass + rootCMass;
}

/////////////////////////////////////////////////////////
void Plant::reproduction()
{
    if(NPP > 0.0) {
        //Sitch et al. 2000
        reproductiveMass += 0.1 * NPP;
        NPP = 0.9 * NPP;
    }
}

/////////////////////////////////////////////////////////
void Plant::nitrogenUptake(int day, double dayLength)
{
    ///Calculation of total N demand (self invented instead of Smith et al. 2014, eqn. C10 (because of strange values)
    /*!< \todo leaf N demand should be dependent on Rubisco as in Smith et al. 2014 but not working */
    uptakenNitrogen = 0.0;
    double leafNDemand = max(0.0, leafCMass / carbToNitroLeaves - leafNMass);
    double rootNDemand = max(0.0, rootCMass / carbToNitroRoots  - rootNMass);
    double sapwoodNDemand = max(0.0, sapwoodCMass / carbToNitroSapwood - sapwoodNMass);

    //Smith et al. 2014, changed eqn. C13 for daily approach (originally yearly time step??)
    double storageNDemand = max(0.0, maxNitrogenStorage() - storageNMass);
    double totNitrogenDemand = leafNDemand + rootNDemand + sapwoodNDemand + storageNDemand; // in kg

    ///Calculation of maximum N uptake
    //Smith et al. 2014, eqn. C15
    double availableNfunc = 0.05 + (avNitrogen / (avNitrogen + kmax * waterSaturationPerDepth() * crownArea)); /*!< \todo not fully clear if waterSaturationPerDepth is correctly calculated? */

    //Smith et al. 2014, eqn. C5
    double soilTempfunc = 0.0326 + 0.00351 * pow(avSoilTemperature, 1.652) - pow((avSoilTemperature / 41.748), 7.19);

    //Smith et al. 2014, eqn. C16
    double NC = (leafNMass + rootNMass) / (leafCMass + rootCMass);
    //von Bloh et al. 2018, eqn. 10
	double CNfunc = max( 0.0, ( NC - (1.0/minLeafCN) ) / ( 2.0/(maxLeafCN + minLeafCN) - 1.0/minLeafCN ) );

    //Smith et al. 2014, eqn. C14
    double maxNUptake = 2.0 * nUptakePerRootC * availableNfunc * soilTempfunc * CNfunc * rootCMass;

    ///Fraction of actual nitrogen demand limited by maximum N uptake
    double fracNDemand = 1.0;
    if(totNitrogenDemand > 0.0) fracNDemand = min(maxNUptake/totNitrogenDemand, 1.0);

    ///Update of N demand for each pool and total demand
    leafNDemand *= fracNDemand;
    rootNDemand *= fracNDemand;
    sapwoodNDemand *= fracNDemand;
    storageNDemand *= fracNDemand;
    totNitrogenDemand = leafNDemand + rootNDemand + sapwoodNDemand + storageNDemand;

    allocationTest[3] = leafNDemand;
    allocationTest[4] = rootNDemand;
    allocationTest[5] = sapwoodNDemand;
    allocationTest[6] = storageNDemand;
    allocationTest[7] = availableNfunc;
    allocationTest[8] = soilTempfunc;
    allocationTest[9] = NC;
    allocationTest[10] = CNfunc;
    allocationTest[11] = avNitrogen;
    allocationTest[12] = 0;

    ///Actual N uptake
    if(totNitrogenDemand > 0.0 && totNitrogenDemand <= avNitrogen) {
        uptakenNitrogen = totNitrogenDemand;
        extractNitrogenFromSoil(totNitrogenDemand);
        leafNMass += leafNDemand;
        rootNMass += rootNDemand;
        sapwoodNMass += sapwoodNDemand;
        storageNMass += storageNDemand;
        nScal = 1.0;
    }

    else {
        ///Calculate fractional demand for each pool given the total N demand without storage demand included
        totNitrogenDemand = leafNDemand + rootNDemand + sapwoodNDemand;
        double fracLeafNDemand = leafNDemand/totNitrogenDemand;
        double fracRootNDemand = rootNDemand/totNitrogenDemand;
        double fracSapwoodNDemand = sapwoodNDemand/totNitrogenDemand;

        if(totNitrogenDemand > 0.0 && totNitrogenDemand <= avNitrogen) {
            uptakenNitrogen = totNitrogenDemand;
            extractNitrogenFromSoil(totNitrogenDemand);
            leafNMass += leafNDemand;
            rootNMass += rootNDemand;
            sapwoodNMass += sapwoodNDemand;
            nScal = 1.0;
        }

        else if(totNitrogenDemand > 0.0 && totNitrogenDemand > avNitrogen) {
            uptakenNitrogen = avNitrogen;
            extractNitrogenFromSoil(avNitrogen);
            leafNMass += fracLeafNDemand * avNitrogen;
            rootNMass += fracRootNDemand * avNitrogen;
            sapwoodNMass += fracSapwoodNDemand * avNitrogen;

            ///Use N storage
            if((totNitrogenDemand - avNitrogen) <= storageNMass) {
                leafNMass += leafNDemand - (fracLeafNDemand * avNitrogen);
                rootNMass += rootNDemand - (fracRootNDemand * avNitrogen);
                sapwoodNMass += sapwoodNDemand - (fracSapwoodNDemand * avNitrogen);
                storageNMass -= leafNDemand - (fracLeafNDemand * avNitrogen) + rootNDemand - (fracRootNDemand * avNitrogen) + sapwoodNDemand - (fracSapwoodNDemand * avNitrogen);
                nScal = 1.0;
            }

            else {
                leafNMass += fracLeafNDemand * storageNMass;
                rootNMass += fracRootNDemand * storageNMass;
                sapwoodNMass += fracSapwoodNDemand * storageNMass;
                nScal = (avNitrogen + storageNMass) / totNitrogenDemand;
                storageNMass = 0.0;
                //Recalculate photosynthesis because of nitrogen limitation
                photosynthesis(day, dayLength, optLambda);
            }
        }
    }
    allocationTest[12] = totNitrogenDemand;
}

/////////////////////////////////////////////////////////
double Plant::maxNitrogenStorage()
{
    double maxNStore = 0.0;
    //Smith et al. 2014, eqn. C12
    maxNStore = k * sapwoodCMass * leafNMass / leafCMass;

    return maxNStore;
}

/////////////////////////////////////////////////////////
void Plant::extractNitrogenFromSoil(double nitrogenAmount)
{
    allocationTest[13] = nitrogenAmount;
    double test = 0;
    if(nitrogenAmount > 0.0 && avNitrogen > 0.0) {
        for(int i = 0; i < intersectedPatches.size(); i++) {
            for(int j = 0; j < input.nSoilLayers; j++) {
                double oldNO3_weight = land->grid[intersectedPatches[i].first.first][intersectedPatches[i].first.second]->nutrientPatch->nutrientLayer[j].NO3_weight;
                double oldNH4_weight = land->grid[intersectedPatches[i].first.first][intersectedPatches[i].first.second]->nutrientPatch->nutrientLayer[j].NH4_weight;

                double fractionNO3 = 0.0;
                double fractionNH4 = 0.0;
                if(oldNO3_weight > 0.0 || oldNH4_weight > 0.0) {
                    fractionNO3 = oldNO3_weight /(oldNO3_weight+oldNH4_weight);
                    fractionNH4 = oldNH4_weight /(oldNO3_weight+oldNH4_weight);
                }

                if(oldNO3_weight > 0.0 && oldNH4_weight > 0.0) {
                    land->grid[intersectedPatches[i].first.first][intersectedPatches[i].first.second]->nutrientPatch->nutrientLayer[j].NO3_weight -= intersectedPatches[i].second * relRootsPerLayer[j] * (land->grid[intersectedPatches[i].first.first][intersectedPatches[i].first.second]->nutrientPatch->nutrientLayer[j].NO3_weight / avNitrogen) * (fractionNO3 * nitrogenAmount);
                    test += intersectedPatches[i].second * relRootsPerLayer[j] * (land->grid[intersectedPatches[i].first.first][intersectedPatches[i].first.second]->nutrientPatch->nutrientLayer[j].NO3_weight / avNitrogen) * (fractionNO3 * nitrogenAmount);
                    land->grid[intersectedPatches[i].first.first][intersectedPatches[i].first.second]->nutrientPatch->nutrientLayer[j].NO3 = land->grid[intersectedPatches[i].first.first][intersectedPatches[i].first.second]->nutrientPatch->nutrientLayer[j].NO3 * (land->grid[intersectedPatches[i].first.first][intersectedPatches[i].first.second]->nutrientPatch->nutrientLayer[j].NO3_weight / oldNO3_weight);
                    land->grid[intersectedPatches[i].first.first][intersectedPatches[i].first.second]->nutrientPatch->nutrientLayer[j].NH4_weight -= intersectedPatches[i].second * relRootsPerLayer[j] * (land->grid[intersectedPatches[i].first.first][intersectedPatches[i].first.second]->nutrientPatch->nutrientLayer[j].NH4_weight / avNitrogen) * (fractionNH4 * nitrogenAmount);
                    land->grid[intersectedPatches[i].first.first][intersectedPatches[i].first.second]->nutrientPatch->nutrientLayer[j].NH4 = land->grid[intersectedPatches[i].first.first][intersectedPatches[i].first.second]->nutrientPatch->nutrientLayer[j].NH4 * (land->grid[intersectedPatches[i].first.first][intersectedPatches[i].first.second]->nutrientPatch->nutrientLayer[j].NH4_weight / oldNH4_weight);
                    test += intersectedPatches[i].second * relRootsPerLayer[j] * (land->grid[intersectedPatches[i].first.first][intersectedPatches[i].first.second]->nutrientPatch->nutrientLayer[j].NH4_weight / avNitrogen) * (fractionNH4 * nitrogenAmount);
                }
                else if(oldNO3_weight <= 0.0 && oldNH4_weight > 0.0) {
                    land->grid[intersectedPatches[i].first.first][intersectedPatches[i].first.second]->nutrientPatch->nutrientLayer[j].NH4_weight -= intersectedPatches[i].second * relRootsPerLayer[j] * (land->grid[intersectedPatches[i].first.first][intersectedPatches[i].first.second]->nutrientPatch->nutrientLayer[j].NH4_weight / avNitrogen) * (fractionNH4 * nitrogenAmount);
                    land->grid[intersectedPatches[i].first.first][intersectedPatches[i].first.second]->nutrientPatch->nutrientLayer[j].NH4 = land->grid[intersectedPatches[i].first.first][intersectedPatches[i].first.second]->nutrientPatch->nutrientLayer[j].NH4 * (land->grid[intersectedPatches[i].first.first][intersectedPatches[i].first.second]->nutrientPatch->nutrientLayer[j].NH4_weight / oldNH4_weight);
                    test += intersectedPatches[i].second * relRootsPerLayer[j] * (land->grid[intersectedPatches[i].first.first][intersectedPatches[i].first.second]->nutrientPatch->nutrientLayer[j].NH4_weight / avNitrogen) * (fractionNH4 * nitrogenAmount);
                }
                else if(oldNO3_weight > 0.0 && oldNH4_weight <= 0.0) {
                    land->grid[intersectedPatches[i].first.first][intersectedPatches[i].first.second]->nutrientPatch->nutrientLayer[j].NO3_weight -= intersectedPatches[i].second * relRootsPerLayer[j] * (land->grid[intersectedPatches[i].first.first][intersectedPatches[i].first.second]->nutrientPatch->nutrientLayer[j].NO3_weight / avNitrogen) * (fractionNO3 * nitrogenAmount);
                    land->grid[intersectedPatches[i].first.first][intersectedPatches[i].first.second]->nutrientPatch->nutrientLayer[j].NO3 = land->grid[intersectedPatches[i].first.first][intersectedPatches[i].first.second]->nutrientPatch->nutrientLayer[j].NO3 * (land->grid[intersectedPatches[i].first.first][intersectedPatches[i].first.second]->nutrientPatch->nutrientLayer[j].NO3_weight / oldNO3_weight);
                    test += intersectedPatches[i].second * relRootsPerLayer[j] * (land->grid[intersectedPatches[i].first.first][intersectedPatches[i].first.second]->nutrientPatch->nutrientLayer[j].NO3_weight / avNitrogen) * (fractionNO3 * nitrogenAmount);
                }
            }
        }
    }
    allocationTest[14] = test;
}

/////////////////////////////////////////////////////////
double Plant::waterSaturationPerDepth()
{
    double waterSat = 0.0;
    double soilDepth = 0.0;

    //assuming that soil depth is the same for each patch (i = 0)
    for(int j = 0; j < input.nSoilLayers; j++) {
        if(relRootsPerLayer[j] > 0.0)
            soilDepth +=  land->grid[intersectedPatches[0].first.first][intersectedPatches[0].first.second]->waterPatch->waterLayer[j].depthLayer/1000.0;
    }

    for(int i = 0; i < intersectedPatches.size(); i++) {
        for(int j = 0; j < input.nSoilLayers; j++) {
            waterSat += intersectedPatches[i].second * relRootsPerLayer[j] * land->grid[intersectedPatches[i].first.first][intersectedPatches[i].first.second]->waterPatch->soilProp.fieldCap;
        }
    }
    return (waterSat * soilDepth);
}

/////////////////////////////////////////////////////////
void Plant::availablePools()
{
    avRelSoilMoisture = 0.0;
    avRelSoilMoistureFC = 0.0;
    avAbsSoilMoistureNew = 0.0;
    avAbsSoilMoisture = 0.0;
    avNitrogen = 0.0;
    avSoilTemperature = 0.0;
    for(int i = 0; i < intersectedPatches.size(); i++) {
        for(int j = 0; j < input.nSoilLayers; j++) {

            // first calculate wilting point for the soil texture in which the plant is rooting here (after van Genuchten, Maidment, p.5.6 / 5.14)
            double alpha = pow(land->grid[intersectedPatches[i].first.first][intersectedPatches[i].first.second]->waterPatch->soilProp.bubPressure, -1);
            double n = land->grid[intersectedPatches[i].first.first][intersectedPatches[i].first.second]->waterPatch->soilProp.poreSize + 1;
            double m = land->grid[intersectedPatches[i].first.first][intersectedPatches[i].first.second]->waterPatch->soilProp.poreSize / n;
            double phi = land->grid[intersectedPatches[i].first.first][intersectedPatches[i].first.second]->waterPatch->soilProp.porosity;
            double rw = land->grid[intersectedPatches[i].first.first][intersectedPatches[i].first.second]->waterPatch->soilProp.residualWater;
            double relWaterWP = pow( 1/ (pow(WP * alpha, n) + 1), m) * (phi-rw)+rw;
            avRelSoilMoisture += intersectedPatches[i].second * relRootsPerLayer[j] * land->grid[intersectedPatches[i].first.first][intersectedPatches[i].first.second]->waterPatch->waterLayer[j].relWater;
            if(land->grid[intersectedPatches[i].first.first][intersectedPatches[i].first.second]->waterPatch->waterLayer[j].relWater > relWaterWP)
                avRelSoilMoistureFC += intersectedPatches[i].second * relRootsPerLayer[j] * (land->grid[intersectedPatches[i].first.first][intersectedPatches[i].first.second]->waterPatch->waterLayer[j].relWater - relWaterWP) / land->grid[intersectedPatches[i].first.first][intersectedPatches[i].first.second]->waterPatch->soilProp.fieldCap; //LPJGuess: water content of soil layers as fraction between wilting point (0) and available water holding capacity (1)
            avAbsSoilMoistureNew += intersectedPatches[i].second * relRootsPerLayer[j] * (land->grid[intersectedPatches[i].first.first][intersectedPatches[i].first.second]->waterPatch->waterLayer[j].relWater - relWaterWP) * land->grid[intersectedPatches[i].first.first][intersectedPatches[i].first.second]->waterPatch->waterLayer[j].depthLayer;
            avAbsSoilMoisture += intersectedPatches[i].second * relRootsPerLayer[j] * land->grid[intersectedPatches[i].first.first][intersectedPatches[i].first.second]->waterPatch->waterLayer[j].absWater;
            avNitrogen += (intersectedPatches[i].second * relRootsPerLayer[j] * land->grid[intersectedPatches[i].first.first][intersectedPatches[i].first.second]->nutrientPatch->nutrientLayer[j].NO3_weight);
            avNitrogen += (intersectedPatches[i].second * relRootsPerLayer[j] * land->grid[intersectedPatches[i].first.first][intersectedPatches[i].first.second]->nutrientPatch->nutrientLayer[j].NH4_weight);
            avSoilTemperature += intersectedPatches[i].second * relRootsPerLayer[j] * land->grid[intersectedPatches[i].first.first][intersectedPatches[i].first.second]->nutrientPatch->nutrientLayer[j].soilTemp;
        }
    }
    relCover = 0.0;
    for(int i = 0; i < intersectedPatches.size(); i++) {
        relCover += intersectedPatches[i].second;
    }
}

/////////////////////////////////////////////////////////
void Plant::dispersal()
{
    //after May et al. 2009, submodel: seed production, dispersal and establishment
    int nSeeds = reproductiveMass / (seedMass / 1000000.0); /*!< \todo only a fraction of the reproductive mass schould go into seed production (the rest was innvested for flowers etc) */
    reproductiveMass = 0.0;

    //transfrom measured mean and sd into sigma and mu needed for the lognormal distribution
    double mu = log( meanDispersalDistance / sqrt( 1.0 + sdDispersalDistance/pow(meanDispersalDistance,2) ) );
    double sigma = sqrt( log( 1.0 + sdDispersalDistance/pow(meanDispersalDistance,2) ) );

    //set (random) seed and random generator depending on current time
    myclock::duration d = myclock::now() - beginning;
    unsigned randomSeed = d.count();
    default_random_engine generator(randomSeed);

    //set lognormal ditribution for distance
    lognormal_distribution<double> logDistance(mu, sigma);
    //set uniform ditribution for direction
    uniform_int_distribution<int> uniDirection(0, 360);
    //set uniform ditribution for germination and establishment
    uniform_real_distribution<double> uniGermEstab(0.0,  1.0);

    double establishmentProb = seedMass * 0.000001; /*!< \todo probabiliy to establish; should be (more) PFT-specific and maybe dependent on environmental variables? */

    for(int seed = 0; seed < nSeeds; seed++) {
        //Only if seeds are lucky to germinate and to establish (only based on probabilities) a new seedling will be saved
        double randomNumber = uniGermEstab(generator);
        if(randomNumber <= germinationProb) {
            randomNumber = uniGermEstab(generator);
            if(randomNumber <= establishmentProb) {
                //each seed is dispersed after an log-normal disperal kernel with mean and sd given by plant traits
                double distance = logDistance(generator);
                if(distance <= 0.0) distance = 0.0;
                double direction = uniDirection(generator) * M_PI/180.0;
                seedList.push_back({distance, direction});
            }
            else {
                abovegroundCLitter += (seedMass / 1000000.0);
            }
        }
        else {
            abovegroundCLitter += (seedMass / 1000000.0);
        }
    }
}

/////////////////////////////////////////////////////////
void Plant::resetOutputVariables()
{
    nitrogenLitter = 0;
    carbonLitter = 0;
}

/////////////////////////////////////////////////////////
void Plant::inCaseOfDeath()
{
    abovegroundCLitter += leafCMass;
    abovegroundCLitter += sapwoodCMass;
    abovegroundCLitter += reproductiveMass;
    belowgroundCLitter += rootCMass;
    carbonLitter += leafCMass + rootCMass + sapwoodCMass;
    abovegroundNLitter += leafNMass;
    abovegroundNLitter += sapwoodNMass;
    belowgroundNLitter += rootNMass;
    nitrogenLitter += leafNMass + rootNMass + sapwoodNMass;

    // necessary to save litter in new patch variable as individuals will be removed before patch variables will be calculated
    // also save dead heartwood biomass as standing dead patch biomass
    for(int i = 0; i < intersectedPatches.size(); i++) {
        land->grid[intersectedPatches[i].first.first][intersectedPatches[i].first.second]->vegetationPatch->deadAboveLitterCMass += intersectedPatches[i].second * abovegroundCLitter;
        land->grid[intersectedPatches[i].first.first][intersectedPatches[i].first.second]->vegetationPatch->deadAboveLitterNMass += intersectedPatches[i].second * abovegroundNLitter;
        land->grid[intersectedPatches[i].first.first][intersectedPatches[i].first.second]->vegetationPatch->standingDeadCMass += intersectedPatches[i].second * heartwoodCMass;

        for(int j = 0; j < input.nSoilLayers; j++) {
            land->grid[intersectedPatches[i].first.first][intersectedPatches[i].first.second]->vegetationPatch->deadBelowLitterCMass[j] += intersectedPatches[i].second * relRootsPerLayer[j] * belowgroundCLitter;
            land->grid[intersectedPatches[i].first.first][intersectedPatches[i].first.second]->vegetationPatch->deadBelowLitterNMass[j] += intersectedPatches[i].second * relRootsPerLayer[j] * belowgroundNLitter;
        }
    }
}

/////////////////////////////////////////////////////////
Plant::~Plant()
{
    //dtor
}
