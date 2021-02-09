#include "VegetationLandscape.h"
#include "Patch.h"
#include "Input.h"
#include "WoodyPlant.h"
#include <iostream>
#include <cmath>
#include <random>
#include <algorithm>

/////////////////////////////////////////////////////////
void VegetationLandscape::distributePlants()
{
    for(int plant = 0; plant < iniVegetation.size(); plant++) {
        plantList.push_back(new WoodyPlant(plant, input, land));
        plantIDs++;
    }
}

/////////////////////////////////////////////////////////
void VegetationLandscape::plantsPerPatch()
{
    for(unsigned int x = 0; x < grid.size(); x++) {
        for(unsigned int y = 0; y < grid.at(x).size(); y++) {
            grid[x][y]->vegetationPatch->deletePlantIDList();
        }
    }

    for(int plant = 0; plant < plantList.size(); plant++) {

        int xLengthUntilMyPatch = min(xCells*input.cellSize-input.cellSize, int(plantList[plant]->xCor / input.cellSize) * input.cellSize);
        int yLengthUntilMyPatch = min(yCells*input.cellSize-input.cellSize, int(plantList[plant]->yCor / input.cellSize) * input.cellSize);
        double xDistanceToPreviousPatch = (plantList[plant]->xCor - xLengthUntilMyPatch);
        double yDistanceToPreviousPatch = (plantList[plant]->yCor - yLengthUntilMyPatch) ;
        int xLengthWithMyPatch = min(xCells*input.cellSize, int(plantList[plant]->xCor / input.cellSize + 1) * input.cellSize);
        int yLengthWithMyPatch = min(yCells*input.cellSize, int(plantList[plant]->yCor / input.cellSize  + 1) * input.cellSize);
        double xDistanceUntilNextPatch = (xLengthWithMyPatch - plantList[plant]->xCor);
        double yDistanceUntilNextPatch = (yLengthWithMyPatch - plantList[plant]->yCor);
        int nXcells = 0;
        int nYcells = 0;

        if(xDistanceToPreviousPatch < plantList[plant]->crownRadius)
        {
            //xPatch-1/yPatch
            nXcells = int((plantList[plant]->crownRadius - xDistanceToPreviousPatch) / input.cellSize) + 1;
            for(int xc = 0; xc < nXcells; xc++) {
                if(plantList[plant]->xPatch - (xc+1) >= 0) {
                    grid[plantList[plant]->xPatch - (xc+1)][plantList[plant]->yPatch]->vegetationPatch->presentPlants.push_back(make_pair(plantList[plant]->myID, 0.0));
                }
                else {
                    grid[xCells - (xc+1)][plantList[plant]->yPatch]->vegetationPatch->presentPlants.push_back(make_pair(plantList[plant]->myID, 1));
                }
            }

            if(yDistanceToPreviousPatch < plantList[plant]->crownRadius)
            {
               //xPatch-1/yPatch-1
                nYcells = int((plantList[plant]->crownRadius - yDistanceToPreviousPatch) / input.cellSize) + 1;
                for(int xc = 0; xc < nXcells; xc++) {
                    for(int yc = 0; yc < nYcells; yc++) {
                        if(plantList[plant]->xPatch - (xc+1) >= 0 && plantList[plant]->yPatch - (yc+1) >= 0) {
                            grid[plantList[plant]->xPatch - (xc+1)][plantList[plant]->yPatch - (yc+1)]->vegetationPatch->presentPlants.push_back(make_pair(plantList[plant]->myID, 0.0));
                        }
                        else if(plantList[plant]->xPatch - (xc+1) < 0 && plantList[plant]->yPatch - (yc+1) >= 0) {
                            grid[xCells - (xc+1)][plantList[plant]->yPatch - (yc+1)]->vegetationPatch->presentPlants.push_back(make_pair(plantList[plant]->myID, 1));
                        }
                        else if(plantList[plant]->xPatch - (xc+1) >= 0 && plantList[plant]->yPatch - (yc+1) < 0){
                            grid[plantList[plant]->xPatch - (xc+1)][yCells - (yc+1)]->vegetationPatch->presentPlants.push_back(make_pair(plantList[plant]->myID, 3));
                        }
                        else if(plantList[plant]->xPatch - (xc+1) < 0 && plantList[plant]->yPatch - (yc+1) < 0){
                            grid[xCells - (xc+1)][yCells - (yc+1)]->vegetationPatch->presentPlants.push_back(make_pair(plantList[plant]->myID, 2));
                        }
                    }
                }
            }

            if(yDistanceUntilNextPatch < plantList[plant]->crownRadius)
            {
               //xPatch-1/yPatch+1
                nYcells = int((plantList[plant]->crownRadius - yDistanceUntilNextPatch) / input.cellSize) + 1;
                for(int xc = 0; xc < nXcells; xc++) {
                    for(int yc = 0; yc < nYcells; yc++) {
                        if(plantList[plant]->xPatch - (xc+1) >= 0 && plantList[plant]->yPatch + (yc+1) < yCells) {
                            grid[plantList[plant]->xPatch - (xc+1)][plantList[plant]->yPatch + (yc+1)]->vegetationPatch->presentPlants.push_back(make_pair(plantList[plant]->myID, 0.0));
                        }
                        else if(plantList[plant]->xPatch - (xc+1) < 0 && plantList[plant]->yPatch + (yc+1) < yCells) {
                            grid[xCells - (xc+1)][plantList[plant]->yPatch + (yc+1)]->vegetationPatch->presentPlants.push_back(make_pair(plantList[plant]->myID, 1));
                        }
                        else if(plantList[plant]->xPatch - (xc+1) >= 0 && plantList[plant]->yPatch + (yc+1) >= yCells){
                            grid[plantList[plant]->xPatch - (xc+1)][yc]->vegetationPatch->presentPlants.push_back(make_pair(plantList[plant]->myID, 7));
                        }
                        else if(plantList[plant]->xPatch - (xc+1) < 0 && plantList[plant]->yPatch + (yc+1) >= yCells){
                            grid[xCells - (xc+1)][yc]->vegetationPatch->presentPlants.push_back(make_pair(plantList[plant]->myID, 8));
                        }
                    }
                }
            }
        }

        if(yDistanceToPreviousPatch < plantList[plant]->crownRadius)
        {
            //xPatch/yPatch-1
            nYcells = int((plantList[plant]->crownRadius - yDistanceToPreviousPatch) / input.cellSize) + 1;
            for(int yc = 0; yc < nYcells; yc++) {
                if(plantList[plant]->yPatch - (yc+1) >= 0) {
                    grid[plantList[plant]->xPatch][plantList[plant]->yPatch-(yc+1)]->vegetationPatch->presentPlants.push_back(make_pair(plantList[plant]->myID, 0.0));
                }
                else {
                    grid[plantList[plant]->xPatch][yCells - (yc+1)]->vegetationPatch->presentPlants.push_back(make_pair(plantList[plant]->myID, 3));
                }
            }
        }

        //xPatch/yPatch
        grid[plantList[plant]->xPatch][plantList[plant]->yPatch]->vegetationPatch->presentPlants.push_back(make_pair(plantList[plant]->myID, 0.0));

        if(yDistanceUntilNextPatch < plantList[plant]->crownRadius)
        {
            //xPatch/yPatch+1
            nYcells = int((plantList[plant]->crownRadius - yDistanceUntilNextPatch) / input.cellSize) + 1;
            for(int yc = 0; yc < nYcells; yc++) {
                if(plantList[plant]->yPatch + (yc+1) < yCells) {
                    grid[plantList[plant]->xPatch][plantList[plant]->yPatch+(yc+1)]->vegetationPatch->presentPlants.push_back(make_pair(plantList[plant]->myID, 0.0));
                }
                else {
                    grid[plantList[plant]->xPatch][yc]->vegetationPatch->presentPlants.push_back(make_pair(plantList[plant]->myID, 7));
                }
            }
        }

        if(xDistanceUntilNextPatch < plantList[plant]->crownRadius)
        {
            //xPatch+1/yPatch
            nXcells = int((plantList[plant]->crownRadius - xDistanceUntilNextPatch) / input.cellSize) + 1;
            for(int xc = 0; xc < nXcells; xc++) {
                if(plantList[plant]->xPatch + (xc+1) < xCells) {
                    grid[plantList[plant]->xPatch + (xc+1)][plantList[plant]->yPatch]->vegetationPatch->presentPlants.push_back(make_pair(plantList[plant]->myID, 0.0));
                }
                else {
                    grid[xc][plantList[plant]->yPatch]->vegetationPatch->presentPlants.push_back(make_pair(plantList[plant]->myID, 5));
                }
            }


            if(yDistanceUntilNextPatch < plantList[plant]->crownRadius)
            {
               //xPatch+1/yPatch+1
                nYcells = int((plantList[plant]->crownRadius - yDistanceUntilNextPatch) / input.cellSize) + 1;
                for(int xc = 0; xc < nXcells; xc++) {
                    for(int yc = 0; yc < nYcells; yc++) {
                        if(plantList[plant]->xPatch + (xc+1) < xCells && plantList[plant]->yPatch + (yc+1) < yCells) {
                            grid[plantList[plant]->xPatch + (xc+1)][plantList[plant]->yPatch + (yc+1)]->vegetationPatch->presentPlants.push_back(make_pair(plantList[plant]->myID, 0.0));
                        }
                        else if(plantList[plant]->xPatch + (xc+1) >= xCells && plantList[plant]->yPatch + (yc+1) < yCells) {
                            grid[xc][plantList[plant]->yPatch + (yc+1)]->vegetationPatch->presentPlants.push_back(make_pair(plantList[plant]->myID, 5));
                        }
                        else if(plantList[plant]->xPatch + (xc+1) < xCells && plantList[plant]->yPatch + (yc+1) >= yCells){
                            grid[plantList[plant]->xPatch + (xc+1)][yc]->vegetationPatch->presentPlants.push_back(make_pair(plantList[plant]->myID, 7));
                        }
                        else if(plantList[plant]->xPatch + (xc+1) >= xCells && plantList[plant]->yPatch + (yc+1) >= yCells){
                            grid[xc][yc]->vegetationPatch->presentPlants.push_back(make_pair(plantList[plant]->myID, 6));
                        }
                    }
                }
            }

            if(yDistanceToPreviousPatch < plantList[plant]->crownRadius)
            {
               //xPatch+1/yPatch-1
                 nYcells = int((plantList[plant]->crownRadius - yDistanceToPreviousPatch) / input.cellSize) + 1;
                 for(int xc = 0; xc < nXcells; xc++) {
                     for(int yc = 0; yc < nYcells; yc++) {
                         if(plantList[plant]->xPatch + (xc+1) < xCells && plantList[plant]->yPatch - (yc+1) >= 0) {
                             grid[plantList[plant]->xPatch + (xc+1)][plantList[plant]->yPatch - (yc+1)]->vegetationPatch->presentPlants.push_back(make_pair(plantList[plant]->myID, 0.0));
                         }
                         else if(plantList[plant]->xPatch + (xc+1) >= xCells && plantList[plant]->yPatch - (yc+1) >= 0) {
                             grid[xc][plantList[plant]->yPatch - (yc+1)]->vegetationPatch->presentPlants.push_back(make_pair(plantList[plant]->myID, 5));
                         }
                         else if(plantList[plant]->xPatch + (xc+1) < xCells && plantList[plant]->yPatch - (yc+1) < 0){
                             grid[plantList[plant]->xPatch + (xc+1)][yCells - (yc+1)]->vegetationPatch->presentPlants.push_back(make_pair(plantList[plant]->myID, 3));
                         }
                         else if(plantList[plant]->xPatch + (xc+1) >= xCells && plantList[plant]->yPatch - (yc+1) < 0){
                             grid[xc][yCells - (yc+1)]->vegetationPatch->presentPlants.push_back(make_pair(plantList[plant]->myID, 4));
                         }
                     }
                 }
             }
        }
    }
}

/////////////////////////////////////////////////////////
void VegetationLandscape::calculatePatchVariables()
{
    plantsPerPatch();
    coverPerPlantAndPatch();
    rootsPerPatch();
    patchesPerPlant();
    variablesPerPatch();
}

/////////////////////////////////////////////////////////
/// Based on implementation of
/// Author: Emanuel Jöbstl <emi@eex-dev.net>
/// Date  : 18.06.2011
/// Link  : http://www.eex-dev.net/index.php?id=100
/////////////////////////////////////////////////////////
void VegetationLandscape::coverPerPlantAndPatch()
{
    //The resolution to use for approximation.
    const double resolution = 0.01;

    for(unsigned int x = 0; x < grid.size(); x++) {
        for(unsigned int y = 0; y < grid.at(x).size(); y++) {
            grid[x][y]->vegetationPatch->absVegCover = 0;
            for(int patchPlant = 0; patchPlant < grid[x][y]->vegetationPatch->presentPlants.size(); patchPlant++) {
                for(int plant = 0; plant < plantList.size(); plant++) {
                    if(plantList[plant]->myID == grid[x][y]->vegetationPatch->presentPlants[patchPlant].first)
                    {
                        //A variable storing the nearest horizontal edge of the patch.
                        double nearestRectangleEdge = 0;
                        int boundaryProblem = 0;
                        if(grid[x][y]->vegetationPatch->presentPlants[patchPlant].second > 0) {
                            boundaryProblem = grid[x][y]->vegetationPatch->presentPlants[patchPlant].second;
                            grid[x][y]->vegetationPatch->presentPlants[patchPlant].second = 0.0;
                        }

                        //Determine what is nearer to the plant crown center - the patch top edge or the patch bottom edge
                        else if(abs(plantList[plant]->yCor - (grid[x][y]->vegetationPatch->yCor + input.cellSize/2.0)) >= abs(plantList[plant]->yCor - (grid[x][y]->vegetationPatch->yCor - input.cellSize/2.0)))
                        {
                            nearestRectangleEdge = (grid[x][y]->vegetationPatch->yCor - input.cellSize/2.0);
                        }

                        else if(abs(plantList[plant]->yCor - (grid[x][y]->vegetationPatch->yCor + input.cellSize/2.0)) < abs(plantList[plant]->yCor - (grid[x][y]->vegetationPatch->yCor - input.cellSize/2.0)))
                        {
                            nearestRectangleEdge = (grid[x][y]->vegetationPatch->yCor + input.cellSize/2.0);
                        }

                        //The bounds of our integration
                        double leftBound = 0;
                        double rightBound = 0;

                        if(boundaryProblem == 1 || boundaryProblem == 2) {
                            leftBound = max(-plantList[plant]->crownRadius + (xCells*input.cellSize) - plantList[plant]->xCor, ((xCells*input.cellSize) + ((xCells*input.cellSize) - (grid[x][y]->vegetationPatch->xCor + input.cellSize/2.0))));
                            rightBound = min(plantList[plant]->crownRadius + (xCells*input.cellSize) - plantList[plant]->xCor, ((xCells*input.cellSize) + ((xCells*input.cellSize) - (grid[x][y]->vegetationPatch->xCor - input.cellSize/2.0))));
                        }

                        else if(boundaryProblem == 4) {
                            leftBound = max(-plantList[plant]->crownRadius + (xCells*input.cellSize) - plantList[plant]->yCor, ((xCells*input.cellSize) + ((xCells*input.cellSize) - ((xCells*input.cellSize) - grid[x][y]->vegetationPatch->xCor + input.cellSize/2.0))));
                            rightBound = min(plantList[plant]->crownRadius + (xCells*input.cellSize) - plantList[plant]->yCor, ((xCells*input.cellSize) + ((xCells*input.cellSize) - ((xCells*input.cellSize) - grid[x][y]->vegetationPatch->xCor - input.cellSize/2.0))));
                        }

                        else if(boundaryProblem == 5) {
                            leftBound = grid[x][y]->vegetationPatch->xCor - input.cellSize/2.0;
                            rightBound = min(plantList[plant]->crownRadius - (xCells*input.cellSize - plantList[plant]->xCor) - (grid[x][y]->vegetationPatch->xCor - input.cellSize/2.0), (grid[x][y]->vegetationPatch->xCor + input.cellSize/2.0));
                        }

                        else if(boundaryProblem == 6) {
                            leftBound = max(-plantList[plant]->crownRadius + (xCells*input.cellSize) - ((xCells*input.cellSize) - plantList[plant]->xCor), ((xCells*input.cellSize) + ((xCells*input.cellSize) - (((xCells*input.cellSize) - grid[x][y]->vegetationPatch->xCor) + input.cellSize/2.0))));
                            rightBound = min(plantList[plant]->crownRadius + (xCells*input.cellSize) - ((xCells*input.cellSize) - plantList[plant]->xCor), ((xCells*input.cellSize) + ((xCells*input.cellSize) - (((xCells*input.cellSize) - grid[x][y]->vegetationPatch->xCor) - input.cellSize/2.0))));
                        }

                        else if(boundaryProblem == 7 || boundaryProblem == 3)
                        {
                            leftBound = max(-plantList[plant]->crownRadius + plantList[plant]->xCor, (grid[x][y]->vegetationPatch->xCor - input.cellSize/2.0));
                            rightBound = min(plantList[plant]->crownRadius + plantList[plant]->xCor, (grid[x][y]->vegetationPatch->xCor + input.cellSize/2.0));
                        }

                        else if(boundaryProblem == 8) {
                            leftBound = max(-plantList[plant]->crownRadius + (xCells*input.cellSize) - ((yCells*input.cellSize) - plantList[plant]->yCor), ((xCells*input.cellSize) + ((xCells*input.cellSize) - (grid[x][y]->vegetationPatch->xCor + input.cellSize/2.0))));
                            rightBound = min(plantList[plant]->crownRadius + (xCells*input.cellSize) - ((yCells*input.cellSize) - plantList[plant]->yCor), ((xCells*input.cellSize) + ((xCells*input.cellSize) - (grid[x][y]->vegetationPatch->xCor - input.cellSize/2.0))));
                        }

                        else if(plantList[plant]->yCor <= (grid[x][y]->vegetationPatch->yCor + input.cellSize/2.0) && plantList[plant]->yCor >= (grid[x][y]->vegetationPatch->yCor - input.cellSize/2.0))
                        {
                            //Take care if the plant crown center lies within the patch.
                            leftBound = max(-plantList[plant]->crownRadius + plantList[plant]->xCor, (grid[x][y]->vegetationPatch->xCor - input.cellSize/2.0));
                            rightBound = min(plantList[plant]->crownRadius + plantList[plant]->xCor, (grid[x][y]->vegetationPatch->xCor + input.cellSize/2.0));
                        }

                        else if(plantList[plant]->crownRadius >= abs(nearestRectangleEdge - plantList[plant]->yCor))
                        {
                            //If the plant crown center lies outside of the patch, we can choose optimal bounds.
                            leftBound = max(-sqrt(pow(plantList[plant]->crownRadius, 2) - abs(pow(nearestRectangleEdge - plantList[plant]->yCor, 2))) + plantList[plant]->xCor, (grid[x][y]->vegetationPatch->xCor - input.cellSize/2.0));
                            rightBound = min(sqrt(pow(plantList[plant]->crownRadius, 2) - abs(pow(nearestRectangleEdge - plantList[plant]->yCor, 2))) + plantList[plant]->xCor, (grid[x][y]->vegetationPatch->xCor + input.cellSize/2.0));
                        }

                        double upperBound;
                        double lowerBound;

                        for(double i = leftBound + resolution; i <= rightBound; i += resolution)
                        {
                            if(boundaryProblem == 1){
                                lowerBound = max((grid[x][y]->vegetationPatch->yCor - input.cellSize/2.0), plantList[plant]->yCor - sqrt(pow(plantList[plant]->crownRadius, 2) - pow((i - resolution / 2.0 - ((xCells*input.cellSize) - plantList[plant]->xCor)), 2)));
                                upperBound = min((grid[x][y]->vegetationPatch->yCor + input.cellSize/2.0), plantList[plant]->yCor + sqrt(pow(plantList[plant]->crownRadius, 2) - pow((i - resolution / 2.0 - ((xCells*input.cellSize) - plantList[plant]->xCor)), 2)));
                            }

                            else if(boundaryProblem == 2){
                                lowerBound = max((yCells*input.cellSize + (yCells*input.cellSize - (grid[x][y]->vegetationPatch->yCor + input.cellSize/2.0))), yCells*input.cellSize - plantList[plant]->yCor - sqrt(pow(plantList[plant]->crownRadius, 2) - pow((i - resolution / 2.0 - (xCells*input.cellSize - plantList[plant]->xCor)), 2)));
                                upperBound = min((yCells*input.cellSize + (yCells*input.cellSize - (grid[x][y]->vegetationPatch->yCor - input.cellSize/2.0))), yCells*input.cellSize - plantList[plant]->yCor + sqrt(pow(plantList[plant]->crownRadius, 2) - pow((i - resolution / 2.0 - (xCells*input.cellSize - plantList[plant]->xCor)), 2)));
                            }

                            else if(boundaryProblem == 3){
                                upperBound = min((grid[x][y]->vegetationPatch->yCor + input.cellSize/2.0), ((yCells*input.cellSize) + plantList[plant]->yCor) + sqrt(pow(plantList[plant]->crownRadius, 2) - pow((i - resolution / 2.0 - plantList[plant]->xCor), 2)));
                                lowerBound = max((grid[x][y]->vegetationPatch->yCor - input.cellSize/2.0), min(upperBound, ((yCells*input.cellSize) + plantList[plant]->yCor) - sqrt(pow(plantList[plant]->crownRadius, 2) - pow((i - resolution / 2.0 - plantList[plant]->xCor), 2))));
                            }

                            else if(boundaryProblem == 4){
                                lowerBound = max((yCells*input.cellSize + (yCells*input.cellSize - (grid[x][y]->vegetationPatch->yCor + input.cellSize/2.0))), yCells*input.cellSize - (xCells*input.cellSize - plantList[plant]->xCor) - sqrt(pow(plantList[plant]->crownRadius, 2) - pow((i - resolution / 2.0 - (xCells*input.cellSize - plantList[plant]->yCor)), 2)));
                                upperBound = min((yCells*input.cellSize + (yCells*input.cellSize - (grid[x][y]->vegetationPatch->yCor - input.cellSize/2.0))), yCells*input.cellSize - (xCells*input.cellSize - plantList[plant]->xCor) + sqrt(pow(plantList[plant]->crownRadius, 2) - pow((i - resolution / 2.0 - (xCells*input.cellSize - plantList[plant]->yCor)), 2)));
                            }

                            else if(boundaryProblem == 5){
                                lowerBound = max((grid[x][y]->vegetationPatch->yCor - input.cellSize/2.0), plantList[plant]->yCor - sqrt(pow(plantList[plant]->crownRadius, 2) - pow((i + (xCells*input.cellSize) - resolution / 2.0 - plantList[plant]->xCor), 2)));
                                upperBound = min((grid[x][y]->vegetationPatch->yCor + input.cellSize/2.0), plantList[plant]->yCor + sqrt(pow(plantList[plant]->crownRadius, 2) - pow((i + (xCells*input.cellSize) - resolution / 2.0 - plantList[plant]->xCor), 2)));
                            }

                            else if(boundaryProblem == 6){
                                lowerBound = max((yCells*input.cellSize + (yCells*input.cellSize - ((yCells*input.cellSize - grid[x][y]->vegetationPatch->yCor) + input.cellSize/2.0))), yCells*input.cellSize - (yCells*input.cellSize - plantList[plant]->yCor) - sqrt(pow(plantList[plant]->crownRadius, 2) - pow((i - resolution / 2.0 - (xCells*input.cellSize - (xCells*input.cellSize - plantList[plant]->xCor))), 2)));
                                upperBound = min((yCells*input.cellSize + (yCells*input.cellSize - ((yCells*input.cellSize - grid[x][y]->vegetationPatch->yCor) - input.cellSize/2.0))), yCells*input.cellSize - (yCells*input.cellSize - plantList[plant]->yCor) + sqrt(pow(plantList[plant]->crownRadius, 2) - pow((i - resolution / 2.0 - (xCells*input.cellSize - (xCells*input.cellSize - plantList[plant]->xCor))), 2)));
                            }

                            else if(boundaryProblem == 7){
                                lowerBound = max((yCells*input.cellSize + grid[x][y]->vegetationPatch->yCor - input.cellSize/2.0), plantList[plant]->yCor - sqrt(pow(plantList[plant]->crownRadius, 2) - pow((i - resolution / 2.0 - plantList[plant]->xCor), 2)));
                                upperBound = min((yCells*input.cellSize + grid[x][y]->vegetationPatch->yCor + input.cellSize/2.0), max(lowerBound, plantList[plant]->yCor + sqrt(pow(plantList[plant]->crownRadius, 2) - pow((i - resolution / 2.0 - plantList[plant]->xCor), 2))));
                            }

                            else if(boundaryProblem == 8){
                                lowerBound = max((yCells*input.cellSize + (yCells*input.cellSize - ((yCells*input.cellSize - grid[x][y]->vegetationPatch->yCor) + input.cellSize/2.0))), yCells*input.cellSize - plantList[plant]->xCor - sqrt(pow(plantList[plant]->crownRadius, 2) - pow((i - resolution / 2.0 - (xCells*input.cellSize - (yCells*input.cellSize - plantList[plant]->yCor))), 2)));
                                upperBound = min((yCells*input.cellSize + (yCells*input.cellSize - ((yCells*input.cellSize - grid[x][y]->vegetationPatch->yCor) - input.cellSize/2.0))), yCells*input.cellSize - plantList[plant]->xCor + sqrt(pow(plantList[plant]->crownRadius, 2) - pow((i - resolution / 2.0 - (xCells*input.cellSize - (yCells*input.cellSize - plantList[plant]->yCor))), 2)));
                            }

                            else {
                                lowerBound = max((grid[x][y]->vegetationPatch->yCor - input.cellSize/2.0), plantList[plant]->yCor - sqrt(pow(plantList[plant]->crownRadius, 2) - pow((i - resolution / 2.0 - plantList[plant]->xCor), 2)));
                                upperBound = min((grid[x][y]->vegetationPatch->yCor + input.cellSize/2.0), plantList[plant]->yCor + sqrt(pow(plantList[plant]->crownRadius, 2) - pow((i - resolution / 2.0 - plantList[plant]->xCor), 2)));
                            }
                            grid[x][y]->vegetationPatch->presentPlants[patchPlant].second += (upperBound - lowerBound) * resolution;
                            grid[x][y]->vegetationPatch->absVegCover += (upperBound - lowerBound) * resolution;
                        }
                        break;
                    }
                }
            }
            grid[x][y]->vegetationPatch->relVegCover = grid[x][y]->vegetationPatch->absVegCover / pow(input.cellSize, 2);
        }
    }
}

/////////////////////////////////////////////////////////
void VegetationLandscape::rootsPerPatch()
{
    for(unsigned int x = 0; x < grid.size(); x++) {
        for(unsigned int y = 0; y < grid.at(x).size(); y++) {
            for(unsigned int layer = 0; layer < input.nSoilLayers; layer++) {
                grid[x][y]->vegetationPatch->relativeRoots[layer] = 0;
                for(int patchPlant = 0; patchPlant < grid[x][y]->vegetationPatch->presentPlants.size(); patchPlant++) {
                    for(int plant = 0; plant < plantList.size(); plant++) {
                        if(plantList[plant]->myID == grid[x][y]->vegetationPatch->presentPlants[patchPlant].first)
                        {
                            //caluclate the relative cover for the plant for this patch given its total plant cover
                            double relativeCover = grid[x][y]->vegetationPatch->presentPlants[patchPlant].second / plantList[plant]->crownArea;
                            //caluclate the relative roots for this patch given the single plant covers and their fractional roots in the respective layers
                            grid[x][y]->vegetationPatch->relativeRoots[layer] += relativeCover * plantList[plant]->relRootsPerLayer[layer];
                        }
                    }
                }
            }
        }
    }
}

/////////////////////////////////////////////////////////
void VegetationLandscape::variablesPerPatch()
{
    for(unsigned int x = 0; x < grid.size(); x++) {
        for(unsigned int y = 0; y < grid.at(x).size(); y++) {

            //Set all patch variables to 0
            grid[x][y]->vegetationPatch->abovegroundCMass = 0.0;
            grid[x][y]->vegetationPatch->abovegroundNMass = 0.0;
            grid[x][y]->vegetationPatch->totalAliveCMass = 0.0;
            grid[x][y]->vegetationPatch->relFPC = 0.0;
            for(unsigned int layer = 0; layer < input.nSoilLayers; layer++) {
                grid[x][y]->vegetationPatch->fracRelWaterWP[layer] = 0.0;
                grid[x][y]->vegetationPatch->belowgroundCMass[layer] = 0.0;
                grid[x][y]->vegetationPatch->belowgroundNMass[layer] = 0.0;
            }

            for(int patchPlant = 0; patchPlant < grid[x][y]->vegetationPatch->presentPlants.size(); patchPlant++) {
                for(int plant = 0; plant < plantList.size(); plant++) {
                    if(plantList[plant]->myID == grid[x][y]->vegetationPatch->presentPlants[patchPlant].first)
                    {
                        //calculate the relative plant cover needed to calculate the proportions for each plant variable for each patch
                        double relativeCover = grid[x][y]->vegetationPatch->presentPlants[patchPlant].second / plantList[plant]->crownArea;
                        //calculate the absolute FPC for this patch
                        grid[x][y]->vegetationPatch->relFPC += relativeCover * plantList[plant]->FPC;
                        //calculate the aboveground C and N mass in patch given the single plant covers
                        grid[x][y]->vegetationPatch->abovegroundCMass += relativeCover * (plantList[plant]->leafCMass + plantList[plant]->sapwoodCMass + plantList[plant]->heartwoodCMass);
                        grid[x][y]->vegetationPatch->abovegroundNMass += relativeCover * (plantList[plant]->leafNMass + plantList[plant]->sapwoodNMass);
                        grid[x][y]->vegetationPatch->totalAliveCMass += relativeCover * (plantList[plant]->leafCMass + plantList[plant]->rootCMass + plantList[plant]->sapwoodCMass + plantList[plant]->heartwoodCMass);
                        //calculate the total dead aboveground C and N mass in patch given the single plant covers
                        grid[x][y]->vegetationPatch->abovegroundLitterCMass += relativeCover * plantList[plant]->abovegroundCLitter;
                        grid[x][y]->vegetationPatch->abovegroundLitterNMass += relativeCover * plantList[plant]->abovegroundNLitter;
                        //add dead aboveground C and N mass in patch from dead plants
                        grid[x][y]->vegetationPatch->abovegroundLitterCMass += grid[x][y]->vegetationPatch->deadAboveLitterCMass;
                        grid[x][y]->vegetationPatch->deadAboveLitterCMass = 0.0;
                        grid[x][y]->vegetationPatch->abovegroundLitterNMass += grid[x][y]->vegetationPatch->deadAboveLitterNMass;
                        grid[x][y]->vegetationPatch->deadAboveLitterNMass = 0.0;
                        for(unsigned int layer = 0; layer < input.nSoilLayers; layer++) {
                            //calculate belowground C and N mass
                            grid[x][y]->vegetationPatch->belowgroundCMass[layer] += relativeCover * plantList[plant]->relRootsPerLayer[layer] * plantList[plant]->rootCMass;
                            grid[x][y]->vegetationPatch->belowgroundNMass[layer] += relativeCover * plantList[plant]->relRootsPerLayer[layer] * plantList[plant]->rootNMass;
                            //calculate dead belowground C and N mass
                            grid[x][y]->vegetationPatch->belowgroundLitterCMass[layer] += relativeCover * plantList[plant]->relRootsPerLayer[layer] * plantList[plant]->belowgroundCLitter;
                            grid[x][y]->vegetationPatch->belowgroundLitterNMass[layer] += relativeCover * plantList[plant]->relRootsPerLayer[layer] * plantList[plant]->belowgroundNLitter;
                            //add dead belowground C and N mass from dead plants
                            grid[x][y]->vegetationPatch->belowgroundLitterCMass[layer] += grid[x][y]->vegetationPatch->deadBelowLitterCMass[layer];
                            grid[x][y]->vegetationPatch->belowgroundLitterNMass[layer] += grid[x][y]->vegetationPatch->deadBelowLitterNMass[layer];
                            //calculate the relative water content at wilting point for each layer of this patch

                            // first calculate wilting point for the soil texture in which the plant is rooting here (after van Genuchten, Maidment, p.5.6 / 5.14)
                            double alpha = pow(grid[x][y]->waterPatch->soilProp.bubPressure, -1);
                            double n = grid[x][y]->waterPatch->soilProp.poreSize + 1;
                            double m = grid[x][y]->waterPatch->soilProp.poreSize / n;
                            double phi = grid[x][y]->waterPatch->soilProp.porosity;
                            double rw = grid[x][y]->waterPatch->soilProp.residualWater;
                            double relWaterWP = pow( 1/ (pow(plantList[plant]->WP * alpha, n) + 1), m) * (phi-rw) + rw;

                            grid[x][y]->vegetationPatch->fracRelWaterWP[layer] += relativeCover * plantList[plant]->relRootsPerLayer[layer] * relWaterWP;
                        }
                    }
                }
            }
            //calculate the relative FPC for this patch depending on the patch size
            grid[x][y]->vegetationPatch->relFPC /= pow(input.cellSize, 2);
        }
    }

    for(int plant = 0; plant < plantList.size(); plant++) {
            plantList[plant]->abovegroundCLitter = 0.0;
            plantList[plant]->abovegroundNLitter = 0.0;
            plantList[plant]->belowgroundCLitter = 0.0;
            plantList[plant]->belowgroundNLitter = 0.0;
    }
}

/////////////////////////////////////////////////////////
void VegetationLandscape::patchesPerPlant()
{
    for(int plant = 0; plant < plantList.size(); plant++) {
        plantList[plant]->intersectedPatches.erase(plantList[plant]->intersectedPatches.begin(), plantList[plant]->intersectedPatches.end());
    }

    for(unsigned int x = 0; x < grid.size(); x++) {
        for(unsigned int y = 0; y < grid.at(x).size(); y++) {
            for(int patchPlant = 0; patchPlant < grid[x][y]->vegetationPatch->presentPlants.size(); patchPlant++) {
                for(int plant = 0; plant < plantList.size(); plant++) {
                    if(plantList[plant]->myID == grid[x][y]->vegetationPatch->presentPlants[patchPlant].first)
                    {
                        double relCover = abs(grid[x][y]->vegetationPatch->presentPlants[patchPlant].second) / plantList[plant]->crownArea;
                        //plantList[plant]->intersectedPatches.push_back({int(x), int(y), relCover});
                        plantList[plant]->intersectedPatches.push_back(make_pair(make_pair(int(x),int(y)), relCover));
                    }
                }
            }
        }
    }
}

/////////////////////////////////////////////////////////
void VegetationLandscape::calculateProcesses(int day, double dayLength)
{
    int motherPlants = plantList.size();

    //To make sure that no indvidual profits just because of always being first called here
    auto rng = default_random_engine {};
    shuffle(begin(plantList), end(plantList), rng);

    for(int plant = 0; plant < motherPlants; plant++) {
        plantList[plant]->resetOutputVariables();
        plantList[plant]->availablePools();
        plantList[plant]->photosynthesis(day, dayLength, plantList[plant]->optLambda); //Photosynthesis without water and nitrogen stress
        plantList[plant]->phenology(day);
        plantList[plant]->transpiration(day, dayLength); //Transpiration (if water limitation, recalculate photosynthesis)
        if(plantList[plant]->phenologyStatus > 0) plantList[plant]->nitrogenUptake(day, dayLength); //Nitrogen uptake (if nitrogen limitation, recalculate photosynthesis)
        plantList[plant]->respiration(day);
        plantList[plant]->reproduction();
        plantList[plant]->allocation();
        dispersal(day, plant);
        //plantList[plant]->tissueTurnover(); //not used anymore because of tissue turnover happening in allocation process
        vegMortality(day, plant);
    }
}

/////////////////////////////////////////////////////////
void VegetationLandscape::vegMortality(int day, int plant)
{
    //Mortality once a year at the end of the year but only for plants older than a year
    if((day+1) % 365 == 0 && plantList[plant]->age >= 365) plantList[plant]->mortality();
//    if(input.latitude >= 0) {
//        if(input.dayOfYear[day] == 0) plantList[plant]->mortality();
//    }
//    else
//        if(input.dayOfYear[day] == 182) plantList[plant]->mortality();

    if(plantList[plant]->dead) {
        deadPlants.push_back(plantList[plant]->myID);
    }
}

/////////////////////////////////////////////////////////
void VegetationLandscape::deleteDeadPlants()
{
    for(int deadPlant = 0; deadPlant < deadPlants.size(); deadPlant++) {
        for(int plant = 0; plant < plantList.size(); plant++) {
            if(deadPlants[deadPlant] == plantList[plant]->myID) {
                delete plantList[plant];
                plantList.erase(plantList.begin() + plant);
            }
        }
    }
    deadPlants.erase(deadPlants.begin(), deadPlants.end());
}

/////////////////////////////////////////////////////////
void VegetationLandscape::dispersal(int day, int plant)
{
    if(input.dayOfYear[day] == 182) { /*!< \todo should be PFT-specific and maybe dependent on environmental variables? */
        plantList[plant]->dispersal();

        for(int seed = 0; seed < plantList[plant]->seedList.size(); seed++) {
            //Calculate the new coordinates of the new plant based on mother plant's position, dispersal distance and direction (according to Polar coordinate system)
            double xCor = plantList[plant]->xCor + ((plantList[plant]->seedList[seed].distance * cos(plantList[plant]->seedList[seed].direction))) / 100.0;
            double yCor = plantList[plant]->yCor + ((plantList[plant]->seedList[seed].distance * sin(plantList[plant]->seedList[seed].direction))) / 100.0;

            //If outside landscape, adapt coordinates for closed boundary conditions
            if(xCor < 0.0) xCor = xLength + xCor;
            else if(xCor > xLength) xCor = xCor - xLength;
            if(yCor < 0.0) yCor = yLength + yCor;
            else if(yCor > yLength) yCor = yCor - yLength;

            int xPatch = int(xCor / input.cellSize);
            int yPatch = int(yCor / input.cellSize);

            //Check if there is space within landscape
            /*!< \todo the following not nicely programmed, should be reprogrammed */
            bool bareSoil = true;
            if(grid[xPatch][yPatch]->vegetationPatch->presentPlants.size() > 0) {
                for(int patchPlant = 0; patchPlant < grid[xPatch][yPatch]->vegetationPatch->presentPlants.size(); patchPlant++) {
                    if(!bareSoil) break;
                    int ID = grid[xPatch][yPatch]->vegetationPatch->presentPlants[patchPlant].first;
                    for(int plantInd = 0; plantInd < plantList.size(); plantInd++) {
                        if(plantList[plantInd]->myID == ID) {

                            double xLower = (plantList[plantInd]->xCor - plantList[plantInd]->crownRadius);
                            double yLower = (plantList[plantInd]->yCor - plantList[plantInd]->crownRadius);
                            double xUpper = (plantList[plantInd]->xCor + plantList[plantInd]->crownRadius);
                            double yUpper = (plantList[plantInd]->yCor + plantList[plantInd]->crownRadius);

                            //Taking care of boundary conditions
                            if(xLower < 0.0) xLower = xLength + xLower;
                            if(yLower < 0.0) yLower = yLength + yLower;
                            if(xUpper > xLength) xUpper = xUpper - xLength;
                            if(yUpper > yLength) yUpper = yUpper - yLength;

                            if(xCor > xLower && xCor < xUpper && yCor > yLower && yCor < yUpper) {
                                bareSoil = false;
                                break;
                            }
                        }
                    }
                }
            }

            if(bareSoil) {
                plantList.push_back(new WoodyPlant(plantIDs+1, plantList[plant]->PFTname, xCor, yCor, input, land)); /*!< \todo what about ID of the plant?*/
                plantIDs++;
                newIndividuals += 1;
            }

            else {  //If seeds have no space to establish, their dead mass will be added to the patch's litter pool
                /*!< \todo for a closed carbon cycle should we add the seeds outside the landscape to the litterpool as well? */
                grid[xPatch][yPatch]->vegetationPatch->abovegroundLitterCMass += (plantList[plant]->seedMass / 1000000.0);
            }
    }
    plantList[plant]->seedList.erase(plantList[plant]->seedList.begin(), plantList[plant]->seedList.end());
    }
}

/////////////////////////////////////////////////////////
void VegetationLandscape::shadingAPAR(double extraRadiation)
{
    //Calculation of intersection area based on http://walter.bislins.ch/blog/index.asp?page=Schnittfl%E4che+zweier+Kreise+berechnen
    /*!< \todo not well programmed: to many repetetions of code that could be put into a function */
    /*!< \todo no explicit consideration of where plant overgrow but might be important if ther eis still an unshaded area left */
    for(int plant = 0; plant < plantList.size(); plant++) {
        plantList[plant]->shadedAPAR = 0.0;
        double totalShadedArea = 0.0;
        for(int neighbor = 0; neighbor < plantList.size(); neighbor++) {
            if(plantList[neighbor]->height > plantList[plant]->height && plantList[neighbor]->myID != plantList[plant]->myID) {
                //Calculate distance between plants
                double a = abs(plantList[plant]->xCor - plantList[neighbor]->xCor);
                double b = abs(plantList[plant]->yCor - plantList[neighbor]->yCor);
                double distance = sqrt(pow(a,2) + pow(b,2));
                double myRadius = plantList[plant]->crownRadius;
                double neighborRadius = plantList[neighbor]->crownRadius;
                if(distance <= abs(myRadius-neighborRadius)) {
                    double r = myRadius;
                    if(myRadius > neighborRadius) r = neighborRadius;
                    else r = myRadius;
                    double shadedArea = r*r * M_PI;
                    //Calculate APAR that will be lost through shading
                    if((plantList[plant]->crownArea - totalShadedArea) > 0.0) {
                        calculateAPAR(extraRadiation, shadedArea, totalShadedArea, plant);
                        totalShadedArea += shadedArea; //to not allow shading many times through different overgrowing individuals
                    }
                    else {
                        break;
                    }
                }
                else if(distance < myRadius + neighborRadius) {
                    double shadedArea = calculateIntersectionArea(myRadius, neighborRadius, distance);
                    //Calculate APAR that will be lost through shading
                    if((plantList[plant]->crownArea - totalShadedArea) > 0.0) {
                        calculateAPAR(extraRadiation, shadedArea, totalShadedArea, plant);
                        totalShadedArea += shadedArea; //to not allow shading many times through different overgrowing individuals
                    }
                    else {
                        break;
                    }
                }
                else {
                    //Considering now neigbors due to boundary conditions
                    double txCor = plantList[neighbor]->xCor;
                    double tyCor = plantList[neighbor]->yCor;
                    bool crossingBoundary = false;
                    if(txCor - plantList[neighbor]->crownRadius < 0.0) {
                        txCor = xLength + txCor;
                        crossingBoundary = true;
                    }
                    else if(txCor + plantList[neighbor]->crownRadius > xLength) {
                        txCor = txCor - xLength;
                        crossingBoundary = true;
                    }
                    if(tyCor - plantList[neighbor]->crownRadius < 0.0) {
                        tyCor = yLength + tyCor;
                        crossingBoundary = true;
                    }
                    else if(tyCor + plantList[neighbor]->crownRadius > yLength) {
                        tyCor = tyCor - yLength;
                        crossingBoundary = true;
                    }

                    if(crossingBoundary) {
                        //Calculate distance between plants
                        double a = abs(plantList[plant]->xCor - txCor);
                        double b = abs(plantList[plant]->yCor - tyCor);
                        double distance = sqrt(pow(a,2) + pow(b,2));
                        double myRadius = plantList[plant]->crownRadius;
                        double neighborRadius = plantList[neighbor]->crownRadius;
                        if(distance <= abs(myRadius-neighborRadius)) {
                            double r = myRadius;
                            if(myRadius > neighborRadius) r = neighborRadius;
                            else r = myRadius;
                            double shadedArea = r*r * M_PI;
                            //Calculate APAR that will be lost through shading
                            if((plantList[plant]->crownArea - totalShadedArea) > 0.0) {
                                calculateAPAR(extraRadiation, shadedArea, totalShadedArea, plant);
                                totalShadedArea += shadedArea; //to not allow shading many times through different overgrowing individuals
                            }
                            else {
                                break;
                            }
                        }
                        else if(distance < myRadius + neighborRadius) {
                            double shadedArea = calculateIntersectionArea(myRadius, neighborRadius, distance);
                            //Calculate APAR that will be lost through shading
                            if((plantList[plant]->crownArea - totalShadedArea) > 0.0) {
                                calculateAPAR(extraRadiation, shadedArea, totalShadedArea, plant);
                                totalShadedArea += shadedArea; //to not allow shading many times through different overgrowing individuals
                            }
                            else {
                                break;
                            }
                        }
                    }
                }
            }
        }

        if((plantList[plant]->crownArea - totalShadedArea) > 0.0) {
            //Taking care of boundary conditions for focal plant
            double txCor = plantList[plant]->xCor;
            double tyCor = plantList[plant]->yCor;
            bool crossingBoundary = false;
            if(txCor - plantList[plant]->crownRadius < 0.0) {
                txCor = xLength + txCor;
                crossingBoundary = true;
            }
            else if(txCor + plantList[plant]->crownRadius > xLength) {
                txCor = txCor - xLength;
                crossingBoundary = true;
            }
            if(tyCor - plantList[plant]->crownRadius < 0.0) {
                tyCor = yLength + tyCor;
                crossingBoundary = true;
            }
            else if(tyCor + plantList[plant]->crownRadius > yLength) {
                tyCor = tyCor - yLength;
                crossingBoundary = true;
            }

            if(crossingBoundary) {
                for(int neighbor = 0; neighbor < plantList.size(); neighbor++) {
                    if(plantList[neighbor]->height > plantList[plant]->height && plantList[neighbor]->myID != plantList[plant]->myID) {
                        //Calculate distance between plants
                        double a = abs(txCor - plantList[neighbor]->xCor);
                        double b = abs(tyCor - plantList[neighbor]->yCor);
                        double distance = sqrt(pow(a,2) + pow(b,2));
                        double myRadius = plantList[plant]->crownRadius;
                        double neighborRadius = plantList[neighbor]->crownRadius;
                        if(distance <= abs(myRadius-neighborRadius)) {
                            double r = myRadius;
                            if(myRadius > neighborRadius) r = neighborRadius;
                            else r = myRadius;
                            double shadedArea = r*r * M_PI;
                            //Calculate APAR that will be lost through shading
                            if((plantList[plant]->crownArea - totalShadedArea) > 0.0) {
                                calculateAPAR(extraRadiation, shadedArea, totalShadedArea, plant);
                                totalShadedArea += shadedArea; //to not allow shading many times through different overgrowing individuals
                            }
                            else {
                                break;
                            }
                        }
                        else if (distance < myRadius + neighborRadius) {
                            double shadedArea = calculateIntersectionArea(myRadius, neighborRadius, distance);
                            //Calculate APAR that will be lost through shading
                            if((plantList[plant]->crownArea - totalShadedArea) > 0.0) {
                                calculateAPAR(extraRadiation, shadedArea, totalShadedArea, plant);
                                totalShadedArea += shadedArea; //to not allow shading many times through different overgrowing individuals
                            }
                            else {
                                break;
                            }
                        }
                    }
                }
            }
        }
    }
}

/////////////////////////////////////////////////////////
double VegetationLandscape::calculateIntersectionArea(double myRadius, double neighborRadius, double distance)
{
    double x = (pow(myRadius,2) - pow(neighborRadius,2) + pow(distance,2)) / (2*distance);
    double As = pow(myRadius,2) * acos( x / myRadius );
    double Ad = x * sqrt( pow(myRadius,2) - pow(x,2) );
    double Aa = As - Ad;
    double y = (pow(neighborRadius,2)- pow(myRadius,2) + pow(distance,2)) / (2*distance);
    double Bs = pow(neighborRadius,2) * acos( y / neighborRadius );
    double Bd = y * sqrt( pow(neighborRadius,2) - pow(y,2) );
    double Ba = Bs - Bd;
    return(Aa + Ba);
}

/////////////////////////////////////////////////////////
void VegetationLandscape::calculateAPAR(double extraRadiation, double shadedArea, double totalShadedArea, int plant)
{
    /*!< \todo is this wohle function necessary. Couldnt we just use the relative shaded CA in the photosynthsis to calcualte the unshadedAPAR */
    double relativeCoverShadingPlant = shadedArea / (plantList[plant]->crownArea - totalShadedArea); //only the cover that has not been shaded already by other overgrowing individuals
    //after Sitch et al. 2000, eqn. 4-8
    double LAI = (relativeCoverShadingPlant*plantList[plant]->leafCMass * plantList[plant]->SLA) / (relativeCoverShadingPlant * plantList[plant]->crownArea);             /*!< \todo LAI stays the same as shadedArea is removed */
    double FPC = (1 - exp(-0.5 * LAI)) * relativeCoverShadingPlant * (plantList[plant]->crownArea);
    //after Schapoff et al., eqn 1, 26
    plantList[plant]->shadedAPAR += 0.5 * (extraRadiation * 2450000.0) * FPC * plantList[plant]->phenologyStatus;
}

/////////////////////////////////////////////////////////
VegetationLandscape::VegetationLandscape(const Input& _input, vector< vector<Patch*> > &_grid, Landscape* _land):
    input(_input), xCells(input.xCells), yCells(input.yCells), iniVegetation(input.iniVegetation), grid(_grid), land(_land)
{
    xLength = input.xCells * input.cellSize;
    yLength = input.xCells * input.cellSize;
    plantIDs = 0;
    distributePlants();
    newIndividuals = 0;
}

/////////////////////////////////////////////////////////
VegetationLandscape::~VegetationLandscape()
{
    for (unsigned int i = 0; i < plantList.size(); i++) {
        delete plantList[i];
        plantList[i] = NULL;
    }
}
