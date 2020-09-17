#include "WaterLandscape.h"
#include "Landscape.h"
#include "Patch.h"
#include "Input.h"
#include <iostream>
#include "cmath"

/////////////////////////////////////////////////////////
void WaterLandscape::precipitation(int day)
{
    for(int x = 0; x < xCells; x++) {
        for(int y = 0; y < yCells; y++) {
            grid[x][y]->waterPatch->surfaceWater += input.prec[day];
        }
    }
}

/////////////////////////////////////////////////////////
void WaterLandscape::setAspectAndInclination()
{
    int x = 0; int y = 0;
    int aspect; //direction of slope inclination

    for(int i = 0; i < xCells; i++) {
        for(int j = 0; j < yCells; j++) {
            if(lowestNeighbour[i][j] != -1) { //if lowest neighour not myself

                //calculate x and y coordinates of lowest neighbouring cell
                x = (lowestNeighbour[i][j] - lowestNeighbour[i][j]%yCells) / yCells;
                y = lowestNeighbour[i][j]%yCells;

                //calculate inclination: arctangens of opposite leg / adjacent leg
                //= atan ( (elevation1-elevation2)/length ) in radian
                grid[i][j]->waterPatch->inclination = atan((grid[i][j]->waterPatch->elevation - grid[x][y]->waterPatch->elevation)/(input.cellSize*1.0));
            }
            else {  //the inclination is set equal to the inclination to the neighbour cell
                if(i==0)
                    grid[i][j]->waterPatch->inclination = atan((grid[1][j]->waterPatch->elevation - grid[0][j]->waterPatch->elevation)/(input.cellSize*1.0));
                else if (i==xCells-1)
                    grid[i][j]->waterPatch->inclination = atan((grid[xCells-2][j]->waterPatch->elevation - grid[xCells-1][j]->waterPatch->elevation)/(input.cellSize*1.0));
                else if (j==0)
                    grid[i][j]->waterPatch->inclination= atan((grid[i][1]->waterPatch->elevation - grid[i][0]->waterPatch->elevation)/(input.cellSize*1.0));
                else if (j==yCells-1)
                    grid[i][j]->waterPatch->inclination = atan((grid[i][yCells-2]->waterPatch->elevation - grid[i][yCells-1]->waterPatch->elevation)/(input.cellSize*1.0));
                else
                    grid[i][j]->waterPatch->inclination = 0;
                if((i==0) && (j==0))
                    grid[i][j]->waterPatch->inclination = atan((grid[1][1]->waterPatch->elevation - grid[i][j]->waterPatch->elevation)/(input.cellSize*1.0));
                if((i==0) && (j==yCells-1))
                    grid[i][j]->waterPatch->inclination = atan((grid[1][yCells-2]->waterPatch->elevation -grid[i][j]->waterPatch->elevation)/(input.cellSize*1.0));
                if((i==xCells-1) && (j==0))
                    grid[i][j]->waterPatch->inclination = atan((grid[xCells-2][1]->waterPatch->elevation - grid[i][j]->waterPatch->elevation)/(input.cellSize*1.0));
                if((i==xCells-1) && (j==yCells-1))
                    grid[i][j]->waterPatch->inclination = atan((grid[xCells-2][yCells-2]->waterPatch->elevation - grid[i][j]->waterPatch->elevation)/(input.cellSize*1.0));
            }

            //if the inclinination degree is small, the cell is assumed to be flat
            if ( (grid[i][j]->waterPatch->inclination * 180.0/M_PI) < 5 ) aspect = FL;

            //it is assumed, that the upper boundary of the MAP is oriented north
            else {
                if (i-x == -1){
                    if (j-y == -1) aspect = NW;
                    if (j-y ==  0) aspect = NN;
                    if (j-y ==  1) aspect = NE;
                }
                else if (i-x == 0) {
                    if (j-y == -1) aspect = WW;
                    if (j-y ==  1) aspect = EE;
                }
                else {
                    if (j-y == -1) aspect = SW;
                    if (j-y ==  0) aspect = SS;
                    if (j-y ==  1) aspect = SE;
                }
            }

            //the aspect factor determines how the radiation is affected by the aspect of the slope
            //the values for this factor were found in
            //Shevenell L., 1999. Regional potential evapotranspiration in arid climates based on temperature, topography and calculated solar radiation. Hydrological Processes 13, 577-596
            switch (aspect) {
                    case FL: grid[i][j]->waterPatch->aspectFactor = 1.0; break;
                    case NN: grid[i][j]->waterPatch->aspectFactor = 0.90; break;
                    case NE: grid[i][j]->waterPatch->aspectFactor = 0.95; break;
                    case EE: grid[i][j]->waterPatch->aspectFactor = 0.98; break;
                    case SE: grid[i][j]->waterPatch->aspectFactor = 1.03; break;
                    case SS: grid[i][j]->waterPatch->aspectFactor = 1.10; break;
                    case SW: grid[i][j]->waterPatch->aspectFactor = 1.05; break;
                    case WW: grid[i][j]->waterPatch->aspectFactor = 1.02; break;
                    case NW: grid[i][j]->waterPatch->aspectFactor = 0.97; break;
            }
        }
    }
}

/////////////////////////////////////////////////////////
void WaterLandscape::calculateLowestNeighbor()
{
    int x_lower;    //help variable to store x-coordinate of so far lower cell
    int y_lower;    //help variable to store y-coordinate of so far lower cell

    for(int i = 0; i < xCells; i++) {
        for(int j = 0; j < yCells; j++) {
            x_lower = i;
            y_lower = j;

            //look at surrounding 8 cells (plus itself)
            for(int k = max(0,i-1); k <= min(i+1,xCells-1); k++) {
                for(int l = max(0,j-1); l <= min(j+1,yCells-1); l++) {
                    if(grid[k][l]->waterPatch->elevation < grid[x_lower][y_lower]->waterPatch->elevation) {
                        x_lower = k;
                        y_lower = l;
                    }
                }
            }
            //direction in which the runoff of a cell flows
            //integer number indicates absolute number of cell,
            //e.g. cell[x][y]->yCells*x + y = lowestNeighbour
            //		y = lowestNeighbour%yCells
            //		x = (lowestNeighbour - y)/yCells
            if(x_lower == i && y_lower == j)
                lowestNeighbour[i][j] = -1; //this cell does not have a lowest neighbor
            else
                lowestNeighbour[i][j] = x_lower*yCells + y_lower;
        }
    }
}

/////////////////////////////////////////////////////////
void WaterLandscape::infiltrationAndDrainage()
{
   //// \todo first call functions according to EcoHyD for testing, later: iterate through layers so that it makes sense
   for(int x = 0; x < xCells; x++) {
        for(int y = 0; y < yCells; y++) {
            if(grid[x][y]->waterPatch->surfaceWater > 0.0) grid[x][y]->waterPatch->fastInfiltration(1);
            if(grid[x][y]->waterPatch->surfaceWater > 0.0) grid[x][y]->waterPatch->slowInfiltration(0);
            grid[x][y]->waterPatch->drainage(0);
            grid[x][y]->waterPatch->drainage(1);
        }
    }
}

/////////////////////////////////////////////////////////
void WaterLandscape::diffusion()
{
    for(int x = 0; x < xCells; x++) {
        for(int y = 0; y < yCells; y++) {
            for(int layer = 0; layer < input.nSoilLayers; layer++) {
            grid[x][y]->waterPatch->diffusion(layer);
            }
        }
    }
}

/////////////////////////////////////////////////////////
void WaterLandscape::evaporation(int day)
{
    //potential evaporation [mm/day] dependent on temperature, extraterrestrial radiation, eqn. 10
    double potEP = 0.0023 * (input.tempMean[day] + 17.8) * pow((input.tempMax[day] - input.tempMin[day]), 0.5) * land->extraterrestrialRadiation(day);

    for(int x = 0; x < xCells; x++) {
        for(int y = 0; y < yCells; y++) {
            grid[x][y]->waterPatch->surfaceEvaporation(potEP);
            grid[x][y]->waterPatch->soilEvaporation(potEP);
            /* not necessary anymore as transpiraion is caluclated in the plant module (transpiration) */
            /*for(int layer = 0; layer < input.nSoilLayers; layer++) {
                grid[x][y]->waterPatch->layerEvapotranspiration(layer, potEP);
            }*/
        }
    }
}

/////////////////////////////////////////////////////////
void WaterLandscape::followRunoff(int i, int j, int addrunon, bool first)
{
        //add the additional runon to the cell
        runonCells[i][j] += addrunon;
        numberOfFlows[i][j] += 1;

        //stop condition: no runoff to other cells
        if (lowestNeighbour[i][j] == -1) return;
        //stop condition: the routine was called from outside with a new cell, but the cell has already distributed its runoff
        if (first && calculated[i][j]) return;

        //if the runoff of this cell was not calculated yet, the runoff of the cell itself is added
        if (!calculated[i][j]) {
                addrunon = runonCells[i][j] + 1; //// \todo why "+ 1" ?
                calculated[i][j] = true;
        }

        //call the method again for the cell where the runoff is distributed to
        followRunoff((lowestNeighbour[i][j]-lowestNeighbour[i][j]%yCells)/yCells, lowestNeighbour[i][j]%yCells, addrunon, false);
}

/////////////////////////////////////////////////////////
void WaterLandscape::followShortRunoffPath(int xPos, int yPos)
{
    helpCell = {xPos, yPos, numberOfFlows[xPos][yPos], -1};
    cellList2.push_back(helpCell);

    calculated[xPos][yPos] = true;
    if (lowestNeighbour[xPos][yPos] != -1
        && numberOfFlows[(lowestNeighbour[xPos][yPos]-lowestNeighbour[xPos][yPos]%yCells)/yCells][lowestNeighbour[xPos][yPos]%yCells] == numberOfFlows[xPos][yPos]
        && calculated[(lowestNeighbour[xPos][yPos]-lowestNeighbour[xPos][yPos]%yCells)/yCells][lowestNeighbour[xPos][yPos]%yCells] == false) {
        followShortRunoffPath((lowestNeighbour[xPos][yPos]-lowestNeighbour[xPos][yPos]%yCells)/yCells,lowestNeighbour[xPos][yPos]%yCells);
    }
}

/////////////////////////////////////////////////////////
void WaterLandscape::calculateRunonCells()
{
    //calculate number of cells that produce runon to a cell
    for(int i = 0; i < xCells; i++) {
        for(int j = 0; j < yCells; j++) {
            followRunoff(i, j, 0, true);
        }
    }

    int helpFlow;
    //make new list in order of number of flows that contribute to the runon
    for(int i = 0; i < xCells; i++) {
        for(int j = 0; j < yCells; j++) {
            helpFlow = 0;
            for(int k = max(0,i-1); k <= min(i+1,xCells-1); k++) {
                for(int l = max(0,j-1); l <= min(j+1,yCells-1); l++) {
                    if(k != i && l != j && ((lowestNeighbour[k][l]-lowestNeighbour[k][l]%yCells)/yCells == i) && (lowestNeighbour[k][l]%yCells == j)) {
                        helpFlow = numberOfFlows[k][l];
                    }
                }
            }
            helpCell = {i , j, numberOfFlows[i][j], helpFlow};
            cellList1.push_back(helpCell);
        }
    }

    for(int i = 0; i < xCells; i++) {
        for(int j = 0; j < yCells; j++) {
            calculated[i][j] = false;
        }
    }

    for(int i = 0; i < cellList1.size(); i++) {
        //// \todo Britta starts from second cell?
        //if the runon to this cell is not calulated yet
        if (calculated[cellList1[i].xCor][cellList1[i].yCor] == false){
                followShortRunoffPath(cellList1[i].xCor, cellList1[i].yCor);
        }
    }
}

/////////////////////////////////////////////////////////
void WaterLandscape::runoff()
{
    //the following route follows the order of the cells in the cellList,
    //first the runon to cells with only one flow contributing to the runon are calculated,
    //afterwards the ones with two flows etc
    //this guaranties, that cells are not double counted and that water can already infiltrate into upslope cells before it reaches downslope cells

    int xpos, ypos;

    //follow sorted list and calculate the infiltrated water in the cell and the water that is passed to the next cell
    //// \todo Britta starts from second cell?
    for(int i = 0; i < cellList2.size(); i++) {
        xpos = cellList2[i].xCor;
        ypos = cellList2[i].yCor;
        runon[xpos][ypos] = 0;
        //it is only necessary to evaluate the runoff that reaches a cell, if there are cells, that distribute runoff to the actual cell
        if(runonCells[xpos][ypos] != 0) {
            //look at the surrounding cells and calculate the total runon
            for(int k = max(0,xpos-1); k <= min(xpos+1,xCells-1); k++) {
                for(int l = max(0,ypos-1); l <= min(ypos+1,yCells-1); l++) {
                    //if the cell provides runoff to this cell and it is not the cell itself and it has runoff
                    if ((lowestNeighbour[k][l] == xpos*yCells + ypos) && !(k==xpos && l==ypos) && (grid[k][l]->waterPatch->totalRunOff > 0)) {
                        //does this make sense?
                        runon[xpos][ypos] += grid[k][l]->waterPatch->totalRunOff;
                    }
                }
            }
            grid[xpos][ypos]->waterPatch->surfaceWater += runon[xpos][ypos];
        }

        //Routine according to Manning-Strickler
        //eqn. 15
        grid[xpos][ypos]->waterPatch->totalRunOff = pow(grid[xpos][ypos]->waterPatch->surfaceWater, 2.0/3.0) * (1 - (grid[xpos][ypos]->vegetationPatch->relVegCover)/2.0) * sqrt(grid[xpos][ypos]->waterPatch->inclination);
        grid[xpos][ypos]->waterPatch->totalRunOff = max(0.0, min(grid[xpos][ypos]->waterPatch->surfaceWater, grid[xpos][ypos]->waterPatch->totalRunOff));
        grid[xpos][ypos]->waterPatch->surfaceWater -= grid[xpos][ypos]->waterPatch->totalRunOff;
    }
}

/////////////////////////////////////////////////////////
void WaterLandscape::meanLandscape()
{
    meanAbsSurfaceWater = 0;
    meanDeepDrainedWater = 0;
    for(int layer = 0; layer < input.nSoilLayers; layer++) {
        meanAbsWater[layer] = 0;
        meanRelWater[layer] = 0;
    }

    for(int x = 0; x < xCells; x++) {
         for(int y = 0; y < yCells; y++) {
             meanAbsSurfaceWater += grid[x][y]->waterPatch->surfaceWater;
             meanDeepDrainedWater += grid[x][y]->waterPatch->deepDrainedWater;
             for(int layer = 0; layer < input.nSoilLayers; layer++) {
                 meanAbsWater[layer] += grid[x][y]->waterPatch->waterLayer[layer].absWater;
                 meanRelWater[layer] += grid[x][y]->waterPatch->waterLayer[layer].relWater;
             }
         }
     }

    meanAbsSurfaceWater /= xCells * yCells;
    meanDeepDrainedWater /= xCells * yCells;
    for(int layer = 0; layer < input.nSoilLayers; layer++) {
        meanAbsWater[layer] /= xCells * yCells;
        meanRelWater[layer] /= xCells * yCells;
    }
}

/////////////////////////////////////////////////////////
WaterLandscape::WaterLandscape(const Input& _input, vector< vector<Patch*> >& _grid, Landscape* _land):
    xCells(_input.xCells), yCells(_input.yCells), input(_input), grid(_grid), land(_land)
{
    meanAbsWater.resize(input.nSoilLayers, 0);
    meanRelWater.resize(input.nSoilLayers, 0);
    lowestNeighbour.resize(xCells, vector<int>(yCells, 0));
    runonCells.resize(xCells, vector<int>(yCells, 0));
    calculated.resize(xCells, vector<bool>(yCells, 0));
    numberOfFlows.resize(xCells, vector<int>(yCells, 0));
    runon.resize(xCells, vector<double>(yCells, 0));

    for(int x = 0; x < xCells; x++) {
        for(int y = 0; y < yCells; y++) {
            grid[x][y]->waterPatch->surfaceWater = 0;
        }
    }

    setAspectAndInclination();
    calculateLowestNeighbor();
    calculateRunonCells();
}

/////////////////////////////////////////////////////////
WaterLandscape::~WaterLandscape()
{

}
