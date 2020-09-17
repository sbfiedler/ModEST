//!  WaterLandscape.
/*!
  This class handles all water processes and parameters relevant on the landscape scale.
*/

#ifndef WATERLANDSCAPE_H
#define WATERLANDSCAPE_H

#include <vector>

using namespace std;

//Forward declaration
class Input;
class Patch;
class Landscape;

//! An enum.
/*! Constants to describe aspect. */
enum Aspects
{
    FL = 0, //!< enum value 1 for flat, no aspect = 0
    NN,     //!< enum value 2 for north = 1
    NE,     //!< enum value 3 for north-east = 2
    EE,     //!< enum value 4 for east = 3
    SE,     //!< enum value 5 for south-east = 4
    SS,     //!< enum value 6 for south = 5
    SW,     //!< enum value 7 for south-west = 6
    WW,     //!< enum value 8 for west = 7
    NW      //!< enum value 9 for north-west = 8
};

class WaterLandscape
{
    public:
        /*!
         * \brief A constructor.
         * \param _input pointer to the Input object
         * \param _grid reference to the gird vector
         * \param _land pointer to the Landscape object
         */
        WaterLandscape(const Input& _input, vector< vector<Patch*> >& _grid, Landscape* _land);

        //! A destructor.
        virtual ~WaterLandscape();

        //Member functions
        /*!
         * \brief Adds actual precipitation to surface water per patch.
         * \param day current day
         */
        void precipitation(int day);

        /*!
         * \brief Infiltration and drainage per patch.
         */
        void infiltrationAndDrainage();

        /*!
         * \brief diffusion
         */
        void diffusion();

        /*!
         * \brief evaporation
         * \param day current day
         */
        void evaporation(int day);

        /*!
         * \brief runoff
         */
        void runoff();

        /*!
         * \brief meanLandscape
         */
        void meanLandscape();

        //Member variables
        const int& xCells;   /*!< Number of x-cells */
        const int& yCells;   /*!< Number of y-cells */
        const Input& input;     /*!< Reference to Input object */
        vector< vector<Patch*> >& grid; /*!< Reference to two-dimensional landscape (grid) of the type "Patch" for each element (grid cell) */
        Landscape* land;     /*!< Pointer to Landscape object */

        double meanAbsSurfaceWater; /*!< Mean absolute surface water [mm] */
        double meanDeepDrainedWater;    /*!< Mean absolute deep drained water [mm] */
        vector<double> meanAbsWater;    /*!< Mean absolute water per layer [mm] */
        vector<double> meanRelWater;    /*!< Mean relative water per layer [vol%] */

    private:
        //Member functions
        /*!    
         * \brief setAspectAndInclination
         */
        void setAspectAndInclination();

        /*!
         * \brief calculateLowestNeighbor
         */
        void calculateLowestNeighbor();

        /*!
         * \brief followRunoff
         * \param i
         * \param j
         * \param addrunon
         * \param first
         */
        void followRunoff(int i, int j, int addrunon, bool first);

        /*!
         * \brief followShortRunoffPath
         * \param xPos
         * \param yPos
         */
        void followShortRunoffPath(int xPos, int yPos);

        /*!
         * \brief runonCells
         */
        void calculateRunonCells();

        //Member variables
        vector< vector<int> > lowestNeighbour;  /*!< which is the lowest neighbor (contents the number of the target cell) */
        vector< vector<int> > runonCells;   /*!< number of patches that provide runon */
        vector< vector<bool> > calculated;  /*!< indicates if the water of this cell has already been distributed */
        vector< vector<int> > numberOfFlows; /*!< insert description here */
        vector< vector<double> > runon; /*!< runon of each patch */

        //// \todo is this needed? Don't see where myRunonCell and myPrevRunonCells will be used at all?
        struct Cell {
            int xCor,    /*!< x-coordiante of the cell */
            yCor,     /*!< y-coordiante of the cell*/
            myRunonCells,     /*!< numberOfFlows: insert description here */
            myPrevRunonCells;    /*!< numberOfFlows: insert description here */
        };
        Cell helpCell;
        vector<Cell> cellList1; /*!< list of cells */
        vector<Cell> cellList2; /*!< list of cells */

};

#endif // WATERLANDSCAPE_H
