//!  NutrientLandscape.
/*!
  This class handles all nutrient processes and parameters relevant on the landscape scale.
*/

#ifndef NUTRIENTLANDSCAPE_H
#define NUTRIENTLANDSCAPE_H

#include <vector>

using namespace std;

//Forward declaration
class Input;
class Patch;

class NutrientLandscape
{
    public:
        /*!
         * \brief A constructor.
         * \param input pointer to the Input object
         */
        NutrientLandscape(const Input& _input, vector< vector<Patch*> > &_grid);

        //! A destructor.
        virtual ~NutrientLandscape();


        //Member variables
        const int& xCells;   /*!< Number of x-cells */
        const int& yCells;   /*!< Number of y-cells */
        const Input& input; /*!< Reference to Input object */
        vector< vector<Patch*> >& grid; /*!< Reference to two-dimensional landscape (grid) of the type "Patch" for each element (grid cell) */
        //double potSolRadiation; /*!< Extra-terrestrial solar radiation [MJ*m-2] */
        //double actSolRad; /*!< Solar radiation reaching the soil [MJ*m-2] */


        //Member functions
        void temperatureSurf(int day);
        void temperatureGrid(int day);
        void residualIncomeGrid();
        void setFactorsNcycleGrid();
        void denitrificationGrid();
        void decompositionGrid();
        void newCNsomGrid();
        void ResHumificGrid();
        void netMineralizationGrid();
        void checkNavailGrid();
        void nitrificationVolatilizationGrid();
        void leachingNO3Grid();
        void processesToPoolsGrid();
        //void potSolarRadiation(int day);
        //void actSolarRadiation(int year, int day, double potSolRadiation);

    private:

};

#endif // NUTRIENTLANDSCAPE_H
