//!  Output.
/*!
  This class handles all the output data.
*/

#ifndef OUTPUT_H
#define OUTPUT_H

#include <fstream>   //to read and write files
#include <vector>
#include "VegetationLandscape.h"

using namespace std;

const bool paper = 1; // for specific paper output

//Forward declaration
class Input;
class WaterLandscape;
class NutrientPatch;
class Patch;

class Output
{
    public:
        /*!
         * \brief Output
         */
        Output(const Input& _input, vector< vector<Patch*> >& _grid);

        /*!
         * \brief ~Output
         */
        virtual ~Output();

        //Member functions
        /*!
         * \brief initialiseWaterFile
         */
        void initialiseWaterFile();

        /*!
         * \brief initialiseNutrientFile
         */
        void initialiseNutrientFile();

        /*!
         * \brief initialisePlantFile
         */
        void initialisePlantFile();

        /*!
         * \brief initialiseVegetationFile
         */
        void initialiseVegetationFile();

        /*!
         * \brief initialiseDebugFile
         */
        void initialiseDebugFile();

        /*!
         * \brief writeWaterFile
         */
        void writeWaterFile(int day, const WaterLandscape& water);

        /*!
         * \brief writeNutrientFile
         */
        void writeNutrientFile(int day, const NutrientPatch& nutrientPatch);

        /*!
         * \brief writePlantFile
         */
        void writePlantFile(int day, const VegetationLandscape& plants);

        /*!
         * \brief writeVegetationFile
         */
        void writeVegetationFile(int day, const VegetationLandscape& plants);

        /*!
         * \brief writeDebugFile
         */
        void writeDebugFile(int day, const NutrientPatch& nutrientPatch);

        void initialiseTestFile();
        void writeTestFile(int day, const VegetationLandscape& plants);

        void initialiseTestFileB();
        void writeTestFileB(int day);
        void initialiseWaterOutputPaper2();
        void writeWaterOutputPaper2(int day, const WaterLandscape& water);
        void initialisePlantOutputPaper2();
        void writePlantOutputPaper2(int day, const VegetationLandscape& plants);
        void initialiseNutrientOutputPaper2();
        void writeNutrientOutputPaper2(int day, const NutrientPatch& nutrientPatch);
        void initialiseVegetationOutputPaper2();
        void writeVegetationOutputPaper2(int day, const VegetationLandscape& plants);

        void writeFiles(int day, const VegetationLandscape& vegetation, const NutrientPatch& nutrientPatch, const WaterLandscape& water);

        //Member variables
        const Input& input;     /*!< Reference to Input object */
        vector< vector<Patch*> >& grid; /*!< Reference to two-dimensional landscape (grid) of the type "Patch" for each element (grid cell) */
        ofstream dailyWaterFile;    /*!< daily mean (average over landscape) water values */
        ofstream dailyWaterFileP2;    /*!< daily mean (average over landscape) water values */
        ofstream dailyNutrientFile;    /*!< daily nutrient values */
        ofstream dailyNutrientFileP2;    /*!< daily nutrient values */
        ofstream dailyPlantFile;    /*!< daily plant properties */
        ofstream dailyPlantFileP2;    /*!< daily plant properties */
        ofstream dailyVegetationFile;    /*!< daily plant properties */
        ofstream dailyVegetationFileP2;    /*!< daily plant properties */
        ofstream dailyDebugFile;    /*!< daily nutrient values */
        ofstream dailyTestFile;    /*!< daily test values */
        ofstream dailyTestFileB;    /*!< daily test values */

    protected:

    private:
};

#endif // OUTPUT_H
