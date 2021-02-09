#ifndef INPUT_H
#define INPUT_H

#include <vector>
#include <string>
#include <chrono>

using namespace std;

const bool exe = 0; // set 1 when creating an exe because paths change in QtCreator

//To define delimiter depending on operating system
const string delimiter =
    #ifdef _WIN32
        "\\";        //for windows
    #else
        "/";         //for linux
    #endif

class Input
{
    public:
        //! A constructor.
        /*!
          A more elaborate description of the constructor.
        */
        Input();

        //! A destructor.
        /*!
          A more elaborate description of the destructor.
        */
        virtual ~Input();

        //Member functions
        //! Reads in the general model parameters.
        void readMaster();

        //! Calculates the number of simulation days.
        void calculateSimDays();

        //! Reads in the weather timeseries data.
        void readWeather();

        //! Calculates mean annual temperature.
        void calculateMAT();

        //! Reads in the elevation of the landscape per patch (grid cell).
        void readTopography();

        //! Reads in the soil type of the landscape per patch.
        void readSoilMap();

        //! Reads in the soil type database including all the soil properties per soil type.
        void readSoilProperties();

        //! Reads in the intial vegetation.
        void readIniVegetation();

        //! Reads in plant traits per PFT.
        void readPlantTraits();

        //Member variables
        ////General model parameters
        int xCells;  /*!< Number of patches in x direction (i.e. the rows) */
        int yCells;  /*!< Number of patches in y direction (i.e. the columns)  */
        int cellSize;  /*!< Size of each patch's edge in m */
        int nSoilLayers;  /*!< Number of soil layers per patch */
        double latitude;    /*!< latitude */
        string simStartDate; /*!< First simulation date */
        string simEndDate; /*!< Last simulation date */
        vector<int> dayOfYear; /*!< The day of the year (vector length: simulations days) */
        vector<int> yearOfDay; /*!< The year for each day (vector length: simulations days) */
        int simDays; /*!< Number of simulation days */
        vector<int> depthLayers;  /*!< Depth of each soil layer per patch */
        double iniSurfaceWater;  /*!< Inital surface water in mm */
        vector<double> iniWaterLayer;  /*!< Inital water content of the layers in vol% */
        bool measuredSolRadiation; /*!< Is there measured soil radiation available to read-in? */
        double iniNH4;  /*!< Initial ammonium per kg soil [mg/kg] */
        double iniNO3;  /*!< Initial nitrate per kg soil [mg/kg] */
        double dailyNdeposition; /*!< Daily nitrogen deposition (sum of NO3 and NH4) [kg/ha] */

        ////Weather data
        vector<string> date;  /*!< Date of the weather data */
        vector<double> prec;  /*!< Daily precipitation [mm] */
        vector<double> tempMean;  /*!< Mean daily temperature [°C] */
        vector<double> tempMin;  /*!< Daily minimum temperature [°C] */
        vector<double> tempMax;  /*!< Daily maximum temperature [°C] */
        vector<double> ambientCO2; /*!< Daily ambient partial pressure of CO2 [ppmv] */
        vector<double> solRadiation;  /*!< Daily solar radiation [W*m-2] */
        vector<double> MAT; /*!< mean annual temperature [degC] */

        ////Elevation
        vector< vector<double> > elevation;  /*!< Elevation of the landscape per patch */

        ////Vegetation
        int PFTs; /*!< number of plant functional types used in this approach */
        bool spatialVegInput;   /*!< Is there spatial vegetation input available? */
        int nPFTs;  /*!< number of plant functional types the landscape should be initiated with */
        vector<string> pftNames;    /*!< List of all PFT names */
        vector<int> nIndividuals;   /*!< Number of individuals per PFT */
        int totalIndividuals;   /*!< Total initial plant individuals */
        vector<double> meanHeight;   /*!< Mean height per PFT */
        vector<double> sdHeight;   /*!< Standard deviation of the height per PFT */
        vector< vector<string> > iniVegetation; /*!< Plants within the landscape (coordinate and initial properties of the plant) */
        vector< vector<string> > plantTraits; /*!< Plant traits per PFT */

        ////Soil properties
        vector< vector<string> > soilType;  /*!< Soil type of the landscape per patch */
        vector< vector<string> > soilProperties;  /*!< A database of all soil types and their properties */

        //for getting a random seed depending on the time
        typedef chrono::high_resolution_clock myclock;
        myclock::time_point begin;
};

#endif // INPUT_H
