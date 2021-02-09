//!  Plant.
/*!
  This class describes all plant processes and parameters.
*/

#ifndef PLANT_H
#define PLANT_H

#include <string>
#include <vector>
#include <chrono>

using namespace std;

//Forward declaration
class Input;
class Landscape;

const double cMass = 12.0107; /*!< Atomic mass of carbon  */
const double atmPa = 1e5;	/*!< Atmospheric pressure [Pa] */
const double ambientO2Pa = 20900.0; /*!< Partial pressure of O2 in [Pa] */
const double CO2ConvFactor = 1.0e-6;  /*!< Conversion factor for CO2 from ppmv to mole fraction */
const double theta = 0.7; /*!< Co-limitation (shape) parameter */
const double alpha = 1.1; /*!< Maximum Priest-Taylor coefficient = 1.391 in Schaphoff et al. 2017 (changed to 1.1 after Monteith 1995) */
const double scalingFactor = 5.0; /*!< Conductance scaling factor from Sitch et al., 2003 */

class Plant
{
    public:

        Plant(int index, const Input& _input, Landscape* _land);
        Plant(int _myID, string _PFTname, double _xCor, double _yCor, const Input& _input, Landscape* _land);

        virtual ~Plant();

        void calculateIniPlantProperties();
        void updatePlantProperties();
        void tissueTurnover();
        void phenology(int day);
        void photosynthesis(int day, double dayLength, double lambda);
        void vmax(double c1, double c2, double dayLength, int day);
        void transpiration(int day, double dayLength);
        void allocation();
        double allocationFunction(double leafMassInc);
        void respiration(int day);
        void reproduction();
        void mortality();
        void inCaseOfDeath();
        void nitrogenUptake(int day, double dayLength);
        void extractNitrogenFromSoil(double nitrogenAmount);

        /*!
         * \brief maximum capacity of nitrogen storage [kg?]
         */
        double maxNitrogenStorage();

        /*!
         * \brief mean water saturation in rooting zones per soil depth
         */
        double waterSaturationPerDepth();
        void availablePools();
        void dispersal();
        void calculatewScal(int day, double phen, double potLeafCMass, double potFPC);

        //Member variables
        int xPatch; /*!< x-patch of the plant [-] */
        int yPatch; /*!< y-patch of the plant [-] */

        //Public member variables
        int myID; /*!< ID of the individual */
        double xCor;     /*!< x-coordinate of the plant [m] */
        double yCor;     /*!< y-coordinate of the plant [m] */

        ///Plant traits
        string PFTname;    /*!< name of the PFT */
        double optLambda; /*!< optimal ci/ca for non-water stressed conditions [-] */
        vector<double> relRootsPerLayer; /*!< relative roots for each soil layer density [-] */
        double seedMass;    /*!< seed mass [mg] */
        double germinationProb;    /*!< germination prob of the seeds [-] */
        double WP;  /*!< wilting point [cm of water] */
        double SLA; /*!< specific leaf area [m2 * kgC-1] */

        ///State variables
        double FPC; /*!< the area of ground covered by foliage above it [m2] */
        double GPP; /*!< daily gross primary productivity [kgC*day-1] */
        double NPP; /*!< daily net primary productivity [kgC*day-1] */
        double leafCMass;    /*!< leaf carbon mass [kg] */
        double rootCMass;    /*!< root carbon mass [kg] */
        double sapwoodCMass;     /*!< sapwood carbon mass [kg] */
        double heartwoodCMass;   /*!< heartwood carbon mass [kg] */
        double totalCMass;   /*!< total plant's carbon mass [kg] */
        double dailyCInc;   /*!< daily plant's carbon mass increment [kg] */
        double leafNMass;    /*!< leaf nitrogen mass [kg] */
        double rootNMass;    /*!< root nitrogen mass [kg] */
        double sapwoodNMass;     /*!< sapwood nitrogen mass [kg] */
        double reproductiveMass;    /*!< reproductive mass [kg] */
        double height;      /*!< plant height [m] */
        double crownRadius; /*!< crown radius [m] */
        double crownArea;   /*!< crown area [m2] */
        double abovegroundCLitter;   /*!< dead carbon litter aboveground [kg] */
        double belowgroundCLitter;   /*!< dead carbon litter belowground [kg] */
        double abovegroundNLitter;   /*!< dead nitrogen litter aboveground [kg] */
        double belowgroundNLitter;   /*!< dead nitrogen litter belowground [kg] */
        bool dead; /*!< information if the plant died in this year or not [no:0 or yes:1] */
        double nScal; /*!< nutrient limitation factor for allocation and rubisco activity [-] \todo Differently calculated as suggested in Smith et al. 2014 */
        double shadedAPAR; /*!< APAR extinct through shading [W*m-2 ??] */
        double phenologyStatus; /*!< leaf phenology status [-] */
        double avNitrogen; /*!< nitrogen availability for the plant [kg] */
        int age; /*!< age of the plant since planting or establishment [days] */
        double waterSupply; /*!< current water supply [mm/day?] */
        double waterDemand; /*!< current water demand [mm/day?] */
        double wScal; /*!< water limitation factor for allocation and phenology [-] */
        double avRelSoilMoistureFC; /*!< relative soil moisture availability as a fraction of field capacity [1] and wilting point [0] for the plant [-] */
        double avAbsSoilMoisture; /*!< absolute soil moisture availability for the plant [-] */

        vector<pair<pair<int, int>, double> > intersectedPatches; /*!< patches that are intersected by this plant with first.first = x-coordinate of the patch [-], first.second = y-coordinate of the patch [-], second = relative cover of the plant [-] */

        struct Seed {
            double distance; /*!< distance of this soon to be dispersed seed from the mother plant [cm] */
            double direction; /*!< direction of this soon to be dispersed seed from the mother plant [rad] */
        };
        vector<Seed> seedList; /*!< seeds that should be dispersed in the landscape [-] */

        ///OUTPUT ONLY
        double transpiredWater; /*!< absolute transpired water for the plant [mm] */
        double uptakenNitrogen; /*!< uptaken nitrogen for the plant [kg] */
        double carbonLitter; /*!< total carbon litter [kg] */
        double nitrogenLitter; /*!< total nitrogen litter [kg] */
        double actLeafRespiration; /*!< leaf respiration [kg] */
        double totalRespiration; /*!< total respiration [kg] */
        vector<double> allocationTest; /*!< test vector for allocation output values */
        void resetOutputVariables();
        double relCover;
        double potLeafCMass;
        double potentialLAI;
        double potentialFPC;
        double wScal_phen;

     private:
        //Member variables
        const Input& input;    /*!< Reference to Input object */
        Landscape* land;    /*!< Pointer to Landscape object */

        ///Plant traits
        string phenologyType;   /*!< phenology type [-] */
        int GDD5;   /*!< number of growing degree days on a 5Â°C base [number of days] */
        double laToSa;   /*!< constant for tranforming leaf area to sapwood area [-] */
        double lmToRm;    /*!< constant for tranforming leaf mass to root mass [-] */
        double allom2;  /*!< first allometry parameter for height diameter relationship [-] */
        double allom3;  /*!< second allometry parameter for height diameter relationship [-]  */
        double allom1;  /*!< allometry parameter for self thinning constraint [-] */
        double reinickeRp;  /*!< parameter for self thinning constraint [-] */
        double maxCrownArea;    /*!< maximum crown area [m2] */
        string photosyntheticPathway;   /*!< photosyntheic pathway */
        double CO2UptakeEfficiency; /*!< intrinsic qunatum efficiency of CO2 uptake */
        double temp1;    /*!< temperature limit 1 [Â°C] */
        double temp2;    /*!< temperature limit 2 [Â°C] */
        double temp3;    /*!< temperature limit 3 [Â°C] */
        double temp4;    /*!< temperature limit 4 [Â°C] */
        double gMin;    /*!< minimum canopy conductance [mm/s] */
        double eMax;    /*!< maximum water transport capacity [-] */
        double b;  /*!< leaf respiration as fraction of Vmax [-]*/
        double respirationCoeff; /*!< respiration coefficient [-] */
        double fracGResp; /*!< fractional growth respiration [-] */
        double carbToNitroLeaves;   /*!< carbon to nitrogen ratio leaves [-] */
        double carbToNitroSapwood;   /*!< carbon to nitrogen ratio sapwood [-] */
        double carbToNitroRoots;   /*!< carbon to nitrogen ratio roots [-] */
        double minLeafCN;   /*!< minimum carbon to nitrogen ratio in the leaves [-] */
        double maxLeafCN;   /*!< maximum carbon to nitrogen ratio in the leaves [-] */
        double k;   /*!< constant in the maximum N storage equation [-] */
        double kmax; /*!< half-saturation concentration for N uptake [kgN*m-3] */
        double nUptakePerRootC; /*!< maximum N uptake per unit fine root biomass [kgN*kgC-1*day-1]  */
        double leafTurnoverTime;    /*!< leaf turnover time [years] */
        double rootTurnoverTime;    /*!< root turnover time [years] */
        double sapwoodTurnoverTime;    /*!< sapwood turnover time [years] */
        double meanDispersalDistance;   /*!< mean dispersal distance [cm] */
        double sdDispersalDistance;   /*!< standard deviation of dispersal distance [cm] */
        double woodDensity; /*!< wood density [kgC * m-3] */
        double mortThreshold; /*!< mortality threshold value [-] */

        ///State variables
        int currentGDD5;    /*!< number of growing degree days on a 5deg°C base for this year [number of days] */
        double sapwoodArea;     /*!< sapwood cross sectional area [m2] */
        double stemDiameter;     /*!< stem diameter [m] */
        double leafArea;    /*!< average leaf area [m2] */
        double LAI; /*!< leaf area index [-] */
        double APAR;    /*!< absorbed photosynthetically active ratiation [W*m-2 ??] */
        double gPot; /*!< potential canopy conductance [mm/s ??] */
        double gAct; /*!< actual canopy conductance [mm/s ??] */
        double potEP; /*!< potential evaporation [mm/day] */
        double debtMass; /*!< current mass for future carbon debts [kg] */
        double leafCMassPreviousYear; /*!< leaf carbon mass of the last year [kg] */
        double rootCMassPreviousYear; /*!< root carbon mass of the last year [kg] */
        double sapwoodCMassPreviousYear; /*!< sapwood carbon mass of the last year [kg] */
        double totalCMassPreviousYear; /*!< total carbon mass of the last year [kg] */
        double ANPP;  /*!< annual net primary productivity [kg] */
        double storageNMass; /*!< nitrogen storage mass [kg] */
        double NPPyesterday; /*!< yesterdays net primary productivity [kgC*day-1] */
        double nitrogenTurnover; /*!< amount of N to be reallocated form turnover of leaves, fine roots and sapwood [kg] */
        double lmToRmScal; /*!< leaf mass to root mass affected by most limitating stress factor (either by nScal or wScal) [-] */
        double avAbsSoilMoistureNew; /*!< absolute soil moisture minus wilting point [mm] */
        double avRelSoilMoisture; /*!< relative soil moisture availability for the plant [-] */
        double avSoilTemperature;  /*!< soil tenperature in the rooting zone of the plant [degC] */
        double Vmax; /*!< non-water-stressed rubisco capacity [?] */
        double NPPinMM; /*!< net daytime photosynthesis [mm*m2-1*day-1] */

        //for getting a random seed depending on the time
        typedef chrono::high_resolution_clock myclock;
        myclock::time_point beginning;

    protected:
};

#endif // PLANT_H
