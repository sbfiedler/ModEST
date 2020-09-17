//!  NutrientPatch.
/*!
  This class handles all nutrient processes and parameters relevant on the patch scale.
*/

#ifndef NUTRIENTPATCH_H
#define NUTRIENTPATCH_H

#include <string>
#include <vector>
#include <cmath>

using namespace std;

//Forward declaration
class Input;
class Patch;

class NutrientPatch
{
    public:
        /*!
         * \brief A constructor.
         * \param input pointer to the Input object
         */
        NutrientPatch(int index, int _xCor, int _yCor, const Input& input, const Patch& _patch);
        virtual ~NutrientPatch();

        //Member variables
        double TsoilSurfBare; /*!< soil surface temperature for bare soil [degree C] */
        double TsoilSurf; /*!< soil surface temperature [degree C] */
        double TsoilCoverF; /*!< soil surface temperature vegetation cover factor [-] */
        double TsoilLag; /*!< soil temperature lag */
        double dampingDepthMaximum; /*!< maximum damping depth [mm] */
        double totalDepth; /*!< total depth of the patch [mm] */
        double mineralizationRateCoef; /*!< mineralization rate coefficient */
        double turnoverRate; /*!< SOM turnover rate */ // 0.055 as default in SWAT, cfdec in hru
        double dailyNO3deposition; /*!< daily NO3 deposition [kg ha-1] */ // in kg ha-1
        double dailyNH4deposition; /*!< daily NH4 deposition [kg ha-1] */ // in kg ha-1

        struct NutrientLayer {
            int myLayer; /*!< the layer number */
            double thicknessLayer; /*!< thickness of the layer [mm] */
            double meanDepthLayer; /*!< depth at the vertical center of the layer [mm] */
            double maxDepthLayer; /*!< depth at the bottom of the layer [mm] */
            double soilTemp; /*!< soil temperature for the vertical depth of the layer [degree C] */
            double soilDepthTF; /*!< soil depth temperature factor */
            double soilWaterScalingFactor;  /*!< scaling factor for soil water */
            double dampingDepth; /*!< damping depth [m] */
            double dampingDepthRatio; /*!< ratio of the depth at the center of the soil layer to the damping depth */

            double hectarWeight; /*!< soil layer weight scaled to 1 hectar [kg] */
            double weight; /*!< soil layer weight [kg] */
            double fc; /*!< field capacity [mm] */
            double soilPorosity; /*!< soil porosity (fraction) */
            double pClay; /*!< percentage of clay [%] */

            double WPcalc; /*!< wilting point (fraction) */
            double AWC; /*!< available water content (fraction) */
            double awc; /*!< available water content [mm] */
            double wp; /*!< wilting point [mm] */
            double voidSoil; /*!< soil void (fraction) */

            double Cres; /*!< residue C content [kg ha-1] */
            double Nres; /*!< residue N content [kg ha-1] */
            double CNres; /*!< residue C:N ratio */
            double CresNew; /*!< C content in newly formed residue [kg ha-1] */
            double NresNew; /*!< N content in newly formed residue [kg ha-1] */
            double CNresNew; /*!< C:N ratio in newly formed residue */

            double Csom; /*!< SOM C content [kg ha-1] */
            double Nsom; /*!< SOM N content [kg ha-1] */
            double CNsomNew; /*!< new SOM C:N content */
            double CsomPerc; /*!< SOM C percentage mass [%] */
            double NsomPerc; /*!< SOM N percentage mass [%] */
            double CNsom; /*!< SOM C:N ratio */
            double NO3; /*!< soil NO3 content [kg ha-1] */
            double NH4; /*!< soil NH4 content [kg ha-1] */
            double NO3_weight; /*!< soil patch NO3 content [kg] */
            double NH4_weight; /*!< soil patch NH4 content [kg] */
            double denit; /*!< total NO3 denitrified in the timestep [kg ha-1] */
            double decompSOC; /*!< total SOC decomposed in the timestep [kg ha-1] */
            double resDecompFraction; /*!< fraction of residue to be decomposed in the timestep, passed to net mineralization */

            double NresComposF; /*!< residue N decomposition rate */
            double NcyclTF; /*!< N cycle temperature factor */
            double NcyclWF; /*!< N cycle water factor */
            double NcyclOF; /*!< N cycle oxygen factor */
            double NcyclEF; /*!< N cycle environment factor */
            double combNcyclTFWF; /*!< N cycle combined temperature and water factor */
            double CsomSat; /*!< SOM C percentage mass that saturates soil */

            double decompRes; /*!< total C residue decomposed in the timestep [kg ha-1] */
            double humificRes; /*!< humified portion of residue decomposition in the timestep [kg ha-1] */
            double humificRate; /*!< residue humification rate [kg ha-1]*/
            double respRes; /*!< C residue respiration rate */
            double minNetRes; /*!< N net mineralization from residue in the timestep [kg ha-1] */
            double minSOM; /*!< N mineralization from SOM in the timestep [kg ha-1] */
            double humificResN; /*!< total N net humification in the timestep [kg ha-1] */
            double BGlitter;
            double CNBGlitter;
            double CfracBGLitter; // To be substituted by Basti's input
            double nitrified; // amount of NH4 nitrified
            double volatilized; // amount of amonium volatilized

            double concAvailNO3; // concentration of NO3 in mobile water, kg N mm-1
            double mobileH2O; // mobile water, leaving the layer mm
            double leaching_out; // vertical movement of NO3 leaving the layer (kg ha-1)
            double leaching_in; // vertical movement of NO3 arriving at the layer (kg ha-1)
            double debitNetMin; // nitrogen used in immobilization withdrawn from NH4 or NO3, to be added to netMin to have the flow output correct.
        };

        vector< vector<double> > nutrientTest; /*!< test vector for each layer */

        const Patch& patch;     /*!< Reference to Patch object */
        const Input& input;     /*!< Reference to Input object */
        vector<NutrientLayer> nutrientLayer;   /*!< List of nutrient layers per patch */

        //Member functions
        void calcSoilSurfTemperature(int day);
        void calcSoilTemperature(int layer, int day);
        void incomeRes(int layer);
        void setFactorsNcycle(int layer);
        void denitrification(int layer);
        void SOCdecomp(int layer);
        void ResDecomp(int layer);
        void newCNsom(int layer);
        void ResHumific(int layer);
        void netMineralization(int layer);
        void checkNavail(int layer);
        void nitrificationVolatilization(int layer);
        void processesToPools(int layer);
        void leachingNO3(int layer);
        void Ndeposition();

    private:
        //Member variables
        int xCor;   /*!< x-coordinate of the patch */
        int yCor;   /*!< y-coordinate of the patch */

        struct soilproperties {
            string soilname;    /*!< name of the soil texture */
            double suction,     /*!< suction (S_f) [mm] */
            satHydrCon,     /*!< saturated hydraulic conductivity [mm/h] */
            fieldCap,    /*!< field capacity in vol% */
            evapoFactor, /*!< evaporation factor to reduce/increase EP from the upper layer */
            residualWater,   /*!< residual water content, the soil cannot dry out below this content [vol%] */
            maxInfRate,  /*!< rate of maximal infiltration (F_LX,frac) */
            maxTotalInf, /*!< maximal total infiltration (F_LX,max) [mm/h] */
            diffSpeed,   /*!< diffusion speed between layers */
            bulkDensity, /*!< bulk density [g*cm-3] */
            clayContent, /*!< clay content [m3*m-3] */
            poreSize,   /*!< pore size distribution lambda [-] */
            bubPressure,   /*!< bubbling pressure hb [cm] */
            porosity; /*!< porosity phi [m3*m-3] */

        } soilProp;

};

#endif // NUTRIENTPATCH_H
