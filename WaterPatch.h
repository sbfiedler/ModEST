//!  WaterPatch.
/*!
  This class handles all water processes and parameters relevant on the patch scale.
*/

#ifndef WATERPATCH_H
#define WATERPATCH_H

#include <string>
#include <vector>

using namespace std;

//Forward declaration
class Input;
class Patch;

class WaterPatch
{
    public:
        /*!
         * \brief A constructor.
         * \param index position of soil type in the soil type vector
         * \param _xCor x-coordinate of the patch
         * \param _yCor y-coordinate of the patch
         * \param input pointer to the Input object
         * \param patch pointer to the Patch object
         */
        WaterPatch(int index, int _xCor, int _yCor, const Input& input, const Patch& _patch);

        //! A destructor.
        virtual ~WaterPatch();

        //Member functions
        /*!
         * \brief Fast, preferential infiltration into layers via macropores.
         * \param layer current water layer
         */
        void fastInfiltration(int layer);

        /*!
         * \brief Slow infiltration into the upper layer vai wetting front.
         * \param layer current water layer
         */
        void slowInfiltration(int layer);

        /*!
         * \brief Water flow into the layers if field capacity is exceeded.
         * \param layer current water layer
         */
        void drainage(int layer);

        /*!
         * \brief Water flow between the layers.
         * \param layer current water layer
         */
        void diffusion(int layer);

        /*!
         * \brief Evaporation from the surface.
         * \param potEP potential evaporation
         */
        void surfaceEvaporation(double potEP);

        /*!
         * \brief Evaporation from the first soil layer.
         * \param potEP potential evaporation
         */
        void soilEvaporation(double potEP);

        /*!
         * \brief layerEvapotranspiration
         * \param layer
         * \param potEP potential evaporation
         */
        void layerEvapotranspiration(int layer, double potEP);

        //Member variables
        struct soilproperties {
            string soilname;    /*!< name of the soil texture */
            double suction,     /*!< effective suction [mm] */
            satHydrCon,     /*!< saturated hydraulic conductivity [mm*day-1] */
            fieldCap,    /*!< field capacity [m3*m-3] */
            evapoFactor, /*!< evaporation factor to reduce/increase EP from the upper layer */
            residualWater,   /*!< residual water content, the soil cannot dry out below this content [m3*m-3] */
            infRate,  /*!< rate of maximal infiltration [-] */
            maxTotalInf, /*!< maximal total infiltration [mm*day-1] */
            diffSpeed,  /*!< diffusion speed between layers [-] */
            bulkDensity, /*!< bulk density [g*cm-3] */
            clayContent, /*!< clay content [m3*m-3] */
            poreSize,   /*!< pore size distribution lambda [-] */
            bubPressure,   /*!< bubbling pressure hb [cm] */
            porosity; /*!< porosity phi [m3*m-3] */

        } soilProp;

        struct WaterLayer {
            int depthLayer; /*!< Depth of the soil layer */
            double absWater;    /*!< Absolute water of the layer [mm] */
            double relWater;    /*!< Relative water of the layer [vol%] */
            double transpiredAbsWater; /*!< Absolute transpired water of the layer [mm] */
            double infWater;    /*!< total infiltrated water coming from this layer [mm*day-1] */
            double drainedWater; /*!< total drained water coming from this layer [mm*day-1] */
            double diffusedWater; /*!< total diffused water coming from this layer [mm*day-1] */
            double infWaterSat;  /*!< infiltrated water until saturation/ponding [mm*day-1] */
            double timeSat; /*!< time until saturation/ponding [days] */
            bool wettingFrontFlag;  /*!< reminder whether ponding already occured in this time step */
        };

        double surfaceWater;    /*!< Absolute surface water of the patch [mm] */
        int elevation;  /*!< Elevation of the patch */
        double aspectFactor; /*!< impact of aspect on radiation */
        double inclination; /*!< inclination of patch to lowest neighbour [radians] */
        int xCor;   /*!< x-coordinate of the patch */
        int yCor;   /*!< y-coordinate of the patch */
        double totalRunOff;  /*!< total run-off [mm*day-1] */

        double unsatHydrCon;     /*!< unsaturated hydraulic conductivity [mm*day-1] */
        double geomHydrCon;     /*!< geometric hydraulic conductivity [mm*day-1] */
        double infRateBare; /*!< determines infiltration rate without vegetation [0, 1] */
        double deepDrainedWater; /*!< Sum of deep drianed water per patch over the simulation [mm] */

        const Patch& patch;     /*!< Reference to Patch object */
        vector<WaterLayer> waterLayer;   /*!< List of water layers per patch */
};

#endif // WATERPATCH_H
