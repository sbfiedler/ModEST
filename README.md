# Introducing ModEST

ModEST (Modelling Ecosystem Functions and Services based on Traits) simulates coupled daily dynamics of nutrients, water, and individual woody plants, from which different ecosystem functions and services can be estimated. The model landscape is subdivided into grid cells (5 x 5 m²), two soil layers, and individual plants characterized by coordinates within the landscape. The model can be run for different environmental settings concerning soil texture, climatic conditions, topography, initial plant composition and their traits, with full descriptions given in Fiedler et al. (unpublished). In the following, we briefly describe the three coupled modules of ModEST.

The nutrient module is based on processes for simulating soil nitrogen and soil carbon described in the model SWAT (Kemanian et al., 2011). Daily dynamics of soil organic matter (SOM), nitrate, and ammonium in two soil layers are driven by nitrogen deposition from the atmosphere, decomposition and humification of plants’ residue to SOM, immobilization, mineralization to ammonium, nitrification to nitrate as well as nutrient losses through volatilization, denitrification and leaching. 

We based the hydrological module on the approach of Tietjen et al. (2009), who simulated surface water as well as soil moisture in two soil layers. Daily water dynamics are driven by precipitation, lateral water redistribution of surface water, infiltration, and vertical fluxes, and by water losses via evaporation and transpiration. For ModEST, we adopted these processes with exception of transpiration implemented after the model LPJ (Sitch et. al., 2003) and LPJmL (Schapoff et al., 2017), to better account for stomatal conductance (see description of the transpiration process in Fiedler et al., unpublished). 

The plant module is mainly based on LPJ and LPJmL (Schaphoff et al., 2017; Sitch et al., 2003; Smith et al., 2014) and local processes as described for an individual-based plant model by May et al. (2009). The module simulates the life cycle of individual woody plants placed in the landscape, their dynamic below- and aboveground carbon and nitrogen pools as well as structural components (e.g. plant height, crown area) based on plant traits and abiotic conditions. We adopted – with some changes – the plant processes photosynthesis, transpiration, respiration, reproduction, and allocation after Sitch et. al. (2003) and Schapoff et al. (2017), nitrogen uptake after Smith et al. (2014), as well as dispersal and establishment after May et al. (2009). We added a simple plant mortality process based on annual plant growth and a species-specific growth threshold below which the individual plant dies. Given these adaptations, we fully describe this module in Fiedler et al. (unpublished).

### References
Fiedler et al. (unpublished). Restoration ecologists might not get what they want: Global change shifts trade-offs among ecosystem functions.

Kemanian, A.R., Julich, S., Manoranjan, V.S., Arnold, J.R., 2011. Integrating soil carbon cycling with that of nitrogen and phosphorus in the watershed model SWAT: Theory and model testing. Ecol. Modell. 222, 1913–1921. https://doi.org/10.1016/j.ecolmodel.2011.03.017

May, F., Grimm, V., Jeltsch, F., 2009. Reversed effects of grazing on plant diversity: The role of below-ground competition and size symmetry. Oikos 118, 1830–1843. https://doi.org/10.1111/j.1600-0706.2009.17724.x

Schaphoff, S., von Bloh, W., Rammig, A., Thonicke, K., Biemans, H., Forkel, M., Gerten, D., Heinke, J., Jägermeyr, J., Knauer, J., Langerwisch, F., Lucht, W., Müller, C., Rolinski, S., Waha, K., 2017. LPJmL4 - a dynamic global vegetation model with managed land: Part I - Model description. Geosci. Model Dev. Discuss. 1–59. https://doi.org/10.5194/gmd-2017-145

Sitch, S., Smith, B., Prentice, I.C., Arneth,  a., Bondeau,  a., Cramer, W., Kaplan, J.O., Levis, S., Lucht, W., Sykes, M.T., Thonicke, K., Venevsky, S., 2003. Evaluation of ecosystem dynamics, plant geography and terrestrial carbon cycling in the LPJ dynamic global vegetation model. Glob. Chang. Biol. 9, 161–185. https://doi.org/10.1046/j.1365-2486.2003.00569.x

Smith, B., Wärlind, D., Arneth, A., Hickler, T., Leadley, P., Siltberg, J., Zaehle, S., 2014. Implications of incorporating N cycling and N limitations on primary production in an individual-based dynamic vegetation model. Biogeosciences 11, 2027–2054. https://doi.org/10.5194/bg-11-2027-2014

Tietjen, B., Zehe, E., Jeltsch, F., 2009. Simulating plant water availability in dry lands under climate change: A generic model of two soil layers. Water Resour. Res. 45, 1–14. https://doi.org/10.1029/2007WR006589
