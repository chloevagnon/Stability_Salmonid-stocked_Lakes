# Stability_Salmonid-stocked_Lakes
 This repository is the repository for the article: 
 Changes in vertical and horizontal diversities mediated by the size structure of introduced fish collectively shape food-web stability \
 Chloé VAGNON1,2, Justin POMERANZ3, Bertrand LOHEAC4, Manuel VALLAT4, Jean GUILLARD1,2, Jean-Claude RAYMOND2,5, Arnaud SENTIS2,6, Victor FROSSARD1,2\

1Univ. Savoie Mont Blanc, INRAE, UMR CARRTEL, 74200 Thonon-les-Bains, France.
2Pôle R&D Ecosystèmes Lacustres (ECLA), OFB-INRAE-USMB, 13182 Aix-en-Provence, France. 
3Colorado Mesa Univ., CO 81501, Colorado, United states.
4Fédération de Savoie pour la Pêche et la Protection du Milieu Aquatique (FDPPMA 73), 73230 Saint-Alban-Leysse, France
5Office Français pour la Biodiversité, Unité Spécialisée Milieux Lacustres, 74200 Thonon-les-Bains, France.
6INRAE, Aix Marseille Univ., UMR RECOVER, 13182 Aix-en-Provence, France.\
\

 This repository is composed of 6 elements: \
    - cv_Functions_Stability_Salmonid_Stocked_Lakes: the functions used in the framework of the Figure 1 \
    - Param_reginvert.Rdata : quantile regressions required to compute niche boundaries of invertebrate consumers (from the aNM; Vagnon et al. 2021) \
    - Param_regvert.Rdata :  quantile regressions required to compute niche boundaries of vertebrate consumers (from the aNM; Vagnon et al. 2021) \
    \
    - cv_Codes_Stability_Salmonid-stocked_Lakes: the code for applying the functions to the simple example of Test_DATA.Rdata \
    - Test_Data: a test dataset to apply quickly the different functions and understand how can be computed stability \
    - BigData: the whole dataset of the 195 nodes included in the metaweb of the article.\
\
In datasets, the clumn "Category" corresponds to the category of the studied node (can be fish, invertebrate, zooplankton, producer,...), "Type" is similar, "Organism" is the type of organism (here identified at the family, genus or species), "BodySize" is the total mean body length for the node expressed in micrometer, "Trophy" is the diet type of the node, "Name" is the name of the node, "Density" is the density of the node in the corresponding lake for which the unit is different depending on the category of the node (See the method in the article for more details).
\
 For any further information, please contact chloe.vagnon@gmail.com
\
\
\
Vagnon, C., Cattanéo, F., Goulon, C., Grimardias, D., Guillard, J. et al. (2021). An allometric niche model for species interactions in temperate freshwater ecosystems. Ecosphere, 12, e03420.


