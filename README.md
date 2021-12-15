# Linking protein structural and functional change to mutation using amino acid networks

by Cristina Sotomayor-Vivas, Enrique Hernández-Lemus, and Rodrigo Dorantes-Gilardi


## Abstract

The function of a protein is strongly dependent on its structure. During evolution, proteins acquire new functions through mutations in the amino-acid sequence. Given the advance in deep mutational scanning, recent findings have found functional change to be position dependent, notwithstanding the chemical properties of mutant and mutated amino acids. This could indicate that structural properties of a given position are potentially responsible for the functional relevance of a mutation. Here, we looked at the relation between structure and function of positions using five proteins with experimental data of functional change available. In order to measure structural change, we modeled mutated proteins via amino-acid networks and quantified the perturbation of each mutation. We found that structural change is position dependent, and strongly related to functional change. Strong changes in protein structure correlate with functional loss, and positions with functional gain due to mutations tend to be structurally robust. Finally, we constructed a computational method to predict functionally sensitive positions to mutations using structural change that performs well on all five proteins with a mean precision of 74.7% and recall of 69.3% of all functional positions.

![loss_predictions](https://github.com/CrisSotomayor/perturbation-networks/blob/main/figures/loss_text.png)

_A) Experimentally obtained functional data from deep mutational scan of VIM-2 protein, with darker values representing higher functional disruption, specifically blue is loss of function while red represents gain of function [2]. B) Standardized data of the number of nodes perturbed by each mutation where each entry is the number of standard deviations from the mean of the distribution. The perturbation network was constructed using a threshold of 9 Å; blue represents highest structural perturbation, and red represents lowest. C) Predictions maximizing precision. X-axis has the sequence positions, Y-axis has the experimentally obtained mean functional value. Blue dots are SSPs---our predictions for FSPs---while shaded blue area contains the 40% of sequence positions with lowest functional scores representing strongest functional loss. Top row shows the functional values experimentally obtained for VIM-2 protein, bottom row the other four proteins studied. D) Predictions maximizing recall. E) Predictions maximizing both measures._


## Data

All the data used can be found under `data`, including the three-dimensional coordinate files for the proteins selected (under `data/pdb`), deep mutational scanning data for the five proteins selected (under `data`, files `functional_{pdb_id}.csv`) and the data obtained from the different perturbation networks (under `data/structure`). The proteins selected and corresponding papers are PSD95^{pdz3} (PDB: 1BE9) [4], phosphatase and tensin homolog (PTEN) (PDB: 1D5R) [5], APH(3')II (PDB: 1ND4) [6], Src kinase catalytic domain (Src CD, PDB: 3DQW) [1], and VIM-2 metallo-$\beta$-lactamase (PDB: 4BZ3) [2].


## Software implementation

Code to obtain the data is found in `.py` files.
- `getdata.py` contains functions to get data from amino acid mutation perturbation networks as CSV files.
- `getmutations.py` contains functions and code to obtain mutations using FoldX [7].
- `iterate.py` contains functions and code to run both `getmutations.py` and `getdata.py`, including the functions from the module `biographs` that were used.
- `main_functions.py` contains functions used in jupyter notebooks for figures and analysis, specifically to read the stored data and obtain predicted positions
- `parallel_getdata.py` contains functions to get data from amino acid mutation perturbation networks to use with multiprocessing.

To use, you will first have to generate the mutations from a specific pdb file, using `getmutations.py`. Modify the main function with the protein you want and the path where protein.pdb file is stored, as well as the path where the FoldX [7] software is installed. To obtain the data from the resulting pdb files, either manually call the `GetData` function from `getdata.py` for each desired protein, or run the `parallel_getdata.py` file, changing the list of proteins to use. In either case, the thresholds list can be modified as desired.


## Requirements

Python 3 is required to run the code, as well as the libraries `Biopyhton, seaborn, matplotlib, pandas, numpy, scipy, biographs, NetworkX` and `scikit-learn`.


## References

1. Ahler, E., Register, A. C., Chakraborty, S., Fang, L., Dieter, E. M., Sitko, K. A., Vidadala, R. S. R., Trevillian, B. M., Golkowski, M., Gelman, H., Stephany, J. J., Rubin, A. F., Merritt, E. A., Fowler, D. M., & Maly, D. J. (2019). A Combined Approach Reveals a Regulatory Mechanism Coupling Src’s Kinase Activity, Localization, and Phosphotransferase-Independent Functions. Molecular Cell. https://doi.org/10.1016/j.molcel.2019.02.003
2. Chen, J. Z., Fowler, D. M., & Tokuriki, N. (2020). Comprehensive exploration of the translocation, stability and substrate recognition requirements in vim-2 lactamase. ELife. https://doi.org/10.7554/eLife.56707
3. Dorantes-Gilardi, R. (2020). Biographs: Amino acid networks in python.
4. McLaughlin, R. N., Poelwijk, F. J., Raman, A., Gosal, W. S., & Ranganathan, R. (2012). The spatial architecture of protein function and adaptation. Nature. https://doi.org/10.1038/nature11500
5. Melnikov, A., Rogov, P., Wang, L., Gnirke, A., & Mikkelsen, T. S. (2014). Comprehensive mutational scanning of a kinase in vivo reveals substrate-dependent fitness landscapes. Nucleic Acids Research. https://doi.org/10.1093/nar/gku511
6. Mighell, T. L., Evans-Dutson, S., & O’Roak, B. J. (2018). A Saturation Mutagenesis Approach to Understanding PTEN Lipid Phosphatase Activity and Genotype-Phenotype Relationships. American Journal of Human Genetics. https://doi.org/10.1016/j.ajhg.2018.03.018
7. Schymkowitz, J., Borg, J., Stricher, F., Nys, R., Rousseau, F., & Serrano, L. (2005). The FoldX web server: An online force field. Nucleic Acids Research. https://doi.org/10.1093/nar/gki387
