# Unsupervised Manifold Mapping of Connectomic Fragility in Glioma: Simulating Global Efficiency Decay and Lobar Vulnerability

> **Note for Reviewers:** This repository hosts the official computational framework and reproducible workflows corresponding to the conference paper prepared for **6to Workshop Chileno sobre Reconocimiento de Patrones (CWPR 2026)**.

This repository contains the official **unsupervised pattern recognition framework** designed to map glioma-induced connectomic fragility, discover latent structural subphenotypes, and simulate cascading network attacks to assess systemic brain resilience.

### Authors & Affiliations

<p align="left">
  <strong>Pamela Franco</strong> <a href="https://orcid.org/0000-0001-7629-3653"><img src="https://img.shields.io/badge/ORCID-0000--0001--7629--3653-A6CE39?logo=orcid&logoColor=white&style=flat-square" height="16"></a><br>
  <small>• Energy Transformation Center, Universidad Andrés Bello, Santiago, Chile.</small>
</p>


<p align="left">
  <strong>Cristian Montalba</strong> <a href="https://orcid.org/0000-0003-3370-0233"><img src="https://img.shields.io/badge/ORCID-0000--0003--3370--0233-A6CE39?logo=orcid&logoColor=white&style=flat-square" height="16"></a><br>
  <small>• Biomedical Imaging Center and Radiology Department, Pontificia Universidad Católica de Chile.<br>• Millennium Institute for Intelligent Healthcare Engineering (iHEALTH).</small>
</p>

<p align="left">
  <strong>Ignacio Espinoza</strong> <a href="https://orcid.org/0000-0003-2400-4498"><img src="https://img.shields.io/badge/ORCID-0000--0003--2400--4498-A6CE39?logo=orcid&logoColor=white&style=flat-square" height="16"></a><br>
  <small>• Institute of Physics, Pontificia Universidad Católica de Chile, Santiago, Chile.</small>
</p>

<p align="left">
  <strong>M. Daniela Cornejo</strong> <a href="https://orcid.org/0009-0003-0425-5721"><img src="https://img.shields.io/badge/ORCID-0009--0003--0425--5721-A6CE39?logo=orcid&logoColor=white&style=flat-square" height="16"></a><br>
  <small>• Institute of Physics / Department of Psychiatry, School of Medicine, Pontificia Universidad Católica de Chile.</small>
</p>

<p align="left">
  <strong>Francisco Torres</strong> <a href="https://orcid.org/0000-0002-0003-2446"><img src="https://img.shields.io/badge/ORCID-0000--0002--0003--2446-A6CE39?logo=orcid&logoColor=white&style=flat-square" height="16"></a><br>
  <small>• Radiology Department, Hospital Carlos van Buren / Universidad de Valparaíso, Valparaíso, Chile.</small>
</p>


<p align="left">
  <strong>Carlos Bennet</strong> <a href="https://orcid.org/0009-0007-1434-273X"><img src="https://img.shields.io/badge/ORCID-0009--0007--1434--273X-A6CE39?logo=orcid&logoColor=white&style=flat-square" height="16"></a><br>
  <small>• Neurosurgery Department, Hospital Carlos van Buren, Valparaíso, Chile.</small>
</p>

<p align="left">
  <strong>Steren Chabert</strong> <a href="https://orcid.org/0000-0002-2890-5077"><img src="https://img.shields.io/badge/ORCID-0000--0002--2890--5077-A6CE39?logo=orcid&logoColor=white&style=flat-square" height="16"></a><br>
  <small>• Biomedical Engineering School, Universidad de Valparaíso.<br>• Millennium Institute for Intelligent Healthcare Engineering (iHEALTH).<br>• Center of Interdisciplinary Biomedical and Engineering Research for Health (MEDING).</small>
</p>

<p align="left">
  <strong>Rodrigo Salas</strong> *(Corresponding Author)* <a href="https://orcid.org/0000-0002-0350-6811"><img src="https://img.shields.io/badge/ORCID-0000--0002--0350--6811-A6CE39?logo=orcid&logoColor=white&style=flat-square" height="16"></a><br>
  <small>• Biomedical Engineering School, Universidad de Valparaíso.<br>• Millennium Institute for Intelligent Healthcare Engineering (iHEALTH).<br>• Center of Interdisciplinary Biomedical and Engineering Research for Health (MEDING).<br> Email: rodrigo.salas@uv.cl</small>
</p>

---

## Abstract
Characterizing structural white matter (WM) disruption induced by primary brain gliomas yields highly complex, high-dimensional topological architectures. While conventional neurooncology paradigms rely heavily on histopathological tumor grading, they inherently fail to capture the latent structural heterogeneity and systemic network vulnerability across patients. 

In this study, we propose an unsupervised pattern recognition framework designed to map glioma-induced connectomic fragility. Leveraging high-dimensional graph-theoretical metrics mapped via non-linear Uniform Manifold Approximation and Projection (UMAP), we identify stable topological subphenotypes ($K=2$) that represent distinct profiles of network resilience. We evaluate these cohorts using simulated targeted node attacks and map exponential decay cascades of network Global Efficiency. Furthermore, we cross-validate the integrity of the core structural network against clinical lobar infiltration zones and demographic parameters. Our results demonstrate that unsupervised manifold mapping successfully segregates patients into distinct risk-resilience profiles, proving that tumor-masked connectomics can model the macroscopic dynamics of network collapse independent of classic grading boundaries.

---

## Framework & Pipeline Overview

The complete pattern recognition pipeline is structured into five core phases:
1. **Multimodal Neuroimaging Preprocessing**: Head motion and eddy-current distortion correction of Diffusion Tensor Imaging (DTI) sequences using a non-parametric Gaussian Process framework to replace signal dropouts.
2. **Tumor-Masked Network Matrix Construction**: To avoid structural misregistration from mass effects, an adaptive threshold (50% maximum tumor intensity) excludes tumor-encroached White Matter regions of interest (ROIs) intersecting the JHU ICBM-DTI-81 Atlas. Fiber orientations are estimated using BEDPOSTX, and probabilistic tractography is propagated using ProbtrackX2.
3. **Multi-Level Connectomic Feature Extraction**: Individual structural matrices are converted into brain graphs across three distinct tiers:
   * **Tier 1 (Standardization)**: Raw matrix parameters (active nodes, edge weights, self-loop removal).
   * **Tier 2 (Global Metrics)**: Edge count, density, streamline sum, and Global Efficiency.
   * **Tier 3 (Local Descriptors)**: Node Strength, Node Degree, PageRank, Closeness, Betweenness Centrality, and weighted Local Clustering Coefficient across a uniform adjusted 307-dimensional master space ($X \in \mathbb{R}^{35 \times 307}$).
4. **Non-linear Manifold Learning & Subphenotypical Clustering**: High-dimensional standardization paired with non-linear **UMAP** dimensionality reduction ($Y \in \mathbb{R}^{35 \times 2}$) optimizing local structural boundaries via fuzzy set cross-entropy loss. Latent coordinates are partitioned via **Spectral Clustering** ($K=2$) to group patients purely by structural network geometry.
5. **Connectomic Vulnerability & Network Attack Simulation**: Evaluation of systemic network resilience modeling targeted node deletions over an exponential decay function:
   $$\kappa(s) = \kappa_{0} \cdot e^{-s \cdot \lambda_{c}}$$
   Unveiling clear divergence in disintegration kinetics between the discovered cohorts ($\lambda_{HGG} = 2.1$ vs. $\lambda_{LGG} = 1.2$). Additionally, deep callosal structural integrity is assessed using a customized *Structural Core Index (SCI)*.

---

## Repository Structure

```

├── Codes/
│   ├── ML_code.py                                     # Main Python code: Z-score standardization, UMAP projection, and Spectral Clustering.
│   └── requirements.txt                               # Python package dependencies.
├── Dataset/
│   ├── dataset_conectomicas_with_patient_details.csv  # Pre-calculated 307 master topological feature dataset from the 35-patient glioma cohort.
│   └── Patients_Project.csv                           # Clinical and demographic data
└── Results/
    ├── manifold_embedding.png                         # 2D Riemannian UMAP embedding visualization separating the subphenotypes.
    │                                                    - Description: A high-resolution scatter plot showing the low-dimensional geometric 
    │                                                      projection of structural connectomes. It visually demonstrates the clear unsupervised 
    │                                                      topological separation between the LGG (Low-Grade Glioma) and HGG (High-Grade Glioma) 
    │                                                      subphenotype clusters.
    │
    ├── resilience_decay.png                           # Exponential Global Efficiency decay curves under simulation.
    │                                                    - Description: A line plot illustrating the simulated connectomic vulnerability 
    │                                                      analysis. It contrasts the network global efficiency ($\kappa$) degradation pathways 
    │                                                      between LGG and HGG cohorts under cumulative targeted node deletion stages (0% to 100%).
    │
    ├── figure_lobar_core_integrity.png                # Boxplot of structural core integrity across infiltrated lobes.
    │                                                    - Description: A hybrid box-and-strip plot that tracks the Structural Core Connection 
    │                                                      Index (derived from corpus callosum metrics) grouped by the primary tumor-infiltrated 
    │                                                      cerebral lobe (Frontal, Temporal, Parietal) to evaluate localized structural attenuation.
    │
    ├── figure_distribution_sex.png                    # Stacked bar chart for demographic gender validation across clusters.
    │                                                    - Description: A stacked frequency bar plot evaluating the distribution of Female and 
    │                                                      Male patients stratified by the unsupervised connectomic subphenotypes (LGG vs. HGG), 
    │                                                      supporting the corresponding Chi-Squared independence tests.
    │
    └── figure_distribution_tumor_hemisphere.png       # Stacked bar chart for hemispheric tumor localization across clusters.
                                                         - Description: A stacked frequency bar plot displaying the distribution of Left, 
                                                           Right, and Bilateral tumor hemispheric localization across the stratified LGG 
                                                           and HGG structural subphenotypes.                                          

```

---

## Reproducibility & Data Availability

- **Code**: The ML pipelines script is hosted in this repository.
- **Data Privacy**: MRI datasets and structural connectivity matrices generated during this study are not publicly available due to patient data privacy restrictions imposed by the Ethics Committee of the Servicio de Salud Valparaíso San Antonio (ORD.001413).
- **Clinical Inquiries**: Anonymized data may be made available upon reasonable request and subject to institutional approval by contacting **Steren Chabert (steren.chabert@uv.cl)**.

---

## How to Run

1. **Clone the repository**:
   ```bash
   git clone [https://github.com/pamelaFranco/unsupervised-glioma-connectomics.git](https://github.com/pamelaFranco/unsupervised-glioma-connectomics.git)
cd unsupervised-glioma-connectomics
   ```

2. **Set up neuroimaging environment**: Ensure you have FSL (v6.0) installed and configured in your environment path.

3. **Install Python dependencies**:
   ```bash
   pip install -r requirements.txt
   ```
4. Execute the connectomics pipeline:  
   ```bash
   # Run machine learning pipeline
   python ML_code.py
   ```
  
---

## Acknowledgements
This work was supported by the National Agency for Research and Development (ANID) of Chile through:

* FONDECYT N°1221938 ("An Explainable Deep Neuro-Fuzzy Inference System for the segmentation of BT in multi-contrast magnetic resonance imaging").

* ANID Millennium Science Initiative Program ICN2021_004 (Millennium Institute for Intelligent Healthcare Engineering - iHEALTH).

* Additionally, this work was funded by the Endowment I+D in Health Competition of the Universidad Andrés Bello (UNAB) 2025, project No. DI-07-25/ICS.

---
## Citation

If you find this pipeline useful for your research, please cite our preliminary work prepared for **Neuroradiology**:

```bibtex
@inproceedings{franco2026unsupervised-glioma-connectomics,
  title={Unsupervised Manifold Mapping of Connectomic Fragility in Glioma: Simulating Global Efficiency Decay and Lobar Vulnerability},
  author={Franco, Pamela and Montalba, Cristian and Espinoza, Ignacio and Cornejo, M. Daniela and Torres, Francisco and Bennett, Carlos and Chabert, Steren and Salas, Rodrigo},
  booktitle={Proceedings of the Chilean Computer Science Days (JCC 2026) - 6th Chilean Workshop on Pattern Recognition (CWPR)},
  publisher={IEEE},
  year={2026},
  note={Submitted for publication / Under review}
}
```

---

## License

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)