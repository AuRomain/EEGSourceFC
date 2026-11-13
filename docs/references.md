# References

What to cite when you publish results obtained with EEGSourceFC:

1. the pipeline itself (README, "How to cite"; `CITATION.cff`);
2. the code it was built from (below);
3. the toolboxes it calls, and the methods of the parts you used and
   report (an inverse solution, a connectivity measure, the statistics).

Every reference below was checked against PubMed, Crossref, DataCite or
arXiv.

## The code it was built from

EEGSourceFC was built from two public repositories of the same research
group, then corrected, generalised and extended:

- **var-EEG-FC** (Sahar Allouch; <https://github.com/sahar-allouch/var-EEG-FC>):
  stage 1, i.e. the inverse solutions (eLORETA, LCMV, wMNE: `+wmne`), the
  connectivity measures (`+connectivity`, `+ROInets`), the conversion of
  Brainstorm's files to FieldTrip's (`BS_to_ft_*`) and `run_all`.
  - Allouch S, Kabbara A, Duprez J, Paban V, Khalil M, Modolo J, Hassan M
    (2023). Effect of analytical variability in estimating EEG-based
    functional connectivity. *bioRxiv* 2023.08.17.553675.
    <https://doi.org/10.1101/2023.08.17.553675>
- **dynCogPD** (Judie Tabbal, <https://github.com/judytabbal/dynCogPD>;
  Joan Duprez's fork, <https://github.com/jduprez/dynCogPD>): the main
  script of stage 1 (`run_pipeline`) and the sliding windows
  (`+connectivity`), which give the window mean here.
  - Duprez J, Tabbal J, Hassan M, Modolo J, Kabbara A, Mheich A, et al.
    (2022). Spatio-temporal dynamics of large-scale electrophysiological
    networks during cognitive action control in healthy controls and
    Parkinson's disease patients. *NeuroImage* 258:119331.
    <https://doi.org/10.1016/j.neuroimage.2022.119331>

## Toolboxes

Called by the pipeline (the versions validated: README, section 1, and `toolbox_reference.json`, which `toolbox_check` compares the installed files with):

| Toolbox | Used for | Cite |
|---|---|---|
| FieldTrip | reading the files, eLORETA, LCMV, the debiased wPLI (`ft_connectivity_wpli`), IRASA (`ft_freqanalysis`); its copy of Brainstorm's specparam | Oostenveld et al. 2011 |
| Brainstorm | band-pass filters, region values (scouts), the ICBM152 head, atlas and caps | Tadel et al. 2011 |
| OpenMEEG | the leadfields shipped in `inputs/`, and that of a new cap | Gramfort et al. 2010; Kybic et al. 2005 |
| Brain Connectivity Toolbox (its definitions, in `gt_window_metrics`; not called, a reference in the tests) | graph metrics | Rubinov & Sporns 2010 |
| NBS 1.2 (included, `3- Analysis/fun/NBS`) | network-based statistic | Zalesky et al. 2010 |
| restingIAF | individual alpha frequency | Corcoran et al. 2018 |
| ROInets (included, `+ROInets`) | symmetric orthogonalisation (leakage correction) | Colclough et al. 2015 |
| *Optional, other head models only* (`tools/make_anatomy.m`): SPM12 (<https://www.fil.ion.ucl.ac.uk/spm/>), CAT12, FreeSurfer | an individual or another template head | Gaser et al. 2024; Fischl 2012 |

Adapted, or used as references in the tests (`tests/test_reference_*.m`):

| Code | What | Cite |
|---|---|---|
| PALM | permutation tests of the general linear model, the z of the Welch t | Winkler et al. 2014 |
| PERMUTOOLS | the max-statistic correction, Hedges' g | Crosse et al. 2024 |
| statsmodels | the HC3 standard error of a regression coefficient (`OLS(...).fit(cov_type='HC3')`) | Seabold & Perktold 2010 |
| wildboottest (0.3.2, github.com/s3alfisc/wildboottest; Python port of fwildclusterboot and boottest) | the restricted wild bootstrap of a regression coefficient | Roodman et al. 2019 |
| BrainSpace; rotate_parcellation (github.com/frantisekvasa/rotate_parcellation) | spin permutations of the regions | Vos de Wael et al. 2020; Váša et al. 2018 |
| FieldTrip `clusterstat`, `findcluster` (adapted in `cluster_perm_test.m`) | the cluster-based permutation test of the PSD analysis | Maris & Oostenveld 2007 |
| restingIAF `restingIAF.m` (adapted in `resting_iaf_pieces.m`) | the IAF of recordings cut into pieces | Corcoran et al. 2018 |
| MNE-Python (`mne/rank.py`, `estimate_rank`) | the rank rule of `connectivity.signal_rank` | Gramfort et al. 2013 |

## Head model and atlas

- **ICBM152 2009 template** (the shipped head): Fonov et al. 2011.
- **Desikan–Killiany atlas** (68 regions): Desikan et al. 2006.
- **Networks of the regions** (`3- Analysis/DK68_networks.csv`): Yeo et
  al. 2011, by overlap with Schaefer et al. 2018.

## Methods

**Stage 1, sources and connectivity**
- eLORETA: Pascual-Marqui 2007; Pascual-Marqui et al. 2011.
- LCMV beamformer: Van Veen et al. 1997.
- Weighted minimum norm (wMNE): Hämäläinen & Ilmoniemi 1994.
- PLV: Lachaux et al. 1999.
- PLI: Stam et al. 2007.
- wPLI and debiased wPLI: Vinck et al. 2011.
- ciPLV: Bruña et al. 2018.
- AEC: Brookes et al. 2011. Orthogonalised AEC: Hipp et al. 2012;
  pairwise regression, Brookes et al. 2012; symmetric orthogonalisation,
  Colclough et al. 2015.
- AEC windows of 6 s every 0.5 s (window mean): O'Neill et al. 2017.
- Individual alpha frequency: Corcoran et al. 2018; IRASA, Wen & Liu
  2016; specparam, Donoghue et al. 2020.
- Data rank (`connectivity.signal_rank`): the rule of MNE-Python,
  Gramfort et al. 2013.

**Stage 2, power spectra**
- Aperiodic and periodic components (specparam): Donoghue et al. 2020.
- Power of each source from the channels' cross-spectra (the source maps):
  Frei et al. 2001.
- Cluster-based permutation test: Maris & Oostenveld 2007.

**Stage 3, statistics**
- Permutation tests, z of the Welch t: Winkler et al. 2014; with
  covariates, Freedman & Lane 1983; the max-statistic over a family:
  Westfall & Young 1993; maxima over tests of different degrees of
  freedom: Alberton et al. 2020.
- Network-based statistic: Zalesky et al. 2010.
- False discovery rate: Benjamini & Hochberg 1995; over families of
  tests: Benjamini & Bogomolov 2014.
- Hedges' g: Hedges 1981.
- Correlations, a studentised permutation test: DiCiccio & Romano 2017;
  the HC3 standard error: MacKinnon & White 1985; Long & Ervin 2000; the
  wild bootstrap: Davidson & Flachaire 2008; MacKinnon 2013; Roodman et
  al. 2019.
- Effect size of an adjusted group difference: What Works Clearinghouse
  2022.
- Spin permutations: Alexander-Bloch et al. 2018; Váša et al. 2018; Vos
  de Wael et al. 2020.
- Graph metrics: Rubinov & Sporns 2010.
- Reliability within recordings: split halves with the Spearman–Brown
  correction, Spearman 1910; Brown 1910; identification of recordings,
  Finn et al. 2015.

## Example data

`examples/settings_openneuro.m` is set up for two public datasets (CC0):

- Salisbury D, Seebold D, Coffman B (2021). EEG: First Episode Psychosis
  vs. Control Resting Task 1. *OpenNeuro* ds003944, version 1.0.0.
  <https://doi.org/10.18112/openneuro.ds003944.v1.0.0>
- Salisbury D, Seebold D, Coffman B (2021). EEG: First Episode Psychosis
  vs. Control Resting Task 2. *OpenNeuro* ds003947, version 1.0.0.
  <https://doi.org/10.18112/openneuro.ds003947.v1.0.0>

## Full references

- Alberton BAV, Nichols TE, Gamba HR, Winkler AM (2020). Multiple testing correction over contrasts for brain imaging. *NeuroImage* 216:116760. <https://doi.org/10.1016/j.neuroimage.2020.116760>
- Alexander-Bloch AF, Shou H, Liu S, Satterthwaite TD, Glahn DC, Shinohara RT, et al. (2018). On testing for spatial correspondence between maps of human brain structure and function. *NeuroImage* 178:540–551. <https://doi.org/10.1016/j.neuroimage.2018.05.070>
- Allouch S, Kabbara A, Duprez J, Paban V, Khalil M, Modolo J, Hassan M (2023). Effect of analytical variability in estimating EEG-based functional connectivity. *bioRxiv* 2023.08.17.553675. <https://doi.org/10.1101/2023.08.17.553675>
- Benjamini Y, Bogomolov M (2014). Selective inference on multiple families of hypotheses. *J R Stat Soc Series B Stat Methodol* 76:297–318. <https://doi.org/10.1111/rssb.12028>
- Benjamini Y, Hochberg Y (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. *J R Stat Soc Series B Stat Methodol* 57:289–300. <https://doi.org/10.1111/j.2517-6161.1995.tb02031.x>
- Brookes MJ, Hale JR, Zumer JM, Stevenson CM, Francis ST, Barnes GR, et al. (2011). Measuring functional connectivity using MEG: methodology and comparison with fcMRI. *NeuroImage* 56:1082–1104. <https://doi.org/10.1016/j.neuroimage.2011.02.054>
- Brookes MJ, Woolrich MW, Barnes GR (2012). Measuring functional connectivity in MEG: a multivariate approach insensitive to linear source leakage. *NeuroImage* 63:910–920. <https://doi.org/10.1016/j.neuroimage.2012.03.048>
- Brown W (1910). Some experimental results in the correlation of mental abilities. *Br J Psychol* 3:296–322. <https://doi.org/10.1111/j.2044-8295.1910.tb00207.x>
- Bruña R, Maestú F, Pereda E (2018). Phase locking value revisited: teaching new tricks to an old dog. *J Neural Eng* 15:056011. <https://doi.org/10.1088/1741-2552/aacfe4>
- Colclough GL, Brookes MJ, Smith SM, Woolrich MW (2015). A symmetric multivariate leakage correction for MEG connectomes. *NeuroImage* 117:439–448. <https://doi.org/10.1016/j.neuroimage.2015.03.071>
- Corcoran AW, Alday PM, Schlesewsky M, Bornkessel-Schlesewsky I (2018). Toward a reliable, automated method of individual alpha frequency (IAF) quantification. *Psychophysiology* 55:e13064. <https://doi.org/10.1111/psyp.13064>
- Crosse MJ, Foxe JJ, Molholm S (2024). PERMUTOOLS: a MATLAB package for multivariate permutation testing. *arXiv* 2401.09401. <https://doi.org/10.48550/arXiv.2401.09401>
- Davidson R, Flachaire E (2008). The wild bootstrap, tamed at last. *J Econom* 146:162–169. <https://doi.org/10.1016/j.jeconom.2008.08.003>
- Desikan RS, Ségonne F, Fischl B, Quinn BT, Dickerson BC, Blacker D, et al. (2006). An automated labeling system for subdividing the human cerebral cortex on MRI scans into gyral based regions of interest. *NeuroImage* 31:968–980. <https://doi.org/10.1016/j.neuroimage.2006.01.021>
- DiCiccio CJ, Romano JP (2017). Robust permutation tests for correlation and regression coefficients. *J Am Stat Assoc* 112:1211–1220. <https://doi.org/10.1080/01621459.2016.1202117>
- Donoghue T, Haller M, Peterson EJ, Varma P, Sebastian P, Gao R, et al. (2020). Parameterizing neural power spectra into periodic and aperiodic components. *Nat Neurosci* 23:1655–1665. <https://doi.org/10.1038/s41593-020-00744-x>
- Duprez J, Tabbal J, Hassan M, Modolo J, Kabbara A, Mheich A, et al. (2022). Spatio-temporal dynamics of large-scale electrophysiological networks during cognitive action control in healthy controls and Parkinson's disease patients. *NeuroImage* 258:119331. <https://doi.org/10.1016/j.neuroimage.2022.119331>
- Finn ES, Shen X, Scheinost D, Rosenberg MD, Huang J, Chun MM, et al. (2015). Functional connectome fingerprinting: identifying individuals using patterns of brain connectivity. *Nat Neurosci* 18:1664–1671. <https://doi.org/10.1038/nn.4135>
- Fischl B (2012). FreeSurfer. *NeuroImage* 62:774–781. <https://doi.org/10.1016/j.neuroimage.2012.01.021>
- Fonov V, Evans AC, Botteron K, Almli CR, McKinstry RC, Collins DL, et al. (2011). Unbiased average age-appropriate atlases for pediatric studies. *NeuroImage* 54:313–327. <https://doi.org/10.1016/j.neuroimage.2010.07.033>
- Freedman D, Lane D (1983). A nonstochastic interpretation of reported significance levels. *J Bus Econ Stat* 1:292–298. <https://doi.org/10.1080/07350015.1983.10509354>
- Frei E, Gamma A, Pascual-Marqui R, Lehmann D, Hell D, Vollenweider FX (2001). Localization of MDMA-induced brain activity in healthy volunteers using low resolution brain electromagnetic tomography (LORETA). *Hum Brain Mapp* 14:152–165. <https://doi.org/10.1002/hbm.1049>
- Gaser C, Dahnke R, Thompson PM, Kurth F, Luders E, the Alzheimer's Disease Neuroimaging Initiative (2024). CAT: a computational anatomy toolbox for the analysis of structural MRI data. *GigaScience* 13:giae049. <https://doi.org/10.1093/gigascience/giae049>
- Gramfort A, Luessi M, Larson E, Engemann DA, Strohmeier D, Brodbeck C, et al. (2013). MEG and EEG data analysis with MNE-Python. *Front Neurosci* 7:267. <https://doi.org/10.3389/fnins.2013.00267>
- Gramfort A, Papadopoulo T, Olivi E, Clerc M (2010). OpenMEEG: opensource software for quasistatic bioelectromagnetics. *Biomed Eng Online* 9:45. <https://doi.org/10.1186/1475-925X-9-45>
- Hämäläinen MS, Ilmoniemi RJ (1994). Interpreting magnetic fields of the brain: minimum norm estimates. *Med Biol Eng Comput* 32:35–42. <https://doi.org/10.1007/BF02512476>
- Hedges LV (1981). Distribution theory for Glass's estimator of effect size and related estimators. *J Educ Stat* 6:107–128. <https://doi.org/10.3102/10769986006002107>
- Hipp JF, Hawellek DJ, Corbetta M, Siegel M, Engel AK (2012). Large-scale cortical correlation structure of spontaneous oscillatory activity. *Nat Neurosci* 15:884–890. <https://doi.org/10.1038/nn.3101>
- Kybic J, Clerc M, Abboud T, Faugeras O, Keriven R, Papadopoulo T (2005). A common formalism for the integral formulations of the forward EEG problem. *IEEE Trans Med Imaging* 24:12–28. <https://doi.org/10.1109/TMI.2004.837363>
- Lachaux JP, Rodriguez E, Martinerie J, Varela FJ (1999). Measuring phase synchrony in brain signals. *Hum Brain Mapp* 8:194–208. <https://doi.org/10.1002/(SICI)1097-0193(1999)8:4%3C194::AID-HBM4%3E3.0.CO;2-C>
- Long JS, Ervin LH (2000). Using heteroscedasticity consistent standard errors in the linear regression model. *Am Stat* 54:217–224. <https://doi.org/10.1080/00031305.2000.10474549>
- MacKinnon JG (2013). Thirty years of heteroskedasticity-robust inference. In: *Recent Advances and Future Directions in Causality, Prediction, and Specification Analysis*, Springer, 437–461. <https://doi.org/10.1007/978-1-4614-1653-1_17>
- MacKinnon JG, White H (1985). Some heteroskedasticity-consistent covariance matrix estimators with improved finite sample properties. *J Econom* 29:305–325. <https://doi.org/10.1016/0304-4076(85)90158-7>
- Maris E, Oostenveld R (2007). Nonparametric statistical testing of EEG- and MEG-data. *J Neurosci Methods* 164:177–190. <https://doi.org/10.1016/j.jneumeth.2007.03.024>
- O'Neill GC, Tewarie PK, Colclough GL, Gascoyne LE, Hunt BAE, Morris PG, Woolrich MW, Brookes MJ (2017). Measurement of dynamic task related functional networks using MEG. *NeuroImage* 146:667–678. <https://doi.org/10.1016/j.neuroimage.2016.08.061>
- Oostenveld R, Fries P, Maris E, Schoffelen JM (2011). FieldTrip: open source software for advanced analysis of MEG, EEG, and invasive electrophysiological data. *Comput Intell Neurosci* 2011:156869. <https://doi.org/10.1155/2011/156869>
- Pascual-Marqui RD (2007). Discrete, 3D distributed, linear imaging methods of electric neuronal activity. Part 1: exact, zero error localization. *arXiv* 0710.3341. <https://doi.org/10.48550/arXiv.0710.3341>
- Pascual-Marqui RD, Lehmann D, Koukkou M, Kochi K, Anderer P, Saletu B, et al. (2011). Assessing interactions in the brain with exact low-resolution electromagnetic tomography. *Philos Trans A Math Phys Eng Sci* 369:3768–3784. <https://doi.org/10.1098/rsta.2011.0081>
- Roodman D, Nielsen MØ, MacKinnon JG, Webb MD (2019). Fast and wild: bootstrap inference in Stata using boottest. *Stata J* 19:4–60. <https://doi.org/10.1177/1536867X19830877>
- Rubinov M, Sporns O (2010). Complex network measures of brain connectivity: uses and interpretations. *NeuroImage* 52:1059–1069. <https://doi.org/10.1016/j.neuroimage.2009.10.003>
- Salisbury D, Seebold D, Coffman B (2021). EEG: First Episode Psychosis vs. Control Resting Task 1. *OpenNeuro* ds003944, version 1.0.0. <https://doi.org/10.18112/openneuro.ds003944.v1.0.0>
- Salisbury D, Seebold D, Coffman B (2021). EEG: First Episode Psychosis vs. Control Resting Task 2. *OpenNeuro* ds003947, version 1.0.0. <https://doi.org/10.18112/openneuro.ds003947.v1.0.0>
- Schaefer A, Kong R, Gordon EM, Laumann TO, Zuo XN, Holmes AJ, et al. (2018). Local-global parcellation of the human cerebral cortex from intrinsic functional connectivity MRI. *Cereb Cortex* 28:3095–3114. <https://doi.org/10.1093/cercor/bhx179>
- Seabold S, Perktold J (2010). Statsmodels: econometric and statistical modeling with Python. *Proc 9th Python in Science Conference*, 92–96. <https://doi.org/10.25080/Majora-92bf1922-011>
- Spearman C (1910). Correlation calculated from faulty data. *Br J Psychol* 3:271–295. <https://doi.org/10.1111/j.2044-8295.1910.tb00206.x>
- Stam CJ, Nolte G, Daffertshofer A (2007). Phase lag index: assessment of functional connectivity from multi channel EEG and MEG with diminished bias from common sources. *Hum Brain Mapp* 28:1178–1193. <https://doi.org/10.1002/hbm.20346>
- Tadel F, Baillet S, Mosher JC, Pantazis D, Leahy RM (2011). Brainstorm: a user-friendly application for MEG/EEG analysis. *Comput Intell Neurosci* 2011:879716. <https://doi.org/10.1155/2011/879716>
- Van Veen BD, van Drongelen W, Yuchtman M, Suzuki A (1997). Localization of brain electrical activity via linearly constrained minimum variance spatial filtering. *IEEE Trans Biomed Eng* 44:867–880. <https://doi.org/10.1109/10.623056>
- Váša F, Seidlitz J, Romero-Garcia R, Whitaker KJ, Rosenthal G, Vértes PE, et al. (2018). Adolescent tuning of association cortex in human structural brain networks. *Cereb Cortex* 28:281–294. <https://doi.org/10.1093/cercor/bhx249>
- Vinck M, Oostenveld R, van Wingerden M, Battaglia F, Pennartz CMA (2011). An improved index of phase-synchronization for electrophysiological data in the presence of volume-conduction, noise and sample-size bias. *NeuroImage* 55:1548–1565. <https://doi.org/10.1016/j.neuroimage.2011.01.055>
- Vos de Wael R, Benkarim O, Paquola C, Lariviere S, Royer J, Tavakol S, et al. (2020). BrainSpace: a toolbox for the analysis of macroscale gradients in neuroimaging and connectomics datasets. *Commun Biol* 3:103. <https://doi.org/10.1038/s42003-020-0794-7>
- Wen H, Liu Z (2016). Separating fractal and oscillatory components in the power spectrum of neurophysiological signal. *Brain Topogr* 29:13–26. <https://doi.org/10.1007/s10548-015-0448-0>
- Westfall PH, Young SS (1993). *Resampling-Based Multiple Testing: Examples and Methods for p-Value Adjustment*. Wiley, New York.
- What Works Clearinghouse (2022). *What Works Clearinghouse Procedures and Standards Handbook, Version 5.0*, appendix E. U.S. Department of Education, Institute of Education Sciences, National Center for Education Evaluation and Regional Assistance. <https://ies.ed.gov/ncee/wwc/Handbooks>
- Winkler AM, Ridgway GR, Webster MA, Smith SM, Nichols TE (2014). Permutation inference for the general linear model. *NeuroImage* 92:381–397. <https://doi.org/10.1016/j.neuroimage.2014.01.060>
- Yeo BT, Krienen FM, Sepulcre J, Sabuncu MR, Lashkari D, Hollinshead M, et al. (2011). The organization of the human cerebral cortex estimated by intrinsic functional connectivity. *J Neurophysiol* 106:1125–1165. <https://doi.org/10.1152/jn.00338.2011>
- Zalesky A, Fornito A, Bullmore ET (2010). Network-based statistic: identifying differences in brain networks. *NeuroImage* 53:1197–1207. <https://doi.org/10.1016/j.neuroimage.2010.06.041>
