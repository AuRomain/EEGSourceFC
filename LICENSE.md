# Licence

SPDX-License-Identifier: GPL-3.0-or-later

This pipeline is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your option)
any later version. It is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
Public License for more details: <https://www.gnu.org/licenses/gpl-3.0.html>.

The full text of the licence is in `LICENSE`: the GNU General Public
License, version 3, of 29 June 2007, as the Free Software Foundation
published it.

## Third-party code and data included here

| Files | Source | Licence |
|---|---|---|
| `1- EEG source connectivity/+ROInets/` (6 files) | ROInets, OHBA, University of Oxford (G. Colclough; <https://github.com/OHBA-analysis/MEG-ROI-nets>), as copied in var-EEG-FC; `cols.m`, `rows.m` and `scale_cols.m` are its small helpers | GPL-3.0-or-later |
| `1- EEG source connectivity/BS_to_ft_tess.m` | Brainstorm (F. Tadel), `out_fieldtrip_tess.m`, renamed and changed (reads the surface files with `load`) | GPL-3.0 |
| `1- EEG source connectivity/resting_iaf_pieces.m` | adapted from `restingIAF.m` of restingIAF (A. W. Corcoran; <https://github.com/corcorana/restingIAF>): each channel's spectrum is averaged over the continuous pieces | GPL-3.0 |
| `3- Analysis/fun/NBS/NBSstats.m`, `get_components.m` | NBS 1.2, Network-Based Statistic (A. Zalesky; Zalesky et al. 2010; `get_components.m` by J. Goñi), unchanged | GPL-3.0-or-later |
| `3- Analysis/fun/cluster_perm_test.m` | adapted from FieldTrip's `private/clusterstat.m` and `findcluster.m` (commit 0beabe6; <https://github.com/fieldtrip/fieldtrip>) | GPL-3.0-or-later |
| `tools/brainstorm_base64_952475b.patch` | the changes of Brainstorm commit 952475b (<https://github.com/brainstorm-tools/brainstorm3>), applied to Brainstorm 3.250107 by `tools/fetch_toolboxes.sh` | GPL-3.0 (as Brainstorm) |
| `1- EEG source connectivity/+wmne/` | var-EEG-FC (S. Yassine, A. Mheich), after Brainstorm's wMNE | as Brainstorm, GPL-3.0 |
| `1- EEG source connectivity/inputs/tess_*.mat`, `scout_Desikan-Killiany_68.mat`, `channel_*_bs.mat` | Brainstorm ICBM152 template anatomy, atlas and caps (MNI ICBM152 2009 non-linear; Desikan et al. 2006) | Brainstorm GPL; ICBM152 template © McConnell Brain Imaging Centre, free use with attribution |
| the other atlases stored in `inputs/tess_cortex_pial_low.mat` (DKT, Destrieux, Brodmann, Brainnetome, HCP-MMP1, PALS-B12, Schaefer 2018, vcAtlas, OASIS), as Brainstorm distributes the file | their authors; HCP-MMP1 (Glasser et al. 2016) under the HCP data use terms | not used by the pipeline (it reads the surface, the normals and the sphere; the atlas comes from the scout file). Removing them would change the file's fingerprint, hence the key of every stage 1 output |
| `1- EEG source connectivity/inputs/ftLeadfield_*.mat` | computed with OpenMEEG (BSD-like CeCILL-B) on the files above | as above |
| `3- Analysis/DK68_networks.csv` | region centroids from the atlas above; Yeo et al. (2011) networks by overlap with Schaefer et al. (2018) | derived data, cite the papers |

The other toolboxes the pipeline calls (FieldTrip, Brainstorm, restingIAF;
README, section 1) are not included; each keeps the licence it is
distributed under.
