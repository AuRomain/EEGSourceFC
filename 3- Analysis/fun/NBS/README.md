# NBS 1.2: the two functions the pipeline calls

`NBSstats.m` and `get_components.m` of the Network-Based Statistic toolbox,
version 1.2 (Zalesky, Fornito & Bullmore 2010, NeuroImage 53:1197–1207;
<https://www.nitrc.org/projects/nbs>), copied unchanged. NBS is free software
under the GNU General Public License, version 3 or later (Copyright (C) 2012,
A. Zalesky; `help/about.txt` of the toolbox). They are included because NBS
is downloaded from NITRC behind a login, so it cannot be fetched
automatically; `toolbox_check` compares them with the validated version
like any toolbox. `compute_NBS_analysis` passes its own test statistics, so
NBSglm and the rest of the toolbox are not needed.
