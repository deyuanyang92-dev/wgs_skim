#!/usr/bin/env bash
# bioconda mitofinder 1.4.1 (py27h9801fc8_1) ships bin/mitofinder without what its install.sh would create:
#   bin/install.sh.ok     - without it: "MitoFinder is not installed ... Aborting" (and exit 0)
#   bin/Mitofinder.config - without it: IOError ... bin/Mitofinder.config
#   bin/lib/gcode         - MitoFinder runs tRNAscan-SE in bin/ with `-g lib/gcode/gcode.invmito`; conda keeps the
#                           codes in lib/tRNAscan-SE/gcode: "FATAL: Could not find lib/gcode/gcode.invmito", which
#                           MitoFinder ignores -> 0 tRNAs (wgs_skim_test pisione: 0, with the link 22 lines of hits)
# The programs it needs are in this environment's bin (MEGAHIT, BLAST, IDBA, metaSPAdes, ARWEN, tRNAscan-SE 2.0);
# MiTFi is not (it needs its own Infernal 1.0.2): with bioconda MitoFinder, mitofinder_trna mitfi -> trnascan
# (rules/common.smk). Snakemake runs this once after creating the environment.
set -euo pipefail
b="${CONDA_PREFIX}/bin/"
touch "${b}install.sh.ok"
mkdir -p "${b}lib"
ln -sfn "${CONDA_PREFIX}/lib/tRNAscan-SE/gcode" "${b}lib/gcode"
cat > "${b}Mitofinder.config" <<CFG
# written by wgs_skim (workflow/envs/mitofinder.post-deploy.sh): the tools of this conda environment
megahitfolder = ${b}
arwenfolder = ${b}
blastfolder = ${b}
idbafolder = ${b}
metaspadesfolder = ${b}
trnascanfolder = ${b}
mitfifolder = default
CFG
