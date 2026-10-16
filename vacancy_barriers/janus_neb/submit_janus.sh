#!/bin/bash
# Stage and submit the Janus NEB array on the free notchpeak-guest partition.
# Needs a ControlMaster to notchpeak: ssh -fNM -o ControlPath=~/.ssh/cm/%r@%h:%p "$CHPC_USERNAME@notchpeak.chpc.utah.edu"
# MPT_JANUS_SUBMIT=1 submits; otherwise only sbatch --test-only.
set -euo pipefail
cd "$(dirname "$0")"
SSHOPT="-o ControlPath=$HOME/.ssh/cm/%r@%h:%p"
G="$CHPC_USERNAME@notchpeak.chpc.utah.edu"
R() { ssh $SSHOPT "$G" "bash -l -s" <<<"$*"; }
SCR=$(R 'echo /scratch/general/vast/$USER')
D="$SCR/mpt_janus"
R "mkdir -p $D/inputs $D/results $D/pseudos && cd $D/pseudos && \
 for f in mo_pbe_v1.uspp.F.UPF w_pbe_v1.2.uspp.F.UPF s_pbe_v1.4.uspp.F.UPF \
  se_pbe_v1.uspp.F.UPF te_pbe_v1.uspp.F.UPF; do \
  [ -s \$f ] || wget -q https://www.physics.rutgers.edu/gbrv/\$f; done; ls *.UPF | wc -l"
scp -q $SSHOPT run_janus.py janus.sbatch "$G:$D/"
scp -q $SSHOPT inputs/manifest.json "$G:$D/inputs/"
ARGS="--account=owner-guest --partition=notchpeak-guest --qos=notchpeak-guest --export=ALL,MPT_JANUS_ROOT=$D"
R "cd $D && python3 -m py_compile run_janus.py && echo COMPILE_OK && sbatch --test-only $ARGS janus.sbatch"
if [ "${MPT_JANUS_SUBMIT:-0}" = "1" ]; then
  R "cd $D && sbatch $ARGS janus.sbatch" | tee .janus_jobid
  R "squeue -u \$USER -n mpt_janus -o '%.18i %.9P %.2t %.10M %R'"
fi
