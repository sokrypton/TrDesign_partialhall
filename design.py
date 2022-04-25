##########################
## SUPRESS ALL WARNINGS ##
##########################
import warnings, logging, os
warnings.filterwarnings('ignore',category=FutureWarning)
logging.disable(logging.WARNING)
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"

import sys, subprocess
from subprocess import DEVNULL
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as plt

####################
## load libraries ##
####################
from utils import *
from models import *
from to_pdb import *

def main(argv):
  ag = parse_args()
  ag.txt("-------------------------------------------------------------------------------------")
  ag.txt("TrRosetta for Design")
  ag.txt("-------------------------------------------------------------------------------------")
  ag.add(["len=",  "l:"], None,  int,   ["set length for unconstrained design"])
  ag.add(["out=",  "o:"], None,  str,   ["filename prefix for output"])
  ag.add(["num=",  "n:"], 1,     int,   ["number of designs"])
  ag.txt("-------------------------------------------------------------------------------------")
  ag.txt("Backbone Design (if PDB provided)")
  ag.txt("-------------------------------------------------------------------------------------")
  ag.add(["pdb=",  "p:"], None,  str,   ["PDB for fixed backbone design"])
  ag.add(["chain=","c:"], None,  str,   ["specify chain to use"])
  ag.add(["mask=", "m:"], None,  str,   ["set positions to rebuild ['start-end:min-max,...']",
                                         "NOTE: the positions are zero-indexed, first position=0",
                                         "use 'start-end' to specify region to rebuild",
                                         "use 'min-max'   to specify length of indel",
                                         "ex: '10-15:0-5' replace pos. 10-15 with variable-loop of length 0-5",
                                         "ex: '10-15:5'   replace pos. 10-15 with fixed-loop of length 5",
                                         "ex: '10-15'     remove pdb cst. from positions 10 to 15",
                                         "ex: '10'        remove pdb cst. from position 10"])
  ag.add(["sam=",  "s:"], 10,    int,   ["number of samples per design run [if --mask]"])
  ag.txt("-------------------------------------------------------------------------------------")
  ag.txt("Extras")
  ag.txt("-------------------------------------------------------------------------------------")
  ag.add(["aa_weight="],  0.0,   float, ["weight for aa loss"])
  ag.add(["rm_aa="],      None,  str,   ["disable specific amino acids from being sampled",
                                         "ex: 'C' or 'W,Y,F'"])
  ag.add(["save_img"],    False, None,  ["save image of contact map"])
  ag.add(["save_npz"],    False, None,  ["save data for PyRosetta"])
  ag.add(["save_pdb"],    False, None,  ["use magic to quickly generate PDB structure"])
  ag.add(["scwrl"],       False, None,  ["use scwrl to add sidechains [if --save_pdb]"])
  ag.txt("-------------------------------------------------------------------------------------")
  ag.txt("Experimental options")
  ag.txt("-------------------------------------------------------------------------------------")
  ag.add(["pssm_design"], False, None,  ["design a PSSM instead of a single sequence"])
  ag.add(["pssm_mask"],   False, None,  ["use PSSM in regions constrained by PDB"])
  ag.add(["msa_design"],  False, None,  ["design a MSA instead of a single sequence"])
  ag.add(["msa_num="],    1000,  int,   ["number of sequences in MSA"])
  ag.add(["feat_drop="],  0,     float, ["dropout rate for features",
                                         "for --msa_design, we recommend 0.8"])
  ag.add(["cce_cutoff="], None,  float, ["filter cce to CB ≤ x"])
  ag.add(["spike="],      0.0,   float, ["initialize design from PDB seq"])
  ag.add(["diag="],       0.4,   float)
  ag.add(["lid="],        0.3,   float)
  ag.add(["lid_scale="],  18.0,  float)
  ag.add(["aa_ref"],      False, None)
  ag.txt("-------------------------------------------------------------------------------------")
  ag.txt("Optimization settings")
  ag.txt("-------------------------------------------------------------------------------------")
  ag.add(["opt_rate="],   1.0,   float)
  ag.add(["opt_iter="],   100,   int,   ["number of iterations"])
  ag.add(["opt_repeat="], 1,     int)
  ag.add(["opt_adam"],    False, None,  ["use ADAM optimizer"])
  ag.add(["opt_decay"],   False, None,  ["use GD+Decay optimizer"])
  ag.add(["opt_sample"],  False, None,  ["sample from PSSM instead of taking argmax of PSSM"])
  ag.add(["serial"],      False, None,  ["enable approx. serial mode"])
  ag.add(["n_models="],   5,     int,   ["number of TrRosetta models to load into memory"])
  ag.txt("-------------------------------------------------------------------------------------")
  o = ag.parse(argv)

  if o.out is None:
    ag.usage(f"ERROR: Output file not defined --out={o.out}")
  if o.pdb is None and o.len is None:
    ag.usage(f"ERROR: --pdb={o.pdb} or --len={o.len} must be defined")

  # configure params
  s_inputs = {"sample":o.opt_sample, "n_models":o.n_models, "feat_drop":o.feat_drop,
              "pssm_design":o.pssm_design, "pssm_mask":o.pssm_mask, "serial":o.serial,
              "DB_DIR":DB_DIR, "lid":o.lid, "lid_scale":o.lid_scale, "diag":o.diag}

  d_inputs = {"opt_iter":o.opt_iter, "opt_rate":o.opt_rate, "opt_repeat":o.opt_repeat}

  if o.pdb is not None:                   s_inputs["add_pdb"] = True
  if o.pdb is None or o.mask is not None: s_inputs["add_bkg"] = True
  if o.mask is not None:                  s_inputs["add_pdb_mask"] = True
  if o.opt_adam:                          d_inputs["opt_method"] = "ADAM"
  if o.opt_decay:                         d_inputs["opt_method"] = "GD_decay"

  if o.aa_weight > 0:
    if o.aa_ref: s_inputs["add_aa_ref"] = True
    else: s_inputs["add_aa_comp"] = True

  if o.msa_design:
    s_inputs["msa_design"] = True
    d_inputs["num"] = o.msa_num
    if o.feat_drop == 0:
      print("WARNING ========================================================")
      print("WARNING for --msa_design we recommmend including --feat_drop=0.8")
      print("WARNING ========================================================")

  if o.rm_aa is not None:
    rm_aa = [aa_1_N[x] for x in o.rm_aa.split(",")]
    d_inputs["rm_aa"] = rm_aa

  # extract pdb features
  if o.pdb is not None:
    print(f"extracting features from pdb={o.pdb}")
    pdb_out = prep_input(o.pdb, chain=o.chain)
    pdb_feat = pdb_out["feat"][None]
    L = pdb_feat.shape[1]
    desired_feat = np.copy(pdb_feat)

    # spike in the pdb sequence
    seq_start = None
    if o.spike > 0:
      print(f"spiking in the sequence from pdb using {o.spike}")
      pdb_seq = np.eye(21)[AA_to_N(pdb_out["seq"])]
      seq_start = o.spike * pdb_seq[None,None,:,:20]
      if o.msa_design:
        seq_start = np.tile(seq_start,[1,o.msa_num,1,1])

    # christoffer norn's CE8 or CE10
    if o.cce_cutoff is not None:
      pdb_dist = pdb_out["dist_ref"][None]
      desired_feat[pdb_dist > o.cce_cutoff] = 0

    # setup moves (insertions/deletions)
    if o.mask is not None:
      moves, pdb_mask, min_L, max_L = mask_to_moves(o.mask,L)
      print(f"regions constrainted by the pdb {pdb_mask}")

  # background distributions
  if "add_bkg" in s_inputs:
    if o.pdb is None: (min_L,max_L) = (o.len,o.len)
    L_range = np.arange(min_L,max_L+1).tolist()
    print(f"computing background distribution for {L_range[0]}-{L_range[-1]}")
    bkg = get_bkg(L_range,DB_DIR=DB_DIR)

  # setup design
  print("setting up design model")
  model = mk_design_model(**s_inputs)

  # start design
  lines = ["#prefix total_loss loss accuracy length sequence"]
  for n in range(o.num):
    if o.pdb is None:
      # unconstrained backbone design
      inputs = {"bkg":bkg[o.len][None]}
      weights = {"aa":o.aa_weight}
      output = model.design(inputs=inputs, weights=weights, **d_inputs)
      lines.append(save_result(output,f"{o.out}_{n}",o))

    elif o.mask is None:
      # fixed backbone design
      inputs = {"pdb":desired_feat,"I":seq_start}
      weights = {"aa":o.aa_weight}
      output = model.design(inputs=inputs, weights=weights, **d_inputs)
      output["acc"] = get_dist_acc(output["feat"], pdb_feat)[0]
      lines.append(save_result(output,f"{o.out}_{n}",o))

    else:
      # combined-loss backbone design
      print(f"initial design: {n}")
      if seq_start is not None: seq_start *= pdb_mask[...,None]
      inputs = {"pdb":      desired_feat,
                "pdb_mask": pdb_mask[None],
                "I":        seq_start,
                "bkg":      bkg[L][None]}
      weights = {"aa":o.aa_weight,"bkg":0}
      output = model.design(inputs=inputs, weights=weights, **d_inputs)
      seq_start = (output["I"] * pdb_mask[...,None])[None]

      # samplin' time
      for s in range(o.sam):
        move,loop_len = sample_move(moves)
        new_len = L + len(move[0]) - len(move[1])

        # apply move to pdb/mask/msa/bkg
        inputs = {"pdb":      apply_move(desired_feat,move,[1,2]),
                  "pdb_mask": apply_move(pdb_mask[None],move,[1]),
                  "I":        apply_move(seq_start,move,[2]),
                  "bkg":      bkg[new_len][None]}
        weights = {"aa":o.aa_weight}
        output = model.design(inputs=inputs, weights=weights, **d_inputs)
        pdb_feat_ = apply_move(pdb_feat, move, [1,2])
        output["acc"] = get_dist_acc(output["feat"], pdb_feat_, inputs["pdb_mask"])[0]
        lines.append(save_result(output, f"{o.out}_{n}_{s}_{loop_len}", o))

  # save the results
  with open(f"{o.out}.txt","w") as file:
    file.write(f'{o.__dict__}\n')
    for line in lines: file.write(f"{line}\n")

#################################################################################
def do_scwrl(inputs,ouputs):
  subprocess.run([f"{SCWRL}/Scwrl4","-i",inputs,"-o",ouputs], stdout=DEVNULL, stderr=DEVNULL)

def save_result(out, pre, o):
  loss = out["loss"]
  seq = N_to_AA(out["I"][0].argmax(-1))[0]
  if "acc" in out: acc = out["acc"]
  else: acc = None
  total_loss = sum(loss.values())
  loss = str(loss).replace(" ","")
  line = f"{pre} {total_loss} {loss} {acc} {len(seq)} {seq}"
  print(line)
  # save PDB
  if o.save_pdb:
    xyz, dm = feat_to_xyz(out["feat"])
    save_PDB(f"{pre}.pdb", xyz, dm, seq)
    if o.scwrl: do_scwrl(f"{pre}.pdb",f"{pre}.scwrl4.pdb")
  # save IMG (contact map)
  if o.save_img:
    plt.figure(figsize=(5,5))
    plt.imshow(split_feat(out["feat"])["dist"].argmax(-1))
    plt.savefig(f"{pre}.png", bbox_inches='tight')
    plt.close()
  # save NPZ (for pyrosetta)
  if o.save_npz:
    with open(f"{pre}.fas",'w') as fas:
      fas.write(f">{pre}\n{seq}\n")
    feats = split_feat(out["feat"])
    np.savez_compressed(f"{pre}.npz",**feats)
  # save MSA
  if len(out["I"]) > 1:
    with open(f"{pre}.msa.fas",'w') as fas:
      for n,seq in enumerate(N_to_AA(out["I"].argmax(-1))):
        fas.write(f">{pre}_{n}\n{seq}\n")
  # save PSSM
  if o.pssm_design or o.pssm_mask:
    np.savetxt(f"{pre}.pssm.txt",out["I_pssm"][0], fmt='%.6f')
  return line
#################################################################################
if __name__ == "__main__":
   main(sys.argv[1:])
