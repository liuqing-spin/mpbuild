# Illustrates the SALIGN multiple structure/sequence alignment
#inpstr = input("input template PDB and chain IDs, seperated by whitespace:\n")
#inplst = inpstr.split()
#template_1 = str(inplst[0])
#chainID_1 = str(inplst[1])
#template_2 = str(inplst[2])
#chainID_2 = str(inplst[3])
#template_3 = str(inplst[4])
#chainID_3 = str(inplst[5])

import sys

template_1 = sys.argv[1]
chainID_1 = sys.argv[2]
template_2 = sys.argv[3]
chainID_2 = sys.argv[4]


from modeller import *

log.verbose()
env = Environ()
env.io.atom_files_directory = ['.', '../atom_files/']

aln = Alignment(env)
#for (code, chain) in ((template_1, chainID_1), (template_2, chainID_2), (template_3, chainID_3)):
for (code, chain) in ((template_1, chainID_1),  (template_2, chainID_2)):
    mdl = Model(env, file=code, model_segment=('FIRST:'+chain, 'LAST:'+chain))
    aln.append_model(mdl, atom_files=code, align_codes=code+chain)

for (weights, write_fit, whole) in (((1., 0., 0., 0., 1., 0.), False, True),
                                    #((1., 0.5, 1., 1., 1., 0.), False, True),
                                    ((1., 1., 1., 1., 1., 0.), True, False)):
    aln.salign(rms_cutoff=3.5, normalize_pp_scores=False,
               rr_file='$(LIB)/as1.sim.mat', overhang=30,
               gap_penalties_1d=(-450, -50),
               gap_penalties_3d=(0, 3), gap_gap_score=0, gap_residue_score=0,
               dendrogram_file='template_salign.tree',
               alignment_type='tree', # If 'progresive', the tree is not
                                      # computed and all structues will be
                                      # aligned sequentially to the first
               feature_weights=weights, # For a multiple sequence alignment only
                                        # the first feature needs to be non-zero
               improve_alignment=True, fit=True, write_fit=write_fit,
               write_whole_pdb=whole, output='ALIGNMENT QUALITY')

aln.write(file='template_salign.pap', alignment_format='PAP')
aln.write(file='template_salign.ali', alignment_format='PIR')

aln.salign(rms_cutoff=1.0, normalize_pp_scores=False,
           rr_file='$(LIB)/as1.sim.mat', overhang=30,
           gap_penalties_1d=(-450, -50), gap_penalties_3d=(0, 3),
           gap_gap_score=0, gap_residue_score=0, dendrogram_file='template_salign_2.tree',
           alignment_type='progressive', feature_weights=[0]*6,
           improve_alignment=False, fit=False, write_fit=True,
           write_whole_pdb=False, output='QUALITY')

print ( "#################################################\n")
print ( "#             Edited By Liu Qing                #\n")
print ( "# University of Science and Technology of China #\n")
print ( "#################################################\n")
