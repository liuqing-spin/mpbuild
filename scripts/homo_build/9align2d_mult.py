from modeller import *
import sys
#inpstr = input("input target_sequence ID:\n")
#inplst = inpstr.split()


log.verbose()
env = Environ()

env.libs.topology.read(file='$(LIB)/top_heav.lib')

# Read aligned structure(s):
aln = Alignment(env)
aln.append(file='template_salign.ali', align_codes='all')
aln_block = len(aln)

target_ID = sys.argv[1]
#target_ID = inplst[0]
target_ali = target_ID + ".ali"
# Read aligned sequence(s):
aln.append(file=target_ali, align_codes=target_ID)

# Structure sensitive variable gap penalty sequence-sequence alignment:
aln.salign(output='', max_gap_length=20,
           gap_function=True,   # to use structure-dependent gap penalty
           alignment_type='PAIRWISE', align_block=aln_block,
           feature_weights=(1., 0., 0., 0., 0., 0.), overhang=0,
           gap_penalties_1d=(-450, 0),
           gap_penalties_2d=(0.35, 1.2, 0.9, 1.2, 0.6, 8.6, 1.2, 0., 0.),
           similarity_flag=True)

out_ali = target_ID + "-mult.ali"
out_pap = target_ID + "-mult.pap"
aln.write(file=out_ali, alignment_format='PIR')
aln.write(file=out_pap, alignment_format='PAP')


print ( "#################################################\n")
print ( "#             Edited By Liu Qing                #\n")
print ( "# University of Science and Technology of China #\n")
print ( "#################################################\n")
