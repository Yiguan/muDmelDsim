# set layout for haplotype PNG files
# merge all images to one page for each species

import os
import numpy as np
from PIL import Image
imgs = [x for x in os.listdir() if x.endswith('_snp.png')]

fams = list(set([x[0:7] for x in imgs]))

species = 'dmeN'
s_fams = [x for x in fams if x.startswith(species)]
s_fams = sorted(s_fams, key = lambda kk:int(kk.split('_')[1]))

fam_array = []

for fam in s_fams:
      print(fam)
      ss = [x for x in imgs if x.startswith(fam)]
      maternal = [x for x in ss if x.endswith('maternal_snp.png')]
      paternal = [x for x in ss if x.endswith('paternal_snp.png')]
      mp_order = [*maternal, *paternal]
      tmp = np.concatenate(
            [np.array(Image.open(x)) for x in mp_order], axis=1
      )
      fam_array.append(tmp)

fam_merge = np.concatenate(fam_array, axis=0)

con = Image.fromarray(fam_merge)
con.save(species + '_haps_sortByFamID.pdf')