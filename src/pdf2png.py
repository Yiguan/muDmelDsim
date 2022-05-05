# convert PDF format to PNF format
# Assume only ONE page in PDF file


from pdf2image import convert_from_path
import os
import re

rr = re.compile("^d..._[0-9]{2}$")

mm = [x for x in os.listdir()]
dd = list(filter(rr.match, mm))

for fam in dd:
      pdfs = os.listdir(fam)
      for pdf in pdfs:
          images = convert_from_path('./' + fam + '/' + pdf)
          outname = fam + '_' + pdf.replace('pdf','png')
          print(outname)
          images[0].save(outname, 'PNG')  
