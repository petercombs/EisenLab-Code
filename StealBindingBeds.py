import urllib
from os import path
from glob import glob
import time

url_base = 'http://bdtnp.lbl.gov/Fly-Net/archives/chipper/Post_Processing/{chip}/{chip}-sym-25_primary_peaks.bed'
file_base = 'Reference/peaks-25-dm2/{chip}-sym-25_peaks.bed'
files = glob('prereqs/BDTNP_in_vivo_binding_Release.2.1/Supplemental_Tables/*.txt')

for file in files:
    chip = path.basename(file).split('-')[0]
    print chip
    fh = urllib.urlopen(url_base.format(chip=chip), )
    fh.readline() # The track line makes liftOver unhappy
    of = open(file_base.format(chip=chip), 'w')
    for line in fh:
        data = line.split()
        data[2] = str(int(data[2]) + 1)
        of.write('\t'.join(data) + '\n')
        #of.write(line.replace('\t',' '))
        #of.write(line)
    time.sleep(1)
