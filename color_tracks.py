import os

dirs = [d for d in dirs if d.startswith('Z')]

for d in dirs:
    os.chdir(d)
    fnames = [f for f in os.listdir('.') if f.startswith('Z')]
    z = float(d.split('Y')[0].replace('Z', ''))
    for fname in fnames:
        os.system('/Volumes/tehom/research/padova_apps/fromHR2mags/fromHR2mags UVbright {} 4 5 1 {}'.format(fname, z))
    os.chdir('..')
