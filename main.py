# Reuses some info from Stephen Bailey shared on [desi-data 3401] "running fiber assignment on a real target catalog"
import os
import subprocess
from astropy.table import Table, join
import numpy as np
from desitarget.targetmask import desi_mask, bgs_mask, mws_mask, obsmask, obsconditions
import fitsio
import glob
from desisim.quickcat import quickcat
import desimodel.io
import argparse

parser = argparse.ArgumentParser(description='Define parameters')
parser.add_argument('--program', type=str, required=True,
                    help='dark or bright')
parser.add_argument('--size', type=str, required=True,
                    help='small or large')
args = parser.parse_args()

program = args.program
size = args.size

#directories
datadir = "./{}_{}/".format(program, size)
fiberdir = "./{}_{}/fiber_output/".format(program, size)
if not os.path.exists(datadir):
    os.makedirs(datadir)  
if not os.path.exists(fiberdir):
    os.makedirs(fiberdir)  

#filenames
paths = {"targets": "/project/projectdirs/desi/target/catalogs/dr7.1/PR372/", 
         "skies": "/project/projectdirs/desi/target/catalogs/dr7.1/0.22.0/", 
}

names = {"targets": "dr7.1-PR372.fits", "skies":"dr7.1-0.22.0.fits"}

mtlfile = os.path.join(datadir, 'mtl_{}.fits'.format(size, program))
starfile = os.path.join(datadir, 'std_{}.fits'.format(size, program))
targetfile = os.path.join(paths["targets"], "targets-{}".format(names["targets"]))
skyfile = os.path.join(paths["skies"], "skies-{}".format(names["skies"]))
tilefile = os.path.join(datadir, "input_tiles.fits")

# tile selection
if program == "bright":
    if not os.path.exists(tilefile):
        tiles = desimodel.io.load_tiles()
        bright = tiles['PROGRAM']=='BRIGHT'

        if size=="small":
            bright = bright & ((tiles['RA']>12) & (tiles['RA']<38) & (tiles['DEC']<13) & (tiles['DEC']>-13))

        if program=="bright":
            Table(tiles[~bright]).write(tilefile)
        else:
            Table(tiles[bright]).write(tilefile)
print("wrote tiles to {}".format(tilefile))

# target selectionza
if (not os.path.exists(mtlfile)) or (not os.path.exists(starfile)):
    columns=['TARGETID','SUBPRIORITY', 'BRICKID', 'BRICK_OBJID', 'REF_ID',
            'PMRA', 'PMDEC', 'PMRA_IVAR', 'PMDEC_IVAR', 'FLUX_G', 'FLUX_R', 'FLUX_Z',
            'FLUX_W1', 'FLUX_W2', 'FLUX_IVAR_G', 'FLUX_IVAR_R', 'FLUX_IVAR_Z',
            'FLUX_IVAR_W1', 'FLUX_IVAR_W2', 'RA_IVAR', 'DEC_IVAR',
            'EBV', 'MORPHTYPE',
            'MW_TRANSMISSION_G', 'MW_TRANSMISSION_R', 'MW_TRANSMISSION_Z',
            'PHOTSYS',
            'FIBERFLUX_G', 'FIBERFLUX_R', 'FIBERFLUX_Z', 
             'FIBERTOTFLUX_G', 'FIBERTOTFLUX_R', 'FIBERTOTFLUX_Z', 'HPXPIXEL',
        'TARGETID', 'RA', 'DEC', 'SUBPRIORITY', 'BRICKNAME',
        'DESI_TARGET', 'BGS_TARGET', 'MWS_TARGET']

    print('Started reading {}'.format(targetfile))
    targetdata = fitsio.read(targetfile, 'TARGETS', columns=columns)
    if size=="small":
        ii = (targetdata['RA']>10) &  (targetdata['RA']<40) & (targetdata['DEC']<15) & (targetdata['DEC']>-15)
        targetdata = targetdata[ii]
    print('Done reading target data to comput mtl + star')

#compute MTL
if not os.path.exists(mtlfile):
    print('computing mtl')
    import desitarget.mtl
    mtl = desitarget.mtl.make_mtl(targetdata)

    # only include BGS and MWS
    isbgsmws = ((mtl['BGS_TARGET']!=0) | (mtl['MWS_TARGET']!=0))
    if program=="bright":
        mtl = mtl[isbgsmws]
    else:
        mtl = mtl[~isbgsmws]


    mtl.meta['EXTNAME'] = 'MTL'
    # rewrite NUMOBS for BGS targets
    mtl.write(mtlfile)
    

    #print some stats
    print('MWS_TARGETS: {}'.format(np.count_nonzero(mtl['MWS_TARGET']!=0)))
    print('BGS_TARGETS: {}'.format(np.count_nonzero(mtl['BGS_TARGET']!=0)))
    print('DESI_TARGETS: {}'.format(np.count_nonzero(mtl['DESI_TARGET']!=0)))
    print('finished computing mtl')

#standards
if not os.path.exists(starfile):
    std_mask = 0
    for name in ['STD', 'STD_FSTAR', 'STD_WD',
             'STD_FAINT', 'STD_FAINT_BEST',
             'STD_BRIGHT', 'STD_BRIGHT_BEST']:
        if name in desi_mask.names():
            std_mask |= desi_mask[name]

    starstd = (targetdata['DESI_TARGET'] & std_mask) != 0
    stardata = targetdata[starstd]

    if program=="bright":
        obscond = np.int_(np.repeat(obsconditions['BRIGHT'], len(stardata)))
    else:
        obscond = np.int_(np.repeat(obsconditions['DARK']|obsconditions['GRAY'], len(stardata))) 

    stardata = np.lib.recfunctions.append_fields(stardata, 'OBSCONDITIONS', obscond)  
        
    fitsio.write(starfile, stardata, extname='STD')
    print('{} dark standards'.format(np.count_nonzero(stardata)))
    print('Finished with standards')
    
# Running fiberassign
cmd = "fiberassign --mtl {} ".format(mtlfile)
cmd += " --sky {} ".format(skyfile)
cmd += " --stdstar {} ".format(starfile)
cmd += " --fibstatusfile ./fiberstatus.ecsv"
cmd += " --footprint {} ".format(tilefile)
cmd += " --outdir {} ".format(fiberdir)

print(cmd)
print('starting fiberassign')
os.system(cmd)
print('finished fiberassign')
