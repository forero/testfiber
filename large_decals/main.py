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
# tile selection
tilefile = "./data/dark_gray_tiles.fits"
if not os.path.exists(tilefile):
    tiles = desimodel.io.load_tiles()
    bright = tiles['PROGRAM']=='BRIGHT'
    Table(tiles[~bright]).write(tilefile)

# target selection
paths = {"targets": "/project/projectdirs/desi/target/catalogs/dr7.1/PR372/", 
         "skies": "/project/projectdirs/desi/target/catalogs/dr7.1/0.22.0/", 
}

names = {"targets": "dr7.1-PR372.fits", "skies":"dr7.1-0.22.0.fits"}

targetfile = os.path.join(paths["targets"], "targets-{}".format(names["targets"]))

mtlfile = './data/mtl.fits'
starfile = './data/std.fits'
if (not os.path.exists(mtlfile)) or (not os.path.exists(starfile)):
    targetdata = fitsio.read(targetfile, 'TARGETS')
    print('Done reading target data to comput mtl + star')

#compute MTL
if not os.path.exists(mtlfile):
    print('computing mtl')
    import desitarget.mtl
    mtl = desitarget.mtl.make_mtl(targetdata)

    # only include BGS and MWS
    isbgsmws = (mtl['BGS_TARGET']!=0) | (mtl['MWS_TARGET']!=0)
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

    fitsio.write(starfile, stardata, extname='STD')
    print('{} dark standards'.format(np.count_nonzero(stardata)))
    print('Finished with standards')

output_bright = 'output/'
if not os.path.exists(output_bright):
    os.makedirs(output_bright)
    
skyfile = '/project/projectdirs/desi/target/catalogs/dr7.1/0.22.0/skies-dr7.1-0.22.0.fits'
cmd = "fiberassign --mtl {} ".format(mtlfile)
cmd += " --sky {} ".format(skyfile)
cmd += " --stdstar {} ".format(stdfile)
cmd += " --fibstatusfile ./data/fiberstatus.ecsv"
cmd += " --footprint {} ".format(tilefile)
cmd += " --outdir output/ "
cmd += " --nstarpetal 20 "
cmd += " --nskypetal 80"
print(cmd)
print('starting fiberassign')
os.system(cmd)
print('finished fiberassign')
