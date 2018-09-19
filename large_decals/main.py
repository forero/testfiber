# Reuses some info from Stephen Bailey shared on [desi-data 3401] "running fiber assignment on a real target catalog"
import os
import subprocess
from astropy.table import Table, join
import numpy as np
from desitarget.targetmask import desi_mask, bgs_mask, mws_mask, obsmask, obsconditions
import fitsio
import glob
from desisim.quickcat import quickcat

# target selection
paths = {"targets": "/project/projectdirs/desi/target/catalogs/dr7.1/PR372/", 
         "skies": "/project/projectdirs/desi/target/catalogs/dr7.1/0.22.0/", 
}

os.environ['DECALS_PATH'] = '/project/projectdirs/desi/target/catalogs/dr7.1/0.23.0/'
targetfile = os.path.join(os.environ['DECALS_PATH'], "targets-dr7.1-0.23.0.fits")

columns = [
    'TARGETID', 'RA', 'DEC', 'SUBPRIORITY', 'BRICKNAME',
    'DESI_TARGET', 'BGS_TARGET', 'MWS_TARGET',
]

std_mask = 0
for name in ['STD', 'STD_FSTAR', 'STD_WD',
             'STD_FAINT', 'STD_FAINT_BEST',
             'STD_BRIGHT', 'STD_BRIGHT_BEST']:
    if name in desi_mask.names():
        std_mask |= desi_mask[name]

#truth file
truthfile = "data/truth.fits"
if not os.path.exists(truthfile):
    import desitarget.mock.mockmaker as mb
    from desitarget.targetmask import desi_mask, bgs_mask, mws_mask

    targets = fitsio.read(targetfile, 'TARGETS', columns=columns)


    colnames = list(targets.dtype.names)
    print(colnames)
    nobj = len(targets)
    truth = mb.empty_truth_table(nobj=nobj)
    print(truth.keys())

    for k in colnames:
        if k in truth.keys():
            print(k)
            truth[k][:] = targets[k][:]

    nothing = '          '
    truth['TEMPLATESUBTYPE'] = np.repeat(nothing, nobj)



    dict_truespectype = {'BGS_ANY':'GALAXY', 'ELG':'GALAXY', 'LRG':'GALAXY', 'QSO':'QSO', 
                         'STD':'STAR'}
    dict_truetemplatetype = {'BGS_ANY':'BGS', 'ELG':'ELG', 'LRG':'LRG', 'QSO':'QSO', 
                        'STD':'STAR'}
    masks = {'BGS_ANY':desi_mask.mask('BGS_ANY'), 'ELG':desi_mask.mask('ELG'), 
             'LRG':desi_mask.mask('LRG'), 'QSO':desi_mask.mask('QSO'), 'STD':std_mask}

    for m in masks.keys():
        istype = (targets['DESI_TARGET'] & masks[m])!=0
        print(np.count_nonzero(istype))
        truth['TRUESPECTYPE'][istype] = np.repeat(dict_truespectype[m], np.count_nonzero(istype))
        truth['TEMPLATETYPE'][istype] = np.repeat(dict_truetemplatetype[m], np.count_nonzero(istype))
        truth['MOCKID'][istype] = targets['TARGETID'][istype]

    del targets
    print('writing truth')
    truth.write(truthfile, overwrite=True)
    print('done truth')



mtlfile = 'data/mtl.fits'
starfile = 'data/std.fits'
if (not os.path.exists(mtlfile)) or (not os.path.exists(starfile)):
    targetdata = fitsio.read(targetfile, 'TARGETS', columns=columns)
    print('Done reading target data to comput mtl + star')

#compute MTL
if not os.path.exists(mtlfile):
    print('computing mtl')
    import desitarget.mtl
    mtl = desitarget.mtl.make_mtl(targetdata)

    # only include BGS and MWS
    isbgsmws = (mtl['BGS_TARGET']!=0) | (mtl['MWS_TARGET']!=0)
    mtl = mtl[isbgsmws]

    mtl.meta['EXTNAME'] = 'MTL'
    # rewrite NUMOBS for BGS targets
    ii = mtl['BGS_TARGET']!=0
    mtl['NUMOBS_MORE'][ii] = 4

    mtl.write(mtlfile)
    

    #print some stats
    print('MWS_TARGETS: {}'.format(np.count_nonzero(mtl['MWS_TARGET']!=0)))
    print('BGS_TARGETS: {}'.format(np.count_nonzero(mtl['BGS_TARGET']!=0)))
    print('DESI_TARGETS: {}'.format(np.count_nonzero(mtl['DESI_TARGET']!=0)))
    print('finished computing mtl')

#standards
if not os.path.exists(starfile):
    starstd = (targetdata['DESI_TARGET'] & std_mask) != 0
    stardata = targetdata[starstd]
    obscond = np.int_(np.repeat(4, len(stardata))) # 4 represents bright time
    stardata = np.lib.recfunctions.append_fields(
    stardata, 'OBSCONDITIONS', obscond)    
    fitsio.write(starfile, stardata, extname='STD')
    print('{} dark standards'.format(np.count_nonzero(stardata)))
    print('Finished with standards')

output_bright = 'output/'
if not os.path.exists(output_bright):
    os.makedirs(output_bright)
    
skyfile = '/project/projectdirs/desi/target/catalogs/dr7.1/0.22.0/skies-dr7.1-0.22.0.fits'
cmd = "fiberassign --mtl data/mtl.fits "
cmd += " --sky {} ".format(skyfile)
cmd += " --stdstar data/std.fits "
cmd += " --fibstatusfile data/fiberstatus.ecsv"
cmd += " --footprint data/bgs_sv.fits " 
cmd += " --outdir output/ "
cmd += " --nstarpetal 20 "
cmd += " --nskypetal 80"
print(cmd)
print('starting fiberassign')
os.system(cmd)
#subprocess.call(cmd.split())
print('finished fiberassign')
