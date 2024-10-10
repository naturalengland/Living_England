"""-------------------------------------------------------------
    Script Name:   snappy_bgz_END_clean.py
    Description:   SNAP processing chain for VV/VH S1 coherence data - step END
    Created By:    Chris Moore
    Date Created:  19.02.21
    Date Modified: 30.05.23
-------------------------------------------------------------"""

# step END (Clean): Delete processing files

import time
import datetime
import sys
import os
import glob
import shutil
import zipfile

# start timer
tic = time.time()
# get inputs
wdir = sys.argv[1]
min_btemp = sys.argv[2]
max_btemp = sys.argv[3]
start_str = sys.argv[4]
end_str = sys.argv[5]
io_base = sys.argv[6]

# get dates
dates = [s1_zip.split("_IW_SLC__1SDV_")[1][:8] for s1_zip in glob.glob(wdir + '\\*.zip')]
# sort unique dates
dates = [*set(dates)]
dates.sort()
# get date pairs within temporal baseline limits
ifgs = []
for d1 in range(len(dates)):
    for d2 in range(d1 + 1, len(dates)):
        btemp = abs((datetime.datetime.strptime(dates[d1], "%Y%m%d") -
                     datetime.datetime.strptime(dates[d2], "%Y%m%d")).days)
        if int(min_btemp) <= btemp <= int(max_btemp):
            ifgs.append(dates[d1] + '_' + dates[d2])

# get steps
clean_steps = range(int(start_str), int(end_str))  # not cleaning last step

# clear local snap cache directories
snap_dir = 'C:\\Users\\ne.user\\.snap\\var'
for sub_dir in ['log', 'cache\\temp', 'quicklooks_cache']:
    for tmp_file in glob.glob(snap_dir + '\\' + sub_dir + '\\*'):
        os.remove(tmp_file)
        shutil.rmtree(tmp_file, ignore_errors=True)
# clear local temp directory
tmp_dir = 'C:\\Users\\ne.user\\AppData\\Local\\Temp\\snap-ne.user'
for tmp_file in glob.glob(tmp_dir + '\\*'):
    shutil.rmtree(tmp_file, ignore_errors=True)

# -------------------- step 0 -------------------- #

if 0 in clean_steps:
    for s1_zip in glob.glob(wdir + '\\*.zip'):
        os.remove(s1_zip)
        # replace with placeholder
        if int(end_str) < 14:
            zipf = zipfile.ZipFile(wdir + '\\' + s1_zip[44::], 'w', zipfile.ZIP_DEFLATED)
            zipf.close()

# -------------------- step 1 -------------------- #

if 1 in clean_steps:
    # loop through dates
    for d1 in dates:
        # get slc
        io_slc = wdir + '\\' + d1 + '_SLC'
        # if have slc
        if os.path.exists(io_slc + '.dim'):
            os.remove(io_slc + '.dim')
            shutil.rmtree(io_slc + '.data', ignore_errors=True)

# -------------------- step 2 -------------------- #

if 2 in clean_steps:
    # loop through dates
    for d1 in dates:
        # loop through subswaths
        for iw in [1, 2, 3]:
            # get slc
            io_IWslc = wdir + '\\' + d1 + '_IW' + str(iw) + '_SLC'
            # if have slc
            if os.path.exists(io_IWslc + '.dim'):
                os.remove(io_IWslc + '.dim')
                shutil.rmtree(io_IWslc + '.data', ignore_errors=True)

# -------------------- step 3 -------------------- #

if 3 in clean_steps:
    # loop through dates
    for d1 in dates:
        # loop through subswaths
        for iw in [1, 2, 3]:
            # get slc
            io_slc_orb = wdir + '\\' + d1 + '_IW' + str(iw) + '_SLC_Orb'
            # if have slc
            if os.path.exists(io_slc_orb + '.dim'):
                os.remove(io_slc_orb + '.dim')
                shutil.rmtree(io_slc_orb + '.data', ignore_errors=True)

# -------------------- step 4 -------------------- #

if 4 in clean_steps:
    # loop through dates
    for ifg in ifgs:
        # loop through subswaths
        for iw in [1, 2, 3]:
            # get rslc
            io_rslc = wdir + '\\' + ifg + '_IW' + str(iw) + '_RSLC'
            # if have rslc
            if os.path.exists(io_rslc + '.dim'):
                os.remove(io_rslc + '.dim')
                shutil.rmtree(io_rslc + '.data', ignore_errors=True)

# -------------------- step 5 -------------------- #

if 5 in clean_steps:
    # loop through dates
    for ifg in ifgs:
        # loop through subswaths
        for iw in [1, 2, 3]:
            # get rslc
            io_rslc_esd = wdir + '\\' + ifg + '_IW' + str(iw) + '_RSLC_esd'
            # if have rslc
            if os.path.exists(io_rslc_esd + '.dim'):
                os.remove(io_rslc_esd + '.dim')
                shutil.rmtree(io_rslc_esd + '.data', ignore_errors=True)

# -------------------- step 6 -------------------- #

if 6 in clean_steps:
    # loop through dates
    for ifg in ifgs:
        # loop through subswaths
        for iw in [1, 2, 3]:
            # get cc
            io_cc = wdir + '\\' + ifg + '_IW' + str(iw) + '_cc'
            # if have cc
            if os.path.exists(io_cc + '.dim'):
                os.remove(io_cc + '.dim')
                shutil.rmtree(io_cc + '.data', ignore_errors=True)

# -------------------- step 7 -------------------- #

if 7 in clean_steps:
    # loop through dates
    for ifg in ifgs:
        # loop through subswaths
        for iw in [1, 2, 3]:
            # get cc
            io_cc_deb = wdir + '\\' + ifg + '_IW' + str(iw) + '_cc_deb'
            # if have cc
            if os.path.exists(io_cc_deb + '.dim'):
                os.remove(io_cc_deb + '.dim')
                shutil.rmtree(io_cc_deb + '.data', ignore_errors=True)

# -------------------- step 8 -------------------- #

if 8 in clean_steps:
    # loop through dates
    for ifg in ifgs:
        # get cc
        io_mrg = wdir + '\\' + ifg + '_cc'
        # if have cc
        if os.path.exists(io_mrg + '.dim'):
            os.remove(io_mrg + '.dim')
            shutil.rmtree(io_mrg + '.data', ignore_errors=True)

# -------------------- step 9 -------------------- #

if 9 in clean_steps:
    # loop through dates
    for ifg in ifgs:
        # get cc
        io_mrg_ml = wdir + '\\' + ifg + '_cc_ML'
        # if have cc
        if os.path.exists(io_mrg_ml + '.dim'):
            os.remove(io_mrg_ml + '.dim')
            shutil.rmtree(io_mrg_ml + '.data', ignore_errors=True)

# -------------------- step 10 -------------------- #

if 10 in clean_steps:
    # loop through dates
    for ifg in ifgs:
        # get cc
        io_mrg_ml_tc = wdir + '\\' + ifg + '_cc_ML_TC'
        # if have cc
        if os.path.exists(io_mrg_ml_tc + '.dim'):
            os.remove(io_mrg_ml_tc + '.dim')
            shutil.rmtree(io_mrg_ml_tc + '.data', ignore_errors=True)

# -------------------- step 11 -------------------- #

# if 11 in clean_steps:
#     # loop through dates
#     for ifg in ifgs:
#         # get cc
#         io_mrg_ml_tc_msk = wdir + '\\' + ifg + '_cc_ML_TC_msk'
#         # if have cc
#         if os.path.exists(io_mrg_ml_tc_msk + '.dim'):
#             os.remove(io_mrg_ml_tc_msk + '.dim')
#             shutil.rmtree(io_mrg_ml_tc_msk + '.data', ignore_errors=True)

# -------------------- step 12 -------------------- #

if 12 in clean_steps:
    # get stack
    io_stack = io_base + '_stack'
    # if have stack
    if os.path.exists(io_stack + '.dim'):
        os.remove(io_stack + '.dim')
        shutil.rmtree(io_stack + '.data', ignore_errors=True)

# -------------------- step 13 -------------------- #

if 13 in clean_steps:
    # get stack
    io_stack = io_base + '_stack'
    # loop through polarisations
    for pol in ['VH', 'VV']:
        # if have stack
        if os.path.exists(io_stack + '_' + pol + '_BM' + '.dim'):
            os.remove(io_stack + '_' + pol + '_BM' + '.dim')
            shutil.rmtree(io_stack + '_' + pol + '_BM' + '.data', ignore_errors=True)

# -------------------- end -------------------- #

toc = time.time()
print(str(datetime.timedelta(seconds=round(toc - tic))))
