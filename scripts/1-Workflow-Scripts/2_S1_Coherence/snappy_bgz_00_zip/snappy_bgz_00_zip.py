"""-------------------------------------------------------------
    Script Name:   snappy_bgz_00_zip.py
    Description:   SNAP processing chain for VV/VH S1 coherence data - step 00b
    Created By:    Chris Moore
    Date Created:  26.07.22
    Date Modified: 24.05.23
-------------------------------------------------------------"""

# step 0 (Zip): Download zip files from ASF urls

import time
import datetime
import sys
import asf_search

# start timer
tic = time.time()
# get inputs
url_name = sys.argv[1]
wdir = sys.argv[2]

# catch error
try:
    # get user credentials
    asf_user = input("enter NASA Earthdata username: ")
    asf_pass = input("enter NASA Earthdata password: ")

    # download data
    asf_search.download_url(url=url_name, path=wdir,
                            session=asf_search.ASFSession().auth_with_creds(asf_user, asf_pass))

except Exception:
    # sys.exc_clear()
    pass

# output
toc = time.time()
print(str(datetime.timedelta(seconds=round(toc - tic))))
