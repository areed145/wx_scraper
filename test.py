# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 21:17:49 2016

@author: areed145
"""

# Include the Dropbox SDK
import os
from dropbox.client import DropboxClient

access_token = '9-dK9AG5fnkAAAAAAAAySl9sqq6pdHP84Fx5jhGn8adduEqGcmI_xvUS33fwhVO5'
local_directory = '/Users/areed145/Dropbox/GitHub/wx_scraper/plots/KCABAKER38'
dropbox_destination = '/plots/KCABAKER38'

client = DropboxClient(access_token)

# enumerate local files recursively
for root, dirs, files in os.walk(local_directory):

    for filename in files:

        # construct the full local path
        local_path = os.path.join(root, filename)

        # construct the full Dropbox path
        relative_path = os.path.relpath(local_path, local_directory)
        dropbox_path = os.path.join(dropbox_destination, relative_path)

        # upload the file
        with open(local_path, 'rb') as f:
            client.put_file(dropbox_path, f, overwrite=True)