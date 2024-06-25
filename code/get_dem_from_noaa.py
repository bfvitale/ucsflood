#!/opt/python-3.11.5/bin/python3

import os
import pathlib
import re
import requests
import subprocess
import sys
from typing import List
import urllib
import zipfile

from ucs_constants import *

def get_and_parse_controllers_js() -> List[str]:
    URL = 'https://coast.noaa.gov/slrdata/js/controllers.js'
    r = requests.get(URL)
    assert r.status_code == 200
    # TODO: Does NOAA have an API we can use? This scraping is not robost.
    # TODO: parse the JS, e.g. esprima? Here we just treat it as text.
    # Evaluating remote code might raise security concerns.
    urls = []
    for line in r.text.split('\n'):
        m = re.match(r'\s*"demurl": "([^"]+)"', line)
        if not m: continue
        url = m.group(1)
        assert url.startswith('//')
        url = 'https:' + url  # https instead of http; noaa.gov uses HSTS
        urls.append(url)
    return urls

def remove_zipfile_if_bad(url: str):
    url_parts = urllib.parse.urlparse(url)
    zfn = pathlib.Path(url_parts.path).name
    print("checking", zfn)
    try:
        with zipfile.ZipFile(zfn) as zf:
            zf.testzip()
    except zipfile.BadZipFile:
        print('removing bad zip file', zfn)
        pathlib.Path(zfn).unlink()

def fetch_zip(url: str):
    """Fetch url and put in a file in current directory named with url's
       basename."""
    # Use wget for its --timestamping option. pyCurl or 'requests' might
    # work if we set the HTTP header 'If-Unmodified-Since'.
    subprocess.run(['/usr/bin/wget', '--timestamping', url], check=True)
    
def main(argv):
    pathlib.Path(DEM_INPUT_DIR).mkdir(parents=True, exist_ok=True)
    os.chdir(DEM_INPUT_DIR)

    DEM_URLS_FN = DEM_INPUT_DIR / 'dem-urls'
    path = pathlib.Path(DEM_URLS_FN)
    if not path.exists():
        urls = get_and_parse_controllers_js()
        with open(DEM_URLS_FN, 'w') as dem_url_file:
            dem_url_file.writelines([url + '\n' for url in urls])
    else:
        urls = open(DEM_URLS_FN).read().splitlines()

    for url in urls:
        remove_zipfile_if_bad(url)
        fetch_zip(url)
        # TODO: call zipfile.extractall(). # risky; e.g. zip "bomb"
        # TODO: mv extracted files from subdirs to cwd

    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv))
