"""
download_real_data.py
=====================

Download the 10x Genomics PBMC 3k dataset for use in the real-data pipeline.

Source
------
- 10x Genomics: 3k PBMCs from a Healthy Donor (Chemistry 1.1)
- Originally from: https://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/
- Cached copy: https://falexwolf.de/data/pbmc3k_raw.h5ad
- Scanpy accessor: ``scanpy.datasets.pbmc3k()``

The downloaded file (~5.6 MB) is saved to ``data/real/pbmc3k_raw.h5ad`` and can
be loaded directly by the main pipeline script.
"""

import os
import urllib.request

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data", "real")
FILENAME = "pbmc3k_raw.h5ad"
URL = "https://falexwolf.de/data/pbmc3k_raw.h5ad"


def main():
    os.makedirs(DATA_DIR, exist_ok=True)
    dest = os.path.join(DATA_DIR, FILENAME)

    if os.path.exists(dest):
        print(f"  File already exists: {dest}")
        return

    print(f"  Downloading PBMC 3k dataset ...")
    print(f"  URL: {URL}")
    urllib.request.urlretrieve(URL, dest)
    size_mb = os.path.getsize(dest) / 1e6
    print(f"  Saved to {dest} ({size_mb:.1f} MB)")


if __name__ == "__main__":
    main()
