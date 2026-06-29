#!/usr/bin/env python3
"""Download raw 10x H5 matrices for HGMM ambient RNA benchmark."""

import shutil
import sys
from pathlib import Path
import requests

sys.path.insert(0, str(Path(__file__).parent / "src"))
from cli import DATASETS, build_parser


def download_file(url, dest):
  print(f"  downloading {url}")
  with requests.get(url, stream=True) as r:
    r.raise_for_status()
    with open(dest, "wb") as f:
      shutil.copyfileobj(r.raw, f)
  print(f"  wrote: {dest} ({dest.stat().st_size / 1e6:.1f} MB)")


def main():
  args = build_parser().parse_args()
  print(f"Full command: {' '.join(sys.argv)}")
  for k in ("output_dir", "name", "dataset_name"):
    print(f"  {k}: {getattr(args, k)}")

  out_dir = Path(args.output_dir)
  out_dir.mkdir(parents=True, exist_ok=True)

  info = DATASETS[args.dataset_name]
  out_path = out_dir / f"{args.name}_raw.h5"
  download_file(info["url"], out_path)


if __name__ == "__main__":
  main()
