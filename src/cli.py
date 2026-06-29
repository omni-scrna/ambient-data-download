import argparse

DATASETS = {
  "hgmm1k": {
    "url": "https://cf.10xgenomics.com/samples/cell-exp/7.0.1/10k_hgmm_3p_gemx/10k_hgmm_3p_gemx_raw_feature_bc_matrix.h5",
  },
  "hgmm5k": {
    "url": "https://cf.10xgenomics.com/samples/cell-exp/3.0.2/hgmm_5k_v3/hgmm_5k_v3_raw_feature_bc_matrix.h5",
  },
  "hgmm6k": {
    "url": "https://cf.10xgenomics.com/samples/cell-exp/3.1.0/5k_hgmm_v3_nextgem/5k_hgmm_v3_nextgem_raw_feature_bc_matrix.h5",
  },
  "hgmm10k": {
    "url": "https://cf.10xgenomics.com/samples/cell-exp/7.0.1/10k_hgmm_3p_gemx/10k_hgmm_3p_gemx_raw_feature_bc_matrix.h5",
  },
  "hgmm12k": {
    "url": "https://cf.10xgenomics.com/samples/cell-exp/3.1.0/hgmm_12k_v3_nextgem/hgmm_12k_v3_nextgem_raw_feature_bc_matrix.h5",
  },
  "hgmm20k": {
    "url": "https://cf.10xgenomics.com/samples/cell-exp/7.0.1/20k_hgmm_3p_gemx/20k_hgmm_3p_gemx_raw_feature_bc_matrix.h5",
  },
}


def build_parser():
  parser = argparse.ArgumentParser(
    description="Download raw 10x H5 matrices for HGMM datasets"
  )
  parser.add_argument(
    "--output_dir", type=str, required=True, help="Output directory"
  )
  parser.add_argument(
    "--name", type=str, required=True, help="Module name/identifier"
  )
  parser.add_argument(
    "--dataset_name",
    type=str,
    required=True,
    choices=list(DATASETS.keys()),
    help="HGMM dataset to download",
  )
  return parser
