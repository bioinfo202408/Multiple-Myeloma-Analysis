#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import subprocess

def auto_install(package, import_name=None):
    try:
        if import_name is None:
            import_name = package
        __import__(import_name)
    except ImportError:
        print(f"Package '{package}' not found. Installing...")
        subprocess.check_call([sys.executable, "-m", "pip", "install", package])
        __import__(import_name)

# Standard libraries do not require automatic installation
import os, sys as _sys, json, argparse, logging, warnings, time
from datetime import datetime

# Third-party libraries requiring automatic installation
auto_install("pyyaml", "yaml")
auto_install("numpy")
auto_install("pandas")
auto_install("joblib")
auto_install("psutil")
auto_install("scikit-learn", "sklearn")

# Standard imports
import yaml
import numpy as np
import pandas as pd
import joblib
import psutil
from sklearn.base import clone
from sklearn.utils import resample
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.svm import LinearSVC
from sklearn.linear_model import LogisticRegression

warnings.filterwarnings("ignore")

def setup_hpc_env(n_jobs: int):
    physical = psutil.cpu_count(logical=False) or os.cpu_count() or 8
    # Assign 1-2 BLAS threads to each parallel task to avoid oversubscription
    blas_threads = "2" if (n_jobs and n_jobs <= physical // 2) else "1"
    for var in ["OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS"]:
        os.environ.setdefault(var, blas_threads)

# ======================= FWSE Core ==========================
class FWSE:
    def __init__(self, filter_estimators, wrapper_estimators,
                 n_bootstraps=10, pruning_factor=0.5, n_jobs=-1, random_state=0):
        self.filter_estimators = filter_estimators
        self.wrapper_estimators = wrapper_estimators
        self.n_bootstraps = n_bootstraps
        self.pruning_factor = pruning_factor
        self.n_jobs = n_jobs
        self.random_state = random_state
        self.ranking_ = None
        
    @staticmethod
    def _get_ranking_from_estimator(estimator, X, y):
        est = clone(estimator)
        est.fit(X, y)
        if hasattr(est, "scores_"):
            imp = est.scores_
            return np.argsort(np.argsort(-imp))
        if hasattr(est, "coef_"):
            imp = np.abs(est.coef_).ravel()
            return np.argsort(np.argsort(-imp))
        if hasattr(est, "feature_importances_"):
            imp = est.feature_importances_
            return np.argsort(np.argsort(-imp))
        if hasattr(est, "ranking_"):
            return est.ranking_ - 1
        raise ValueError("Estimator must contain 'scores_', 'coef_', 'feature_importances_', or 'ranking_'.")

    @staticmethod
    def _aggregate_ranks(ranks):
        return np.argsort(np.argsort(np.sum(ranks, axis=0)))

    def _run_bootstrap_stage(self, X, y, estimators, rng):
        def fit_single_task(estimator, seed):
            # Note: Xb and yb must be purely numeric
            Xb, yb = resample(X, y, random_state=seed)
            return self._get_ranking_from_estimator(estimator, Xb, yb)

        tasks = []
        for est in estimators:
            seeds = rng.randint(np.iinfo(np.int32).max, size=self.n_bootstraps)
            for s in seeds:
                tasks.append(joblib.delayed(fit_single_task)(est, s))

        all_ranks = joblib.Parallel(n_jobs=self.n_jobs, backend="threading")(tasks)

        agg_per_est = []
        m = len(estimators)
        for i in range(m):
            beg = i * self.n_bootstraps
            end = beg + self.n_bootstraps
            agg_per_est.append(self._aggregate_ranks(all_ranks[beg:end]))
            logging.info(f" > Aggregated: {estimators[i].__class__.__name__}")
        return self._aggregate_ranks(agg_per_est)

    def fit(self, X, y):
        rng = np.random.RandomState(self.random_state)
        n_features = X.shape[1]

        logging.info("--- [1/3] Filter stage ---")
        agg_filter = self._run_bootstrap_stage(X, y, self.filter_estimators, rng)

        logging.info("--- [2/3] Pruning ---")
        keep = max(1, int(n_features * (1 - self.pruning_factor)))
        idx_keep = np.argsort(agg_filter)[:keep]
        X_kept = X[:, idx_keep]
        logging.info(f"Retained {keep}/{n_features} ({(keep/n_features)*100:.3f}%) features after pruning for the Wrapper stage.")

        logging.info("--- [3/3] Wrapper stage ---")
        agg_wrap = self._run_bootstrap_stage(X_kept, y, self.wrapper_estimators, rng)

        logging.info("--- Generating final ranking ---")
        final = np.full(n_features, n_features, dtype=int)
        final[idx_keep] = agg_wrap
        self.ranking_ = final
        return self

# ======================= Helper Functions ==========================
def snr_score_func(X, y):
    mu1 = X[y==1].mean(axis=0); mu0 = X[y==0].mean(axis=0)
    sd1 = X[y==1].std(axis=0);  sd0 = X[y==0].std(axis=0)
    return np.nan_to_num(np.abs(mu1 - mu0) / (sd1 + sd0 + 1e-9))

class SelectKBestSNR(SelectKBest):
    def __init__(self): super().__init__(score_func=snr_score_func, k="all")

def _resolve_brain_region(meta, meta_cfg, cfg, brain_region):
    if not brain_region: return meta
    col = meta_cfg["brain_region_column"]
    mask = (meta[col] == brain_region)
    if mask.any(): return meta[mask]
    abr_map = cfg.get("brain_region_abbreviations", {})
    rev = {v: k for k, v in abr_map.items()}
    full = rev.get(brain_region, None)
    if full is None:
        target_abr = abr_map.get(brain_region, None)
        if target_abr is None: raise ValueError(f"Brain region '{brain_region}' was not found in the metadata or abbreviation table.")
        mask2 = (meta[col] == brain_region) | (meta[col] == target_abr)
        if not mask2.any(): raise ValueError(f"No samples found for brain region '{brain_region}'.")
        return meta[mask2]
    mask3 = (meta[col] == full)
    if not mask3.any(): raise ValueError(f"No samples found for brain region '{brain_region}' (full name: {full}).")
    return meta[mask3]

def load_data_from_config(config, data_type_key="mRNA"):
    pheno_cfg = config["phenotypes"]
    expr_cfg = config["expression_data"][data_type_key]

    # Read metadata
    meta_cfg = config["metadata"]  # Define meta_cfg by obtaining the metadata dictionary from the configuration
    meta = pd.read_csv(meta_cfg["file"], sep=",", dtype=str, encoding="utf-8-sig")
    meta = meta.set_index(meta_cfg["sample_id_column"])
    meta.index = meta.index.str.strip()
    logging.info(f"Number of metadata samples: {meta.shape[0]}, examples: {meta.index.tolist()[:5]}")

    print("First five rows of metadata:")
    print(meta.head())         # Directly print the first five rows
    logging.info(f"\nFirst five rows of metadata:\n{meta.head()}")  # Also write to the log

    meta["y"] = (meta["group"].str.strip() == pheno_cfg["disease"]).astype(int)
    logging.info(f"Metadata loaded. Case='{pheno_cfg['disease']}', control='{pheno_cfg['control']}'")

    # Read expression matrix
    expr = pd.read_csv(expr_cfg["file"], sep=",", dtype=str)
    logging.info(f"Initial expression matrix shape: {expr.shape}")
    logging.info(f"Example expression matrix column names: {expr.columns.tolist()[:10]}")
    logging.info(f"First rows of expression matrix:\n{expr.head()}")

    expr.index = expr["gene_name"]
    expr = expr.drop(columns=["gene_name"])
    expr.columns = expr.columns.str.strip().str.rstrip('/') 
    logging.info(f"Example column names after removing gene_name from expression matrix: {expr.columns.tolist()[:10]}")
    logging.info(f"Example expression matrix indices: {expr.index.tolist()[:5]}")

    expr = expr.astype(np.float32).T
    expr.index.name = meta_cfg["sample_id_column"]
    expr.index = expr.index.str.strip()
    logging.info(f"Example indices after transposing expression matrix: {expr.index.tolist()[:5]}")
    logging.info(f"First rows of expression matrix:\n{expr.head()}")
    
    logging.info(f"Expression matrix shape after transposition: {expr.shape}")

    # ==== Troubleshoot sample ID alignment ====
    logging.info("Troubleshooting sample alignment:")
    meta_ids = set(meta.index)
    expr_ids = set(expr.index)
    intersect = meta_ids & expr_ids
    only_in_meta = meta_ids - expr_ids
    only_in_expr = expr_ids - meta_ids

    logging.info(f"Number of metadata samples: {len(meta_ids)}")
    logging.info(f"Number of expression matrix samples: {len(expr_ids)}")
    logging.info(f"Number of shared samples: {len(intersect)}")
    logging.info(f"Samples only present in metadata (examples): {list(only_in_meta)[:5]}")
    logging.info(f"Samples only present in expression matrix (examples): {list(only_in_expr)[:5]}")
    logging.info(f"Shared samples (examples): {list(intersect)[:5]}")

    # Print repr to help detect hidden character issues
    logging.info("Metadata sample ID repr (first 5): %s", [repr(x) for x in list(meta.index)[:5]])
    logging.info("Expression matrix sample ID repr (first 5): %s", [repr(x) for x in list(expr.index)[:5]])

    df = expr.join(meta, how="inner")
    logging.info(f"Number of samples after join: {df.shape[0]}, examples: {df.index.tolist()[:5]}")

    if df.empty:
        raise ValueError("The joined table between the expression matrix and metadata is empty. Please check Sample ID alignment.")
    logging.info(f"Loaded {data_type_key}: {df.shape[0]} samples × {expr.shape[1]} features")

    return df, expr.columns.values

def main():
    # Keep the argument parsing section unchanged
    ap = argparse.ArgumentParser(description="Run FWSE on mouse RNA-seq (ncRNA/mRNA).", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--config", default="/home/yjliu/mmProj/scripts/Machine_learning/feature_Select/human/miRNA/config_miRNA.yaml")
    ap.add_argument("--out-root", default="/home/yjliu/mmProj/data_process/Human/Feature_select/miRNA/tuning_groupkfold")
    ap.add_argument("--data-type", default="mRNA", choices=["ncRNA","mRNA"])
    ap.add_argument("--brain-region", default=None, help="Full name or abbreviation can be used, e.g., PFC/NAc/HPC/AMY")
    ap.add_argument("--top-k", type=int, default=500)
    ap.add_argument("--n-bootstraps", type=int, default=10)
    ap.add_argument("--pruning-factor", type=float, default=0.999)
    ap.add_argument("--n-jobs", type=int, default=48)
    ap.add_argument("--random-seed", type=int, default=42)
    args = ap.parse_args()

    # Environment and logging setup sections remain unchanged
    setup_hpc_env(args.n_jobs)
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    region_tag = f"_{args.brain_region}" if args.brain_region else "_all_regions"
    run_name = f"{args.data_type}{region_tag}_pf{args.pruning_factor}_{ts}"
    out_dir = os.path.join(args.out_root, run_name)
    os.makedirs(out_dir, exist_ok=True)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)-7s] %(message)s", handlers=[logging.FileHandler(os.path.join(out_dir,"run.log")), logging.StreamHandler(sys.stdout)])
    logging.info("="*80)
    logging.info(f"▶️ FWSE Run on {args.data_type}")
    logging.info(f"- Config: {args.config}")
    logging.info(f"- OutDir: {out_dir}")
    logging.info(f"- Params: {vars(args)}")
    logging.info("="*80)

    with open(args.config,"r") as f:
        config = yaml.safe_load(f)

    df, feature_names = load_data_from_config(config, args.data_type)


    # Directly use the clean feature name list returned during data loading to construct X, ensuring that X contains only numeric values
    X_train = df[feature_names].values
    y_train = df["y"].values
    # ==========================================================================

    logging.info(f"Matrix used for FWSE: {X_train.shape[0]} samples × {X_train.shape[1]} features")

    # Filter and Wrapper estimator definitions remain unchanged
    filter_estimators = [SelectKBest(f_classif, k="all"), SelectKBestSNR()]
    wrapper_estimators = [
        LinearSVC(penalty="l2", dual=True, class_weight="balanced", max_iter=4000, random_state=args.random_seed),
        LogisticRegression(penalty="l1", solver="liblinear", class_weight="balanced", max_iter=4000, random_state=args.random_seed)
    ]

    # FWSE execution and result saving sections remain unchanged
    fwse = FWSE(filter_estimators, wrapper_estimators,
                n_bootstraps=args.n_bootstraps,
                pruning_factor=args.pruning_factor,
                n_jobs=args.n_jobs, random_state=args.random_seed)

    fwse.fit(X_train, y_train)

    logging.info("Saving results...")
    final_rank = fwse.ranking_ + 1
    print(f"len(feature_names): {len(feature_names)}")
    print(f"len(final_rank): {len(final_rank)}")
    feature_rank_map = dict(zip(feature_names, final_rank))
    # Assume expr_srr is the required SRR order
    rank_df = pd.DataFrame({
    "feature_id": [f for f in feature_names if f in feature_rank_map],
    "rank": [feature_rank_map[f] for f in feature_names if f in feature_rank_map]
    }).sort_values("rank")
    path_rank = os.path.join(out_dir, "fwse_full_ranking.csv")
    rank_df.to_csv(path_rank, index=False)
    topk = rank_df.head(args.top_k)["feature_id"].tolist()
    path_topk = os.path.join(out_dir, f"top_{args.top_k}_features.txt")
    with open(path_topk,"w") as f: f.write("\n".join(topk))
    model_path = os.path.join(out_dir, "fwse_model.joblib")
    joblib.dump(fwse, model_path)
    meta = {
        "run_parameters": vars(args),
        "execution_info": {"start": ts, "host": os.uname().nodename},
        "output_files": {"ranking_csv": os.path.basename(path_rank),
                         "top_k_list": os.path.basename(path_topk),
                         "model_object": os.path.basename(model_path),
                         "log_file": "run.log"},
        "config_path": args.config
    }
    with open(os.path.join(out_dir,"run_metadata.json"),"w") as f:
        json.dump(meta, f, indent=2)

    logging.info("="*80)
    logging.info("✅ FWSE finished. Outputs are available in the output directory.")
    logging.info("="*80)

if __name__ == "__main__":
    main()

# python run_fwse.py --data-type mRNA --pruning-factor VALUE --top-k VALUE --n-bootstraps 10 --n-jobs 48
# python run_fwse.py --data-type ncRNA --pruning-factor 0.7 --top-k 100 --n-bootstraps 10 --n-jobs 48