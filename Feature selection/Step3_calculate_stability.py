#!/usr/bin/env python3
# -*- coding: utf-8 -*-



import os
import sys
import yaml
import json
import argparse
import logging
import warnings
import time
from collections import Counter
import numpy as np
import pandas as pd
import joblib
from sklearn.model_selection import StratifiedKFold
from sklearn.feature_selection import SelectKBest


# Import other components from our main script
from run_fwse import FWSE, SelectKBestSNR, load_data_from_config, setup_hpc_env, f_classif, LinearSVC, LogisticRegression

warnings.filterwarnings("ignore")

def run_fwse_on_fold(train_indices, X_full, y_full, fwse_params, feature_names):
    """Runs FWSE on a single CV fold and returns the set of top-k feature names."""
    X_fold, y_fold = X_full[train_indices], y_full[train_indices]
    
    # Estimator definitions are now correct and consistent
    filter_estimators = [SelectKBest(f_classif, k="all"), SelectKBestSNR()]
    wrapper_estimators = [
        LinearSVC(penalty="l2", dual=True, class_weight="balanced",
                  max_iter=4000, random_state=fwse_params["random_seed"]),
        LogisticRegression(penalty="l1", solver="liblinear",
                           class_weight="balanced", max_iter=4000,
                           random_state=fwse_params["random_seed"])
    ]

    fwse = FWSE(
        filter_estimators=filter_estimators, 
        wrapper_estimators=wrapper_estimators,
        n_bootstraps=fwse_params['n_bootstraps'],
        pruning_factor=fwse_params['pruning_factor'],
        n_jobs=1,
        random_state=fwse_params['random_seed']
    )
    
    fwse.fit(X_fold, y_fold)
    
    top_k_indices = np.argsort(fwse.ranking_)[:fwse_params['top_k']]
    return set(feature_names[top_k_indices])

def main():
    parser = argparse.ArgumentParser(description="Perform stability selection using FWSE across CV-folds.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--run-dir", required=True, help="Path to a completed FWSE run directory (must contain run_metadata.json).")
    parser.add_argument("--n-folds", type=int, default=10, help="Number of folds for cross-validation.")
    parser.add_argument("--freq-threshold", type=int, default=8, help="Minimum frequency for a feature to be considered stable (e.g., 8 out of 10 folds).")
    parser.add_argument("--n-jobs", type=int, default=24, help="Number of parallel jobs for running folds.")
    args = parser.parse_args()

    setup_hpc_env(args.n_jobs)

    log_path = os.path.join(args.run_dir, "stability.log")
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)-7s] %(message)s", handlers=[logging.FileHandler(log_path, mode='w'), logging.StreamHandler(sys.stdout)])

    logging.info("--- 1. Loading metadata and data from previous run ---")
    metadata_path = os.path.join(args.run_dir, "run_metadata.json")
    if not os.path.exists(metadata_path):
        logging.error(f"Error: 'run_metadata.json' not found in {args.run_dir}")
        sys.exit(1)
    
    with open(metadata_path, 'r') as f:
        run_metadata = json.load(f)
    run_params = run_metadata['run_parameters']
    logging.info(f"Using parameters from previous run: {run_params}")
    
    config_path = run_metadata.get("config_path", run_params['config'])
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
        
    df, feature_names = load_data_from_config(config, run_params['data_type'])
    X_full = df[feature_names].values
    y_full = df["y"].values

    logging.info(f"Stability test will be based on: {X_full.shape[0]} samples × {X_full.shape[1]} features")

    fwse_params = {
        "n_bootstraps": run_params.get("n_bootstraps", 10),
        "pruning_factor": run_params.get("pruning_factor", 0.6),
        "top_k": run_params.get("top_k", 300),
        "random_seed": run_params.get("random_seed", 42)
    }

    logging.info(f"--- 2. Running FWSE in parallel across {args.n_folds} CV-folds (using 'threading' backend) ---")
    cv = StratifiedKFold(n_splits=args.n_folds, shuffle=True, random_state=fwse_params['random_seed'])
    

    fold_results = joblib.Parallel(n_jobs=args.n_jobs, backend="threading")(
        joblib.delayed(run_fwse_on_fold)(train_idx, X_full, y_full, fwse_params, feature_names)
        for train_idx, _ in cv.split(X_full, y_full)
    )
    
    logging.info("--- 3. Aggregating results and calculating frequencies ---")
    feature_counts = Counter(gene for fold_set in fold_results for gene in fold_set)
    
    freq_df = pd.DataFrame(feature_counts.items(), columns=['feature_id', 'frequency']).sort_values(by='frequency', ascending=False)
    
    freq_path = os.path.join(args.run_dir, "stability_selection_frequencies.csv")
    freq_df.to_csv(freq_path, index=False)
    logging.info(f"Full feature frequency statistics saved to: {freq_path}")
    
    stable_features = freq_df[freq_df['frequency'] >= args.freq_threshold]
    
    stable_list_path = os.path.join(args.run_dir, f"stable_features_freq_ge{args.freq_threshold}.txt")
    stable_features['feature_id'].to_csv(stable_list_path, index=False, header=False)
    
    logging.info("="*60)
    logging.info(f"Stability Selection Results:")
    logging.info(f" - Frequency Threshold: >= {args.freq_threshold} / {args.n_folds} folds")
    logging.info(f" - Number of stable features found: {len(stable_features)}")
    logging.info(f" - Stable biomarker list saved to: {stable_list_path}")
    logging.info("="*60)

if __name__ == "__main__":
    main()

# python calculate_stability.py --run-dir /home/yjliu/mmProj/Machine_learning/Feature_Select/miRNA/tuning_groupkfold/ncRNA_all_regions_pf0.7_20250913_192933/ --n-folds 10 --freq-threshold 9 --n-jobs 24