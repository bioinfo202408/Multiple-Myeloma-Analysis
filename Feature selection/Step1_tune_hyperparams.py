#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import os, sys, yaml, argparse, logging, warnings, time
from datetime import datetime
import numpy as np, pandas as pd, joblib
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import GroupKFold
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.metrics import roc_auc_score
from lightgbm import LGBMClassifier
from sklearn.feature_selection import SelectKBest

from run_fwse import FWSE, SelectKBestSNR, load_data_from_config, setup_hpc_env, f_classif, LinearSVC, LogisticRegression

warnings.filterwarnings("ignore")

def run_groupkfold_for_params(X, y, groups, pf, k, n_bootstraps, n_jobs, random_seed):
    # When using GroupKFold, shuffling is generally not required because the split is performed by groups
    cv = GroupKFold(n_splits=10)
    fold_aucs = []
    
   
    for fold_idx, (train_indices, test_indices) in enumerate(cv.split(X, y, groups), 1):
        X_train, y_train = X[train_indices], y[train_indices]
        X_test, y_test = X[test_indices], y[test_indices]

        # If the test set contains only one class, skip this fold because AUC cannot be calculated
        if len(np.unique(y_test)) < 2:
            logging.warning(f"  > Fold {fold_idx}/5: Skipped, test set has only one class.")
            continue
            
        logging.info(f"  > Fold {fold_idx}/5: Running FWSE on {len(X_train)} training samples...")
        
        fwse = FWSE(
            filter_estimators=[SelectKBest(f_classif, k="all"), SelectKBestSNR()],
            wrapper_estimators=[
                LinearSVC(penalty="l2", dual=True, class_weight="balanced", max_iter=4000, random_state=random_seed),
                LogisticRegression(penalty="l1", solver="liblinear", class_weight="balanced", max_iter=4000, random_state=random_seed)
            ],
            n_bootstraps=n_bootstraps, pruning_factor=pf, n_jobs=n_jobs, random_state=random_seed
        )
        fwse.fit(X_train, y_train)
        
        ranking_indices = np.argsort(fwse.ranking_)
        top_k_indices = ranking_indices[:k]
        
        X_train_k = X_train[:, top_k_indices]
        X_test_k = X_test[:, top_k_indices]
        
        proxy_model = Pipeline([
            ("scale", StandardScaler(with_mean=False)),
            ("clf", LGBMClassifier(n_estimators=100, n_jobs=1, class_weight="balanced", random_state=random_seed))
        ])
        proxy_model.fit(X_train_k, y_train)
        y_pred_proba = proxy_model.predict_proba(X_test_k)[:, 1]
        
        fold_aucs.append(roc_auc_score(y_test, y_pred_proba))

    if not fold_aucs:
        logging.error(f"For (pf={pf}, k={k}), no valid folds found to calculate AUC. Check group/label distribution.")
        return {"pruning_factor": pf, "k": k, "mean_auc": 0.0}

    mean_auc = np.mean(fold_aucs)
    logging.info(f"--- Finished (pf={pf}, k={k}). Mean AUC = {mean_auc:.4f} ---")
    return {"pruning_factor": pf, "k": k, "mean_auc": mean_auc}

def main():
    ap = argparse.ArgumentParser(description="Tune FWSE with GroupKFold CV.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--config", default="/home/yjliu/mmProj/scripts/Machine_learning/feature_Select/human/miRNA/config_miRNA.yaml")
    ap.add_argument("--out-root", default="/home/yjliu/mmProj/data_process/Human/Feature_select/miRNA/tuning_groupkfold")
    ap.add_argument("--data-type", required=True, choices=["ncRNA","mRNA"])
    ap.add_argument("--brain-region", default=None)
    ap.add_argument("--pruning-factors", nargs="+", type=float, default=[0.995, 0.999, 0.9995])
    ap.add_argument("--k-values", nargs="+", type=int, default=[100, 500, 1000])
    ap.add_argument("--n-jobs", type=int, default=32) 
    ap.add_argument("--random-seed", type=int, default=42)
    ap.add_argument("--n-bootstraps", type=int, default=3, help="Bootstrap iterations for tuning (use a small number).")
    args = ap.parse_args()

    setup_hpc_env(args.n_jobs)
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = os.path.join(args.out_root, f"tuning_{args.data_type}_{ts}")
    os.makedirs(out_dir, exist_ok=True)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)-7s] %(message)s", handlers=[logging.FileHandler(os.path.join(out_dir,"tuning.log")), logging.StreamHandler(sys.stdout)])
    
    logging.info("="*80)
    logging.info(f"▶️ Starting GroupKFold Tuning for {args.data_type} (Most Rigorous)")
    logging.info("="*80)

    with open(args.config,"r") as f: config = yaml.safe_load(f)
    df, feat_names = load_data_from_config(config, args.data_type)
    X = df[feat_names].values
    y = df["y"].values
    
    # Get dataset IDs for grouping
    groups = df[config["metadata"]["dataset_column"]].values
    
    all_results = []
    for pf in args.pruning_factors:
        for k in args.k_values:
            result = run_groupkfold_for_params(X, y, groups, pf, k, args.n_bootstraps, args.n_jobs, args.random_seed)
            all_results.append(result)

    dfres = pd.DataFrame(all_results)
    dfres.to_csv(os.path.join(out_dir,"tuning_results_full.csv"), index=False)
    plt.style.use('seaborn-v0_8-whitegrid')
    pv = dfres.pivot(index="pruning_factor", columns="k", values="mean_auc")
    fig, ax = plt.subplots(figsize=(12,8))
    sns.heatmap(pv, annot=True, fmt=".4f", cmap="viridis", ax=ax, linewidths=.5)
    ax.set_title(f"GroupKFold Nested CV AUC Heatmap ({args.data_type})")
    ax.set_xlabel("Top-K Features"); ax.set_ylabel("Pruning Factor")
    plt.tight_layout()
    fig.savefig(os.path.join(out_dir,"tuning_heatmap.png"), dpi=300)
    logging.info("✅ Rigorous tuning complete.")

if __name__ == "__main__":
    main()

# Output results
# 2025-09-02 17:15:01,942 [INFO   ] --- Finished (pf=0.999, k=2000). Mean AUC = 1.0000 ---
# 2025-09-02 17:15:02,744 [INFO   ] ✅ Rigorous tuning complete.
    
# python tune_hyperparams.py --data-type mRNA --pruning-factors 0.4 0.5 0.6 0.7 0.8 0.85 0.9 --k-values 200 500 1000 1500 2000 --n-bootstraps 10 --n-jobs 20
# python tune_hyperparams_eRNA.py --data-type ncRNA --pruning-factors 0.6 0.7 0.75 0.8 0.85 0.9 0.95 --k-values 200 500 750 1000 1500 --n-bootstraps 10 --n-jobs 20