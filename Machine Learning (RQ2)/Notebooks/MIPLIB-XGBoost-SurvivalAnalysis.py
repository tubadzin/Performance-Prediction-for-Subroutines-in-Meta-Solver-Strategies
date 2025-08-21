#%%
import xgboost as xgb
import json, logging
import numpy as np
import pandas as pd
import optuna
import joblib
from pathlib import Path
from datetime import datetime
from sklearn.model_selection import train_test_split
from sklearn.model_selection import KFold
#%%
NOTEBOOK_NAME = "MIPLIB-XGBoost-SurvivalAnalysis"

timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
RUN_DIR = Path("../results") / f"{timestamp}_{NOTEBOOK_NAME}"
SUBDIRS = ["model", "optuna", "metrics", "logs"]

for sd in SUBDIRS:
    (RUN_DIR / sd).mkdir(parents=True, exist_ok=True)

print(f"All outputs will be saved under: {RUN_DIR.resolve()}")

def save_model(model, name="xgb_runtime_model"):
    model_path = RUN_DIR / "model" / f"{name}.json"
    model.save_model(model_path)
    print(f"Model saved to {model_path}")


def save_optuna(study):
    joblib.dump(study, RUN_DIR / "optuna" / "study.pkl")
    trials_df = study.trials_dataframe()
    trials_df.to_csv(RUN_DIR / "optuna" / "trials.csv", index=False)
    with open(RUN_DIR / "optuna" / "best_params.json", "w") as f:
        json.dump(study.best_params, f, indent=2)
    print("Optuna study & best params stored")


def save_metrics(name: str, value):
    with open(RUN_DIR / "metrics" / f"{name}.json", "w") as f:
        json.dump(value, f, indent=2)

logfile = RUN_DIR / "logs" / "run.log"
logging.basicConfig(filename=logfile, filemode="w",
                    level=logging.INFO,
                    format="%(asctime)s %(levelname)s | %(message)s")
logging.info("Run directory initialised")
#%%
def plot_calibration_on_solved_cases(
    true_log, pred_log,
    save_path=None, show=False,
    lower=1e-4, upper=1e4,
    figsize=(6, 6), dpi=150,
    title="True vs. Predicted runtime"
):
    true_log = np.asarray(true_log, float)
    pred_log = np.asarray(pred_log, float)
    msk = np.isfinite(true_log) & np.isfinite(pred_log)
    t_log, p_log = true_log[msk], pred_log[msk]

    t_sec = np.power(10.0, t_log)
    p_sec = np.power(10.0, p_log)

    # inside/outside plotting range
    inside = (p_sec >= lower) & (p_sec <= upper)
    p_sec_clip = np.clip(p_sec, lower, upper)

    # figure
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    ax.set_xscale("log"); ax.set_yscale("log")

    # points
    ax.scatter(t_sec[inside],  p_sec[inside],  s=18, alpha=0.6, edgecolor="none")
    ax.scatter(t_sec[~inside], p_sec_clip[~inside], s=18, alpha=0.6, marker="x")

    # red dashed diagonal
    ax.plot([lower, upper], [lower, upper], linestyle="--", linewidth=1.0, color="red")

    # axes, ticks, labels
    ax.set_xlim(lower, upper); ax.set_ylim(lower, upper)
    tick_vals = [1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4]
    tick_lbls = [fr"$10^{{{int(np.log10(t))}}}$" for t in tick_vals]
    ax.set_xticks(tick_vals, tick_lbls); ax.set_yticks(tick_vals, tick_lbls)

    ax.set_xlabel("true runtime (s)")
    ax.set_ylabel("predicted runtime (s)")
    ax.set_title(title)
    ax.grid(which="both", linestyle="--", linewidth=0.5, alpha=0.6)
    plt.tight_layout()

    # save/close
    if save_path is not None:
        save_path = Path(save_path)
        save_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, dpi=dpi)
    if show:
        plt.show(block=False)
    else:
        plt.close(fig)

#%%
def plot_aggregated_learning_curves(nll_histories, cidx_histories,
                                    save_path=None, show=False, dpi=150):
    """
    Aggregates live curves across folds.
    Plots mean (solid) and ±1 SD (shaded) over boosting iterations for:
      - aft-nloglik on the left
      - validation C-index on the right
    """
    def _pad_to_matrix(list_of_lists):
        max_len = max(len(h) for h in list_of_lists)
        M = np.full((len(list_of_lists), max_len), np.nan, dtype=float)
        for i, h in enumerate(list_of_lists):
            M[i, :len(h)] = np.asarray(h, float)
        return M

    nll_M  = _pad_to_matrix(nll_histories)
    cidx_M = _pad_to_matrix(cidx_histories)

    x      = np.arange(1, nll_M.shape[1] + 1)
    nll_mu = np.nanmean(nll_M, axis=0)
    nll_sd = np.nanstd(nll_M,  axis=0)

    x2      = np.arange(1, cidx_M.shape[1] + 1)
    cidx_mu = np.nanmean(cidx_M, axis=0)
    cidx_sd = np.nanstd(cidx_M,  axis=0)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 3), dpi=dpi)

    ax1.plot(x, nll_mu, label="aft-nloglik")
    ax1.fill_between(x, nll_mu - nll_sd, nll_mu + nll_sd, alpha=0.2)
    ax1.set_xlabel("# Boosting Iterations")
    ax1.set_ylabel("aft-nloglik")
    ax1.legend(loc="best")

    ax2.plot(x2, cidx_mu, label="C-index (valid)")
    ax2.fill_between(x2, cidx_mu - cidx_sd, cidx_mu + cidx_sd, alpha=0.2)
    ax2.set_xlabel("# Boosting Iterations")
    ax2.set_ylabel("C-index (valid)")
    ax2.legend(loc="best")

    fig.tight_layout()

    if save_path is not None:
        save_path = Path(save_path)
        save_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, dpi=dpi)

    if show:
        plt.show(block=False)
    else:
        plt.close(fig)
#%%
from lifelines.utils import concordance_index

class LiveLearningCurves(xgb.callback.TrainingCallback):
    """
    Live learning curves:
      - left: aft-nloglik (valid if available else train)
      - right: validation C-index (lifelines), updated every `plot_every` rounds
    """
    def __init__(self, dvalid: xgb.DMatrix, ylb_valid: np.ndarray, yub_valid: np.ndarray,
                 plot_every: int = 10, live: bool = True):
        super().__init__()
        self.dvalid = dvalid
        self.ylb_valid = np.asarray(ylb_valid, float)              # durations (sec)
        self.yub_valid = np.asarray(yub_valid, float)              # inf => censored
        self.event_valid = np.isfinite(self.yub_valid).astype(int) # 1=observed, 0=censored
        self.plot_every = plot_every
        self.live = live

        self.nll_hist:  list[float] = []
        self.cidx_hist: list[float] = []

        # figure
        self.fig, (self.ax1, self.ax2) = plt.subplots(1, 2, figsize=(8, 3))
        (self.nll_line,)  = self.ax1.plot([], [], "o-", label="aft-nloglik")
        (self.cidx_line,) = self.ax2.plot([], [], "o-", label="C-index (valid)")
        self.ax1.set_xlabel("# Boosting Iterations")
        self.ax2.set_xlabel("# Boosting Iterations")
        self.ax1.legend(loc="best")
        self.ax2.legend(loc="best")
        self.fig.tight_layout()
        if self.live:
            plt.ion()
            self.fig.show()

    def after_iteration(self, model: xgb.Booster, epoch: int,
                        evals_log: xgb.callback.TrainingCallback.EvalsLog) -> bool:
        if "valid" in evals_log and "aft-nloglik" in evals_log["valid"]:
            self.nll_hist.append(evals_log["valid"]["aft-nloglik"][-1])
        elif "train" in evals_log and "aft-nloglik" in evals_log["train"]:
            self.nll_hist.append(evals_log["train"]["aft-nloglik"][-1])

        y_pred = model.predict(self.dvalid)  # predicted runtimes (sec)
        try:
            cidx = concordance_index(self.ylb_valid, y_pred, event_observed=self.event_valid)
        except Exception:
            cidx = np.nan
        self.cidx_hist.append(float(cidx))

        if (epoch % self.plot_every) == 0:
            x = np.arange(1, len(self.nll_hist) + 1)
            self.nll_line.set_data(x, self.nll_hist)
            self.ax1.set_xlim(1, max(2, len(self.nll_hist)))

            x2 = np.arange(1, len(self.cidx_hist) + 1)
            self.cidx_line.set_data(x2, self.cidx_hist)
            self.ax2.set_xlim(1, max(2, len(self.cidx_hist)))

            # autoscale y
            for ax, y in [(self.ax1, self.nll_hist), (self.ax2, self.cidx_hist)]:
                if len(y) >= 2:
                    ymin, ymax = np.nanmin(y), np.nanmax(y)
                    pad = 0.05 * (ymax - ymin if ymax > ymin else (ymax or 1.0))
                    ax.set_ylim(ymin - pad, ymax + pad)

            if self.live:
                self.fig.canvas.draw()
                self.fig.canvas.flush_events()
        return False

    def save(self, path):
        self.fig.savefig(path, dpi=130)

    def close(self, path=None):
        if path:
            self.save(path)
        plt.close(self.fig)

#%%
NOTEBOOK_NAME = "MIPLIB-XGBoost-SurvivalAnalysis"

timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
RUN_DIR = Path("../results") / f"{timestamp}_{NOTEBOOK_NAME}"
SUBDIRS = ["model", "optuna", "metrics", "logs"]

for sd in SUBDIRS:
    (RUN_DIR / sd).mkdir(parents=True, exist_ok=True)

print(f"All outputs will be saved under: {RUN_DIR.resolve()}")

def save_model(model, name="xgb_runtime_model"):
    model_path = RUN_DIR / "model" / f"{name}.json"
    model.save_model(model_path)
    print(f"Model saved to {model_path}")


def save_optuna(study):
    joblib.dump(study, RUN_DIR / "optuna" / "study.pkl")
    trials_df = study.trials_dataframe()
    trials_df.to_csv(RUN_DIR / "optuna" / "trials.csv", index=False)
    with open(RUN_DIR / "optuna" / "best_params.json", "w") as f:
        json.dump(study.best_params, f, indent=2)
    print("Optuna study & best params stored")


def save_metrics(name: str, value):
    with open(RUN_DIR / "metrics" / f"{name}.json", "w") as f:
        json.dump(value, f, indent=2)

logfile = RUN_DIR / "logs" / "run.log"
logging.basicConfig(filename=logfile, filemode="w",
                    level=logging.INFO,
                    format="%(asctime)s %(levelname)s | %(message)s")
logging.info("Run directory initialised")
#%%
df_features = pd.read_csv('data/miplib_hutter_features.csv', header=0)
df_cplex_runtimes = pd.read_csv('/data/miplib_cplex_runtimes.csv')
df_cplex_runtimes["runtime_sec"] = pd.to_numeric(df_cplex_runtimes["runtime_sec"])
df_features.head()
#%%
df_combined = (
    df_features                         # left table
      .merge(                           # join
          df_cplex_runtimes[["instance", "runtime_sec"]],  # right table (only needed cols)
          on="instance",                # key column
          how="inner"                   # inner = intersection
      )
)
df_combined.head()
#%%
print(df_combined.shape)
#%%
df_features = df_combined.drop(columns=["runtime_sec"])
TIME_LIMIT = 3600.0
df_cplex_runtimes = df_combined["runtime_sec"].astype(float)

solved = df_cplex_runtimes < TIME_LIMIT
y_lb = np.clip(df_cplex_runtimes, 1e-9, None)
y_ub = np.where(solved, y_lb, np.inf)
df_features.set_index("instance", inplace=True)
df_cplex_runtimes.index = df_features.index
#%%
nonconst_mask = df_features.apply(lambda c: c.dropna().nunique() > 1)
df_features = df_features.loc[:, nonconst_mask]
#%%
df_features.replace(-512, np.nan, inplace=True)
#%%
X = df_features
runtimes = df_cplex_runtimes.clip(lower=0.005)
#%%
RANDOM_STATE    = 1234
#%%
outer = KFold(n_splits=10, shuffle=True, random_state=RANDOM_STATE)
nll_outer = []
c_index_outer = []
c_index_outer_uncensored = []
spearman_outer_uncensored = []
rmse_log10_outer = []
share_censored_predicted_lower_lb_outer = []
precision_outer = []
recall_outer = []

oof_true_log = []
oof_pred_log = []

nll_histories_all = []
cidx_histories_all = []



best_study        = None          # will hold the Optuna study with the lowest fold RMSE
best_fold_idx     = None
best_nll = np.inf
best_params = None
best_num_round = None

for fold, (train_idx, test_idx) in enumerate(outer.split(X), 1):
    print(f"Outer fold {fold}/10")
    X_tr, X_te = X.iloc[train_idx], X.iloc[test_idx]
    y_tr, y_te = runtimes.iloc[train_idx], runtimes.iloc[test_idx]

    def objective(trial):
        num_round = trial.suggest_int("num_boost_round", 750, 1500, step=50)
        params = {
            "eval_metric": "aft-nloglik",
            #"aft_loss_distribution": trial.suggest_categorical("aft_loss_distribution", ["normal", "logistic", "extreme"]),
            "aft_loss_distribution":    "extreme",
            # "aft_loss_distribution_scale": trial.suggest_float("aft_loss_distribution_scale", 0.1, 10, log=True),
            "aft_loss_distribution_scale": 1,
            "learning_rate":     trial.suggest_float("learning_rate", 0.02, 0.2, log=True),
            "max_depth":         trial.suggest_int("max_depth", 5, 10),
            "min_child_weight":  trial.suggest_float("min_child_weight", 1e-2, 3, log=True),
            "subsample":         trial.suggest_float("subsample", .2, .7),
            "colsample_bytree":  trial.suggest_float("colsample_bytree", .7, .9),
            "reg_alpha":         trial.suggest_float("reg_alpha", 1e-9, 1e-3, log=True),
            "reg_lambda":        trial.suggest_float("reg_lambda", 1e-3, 1, log=True),
            # fixed:
            "objective":         "survival:aft",
            "tree_method":       "hist",
            "seed":              RANDOM_STATE,
            "gamma":             0,
        }

        dtrain = xgb.DMatrix(X_tr)
        dtrain.set_float_info("label_lower_bound", y_lb[train_idx])
        dtrain.set_float_info("label_upper_bound", y_ub[train_idx])

        cv_results = xgb.cv(
            params,
            dtrain,
            num_boost_round       = num_round,
            nfold           = 2,
            metrics         = "aft-nloglik",
            seed            = RANDOM_STATE,
            stratified      = False,
            shuffle         = True,
            early_stopping_rounds = 50,
            verbose_eval    = False,
        )
        best = cv_results["test-aft-nloglik-mean"].iloc[-1]
        return best


    study = optuna.create_study(direction="minimize", pruner=optuna.pruners.MedianPruner(n_warmup_steps=5), sampler=optuna.samplers.TPESampler(seed=RANDOM_STATE))
    study.optimize(objective, n_trials=30, show_progress_bar=True, n_jobs=1) # TODO: here 40

    fold_params = study.best_params.copy()
    num_round   = fold_params.pop("num_boost_round")

    fold_params |= {
        "objective"   : "survival:aft",
        "eval_metric":  "aft-nloglik",
        "tree_method" : "hist",
        "seed"        : RANDOM_STATE,
    }

    # 10% validation slice from the TRAIN fold for ES + live curves
    tr_idx = np.array(train_idx)
    te_idx = np.array(test_idx)

    X_tr, X_te = X.iloc[tr_idx], X.iloc[te_idx]

    tr_fit_idx, tr_vis_idx = train_test_split(
        tr_idx,
        test_size=0.10,
        random_state=RANDOM_STATE,
        stratify=solved.iloc[tr_idx].astype(int)
    )

    X_fit = X.iloc[tr_fit_idx]
    X_vis = X.iloc[tr_vis_idx]

    dtrain = xgb.DMatrix(X_fit)
    dtrain.set_float_info("label_lower_bound", y_lb[tr_fit_idx])
    dtrain.set_float_info("label_upper_bound", y_ub[tr_fit_idx])

    dvalid = xgb.DMatrix(X_vis)
    dvalid.set_float_info("label_lower_bound", y_lb[tr_vis_idx])
    dvalid.set_float_info("label_upper_bound", y_ub[tr_vis_idx])

    dtest  = xgb.DMatrix(X_te)
    dtest.set_float_info("label_lower_bound", y_lb[te_idx])
    dtest.set_float_info("label_upper_bound", y_ub[te_idx])

    fig_dir = RUN_DIR / "figs" / f"fold_{fold}"
    fig_dir.mkdir(parents=True, exist_ok=True)

    live_cb = LiveLearningCurves(
        dvalid=dvalid,
        ylb_valid=y_lb[tr_vis_idx],
        yub_valid=y_ub[tr_vis_idx],
        plot_every=10,
        live=False
    )

    booster = xgb.train(
        params=fold_params,
        dtrain=dtrain,
        num_boost_round=5000,
        evals=[(dtrain, "train"), (dvalid, "valid")],
        early_stopping_rounds=200,
        callbacks=[live_cb],
        verbose_eval=False,
    )

    nll_histories_all.append(live_cb.nll_hist)
    cidx_histories_all.append(live_cb.cidx_hist)
    live_cb.close()

    best_iter = getattr(booster, "best_iteration", None)
    if best_iter is not None:
        pred_sec = booster.predict(dtest, iteration_range=(0, best_iter + 1))
        eval_str = booster.eval(dtest, name="test", iteration=best_iter)
    else:
        pred_sec = booster.predict(dtest)
        eval_str = booster.eval(dtest, name="test")


    fold_nll = float(eval_str.split(":")[1])
    nll_outer.append(fold_nll)


    dur_test = runtimes.iloc[te_idx].to_numpy()
    evt_test = solved.iloc[te_idx].to_numpy().astype(int)

    risk_scores = pred_sec
    c_idx = concordance_index(dur_test, risk_scores, event_observed=evt_test)
    c_index_outer.append(c_idx)


    # Uncensored-only C-index — removes censoring effects
    mask = solved.iloc[te_idx].to_numpy().astype(bool)
    c_unc = concordance_index(dur_test[mask], pred_sec[mask], event_observed=np.ones(mask.sum(), dtype=bool))
    c_index_outer_uncensored.append(c_unc)
    print("C-index (uncensored only):", c_unc)

    # Monotonicity on uncensored – should be > 0 if direction is right
    from scipy.stats import spearmanr
    rho, _ = spearmanr(dur_test[mask], pred_sec[mask])
    spearman_outer_uncensored.append(rho)
    print("Spearman rho (uncensored):", rho)

    # How many censored are predicted below their LB? (this is the killer)
    mask_cens = evt_test == 0
    share_under = float((pred_sec[mask_cens] < y_lb[te_idx][mask_cens]).mean())
    share_censored_predicted_lower_lb_outer.append(share_under)
    print("Share censored with y_pred < lower_bound:", share_under)

    pred_timeout = pred_sec >= TIME_LIMIT
    true_timeout = evt_test == 0

    tp = int((pred_timeout & true_timeout).sum())   # correct timeout call
    fp = int((pred_timeout & ~true_timeout).sum())  # predicted timeout but solved
    fn = int((~pred_timeout & true_timeout).sum())  # missed timeout
    tn = int((~pred_timeout & ~true_timeout).sum())

    precision = tp / (tp + fp + 1e-9)
    recall    = tp / (tp + fn + 1e-9)
    precision_outer.append(precision)
    recall_outer.append(recall)
    print({"timeout_precision":precision, "timeout_recall":recall})



    uncensored_mask = solved.iloc[test_idx].to_numpy()
    true_seconds = dur_test[uncensored_mask]
    pred_sec_unc = np.clip(pred_sec[uncensored_mask], 1e-9, None)   # avoid log of <=0
    log_true = np.log10(true_seconds)
    log_pred = np.log10(pred_sec_unc)
    oof_true_log.append(log_true)
    oof_pred_log.append(log_pred)


    rmse_unc_log10 = float(np.sqrt(np.mean((log_pred - log_true)**2)))
    rmse_log10_outer.append(rmse_unc_log10)

    if fold_nll < best_nll:
        best_nll = fold_nll
        best_params = fold_params.copy()
        best_num_round = num_round
        best_study = study


    print(f"  NLL = {fold_nll:.4f}")
    print(f"  C-INDEX = {c_idx:.4f}")
    print(f"  RMSE(log10 sec, uncens.) = {rmse_unc_log10:.4f}")
#%%
from optuna.visualization.matplotlib import plot_slice
import matplotlib.pyplot as plt

ax  = plot_slice(best_study)
plt.show()
#%%
print(f"\n10-fold nested-CV NLL: {np.mean(nll_outer):.4f} ± {np.std(nll_outer):.4f}")
print(f"10-fold nested-CV Full C-index: {np.nanmean(c_index_outer):.4f} ± {np.nanstd(c_index_outer):.4f}")
print(f"10-fold nested-CV RMSE: {np.nanmean(rmse_log10_outer):.4f} ± {np.nanstd(rmse_log10_outer):.4f}")

print(f"10-fold nested-CV Uncensored C-index: {np.nanmean(c_index_outer_uncensored):.4f} ± {np.nanstd(c_index_outer_uncensored):.4f}")
print(f"10-fold nested-CV Uncensored Spearman P-score: {np.nanmean(spearman_outer_uncensored):.4f} ± {np.nanstd(spearman_outer_uncensored):.4f}")
print(f"10-fold nested-CV Share Of Censored Which Predict < Cap of 3600s: {np.nanmean(share_censored_predicted_lower_lb_outer):.4f} ± {np.nanstd(share_censored_predicted_lower_lb_outer):.4f}")

print(f"10-fold nested-CV Precision: {np.nanmean(precision_outer):.4f} ± {np.nanstd(precision_outer):.4f}")
print(f"10-fold nested-CV Recall: {np.nanmean(recall_outer):.4f} ± {np.nanstd(recall_outer):.4f}")


print("Best params")
for k, v in best_study.best_params.items():
    print(f"  {k}: {v:.4}" if isinstance(v, float) else f"  {k}: {v}")
#%%
oof_true_log = np.concatenate(oof_true_log)
oof_pred_log = np.concatenate(oof_pred_log)

fig_dir = RUN_DIR / "figs"
plot_calibration_on_solved_cases(
    oof_true_log, oof_pred_log,
    save_path=fig_dir / "calibration_scatter_oof.png", # oof means out of fold
    show=False,
    lower=1e-4, upper=1e4,
    title="XGBoost ATF"
)

plot_aggregated_learning_curves(
    nll_histories_all,
    cidx_histories_all,
    save_path=fig_dir / "learning_curves_aggregated.png",
    show=False
)
#%%
dall = xgb.DMatrix(X)
dall.set_float_info("label_lower_bound", y_lb)
dall.set_float_info("label_upper_bound", y_ub)

final_params = best_params | {
    "objective": "survival:aft",
    "eval_metric": "aft-nloglik",
}

num_round = best_num_round
best_booster = xgb.train(final_params, dall, num_boost_round=num_round)
best_booster.save_model(RUN_DIR / "model" / "runtime_xgb_aft_full.json")

metrics_summary = {
    "nll_outer_mean":                      nll_outer,
    "nll_outer_std":                       nll_outer,
    "c_index_outer_mean":                  c_index_outer,
    "c_index_outer_std":                   c_index_outer,
    "c_index_outer_uncensored_mean":       c_index_outer_uncensored,
    "c_index_outer_uncensored_std":        c_index_outer_uncensored,
    "spearman_outer_uncensored_mean":      spearman_outer_uncensored,
    "spearman_outer_uncensored_std":       spearman_outer_uncensored,
    "rmse_log10_outer_mean":               rmse_log10_outer,
    "rmse_log10_outer_std":                rmse_log10_outer,
    "share_cens_pred_below_lb_mean":       share_censored_predicted_lower_lb_outer,
    "share_cens_pred_below_lb_std":        share_censored_predicted_lower_lb_outer,
    "precision_timeout_mean":              precision_outer,
    "precision_timeout_std":               precision_outer,
    "recall_timeout_mean":                 recall_outer,
    "recall_timeout_std":                  recall_outer,
    # per-fold detail
    "best_fold_idx":                       int(best_fold_idx) if best_fold_idx is not None else None,
    "best_nll":                            float(best_nll) if np.isfinite(best_nll) else None,
    "best_num_round":                      int(best_num_round) if best_num_round is not None else None,
}
save_metrics("survival_nested_cv_summary", metrics_summary)
save_optuna(best_booster)