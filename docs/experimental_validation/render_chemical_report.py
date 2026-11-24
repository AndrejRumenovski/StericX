"""Create readable tables and standalone scientific figures from saved predictions."""
from pathlib import Path
import json,shutil,numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[2];DOC=Path(__file__).parent;BASE=ROOT/'.stericx/experimental_validation';SRC=BASE/'chemical_analysis_v2';OUT=DOC/'results';FIG=DOC/'figures'
OUT.mkdir(exist_ok=True);FIG.mkdir(exist_ok=True)
for name in ['experimental_observations.csv','conformer_descriptors.csv','ensemble_statistics.json','numerical_comparison.csv','descriptive_correlations.csv','primary_design.csv','statistical_results.json','post_stereo_sensitivity.json','certified_model_sensitivity.json']:shutil.copyfile(SRC/name,OUT/name)
shutil.copyfile(BASE/'conformer_sensitivity_v2/summary.json',OUT/'conformer_sensitivity.json');shutil.copyfile(BASE/'stereochemistry_v1/audit.json',OUT/'source_stereochemistry_audit.json');shutil.copyfile(BASE/'sources/correction_table_comparison.json',OUT/'correction_table_comparison.json')
r=json.loads((SRC/'statistical_results.json').read_text())['results'];primary=r['primary_as_distributed'];rows=[]
labels={'A0_train_mean':'A0: training mean','A_electronic':'A: electronic baseline','B_single':'B: first conformer','C_ensemble_min':'C: ensemble minimum','D_published':'D: published reference','E_combined':'E: reference + ensemble'}
for k,label in labels.items():
 v=primary[k];m=v['metrics'];rows.append(f"| {label} | {v['n_features']} | {m['mae']:.3f} | {m['rmse']:.3f} | {m['r2']:.3f} | {m['bias']:.3f} | {v['fold_rmse_sd']:.3f} |")
table='\n'.join(rows);u=primary['paired_uncertainty'];corr=pd.read_csv(SRC/'descriptive_correlations.csv').set_index('feature');cc=corr.loc[['A_electronic','B_single','C_ensemble_min','D_published'],['pearson_r','spearman_rho']]
correlation='\n'.join(f"| {labels[k]} | {v.pearson_r:.3f} | {v.spearman_rho:.3f} |" for k,v in cc.iterrows())
text=f'''# Experimental statistical results

The fixed comparison favors the published-library minimum in its point estimate,
but **does not establish the primary hypothesis**. There are eleven experimental
ligands in four chemical families, source stereochemistry problems, no verified
published per-state energies, and substantial uncertainty. Outcomes were already
known historically; this is an outcome-aware retrospective analysis with a protocol
locked before these new predictions, not a prospective preregistered discovery.

## Response and validation design

Response is apparent |ΔΔG|=RT ln[(100+ee)/(100−ee)], in kcal/mol, at the reported
353.15 K, R=0.0019872042586408316 kcal mol⁻¹ K⁻¹. It represents a competing-pathway
kinetic interpretation of ee; it is not a measured activation barrier, proof of
irreversibility, or a prediction of product-enantiomer sign. Integer ee and original
author ΔΔG are both retained. The source's 2064 ee/ΔΔG inconsistency is not silently
repaired by choosing the more favorable endpoint. The later corrected SI has the
same Table S3.

Primary outer splits leave out one family: cage (401), spiro (498), benzoxaphosphole
(723,724,1057,1058), or BINOL phosphoramidite (785,2062,2063,2064,2067). All models
use identical outer memberships. Every ligand is predicted exactly once; no
conformer is an experimental row. Model fitting excludes its outer test responses.
Nested ridge scaling and alpha selection use only inner training folds. All fold
predictions, coefficients, ranks, condition numbers, hyperparameters and failures
are retained. No primary or secondary fold failed.

The primary B/C descriptor was selected from the previously published scientific
model. That historical choice was informed by these same outcomes; it affects B/C
as well as the published reference D. Locking it now prevents additional adaptive
selection but cannot create an untouched external test dataset. Reported scores
are retrospective cross-validation, with this selection limitation. Current model
fitting/tuning does not leak held-out responses.

## Primary as-distributed library comparison

B is the maximum adjacent-quadrant occupied-volume difference from the lowest
numeric conformer ID; C is its minimum across the available conformers. Both are
one-slope OLS models with intercept, identical descriptor conventions and equal
capacity. A uses published P NBO charge; D uses the published Kraken steric minimum;
E uses D plus C. Conditions are constant and cannot be fitted as independent controls.

| Model | Predictors | MAE kcal/mol | RMSE kcal/mol | R² | Bias kcal/mol | SD of family-fold RMSE |
|---|---:|---:|---:|---:|---:|---:|
{table}

Pooled ΔRMSE(C−B) = **{u['delta_rmse_c_minus_b']:.3f} kcal/mol**. A 10,000-draw
paired family-cluster bootstrap gives the descriptive 95% interval
**[{u['family_bootstrap_95_interval'][0]:.3f}, {u['family_bootstrap_95_interval'][1]:.3f}]**.
The exact family-loss sign-flip one-sided p is **{u['exact_signflip_one_sided_p']:.2f}**
(only sixteen sign patterns). The interval includes no benefit; four clusters
cannot support a strong population inference. Bootstrap resampling conditions on
the saved cross-fitted predictions, does not refit models, and inherits dependence
between overlapping training folds. Ligand-bootstrap intervals are also saved but
do not resolve family dependence. Fold SD is descriptive variability, not an SE.

![Family-held-out predictions](figures/held_out_predictions.svg)

Asterisks mark the two ligands with source/export P-configuration disagreements.

The combined model does not improve on D. Near-redundancy of two steric predictors,
small training samples and extrapolation limit coefficient interpretation. Coefficients
and simple correlations describe association, not independent causal contributions.

## Correlation and uncertainty

| Predictor | Pearson r | Spearman ρ |
|---|---:|---:|
{correlation}

These use all eleven observations and are descriptive, with no multiplicity-adjusted
confirmation. The full eleven-feature correlation table, including constant features,
is retained. No feature was selected by these new correlations.

Conditional OLS t prediction intervals assume a fixed linear model with independent,
homoscedastic Gaussian errors. Held-out coverage is 8/11 for A0/A, 10/11 for B/C/D,
and 11/11 for E. Mean interval widths are respectively 3.672, 4.529, 2.838, 3.008,
2.968 and 3.723 kcal/mol. These wide, model-conditional intervals are not guaranteed
new-family prediction coverage; measurement error, model selection and chemical
domain shift are not included. Experimental ee repeatability was not available.

## Sensitivities and unsuccessful models

| Analysis | N | B RMSE | C RMSE | Interpretation |
|---|---:|---:|---:|---|
| Leave one ligand out | 11 | 0.555 | 0.363 | Easier interpolation; not a new independent dataset |
| Exclude P-discrepant 401/1058 | 9 | 0.982 | 0.484 | Three families; bootstrap ΔRMSE interval still crosses zero |
| Filter 1058, exclude 401 | 10 | 0.962 | 0.465 | P-compatible subset; published D still reflects its original construction |
| Also exclude axial-mixed 2063 | 8 | 1.054 | 0.569 | Conservative post-primary identity sensitivity |
| Filter P/relative-axis records, exclude 401 | 10 | 0.962 | 0.465 | Relative geometry control, not verified absolute atropisomer identities |
| Certified continuous regional midpoint | 11 | 0.850 | 0.395 | Same split/structures, explicit numerical method change |

Author ΔΔG and exclusion of 2064 are fully scored in the JSON. Neither changes the
qualitative conclusion. The certified all-lower/all-upper descriptor endpoint fits
are sensitivity scenarios, **not** certified extrema of a fitted statistical model.
The primary improvement is far larger than this numerical-method effect, while
statistical and structural uncertainty remain much larger still.

Nested ridge (inner leave-one-family-out, training-only scaling, alpha in
0.01/0.1/1/10/100/1000) gives RMSE **1.039** for eleven single-state features,
**1.053** for eleven ensemble means, and **1.093** for 33 minimum/maximum/mean
features. Ensemble averaging and more ensemble features do not improve this test.
Fixed random forests/gradient boosting were pre-permitted to be omitted if sample
size could not support their complexity; four clusters and training sets as small
as six do not justify that additional flexible-model comparison. No favorable
nonlinear endpoint is substituted.

Generated-geometry models are separate exploratory preparation experiments. The
corrected MMFF-weighted descriptor gives 353.15-K RMSE **1.040/1.157** across two
seeds; generated minima give **1.214/1.047**. The first, inadequately stereo-controlled
preparation and all its models are preserved too. Neither source is silently mixed
with the published conformers or presented as measured solution populations.

## Scientific conclusion

One particular ensemble minimum has a promising retrospective association and lower
cross-validated error than an arbitrary documented single conformer. The uncertainty
interval includes no improvement; multivariate ensemble means and generated energy
weights do not reproduce the apparent benefit. Native sterics add no demonstrated
predictive information beyond the published steric reference in this small sample.
The broad conformer-aware hypothesis therefore remains **insufficiently supported**.
The study establishes experimental association in a defined dataset, not mechanism,
causation, universal catalyst-ranking utility or an independent prospective validation.
'''
(DOC/'STATISTICAL_RESULTS.md').write_text(text)
plt.rcParams.update({'font.size':10,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
colors={'cage':'#ad574b','spiro':'#997335','benzoxaphosphole':'#176b86','BINOL':'#487c58'}
fig,axes=plt.subplots(1,2,figsize=(10,5.1),sharex=True,sharey=True,layout='constrained')
allpred=[v['predicted'] for k in ['B_single','C_ensemble_min'] for v in primary[k]['predictions']];low=min(-.1,min(allpred)-.18);high=max(2.3,max(allpred)+.18)
for ax,k,title in zip(axes,['B_single','C_ensemble_min'],['Single conformer · RMSE 0.850','Ensemble minimum · RMSE 0.389']):
 ax.plot([low,high],[low,high],color='#a1a7ad',ls='--',lw=1)
 for v in primary[k]['predictions']:
  ax.scatter(v['observed'],v['predicted'],color=colors[v['family']],s=42,zorder=3);ax.annotate(str(v['ligand_id'])+('*' if v['ligand_id'] in [401,1058] else ''),(v['observed'],v['predicted']),xytext=({498:(5,-13),2063:(6,9),2062:(5,10)}.get(v['ligand_id'],(4,5))),textcoords='offset points',fontsize=8)
 ax.set(title=title,xlabel='Observed apparent |ΔΔG| (kcal/mol)',xlim=(low,high),ylim=(low,high));ax.grid(alpha=.16)
axes[0].set_ylabel('Family-held-out prediction (kcal/mol)');fig.suptitle('Ni-hDA · 11 ligands · same four family splits',fontsize=13)
from matplotlib.lines import Line2D
fig.legend([Line2D([0],[0],marker='o',color='none',markerfacecolor=c,markeredgecolor=c,markersize=6) for c in colors.values()],list(colors),loc='outside lower center',ncol=4,frameon=False)
fig.savefig(FIG/'held_out_predictions.svg');fig.savefig(FIG/'held_out_predictions.png',dpi=160);plt.close(fig)
ensembles=json.loads((SRC/'ensemble_statistics.json').read_text());pub=[r for r in ensembles if r['lane']=='published_dft_library'];precision=pd.read_csv(SRC/'numerical_comparison.csv');fig,ax=plt.subplots(figsize=(8,5.7),layout='constrained');xx=np.arange(len(pub));ranges=[r['features']['max_delta_qvbur']['range'] for r in pub];errors=[precision[precision.ligand_id==r['ligand_id']].delta_grid_minus_continuous.abs().max() for r in pub];ax.scatter(ranges,xx,label='Across-library range',s=42,color='#176b86');ax.scatter(errors,xx,label='Max grid–continuous discrepancy',s=38,marker='x',color='#b45f39');ax.set_yticks(xx,[str(r['ligand_id'])+('*' if r['ligand_id'] in [401,1058] else '') for r in pub]);ax.invert_yaxis();ax.set(xscale='log',xlabel='Adjacent-quadrant volume difference (Å³, logarithmic scale)',ylabel='Kraken ligand ID',title='Conformer variation versus observed numerical discrepancy');ax.grid(axis='x',alpha=.2);ax.legend(fontsize=9);fig.savefig(FIG/'conformer_vs_numerics.svg');fig.savefig(FIG/'conformer_vs_numerics.png',dpi=160);plt.close(fig)
print('Created statistical report, tables and scientific figures.')
