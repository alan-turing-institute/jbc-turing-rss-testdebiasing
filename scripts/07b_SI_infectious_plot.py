import math
import matplotlib

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from matplotlib import gridspec

sns.set_style("whitegrid")

plt.rc("text", usetex=True)
plt.rc("font", family="serif")
matplotlib.rcParams["text.latex.preamble"] = [r"\usepackage{amsmath}", r"\usepackage{amsfonts}"]

#=================================================
# DATA FILES
#=================================================
# 1. The weekly positive counts by PHE region
weekly_counts_path = "data/region.csv"
# 2. The infectious estimates - produced by the infectiousness_estimates.R script
infectious_estimates_path = "data/moment_match_infectious.csv"
# 3. samples for P(PCR+ | time since infection)
# from https://github.com/cmmid/pcr-profile
samples_path = "data/samples_pcr.csv"

#=================================================
# LOAD DATA AND PLOT
#=================================================
df_data = pd.read_csv(weekly_counts_path)
df_data = df_data[['phe_region', 'nt', 'mid_week']]
df_data['mid_week'] = pd.to_datetime(df_data['mid_week'])
dates = df_data['mid_week'].unique()

inf_df = pd.read_csv(infectious_estimates_path)
inf_df = inf_df.loc[inf_df['mid_week'] > "2020-06-14"]

# Get mean and lower/upper CIs for uniform incidence
def get_pcr_plus(path=samples_path):

    df = pd.read_csv(path)
    df = df.drop(['Unnamed: 0'], axis=1)
    arr = df.to_numpy()

    res = []

    for k in range(4):
        start_week = k * 70
        end_week = 70 + k * 70
        mean_samples = arr[:, start_week:end_week].mean(axis=1)
        res.append(mean_samples)

    return np.array(res)

pcr_plus = get_pcr_plus()
multipliers = np.array([6/7, 5/7, 0, 0])

uniform = (pcr_plus * multipliers[:, None]).sum(axis=0)/pcr_plus.sum(axis=0)
lower_uniform, upper_uniform = np.quantile(uniform, 0.025), np.quantile(uniform, 0.975)
BASE_PI = uniform.mean()

def plot_posterior(
    ax,
    posterior_means,
    posterior_confidence_intervals=None,
    title=None,
    legend=False,
    posterior_color="xkcd:orange",
    posterior_label="Posterior mean",
    confidence_interval_color="xkcd:sky blue",
    confidence_interval_label=None,
    obs_label=None,
    tex_plot=False,
    start_idx=0,
    alpha=0.2
):
    if tex_plot:
        plt.rc("text", usetex=True)
        plt.rc("font", family="serif")
        matplotlib.rcParams["text.latex.preamble"] = [r"\usepackage{amsmath}"]
    support = np.arange(len(posterior_means)) + start_idx
    ax.plot(support, posterior_means, "r", label=posterior_label, color=posterior_color)

    if posterior_confidence_intervals is not None:
        upper, lower = (
            posterior_confidence_intervals[:, 0],
            posterior_confidence_intervals[:, 1],
        )
        ax.fill_between(
            support,
            posterior_means,
            upper,
            alpha=alpha,
            color=confidence_interval_color,
            label=confidence_interval_label,
        )
        ax.fill_between(
            support, posterior_means, lower, alpha=alpha, color=confidence_interval_color
        )
    if title is not None:
        ax.set_title(title, size=16)
    if legend:
        ax.legend(loc=0, fontsize=16)

    plt.setp(ax.xaxis.get_majorticklabels(), rotation=70)

fig = plt.figure(figsize=(12, 25))
gs = gridspec.GridSpec(10, 1,height_ratios=[3]+[2]*9)
ax = [plt.subplot(gs[i]) for i in range(10)]

for name, area_df in df_data[df_data['mid_week']>='2020-06-21'].groupby('phe_region'):
    ax[0].plot(np.arange(53), area_df['nt'], label=name)
    ax[0].set_yscale('log')
ax[0].legend(loc=2, fontsize=10)
ax[0].set_xticklabels([])
ax[0].set_xticks(np.arange(53))

for i, name in enumerate(inf_df['phe_region'].unique()):
    area_df = inf_df[inf_df['phe_region']==name]
    mus = area_df['mean'].values
    up = area_df['upper'].values
    low = area_df['lower'].values
    cis = np.array([(u,l) for u,l in zip(up,low)])

    plot_posterior(ax[i+1], mus, cis, posterior_label='Varying Incidence' if i==0 else '')
    ax[i+1].plot([0,52], [BASE_PI, BASE_PI], 'r--', alpha=0.5, label='Uniform Incidence' if i==0 else '')
    ax[i+1].set_xticklabels([])
    ax[i+1].set_xticks(np.arange(53))
    ax[i+1].set_ylim(0.4, 0.8)
    eb = ax[i+1].errorbar([0], [BASE_PI], np.array([upper_uniform - BASE_PI, BASE_PI - lower_uniform])[:, None], alpha=0.5,
                   fmt='r--', capsize=8)
    eb[-1][0].set_linestyle('--')

    txt_anno = ax[i+1].text(1,0.75, name, fontsize=12)
    txt_anno.set_bbox(dict(alpha=0.6, facecolor='white'))

    plt.xticks(np.arange(53), pd.DatetimeIndex(dates[3:]).strftime("%m-%d"))

    plt.setp(ax[i+1].get_xticklabels(), fontsize=12)
    plt.setp(ax[i+1].get_yticklabels(), fontsize=12)

plt.setp(ax[0].get_xticklabels(), fontsize=12)
plt.setp(ax[0].get_yticklabels(), fontsize=12)
ax[0].set_title('Pillar 1+2 raw weekly incidence', fontsize=16)
ax[0].set_ylabel('Positive test counts (log scale)', fontsize=14)
ax[1].set_title(r'$\mathbb{P}(\text{Infectious} \mid \text{PCR+})$', fontsize=16)
ax[1].legend(loc=1, fontsize=12, title='Sampling Model', title_fontsize=12)
ax[-1].set_xlabel('Week', fontsize=16)
ax[5].set_ylabel(r'$\mathbb{P}(\text{Infectious} \mid \text{PCR+})$', fontsize=16)
plt.tight_layout()

plt.savefig('plots/SI_pcrpos_to_infectious', dpi=150)
