import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
    '''
    After running steps 01 and 02, we should have output light curve text files to be plotted.
    This code plots both together on a single plot in different colors.
    '''

    df_clean_b = pd.read_csv("output_lc_b.csv")
    df_clean_r = pd.read_csv("output_lc_r.csv")

    fig = plt.figure(figsize=(10,5))
    ax = fig.add_subplot(111)

    ax2 = ax.secondary_xaxis('top',functions=(lambda x: 24*(x-min(df_clean_b['bjd_tdb'].values)), lambda x: x/24+min(df_clean_b['bjd_tdb'].values)))
    ax2.set_xlabel('Time (hours)', fontsize='large')

    ax.errorbar(df_clean_b["bjd_tdb"], df_clean_b["rel_flux"], yerr=df_clean_b["rel_flux_err"], color='green', marker=".", markersize=3, elinewidth=0.5, capsize=1, linestyle='None', alpha=0.8, label="g-band")
    ax.set_xlabel("JD (days)", fontsize="large")
    ax.set_ylabel("Relative Flux", fontsize="large")
    plt.grid(visible=True,which='both',linestyle="--", alpha=0.5, color='grey', linewidth=0.5)

    ax.errorbar(df_clean_r["bjd_tdb"], df_clean_r["rel_flux"], yerr=df_clean_r["rel_flux_err"], color='red', marker=".", markersize=3, elinewidth=0.5, capsize=1, linestyle='None', alpha=0.8, label="i-band")
    ax.set_ylabel("Relative Flux", fontsize="large")
    plt.legend()
    plt.show()
