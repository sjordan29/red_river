'''
Script to plot historical annual cycles of tair and prcp for 
APHRODITE and CMIP5 realizations.
'''

import sys
sys.path.append('/scratch/smj5vup/CMIP6/omoCMIP/')
sys.path.append('/scratch/smj5vup/CMIP6/omoCMIP/esd/')
sys.path.append('/scratch/smj5vup/CMIP6/omoCMIP/esd/util/')
sys.path.append('/home/smj5vup/.local/bin')

import esd
import matplotlib.pyplot as plt
import os
import pandas as pd
import seaborn as sns
sns.set_context('poster',font_scale=1.3, rc={"lines.linewidth": 3})
sns.set_style('whitegrid')

if __name__ == '__main__':
    
    # Historical annual cycles for aphrodite 1deg and CMIP5 realizations that
    # have NOT been biased corrected and downscaled
    fpath_prcp_csv = os.path.join(esd.cfg.path_cmip6_trends, 'hist_ann_cycle_prcp.csv')
    fpath_tmin_csv = os.path.join(esd.cfg.path_cmip6_trends, 'hist_ann_cycle_tmin.csv')
    fpath_tmax_csv = os.path.join(esd.cfg.path_cmip6_trends, 'hist_ann_cycle_tmax.csv')

    fpath_out_fig_prcp = os.path.join(esd.cfg.path_cmip6_trends,'figures','annual_cycle_prcp.png')
    fpath_out_fig_tmin = os.path.join(esd.cfg.path_cmip6_trends,'figures','annual_cycle_tmin.png')
    fpath_out_fig_tmax = os.path.join(esd.cfg.path_cmip6_trends,'figures','annual_cycle_tmax.png')

    # # Historical annual cycles for original aphrodite .25deg and CMIP5 realizations that
    # # have been biased corrected and downscaled
    # fpath_prcp_csv = os.path.join(esd.cfg.path_cmip5_trends, 'hist_ann_cycle_prcp_p25deg.csv')
    # fpath_tair_csv = os.path.join(esd.cfg.path_cmip5_trends, 'hist_ann_cycle_tair_p25deg.csv')
    # fpath_out_fig_prcp = os.path.join(esd.cfg.path_cmip5_trends,'figures','annual_cycle_prcp_p25deg.png')
    # fpath_out_fig_tair = os.path.join(esd.cfg.path_cmip5_trends,'figures','annual_cycle_tair_p25deg.png')
    
    ############################################################################
    
    # Precipitation
    
    df_prcp = pd.read_csv(fpath_prcp_csv,index_col=0)
    df_prcp.iloc[:,0:-1].plot(legend=False, color='#bdbdbd')
    plt.plot(df_prcp.chirps05, color='#d95f02', label='CHIRPS05',lw=6)
    plt.plot(df_prcp.iloc[:,0:-1].mean(axis=1), color='#7570b3', label='ensemble mean',lw=6)
    
    ax = plt.gca()
    handles, labels = ax.get_legend_handles_labels()
    
    plt.legend(handles[-3:],['ensemble member']+labels[-2:],bbox_to_anchor=(0, 1), loc='upper left', fontsize =10)
    plt.ylabel('Precipitation (mm)', fontsize=15)

    plt.xlabel("month", fontsize=15)

    plt.tight_layout()
    plt.savefig(fpath_out_fig_prcp, dpi=300, bbox_inches='tight')
    
    # Temperature 

    # min
    for csv_name,y_name, output_path in zip([fpath_tmin_csv, fpath_tmax_csv], ["tmin", "tmax"], [fpath_out_fig_tmin, fpath_out_fig_tmax]):
    
        df_tair = pd.read_csv(csv_name,index_col=0)
        df_tair.iloc[:,0:-1].plot(legend=False, color='#bdbdbd')
        plt.plot(df_tair.era5, color='#d95f02', label='ERA5',lw=10)
        plt.plot(df_tair.iloc[:,0:-1].mean(axis=1), color='#7570b3', label='ensemble mean',lw=6)
        
        ax = plt.gca()
        handles, labels = ax.get_legend_handles_labels()
        
        plt.legend(handles[-3:],['ensemble member']+labels[-2:],bbox_to_anchor=(0, 1), loc='upper left', fontsize =10)
        if y_name =="tmin":
            plt.ylabel(u'Minimum Temperature (\u00b0C)', fontsize=15)
        else:
            plt.ylabel(u'Maximum Temperature (\u00b0C)', fontsize=15)

        plt.xlabel("month", fontsize=15)

        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        