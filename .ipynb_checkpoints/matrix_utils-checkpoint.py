import glob
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt

def targets_ykr_ids(grid, targets):    
    # print column names
    # print('grid cols: ', grid.columns)
    # print('targets cols: ', targets.columns)
    
    # check if CRS of layers match (-> False)
    # print('CRS match:', targets.crs == grid.crs)
    
    # reproject targets to grid CRS
    targets = targets.to_crs(grid.crs)
    # CRS should now match
    print('CRS match:', targets.crs == grid.crs)
    
    # join ykr grid info to targets
    targets_ykr = gpd.sjoin(targets, grid, how="inner", op="within")
    # get the ids as list
    ykr_ids = targets_ykr.YKR_ID.tolist()
    print('YKR ids:', ykr_ids)
    
    timecols = []
    names = []
    for idx, target in targets.iterrows():
        columnn = 'pt_r_t_'+str(idx)
        # print('target:', target)
        names.append(target['name'])
        timecols.append(columnn)
        
    return {'ykr_ids': ykr_ids, 'timecols': timecols, 'names': names}

def get_travel_times_to_targets(grid, ykr_ids, folder):
    filepaths = []
    for id in ykr_ids:
        subfolder = str(id)[:4]+'xxx/'
        filename = 'travel_times_to_ '+str(id)+'.txt'
        filepaths.append(folder + subfolder + filename)
        
    # read files from filepaths
    tt_dfs = []
    for idx, fp in enumerate(filepaths):
        data = pd.read_csv(fp, sep=';')
        data = data[['from_id', 'pt_r_t']]
        tt_dfs.append(data)

    # rename travel time columns in dataframes
    # pt_r_t: public transport rush hour traffic
    t_columns = []
    for idx, tt_df in enumerate(tt_dfs):
        columnn = 'pt_r_t_'+str(idx)
        tt_df.columns = ['from_id', columnn]
        t_columns.append(columnn)

    # combine travel times
    ttimes = grid    
    for tt_df in tt_dfs:
        ttimes = pd.merge(ttimes, tt_df, how='inner', left_on='YKR_ID', right_on='from_id')
   
    # select only relevant columns from the data
    columns = ['YKR_ID', 'geometry']
    ttimes = ttimes[columns + t_columns]
    
    # replace -1 with No Data
    ttimes = ttimes.replace(-1, np.nan)
    
    return ttimes
    
def add_min_travel_times_to_df(ttimes):
    # get time column names as list
    timecols = []
    for col in list(ttimes):
        if 'pt_r_t_' in col:
            timecols.append(col)
            
    # calculate minimum travel time to column 'min_t'
    ttimes['min_t'] = ttimes.apply(lambda row: row[timecols].min(), axis=1)
    
    # add index of the smallest travel time to the df
    ttimes['min_idx'] = ttimes.apply(lambda row: row[timecols].astype(float).idxmin(), axis=1)
    
    # drop rows having nan values in 'min_t' or 'min_idx' columns
    count_withnan = len(ttimes)
    ttimes = ttimes.dropna(subset=['min_t', 'min_idx'])
    count_withoutnan = len(ttimes)

    print('dropped', count_withnan - count_withoutnan, 'rows with na values')
    
    # float to int
    ttimes['min_t'] = [int(value) for value in ttimes['min_t']]
        
    return ttimes

def add_population_to_travel_times(ttimes, ykr_pop):
    # join population info (from points) to grid cells
    ttimes_pop = pd.merge(ttimes, ykr_pop, how='inner', on='YKR_ID')
    return ttimes_pop


def calculate_cumulative_pop(ttimes_pop, timecol):
    grouped = ttimes_pop.groupby([timecol])
    times = []
    pops = []
    cumpops = []
    totpop = 0
    for idx, values in grouped:
        times.append(idx)
        grouppop = values['ASUKKAITA'].sum()
        cumpop = sum(pops) + grouppop
        cumpops.append(cumpop)
        pops.append(grouppop)

    cum_pops = gpd.GeoDataFrame(data={'time': times, 'population': pops, 'cumpopulation': cumpops})  
    return cum_pops

def plot_cum_pops(cum_pops_15, cum_pops_18, idx, target_info):
    # prepare fig & ax for plotting
    fig, ax = plt.subplots(figsize=(12,7))

    # plot data
    ax.plot(cum_pops_18['time'], cum_pops_18['cumpopulation'], c='red', label='2018')
    ax.plot(cum_pops_15['time'], cum_pops_15['cumpopulation'], c='black', label='2015')

    # set labels
    targetname = target_info['names'][idx]
    ax.set(xlabel='Travel time (PT)', ylabel='Population reached')
    ax.set_title(target_info['names'][idx], fontsize=16)

    # set axis & ticks
    ax.set_xlim([0,71])
    ax.grid()
    # ticks every 10
    x_ticks = np.arange(0, 71, 10)
    #y_ticks = np.arange(0, 1110000, 100000)
    ax.set_xticks(x_ticks)
    #ax.set_yticks(y_ticks)

    # set font size
    plt.rcParams.update({'font.size': 12})
    ax.xaxis.label.set_size(14)
    ax.yaxis.label.set_size(14)

    # add legend
    ax.legend()
    
    filename = 'cum_pop_'+ str(idx) +'_'+ targetname

    # save plot
    fig.savefig('plots/'+filename, dpi=130)
    # fig.savefig('plots/pop_curve_test.eps', format='eps', dpi=1000)

    # show plot
    plt.show()