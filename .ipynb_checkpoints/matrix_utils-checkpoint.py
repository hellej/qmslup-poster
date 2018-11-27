import glob
import numpy as np
import pandas as pd
import geopandas as gpd

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

    return ykr_ids

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
    
    ykrcolumns = ['YKR_ID', 'geometry']
    ttcolumns = ['min_t', 'min_idx']
    
    return ttimes[ykrcolumns + ttcolumns]


def add_population_to_travel_times(ttimes, ykr_pop):
    # join population info (from points) to grid cells
    ttimes_pop = pd.merge(ttimes, ykr_pop, how='inner', on='YKR_ID')
    return ttimes_pop


def calculate_cumulative_pop(ttimes_pop):
    grouped = ttimes_pop.groupby(['min_t'])
    times = []
    pops = []
    cumpops = []
    totpop = 0
    for idx, values in grouped:
        times.append(idx)
        pops.append(values['ASUKKAITA'].sum())
        totpop = sum(pops) + values['ASUKKAITA'].sum()
        cumpops.append(totpop)

    cum_pops = gpd.GeoDataFrame(data={'time': times, 'population': pops, 'cumpopulation': cumpops})  
    return cum_pops
