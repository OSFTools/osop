from reformat import reformat
from c3s_get_data import get_obs, get_hcast
from c3s_compute_scores import *
from play_scores import *
from c3s_compute_products import *


C3S_DIR = os.path.join('/data/users/edyer', 'seafoam')
SCOREDIR = os.path.join(C3S_DIR, 'scores', 'c3s')

def set_verification_config():
    config = dict(
        list_vars = ['2m_temperature', ],#['Total precipitation',],#['2m_temperature', ],
        var = 't2m',#'tprate',#'t2m',##TODO make this work if above is a list
        aggr = '3m',
        hcstarty = 1993,
        hcendy = 2016,
        start_month = 5, ## Forecast initialised in month 5 (may)
        valid_month = 6, ## Forecast valid for month 6 (june)
        area = "45,-30,-2.5,60",
        scores = ['rocss', 'rps', 'bs', 'roc', 'rel'],
        leads = '2,3,4', ## TODO edit to work for custom length periods (not just 1m or 3m aggregations)
        leads_p = [2,3,4] ,## TODO edit so don't need both of these
        leads_obs = [1,2,3],
        isLagged = True
    )
    ## hindcast info
    origin = 'ukmo.s602'#'cmcc.s35'#'cmcc.s35'# ## find options here https://confluence.ecmwf.int/display/CKB/Description+of+the+C3S+seasonal+multi-system#DescriptionoftheC3Sseasonalmultisystem-system_keyword
    config['origin'], config['system'] = origin.split('.s')

    return config

def main(config):
    ## Get hindcast
    get_hcast(C3S_DIR, config)

    ## make proucts
    # For the re-shaping of time coordinates in xarray.Dataset we need to select the right one 
    #  -> burst mode ensembles (e.g. ECMWF SEAS5) use "time". This is the default option
    #  -> lagged start ensembles (e.g. MetOffice GloSea6) use "indexing_time" (see CDS documentation about nominal start date)
    st_dim_name = 'time' if not config.get('isLagged',False) else 'indexing_time'
 
    hcst_bname = '{origin}_s{system}_stmonth{start_month:02d}_hindcast{hcstarty}-{hcendy}_monthly'.format(**config)
    hcst_fname = f'{C3S_DIR}/{hcst_bname}.grib'

    ## calc anoms
    hcst, hcst_3m = calc_anoms(hcst_fname, hcst_bname, config, st_dim_name, C3S_DIR)
    ## calc terc probs
    prob_terc(hcst_bname, hcst, hcst_3m, C3S_DIR)

    ## get observations - check correct interpretation of naming conventions
    ## need to add an if file exists skip
    obs_fname = '{fpath}/era5_monthly_stmonth{start_month:02d}_{hcstarty}-{hcendy}.grib'.format(fpath=C3S_DIR,**config)
    #if not os.path.exists(obs_fname):
    get_obs(obs_fname, config)

    ## read obs
    era5_1deg, era5_1deg_3m = read_obs(obs_fname, config)

    ## calc scores
    scores_prblstc(era5_1deg, era5_1deg_3m, hcst_bname, config, C3S_DIR)

    ## plot scores
    plot_scores(hcst_bname, config)

def plot_scores(hcst_rf, config):
    titles =prep_titles(config)
    for score in config['scores']:
        print(score)
        config['score'] = score
    
        c3s_score = '{origin}_s{system}_stmonth{start_month:02d}_hindcast{hcstarty}-{hcendy}_monthly.{aggr}.{score}.nc'.format(**config)
        bs_c3s = read_score(os.path.join(C3S_DIR, 'scores', c3s_score))

        if score == 'rps':
            plot_check('c3s', bs_c3s, None, config, score, titles, C3S_DIR)
    
        elif score == 'bs':
            for cat in [0,1,2]:
                plot_check('c3s', bs_c3s, cat, config, score, titles, C3S_DIR)

        elif score == 'rel':
            plot_rel('c3s', bs_c3s, None, config, score, C3S_DIR)

        else:
            for cat in [0,1,2]:
                #config['category'] = cat
                plot_check('c3s', bs_c3s, cat, config, score, titles, C3S_DIR)

if __name__ == "__main__":
    config = set_verification_config()
    main(config)
