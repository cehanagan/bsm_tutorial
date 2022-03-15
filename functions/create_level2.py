### This file contains the imports and function to create a level 2 bsm tsv
### and regionally transformed time series.
## Level 2 data include:
# Metadata, linearized gauge strain with modelled tides,
# barometric pressure correction, linear trend
# Areal and shear strains with the same corrections

## Output files:
# SeedCode.DataRate.StartDateTime.EndDateTime.Level2.txt
# SeedCode.DataRate.StartDateTime.EndDateTime.RegionalStrains.txt
# Metadata is above the headers

import bsm_functions


def input(sta_code, rate, network, start, end):
    odf = bsm_functions.get_data(sta_code, rate, network, start, end)
    df, comment = bsm_functions.linearize(sta_code, odf)
    sta_meta = bsm_functions.meta_get(sta_code)
    df, baro_channel, sf, baro_resp, baro_mean = bsm_functions.baro_corr(network, sta_code, start, end, sta_meta, df)
    df = bsm_functions.lin_corr(df)
    df, tide_in = bsm_functions.spotl_tides_corr(df, sta_code, rate, start, end, sta_meta)
    # Write to file
    fname = sta_code + '.' + rate + '.' + start + '.' + end + '.' + 'Level2.txt'
    with open(fname, 'w') as f:
        f.write(comment)
        f.write('Barometric Pressure Correction: Channel and scale factor (counts to millibars):' + \
                baro_channel + ', ' + str(sf) + '\n gauge 0-3 response: ' + str(baro_resp[0]) + ', ' + str(
            baro_resp[1]) + ', ' + str(baro_resp[2]) + ', ' + str(baro_resp[3]) \
                + 'Mean removed: ' + str(baro_mean) + ', correction time series linearly interpolated.' + '\n')
        f.write('Linear trend is from scipy detrend function \n')
        f.write('SPOTL tides computed with standard input for the hartid function: \n' + \
                str(tide_in) + '\n')
    f.close()
    df.to_csv(fname, sep='\t', index=False, mode='a')
    # Regionally transformed df
    a_inv = bsm_functions.orient_mat(sta_meta)
    tdf = bsm_functions.transform(df, a_inv)
    tfname = sta_code + '.' + rate + '.' + start + '.' + end + '.' + 'RegionalStrains.txt'
    with open(tfname, 'w') as f:
        f.write('Transformed gauge strains and corrections. \n')
        f.write('Transformation matrix from metadata: \n' + \
                str(a_inv) + '\n')
    f.close()
    tdf.to_csv(tfname, sep='\t', index=False, mode='a')
    return df, tdf
