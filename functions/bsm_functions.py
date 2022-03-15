### This file contains functions for downloading and correcting
### Strain data and metadata from IRIS and UNAVCO

import numpy as np
import pandas as pd
from math import ceil
import os
from obspy import read, UTCDateTime, read_inventory
from obspy.clients.iris import Client
client = Client()
from obspy.clients.fdsn import Client
inv_client = Client('IRIS')
import shutil
import urllib.request as request
from contextlib import closing
import xmltodict
from scipy import signal, linalg

# Function to get data, with obspy
def get_data(sta_code, rate, network, start, end):
    '''
    Input for obspy client to get data from IRIS.
        4 character site code (e.g. B921)
        Network code (e.g. PB)
        2 charater code channel for data rate
        (1 sps: LS, 20 sps: BS, 1/600 sps: RS)
        datetime start and end (format YYYY-MM-DDThh:mm:ss.ssssssZ
    '''
    # create "dictionaries" for the gauge channels
    # ch will store the data stream for each gauge
    # (IRIS data is in multiple "stream" objects that need to be combined)
    # ts will store the time and data for each gauge
    ch = {'0': [], '1': [], '2': [], '3': []}

    # This first for loop brings in the data from iris using
    # the ObsPy client
    for channel in ch:
        # assign the gauge channel number
        n = str(int(channel) + 1)
        # Use obspy client to get the timeseries from iris
        ch[channel] = client.timeseries(network, sta_code, 'T0', rate + n, start, end)

    # Establish a dataframe with the raw data to be linearized
    # Replace bad data values flagged with 999999 to None
    tmpd = {'0': None, '1': None, '2': None, '3': None}

    for channel in ch:
        empty = np.array([])
        for i in range(0, len(ch[channel])):
            b = np.append(empty, ch[channel][i].data)
        tmpd[channel] = b

    for i in range(0, len(ch['0'])):
        empty = np.array([])
        tmpt = np.append(empty, ch['0'][i].times('timestamp'))
    odf = pd.DataFrame({'ch0': tmpd['0'], 'ch1': tmpd['1'], 'ch2': tmpd['2'], 'ch3': tmpd['3']}, tmpt).replace(999999,
                                                                                                               None)

    return odf


# Function to grab metadata and linearize gauge strain
def linearize(sta_code, odf):
    '''
    Linearizes the raw gauge strain with info from the
    UNAVCO metadata.
    '''
    metalink = 'https://www.unavco.org/data/strain-seismic/bsm-data/lib/docs/bsm_metadata.txt'
    metadata = pd.read_csv(metalink, sep='\s+', index_col=False, header=0,
                           names=['BNUM', 'NAME', 'LAT', 'LONG', 'ELEV(m)', 'INSTALL_DATE', 'CH0(EofN)', 'BSM_Depth(m)',
                                  'SEISMOMETER_Depth(m)', 'PORE_DEPTH(m)', 'DATA_START', ' ', 'DATA_END', 'GAP(m)',
                                  'L_DATE', 'L0(cnts)', 'L1(cnts)', 'L2(cnts)', 'L3(cnts)', 'REGION'])

    metadata.index = metadata['BNUM']
    metadata = metadata.loc[sta_code, :]
    ut = {'0': None, '1': None, '2': None, '3': None}
    R = metadata['GAP(m)']
    M = 0.087  # in meters
    comment = f'# Linearization info: \n \
               # Reference Gap [m] {R}\n \
               # Gauge diameter [m] {M}\n \
               # The subtracted terms for CH 0-3 [ms]'
    time = []
    for t in range(0, len(odf)):
        time.append(UTCDateTime(odf.index[t]))
    print(comment)
    for i in ut:
        empty = []
        d0 = odf.iat[0, int(i)] / 1e8
        comment = comment + ', ' + str(d0 * R / M * 1e6)
        print(str(d0 * R / M * 1e6))
        for j in range(0, len(odf)):
            dt = odf.iat[j, int(i)] / 1e8
            empty.append(((dt / (1 - dt)) - (d0 / (1 - d0))) * R / M * 1e6)
            ut[i] = empty
    comment = comment + '\n'
    df = pd.DataFrame(
        {'datetime': time, 'ch0 [ms]': ut['0'], 'ch1 [ms]': ut['1'], 'ch2 [ms]': ut['2'], 'ch3 [ms]': ut['3']})
    return df, comment


# Function to call relevant values from the metadata
# Downloads the metadata file. Which sometimes has formatting issues fyi
def meta_get(sta_code):
    xml_file = sta_code + '.xml'
    xml_path = 'ftp://bsm.unavco.org/pub/bsm/level2/' + sta_code + '/'
    with closing(request.urlopen(xml_path + xml_file)) as r:
        with open(xml_file, 'wb') as f:
            shutil.copyfileobj(r, f)
    xmldict = xmltodict.parse(open(xml_file).read(), process_namespaces=True)
    return xmldict

# Function to create pressure correction ts
def baro_corr(network, sta_code, start, end, xmldict, df):
    # We need to convert the raw pressure data to millibars
    # LDO is 1sps, RDO is 30 min. Automatically grab LDO if available.
    # Let's get the pressure data, take an initial look, and get the conversion to geophysical units
    inv = inv_client.get_stations(network=network, station=sta_code, channel='*DO', level='response',
                                  matchtimeseries=True, starttime=start, endtime=end)
    if len(inv.get_contents()['channels'][0]) == 2:
        baro_channel = 'LDO'
    else:
        baro_channel = inv.get_contents()['channels'][0][-3:]
    baro_meta = inv_client.get_stations(network=network, station=sta_code, channel=baro_channel, level='response')
    if baro_channel == 'LDO':
        baro_loc = ''
        # Technically we need to add 1000 to get real units of pressure for LDO
        # but we demean the data so just applying the factor is sufficient
        sf = baro_meta[0][0][0].response.instrument_sensitivity.value
    else:
        # The scale factor is a conversion to KPa for RDO, we need millibars
        baro_loc = 'TS'
        sf = baro_meta[0][0][0].response.instrument_sensitivity.value * 10
    baro = client.timeseries(network, sta_code, baro_loc, baro_channel, start, end)
    # Start on the linearly interpolated correction
    fillme = np.array([])
    for i in range(0,len(baro)):
        p_array = np.append(fillme,baro[i].data)
    # first use sf in instrument response info to get milibar
    # (check that the above output is in hectopascals, hPa=milibar)
    # Note that the setra data requires 1000 to be added after unit conversion for real units,
    # but we remove the mean anyway so it doesn't matter
    millibar = p_array * sf
    # Get the UTC datetimes
    empty = np.array([])
    for i in range(0,len(baro)):
        t = baro[i].times('utcdatetime')
        tmpt = np.append(empty, t)
    # Demean the data
    baro_df = pd.DataFrame(millibar-millibar.mean(),tmpt,columns=['millibar'])
    print('Linearly interpolating '+baro_channel+' data to gauge data.')
    baro_data = interp(baro_df,df.datetime).values

    # compute the barometric pressure correction
    baro_ch = {'0':np.array([]),'1':np.array([]),'2':np.array([]),'3':np.array([])}
    # get the channel response from the xml data
    baro_resp = []
    for channel in baro_ch:
        ch_resp = xmldict['strain_xml']['inst_info']['processing']['bsm_processing_history'][-1]['bsm_processing']['atm_pressure']['apc_g'+channel]
        baro_resp.append(ch_resp)
        baro_ch[channel] = np.array(baro_data*float(ch_resp))
        df['baro_ch'+channel] = baro_ch[channel]
    return df, baro_channel, sf, baro_resp, millibar.mean()

def lin_corr(df):
    print('Computing linear trend...')

    trend_ch = {'0':np.array([]),'1':np.array([]),'2':np.array([]),'3':np.array([])}

    # get the least squares linear trend of the data using the scipy 'detrend' function
    for channel in trend_ch:
        trend = np.array(df['ch'+channel+' [ms]']-df['baro_ch'+channel]) - signal.detrend(np.array(df['ch'+channel+' [ms]']-df['baro_ch'+channel]),type='linear')
        trend_ch[channel] = trend
        df['trend_ch'+channel] = trend_ch[channel]
    return df

# Function to run spotl and create tidal ts
# Error message if spotl doesn't run
# Also loads spotl tides to a dataframe
def spotl_tides_corr(df, sta_code, rate, start, end, xmldict):
    # Directory of SPOTL and directory to store the tidal time series
    print('Accessing SPOTL to compute tides')
    spotl = 'hartid'
    dname = os.getcwd()+'/SpotlTides/'
    os.makedirs(dname, exist_ok=True)
    # xml metadata starting point
    h_stit = xmldict['strain_xml']['inst_info']['processing']['bsm_processing_history'][-1]['bsm_processing']['tidal_parameters']['tide']
    start_dt = str(start)[0:4] + ' ' +str(start)[5:7] + ' ' +str(start)[8:10] + ' ' +str(start)[11:13] + ' ' +str(start)[14:16] + ' ' +str(start)[17:19]
    for ch in ['gauge0','gauge1','gauge2','gauge3']:
        # Start a formatted print statement to run from the command line
        lon = xmldict['strain_xml']['inst_info']['station_information']['coordinate']['long']
        print_str = f'printf \'l\\n{str(lon)}'
        # loop through each constituent present in the metadata
        for i in range(0,len(h_stit)):
            gauge = h_stit[i]['phz']['@kind']
            tide = h_stit[i]['@name']
            if ch == gauge and (tide == 'O1' or tide == 'M2' or tide == 'K1' or tide == 'S2' or tide == 'N2' or tide == 'P1'):
                # Separated Cartwright codes for each harmonic constituent in the metadata
                dood = str.split(h_stit[i]['@doodson'])
                one = dood[0]
                two = dood[1]
                three = dood[2]
                four = dood[3]
                five = dood[4]
                six = dood[5]
                # Amplitude and phase
                amp = h_stit[i]['amp']['#text']
                phz = h_stit[i]['phz']['#text']
                # Create strings with the correct code, amplitude, and phase
                if tide == 'O1':
                    O1 = f'{one:>2}{two:>2}{three:>2}{four:>2}{five:>2}{six:>2} {amp:0<7} {phz:<8}'
                    print_str = print_str + f'\\n{O1}'
                if tide == 'M2':
                    M2 = f'{one:>2}{two:>2}{three:>2}{four:>2}{five:>2}{six:>2} {amp:0<7} {phz:<8}'
                    print_str = print_str + f'\\n{M2}'
                if tide == 'K1':
                    K1 = f'{one:>2}{two:>2}{three:>2}{four:>2}{five:>2}{six:>2} {amp:0<7} {phz:<8}'
                    print_str = print_str + f'\\n{K1}'
                if tide == 'S2':
                    S2 = f'{one:>2}{two:>2}{three:>2}{four:>2}{five:>2}{six:>2} {amp:0<7} {phz:<8}'
                    print_str = print_str + f'\\n{S2}'
                if tide == 'N2':
                    N2 = f'{one:>2}{two:>2}{three:>2}{four:>2}{five:>2}{six:>2} {amp:0<7} {phz:<8}'
                    print_str = print_str + f'\\n{N2}'
                if tide == 'P1':
                    P1 = f'{one:>2}{two:>2}{three:>2}{four:>2}{five:>2}{six:>2} {amp:0<7} {phz:<8}'
                    print_str = print_str + f'\\n{P1}'
        fname = f'{sta_code}.{rate}.{ch}tides{str(start)[0:4]}{str(start)[5:7]}{str(start)[8:10]}.txt'
        # SPOTL will only allow up to 999999 values to be calculated
        # if the dataframe is longer than this, we will have to modify the number of
        # samples and sample interval
        if len(df) > 999999:
            nterms = 999999
            print('Warning, your data has more than 999999 samples, so tides will be interpolated.')
        else:
            nterms = len(df)
        # Finish the print statement
        tot_t = UTCDateTime(end)-UTCDateTime(start)
        samp = ceil(tot_t/nterms)
        print_str = print_str + f'\\n-1\' | docker run -i spotl {spotl} {start_dt} {nterms} {samp} > {dname}{fname} \n'
        # Run in SPOTL
        os.system("%s" % (print_str))
        # Make sure the files wrote
        if int(os.popen('wc -l '+dname+fname).read().split(  )[0]) > 0:
            print(ch+' tides saved to '+dname)
        else:
            print('The '+ch+' tides were not computed.')
    print('Loading tides to dataframe')
    tide = {'ch0':None,'ch1':None,'ch2':None,'ch3':None}
    for file in os.listdir('./SpotlTides'):
        if file.startswith(f'{sta_code}.{rate}.gauge') and file.endswith(f'{str(start)[0:4]}{str(start)[5:7]}{str(start)[8:10]}.txt'):
            n = file[13:14]
            tmpdf = pd.read_csv('./SpotlTides/'+file,header=None,names=['column'])
            # if the dataframe is longer than 999999, interpolate the tides to the dataframe times
            if len(df) > 999999:
                ind = []
                for i in range(0,999999):
                    ind.append(UTCDateTime(start) + samp*i)
                tmpdf.index = ind
                tide['ch'+n] = interp(tmpdf,df.datetime).values
            else:
                print(tmpdf)
                tide['ch'+n] = tmpdf.values
    df['tide_ch0'] = tide['ch0']/1000
    df['tide_ch1'] = tide['ch1']/1000
    df['tide_ch2'] = tide['ch2']/1000
    df['tide_ch3'] = tide['ch3']/1000

    return df, print_str

# Function to get the orientation matrix
def orient_mat(sta_meta):
    '''
    Computed the Moore-Penrose pseudo inverse from the orientation matrix in metadata
    '''
    a_loc = sta_meta['strain_xml']['inst_info']['processing']['bsm_processing_history'][-1]['bsm_processing']['orientation_matrix']
    a1 = np.array([float(a_loc['o11']),float(a_loc['o12']),float(a_loc['o13'])])
    a2 = np.array([float(a_loc['o21']),float(a_loc['o22']),float(a_loc['o23'])])
    a3 = np.array([float(a_loc['o31']),float(a_loc['o32']),float(a_loc['o33'])])
    a4 = np.array([float(a_loc['o41']),float(a_loc['o42']),float(a_loc['o43'])])
    a_mat = np.array([a1,a2,a3,a4])
    # Compute the Moore-Penrose pseudo inverse matrix with SciPy linalg module
    a_inv = linalg.pinv(a_mat)

    return a_inv

def transform(df,a_inv):
    '''Transforms strains with orientation matrix and '''
    gauge = regional_s(df['ch0 [ms]'],df['ch1 [ms]'],df['ch2 [ms]'],df['ch3 [ms]'],df,a_inv)
    gaugeEA = gauge[0]; gaugeES = gauge[1]; gaugeED = gauge[2]
    baros = regional_s(df['baro_ch0'],df['baro_ch1'],df['baro_ch2'],df['baro_ch3'],df,a_inv)
    baroEA = baros[0]; baroES = baros[1]; baroED = baros[2]
    tide = regional_s(df['tide_ch0'],df['tide_ch1'],df['tide_ch2'],df['tide_ch3'],df,a_inv)
    tideEA = tide[0]; tideES = tide[1]; tideED = tide[2]
    trend = regional_s(df['trend_ch0'],df['trend_ch1'],df['trend_ch2'],df['trend_ch3'],df,a_inv)
    trendEA = trend[0]; trendES = trend[1]; trendED = trend[2]
    cols = np.column_stack([df.datetime,gaugeEA,gaugeED,gaugeES,tideEA,tideED,tideES,baroEA,baroED,baroES,trendEA,trendED,trendES])
    moddf = pd.DataFrame(cols,columns=['datetime','gaugeEA','gaugeED','gaugeES','tideEA','tideED','tideES','baroEA','baroED','baroES','trendEA','trendED','trendES'])
    return moddf

# Function to transform ts by calibration matrix
def regional_s(z,y,x,w,df,a_inv):
    '''
    A function that takes 4 dataframe type columns (strain) and applies the orientation matrix to produce areal and shear strains.
    '''
    EA = []; ED = []; ES = []
    for i in range(0,len(df)):
        gauge_strain = np.array([z.iloc[i],y.iloc[i],x.iloc[i],w.iloc[i]])
        Ei = np.matmul(a_inv,gauge_strain)
        EA.append(Ei[0])
        ED.append(Ei[1])
        ES.append(Ei[2])
    return EA, ED, ES

def interp(df, new_index):
    """Return a new DataFrame with all columns values interpolated
    to the new_index values."""
    df_out = pd.DataFrame(index=new_index)
    df_out.index.name = df.index.name

    for colname, col in df.iteritems():
        df_out[colname] = np.interp(new_index, df.index, col) # default is a linear interpolation

    return df_out
