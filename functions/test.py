import bsm_functions as bsm
import numpy as np
from obspy.clients.iris import Client
client = Client()
from obspy.clients.fdsn import Client
inv_client = Client('IRIS')

odf = bsm.get_data('B917','RS','PB',
                  '2019-06-19T00:00:00.000Z',
                  '2019-07-02T23:50:00.000Z')

df, comment = bsm.linearize('B917',odf)

def baro_data(network, sta_code, start, end):
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
        p_array = np.append(fillme,baro[i].data*sf)
    with open(sta_code+'.barodata.txt','w') as f:
        f.write('Atmospheric Pressure\n')
        for i in range(0,len(p_array)):
            f.write(str(p_array[i]-np.mean(p_array[i]))+'\n')
    f.close()
    return f'Done! {sta_code}.barodata.txt written.'

baro_data('PB', 'B917','2019-06-19T00:00:00.000Z',
                '2019-07-02T23:50:00.000Z')
