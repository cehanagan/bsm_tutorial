import create_level2
from obspy import UTCDateTime
# Parameters for strainmeter data
foredt = '2019-07-04T17:33:49.040000Z'
maindt = '2019-07-06T03:19:53.040000Z'
start = foredt; end = str(UTCDateTime(UTCDateTime(foredt)+11*24*60*60))
rate = 'LS' # LS is 1 Hz data
sta = 'B921'

df,tdf = create_level2.input(sta, rate, 'PB',
                              start,
                              end)

print(df.head())
