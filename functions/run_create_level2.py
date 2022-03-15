import create_level2

df, tdf = create_level2.input('B917', 'BS', 'PB',
                              '2019-06-19T00:00:00.000Z',
                              '2019-07-03T00:00:00.000Z')
print(df.head())
