#
# CODE DUMP -> LEAVE COMPLETED RUN STATEMENTS OR WORKING CODE TO BE SAVED 
#


'''
returns = Parallel(n_jobs= 20)(delayed(subhalo_slope_analysis)(i) for i in valid_id)
df3=pd.DataFrame(returns,columns=['linslope','brokenslope'])

linearslopes = list(df3['linslope'])
brokenslopes = list(df3['brokenslope'])

print("Linear min: {}   MAX: {}".format(min(linearslopes),max(linearslopes)))
print("Broken min: {}   MAX: {}".format(min(brokenslopes),max(brokenslopes)))


BCUTILS.MSfilter(dfin,df2,'csv/tng99MAIN.csv')
Key Numbers: 

'''
