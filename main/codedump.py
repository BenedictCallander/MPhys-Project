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

'''
returns = Parallel(n_jobs= 20)(delayed(subhalo_analysis)(i) for i in valid_id)
df2=pd.DataFrame(returns,columns=['slope','met','id','sfr','inside','outside'])
df2.insert(5,'mass', dfin['mass'],True)
df2.dropna()
df2.to_csv("csv/tng33MSslopes.csv")
'''
    def broken_fit(self,dfin,breakpoint,pc):
        r'''
        Pseudocode
        Inputs:

        Dataframe (can be inherited from subhalo object)
        breakpoint (the point at which break in linear fit is placed)
        '''
        dfin.sort_values(by='rad',inplace = True)
        med_data1 = medfilt(12+np.log10(dfin['met2']), kernel_size=21)
        x0 = np.array([min(dfin['rad']), breakpoint, max(dfin['rad'])])
        break1 = list(dfin[dfin['rad']<breakpoint])
        break2 = list(dfin[dfin['rad']>=breakpoint])
        my_pwlf = pwlf.PiecewiseLinFit(dfin['rad'], 12+np.log10(dfin['met2']),weights=1/dfin['sfr'])
        my_pwlf.fit_with_breaks(x0)
        slope1 = my_pwlf.slopes[0]
        slope2 = my_pwlf.slopes[1]
        
        xHat = np.linspace(min(dfin['rad']), max(dfin['rad']), num=10000)
        yHat = my_pwlf.predict(xHat)
        '''
        plt.figure(figsize=(20,12))
        plt.plot(dfin['rad'], med_data1, 'b--')
        plt.plot(xHat,yHat, 'g-')
        plt.xlabel("Radius (Normalised Code Units)")
        plt.ylabel("12+$log_{10}$ $(O/H)$")
        filename = 'brfit/sub_{}_break_snapshot_{}.png'.format(self.subID, self.snapID)
        plt.savefig(filename)
        plt.close()
        '''
        return(slope1,slope2)