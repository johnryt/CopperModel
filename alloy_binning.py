import pandas as pd
import numpy as np

def alloy_binning(lower_mul = 0.95, upper_mul = 1.05, convergence_threshold = 0.005):
    
    # Primarily using prod_spec, with influences from prod_shape
    prod_spec1  = pd.read_excel("Specs\\Prod_spec_20190212.xlsx")
    cu_comp1 = pd.read_csv("Specs\\alloy_shape.csv").iloc[:,13:18].T 
    
    # Deleting some alloys that produce infeasible solutions
    for i in prod_spec1.index:
        if prod_spec1.loc[i,'Primary code'] == 'CC381H' or prod_spec1.loc[i,'Primary code'] == 'CW301G' or prod_spec1.loc[i,'Primary code'] == 'CW302G':
            prod_spec1.drop(index=i,inplace=True)
            
    # Computing mean compositions and storing them on a new column
    prod_spec1['Mean_Cu'] = prod_spec1.iloc[:,[3,4]].mean(axis=1)
    prod_spec1['Mean_Zn'] = prod_spec1.iloc[:,[5,6]].mean(axis=1)
    prod_spec1['Mean_Pb'] = prod_spec1.iloc[:,[7,8]].mean(axis=1)
    prod_spec1['Mean_Sn'] = prod_spec1.iloc[:,[9,10]].mean(axis=1)
    prod_spec1['Mean_Ni'] = prod_spec1.iloc[:,[11,12]].mean(axis=1)
    prod_spec1['Mean_Al'] = prod_spec1.iloc[:,[13,14]].mean(axis=1)
    prod_spec1['Mean_Mn'] = prod_spec1.iloc[:,[15,16]].mean(axis=1)
    prod_spec1['Mean_Fe'] = prod_spec1.iloc[:,[17,18]].mean(axis=1)
    
    # Initializing new DataFrame for recording compositions 
    shape_comp1 = pd.DataFrame(columns = ['Cu','Zn','Pb','Sn','Ni','Al','Mn','Fe'], index = ['Tube','RBS','PSS','Wire','Cast'])
    
    # Formatting the data for alloy production. 2010 is likely the most accurate year
    cu_comp1.columns = [2006,2007,2008,2009,2010]
    cu_comp1.index = shape_comp1.index
    
    # Averaging content by shape and saving to shape_comp
    # Pre-weight averaging for comparison
    for i in shape_comp1.index:
        for k in shape_comp1.columns:
            s = 'Mean_'
            shape_comp1.loc[i,k] = prod_spec1[prod_spec1['Category'] == i].loc[:,s+k].mean()/100

    # Creating weights based on the number of alloy suppliers
    prod_spec1.loc[:,'Weight 1'] = 0.001
    prod_spec1['ComplexLabel'] = prod_spec1['Primary code'] + ' ' + prod_spec1['Category']
    prod_spec1.set_index('ComplexLabel',inplace = True)
    prod = prod_spec1.copy()
    prod[['PSS','RBS','Wire','Tube','Cast']]

    for shape in shape_comp1.index:
        large = 0.7
        med = 0.2
        small = 0.1
        var = prod[prod['Category'] == shape]
        std = var[shape].std()
        mean = var[shape].mean()
        upper_bound = mean + 2*std
        lower_bound = mean + 0.5*std
        if shape == 'Tube':
            large = 0.3
            med = 0.4
            small = 0.3
            lower_bound = mean
        if shape == 'Wire':
            large = 0.3
            med = 0.4
            small = 0.3
            lower_bound = mean
        if shape == 'PSS':
            large = 0.5
            med = 0.3
            small = 0.2
    
        total0 = 0
        total1 = 0
        total2 = 0
    
        for i in prod.index:
            if prod.loc[i,shape] > upper_bound:
                total2 += prod.loc[i,shape]
            elif prod.loc[i,shape] < upper_bound and prod.loc[i,shape] > lower_bound:
                total1 += prod.loc[i,shape]
            elif prod.loc[i,shape] < lower_bound:
                total0 += prod.loc[i,shape]
    
        for i in prod.index:
            if prod.loc[i,shape] > upper_bound:
                prod.loc[i,'Weight 1'] = large * prod.loc[i,shape] / total2
            elif prod.loc[i,shape] < upper_bound and prod.loc[i,shape] > lower_bound:
                prod.loc[i,'Weight 1'] = med * prod.loc[i,shape] / total1
            elif prod.loc[i,shape] < lower_bound:
                if total0 != 0:
                    prod.loc[i,'Weight 1'] = small * prod.loc[i,shape] / total0
                else:
                    prod.loc[i,'Weight 1'] = 0.005
          
    for i in prod.index:
        if prod.loc[i,'Weight 1'] == 0:
            prod.loc[i,'Weight 1'] = 0.005
    prod.reset_index(inplace = True, drop=True)

    # Applying weights and computing new compositions        shape_comp = current value   cu_comp = target    prod_spec = alloy value
    prod_spec = prod.copy()
    shape_comp = shape_comp1.copy()
    cu_comp = cu_comp1.copy()
    for i in shape_comp.index:
        while abs(shape_comp.loc[i,'Cu'] - cu_comp.loc[i,2010]) > convergence_threshold:
            for j in prod_spec.index: 
                if prod_spec.loc[j,'Category'] == i and prod_spec.loc[j,'Mean_Cu']/100 > shape_comp.loc[i,'Cu'] and shape_comp.loc[i,'Cu'] > cu_comp.loc[i,2010]:
                    prod_spec.loc[j,'Weight 1'] *= lower_mul       
                elif prod_spec.loc[j,'Category'] == i and prod_spec.loc[j,'Mean_Cu']/100 > shape_comp.loc[i,'Cu'] and shape_comp.loc[i,'Cu'] < cu_comp.loc[i,2010]:
                    prod_spec.loc[j,'Weight 1'] *= upper_mul
                elif prod_spec.loc[j,'Category'] == i and prod_spec.loc[j,'Mean_Cu']/100 < shape_comp.loc[i,'Cu'] and shape_comp.loc[i,'Cu'] > cu_comp.loc[i,2010]:
                    prod_spec.loc[j,'Weight 1'] *= upper_mul
                elif prod_spec.loc[j,'Category'] == i and prod_spec.loc[j,'Mean_Cu']/100 < shape_comp.loc[i,'Cu'] and shape_comp.loc[i,'Cu'] < cu_comp.loc[i,2010]:
                    prod_spec.loc[j,'Weight 1'] *= lower_mul
            for k in shape_comp.columns:
                s = 'Mean_'
                shape_comp.loc[i,k] = (prod_spec[prod_spec['Category'] == i].loc[:,s+k] * prod_spec[prod_spec['Category'] == i].loc[:,'Weight 1'] / prod_spec[prod_spec['Category'] == i].loc[:,'Weight 1'].sum()/100).sum()
                

    # Normalizing weights
    for s in shape_comp.index:
        sum1 = prod_spec[prod_spec['Category']==s]['Weight 1'].sum()
        for i in prod_spec.index:
            if prod_spec.loc[i,'Category'] == s:
                prod_spec.loc[i,'Weight 1'] = prod_spec.loc[i,'Weight 1']/sum1   

    # Creating dataframe for UNS alloy names

    uns = prod_spec.copy().sort_values('Primary code', axis=0, ascending=True, inplace=False)
    uns.loc[:,'New'] = np.zeros(len(uns['UNS']))
    counter = 1
    
    # finding all the duplicate UNS alloys and giving them unique identifiers so I can sum their weights and average their compositions
    for i in uns.index:
        for j in uns.index:    
            if i != j and uns.loc[i,'UNS'] == uns.loc[j,'UNS'] and uns.loc[i,'Category'] == uns.loc[j,'Category'] and uns.loc[i,'New'] == 0:
                uns.loc[i,'New'] = counter
                uns.loc[j,'New'] = counter
                counter += 1 
            elif i != j and uns.loc[i,'UNS'] == uns.loc[j,'UNS'] and uns.loc[i,'Category'] == uns.loc[j,'Category'] and uns.loc[i,'New'] != 0:
                uns.loc[j,'New'] = uns.loc[i,'New']

    # summing weights and weighted averaging  the compositions of the duplicates. For the most part they duplicates have the
    # same weights, but for the few cases where one is weighted more, the UNS composition is weighted in its favor
    for i in range(1,counter):
        temp = uns[uns['New'] == i]
        x = temp.index
        uns.loc[x[0],'Weight 1'] = temp['Weight 1'].sum()
        for s in shape_comp.columns:
            uns.loc[x[0],'Mean_'+s] = (temp['Mean_'+s]*temp['Weight 1']).sum() / temp['Weight 1'].sum() #temp['Mean_'+s].sum()

        check = np.zeros(counter-1)
        temp2 = uns.copy()
        for i in uns.index:
            if pd.isnull(uns.loc[i,'UNS']):
                temp2.drop(index=i,inplace=True)
            for j in range(1,counter):
                if uns.loc[i,'New'] == j and check[j-1] == 0:
                    check[j-1] += 1
                elif uns.loc[i,'New'] == j and check[j-1] != 0:
                    temp2.drop(index=i,inplace=True)
    uns = temp2.copy()
    
    # Normalizing UNS weights
    for s in shape_comp.index:
        sum1 = uns[uns['Category']==s]['Weight 1'].sum()
        for i in uns.index:
            if uns.loc[i,'Category'] == s:
                uns.loc[i,'Weight 1'] = uns.loc[i,'Weight 1']/sum1

    # Cleaning up CEN dataframe
    new_cen = prod_spec.copy()
    #new_cen['Weight 1'] *= 100
    new_cen.rename(columns = {'Weight 1':'Weight'},inplace = True)
    for col in new_cen.columns:
        if 'Unnamed' in col:
            del new_cen[col]
            
    # Cleaning up UNS dataframe
    new_uns = uns.copy()
    #new_uns['Weight 1'] *= 100
    new_uns.rename(columns = {'Weight 1':'Weight'},inplace = True)
    del new_uns['New']
    for col in new_uns.columns:
        if 'Unnamed' in col:
            del new_uns[col]
    
    
    return(new_cen, new_uns)
    


