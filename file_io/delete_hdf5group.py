import h5py
import os 
import sys
import pandas as pd

#

def delete_groups(fn, group_list):
    #first print current datasets
    print('current datasets present')
    with h5py.File(fn,'r') as hf:
        dataset_names = list(hf.keys())
    print(dataset_names)
    
    #next delete groups 
    for g in group_list:
        with h5py.File(fn, "a") as myfile:
            if myfile[g]:
                print(g + " exists, deleting")
                del myfile[g]
            try:
                myfile[g].value
            except KeyError as err:
                print(err)
                pass
        
    #check deletion
    print('datasets after deletion')
    with h5py.File(fn,'r') as hf:
        dataset_names = list(hf.keys())
    print(dataset_names)
    

        
def print_groups(fn):
    #print current datasets
    print('file ' + fn + ' \n')
    with h5py.File(fn,'r') as hf:
        dataset_names = list(hf.keys())
        
        
        ds = []
        allsizes = []
        for d in dataset_names:
            names = []
            hf[d].visit(names.append)
            
            size = 0
            dataset_count = 0
           
            for n in names:
                if isinstance(hf[d+'/'+n], h5py.Dataset):
                    size += hf[d+'/'+n].size
                    dataset_count += 1
            
            allsizes.append(size)
            ds.append(dataset_count)
            
    info = pd.DataFrame({'group name': dataset_names,'size (bytes)': allsizes,'datasets': ds})
            
    print(info)
    print('\n')
    print("%i MB in %i datasets out of %i items in hdf5 file." % (sum(allsizes)*1e-6, sum(ds), len(names)))
  

def rename_group(fn,old_name,new_name):
    with h5py.File(fn,'r') as hf:
        group = hf[old_group]
        if group:
            print(old_name +" is accessible, changing link to " + new_name +'\n')
        else:
            print("group is inaccessible")
        try:   
            names = []
            hf[old_name].visit(names.append)
            
        
            for n in names:
                hf[old_name+'/'+n] = hf[new_name+'/'+n]
                
            print('datasets after replacement')
            dataset_names = list(hf.keys())
            print(dataset_names)
            
        except:
            print('copy error, closing file')
            
        #iterate through all subgroup names and replace
        
        
# EXAMPLE USE
old_group =  '/Mueller4'
new_group = '/Mueller0'
fn = '/Volumes/StepperData/Measurements/Waveguides/Option2/SpaceMeasurements/Run3_09102021/OnAxis_Option2_530nm-10-Sep-2021.h5'

#rename_group(fn,old_group,new_group)
print_groups(fn)