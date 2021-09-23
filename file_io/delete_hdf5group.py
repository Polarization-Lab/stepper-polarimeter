import h5py 

def delete_groups(fn, group_list):
    #first print current datasets
    print('current datasets present')
    with h5py.File(fn,'r') as hf:
        dataset_names = list(hf.keys())
    print(dataset_names)
    
    #next delete groups 
    for g in group_list:
        with h5py.File(fn, "a") as myfile:
            del myfile[g]
            try:
                myfile[g].value
            except KeyError as err:
            # print(err)
                pass
        
    #check deletion
    print('datasets after deletion')
    with h5py.File(fn,'r') as hf:
        dataset_names = list(hf.keys())
    print(dataset_names)
    
# EXAMPLE USE
#group_list =  [ '/Mueller1','/Mueller2','/MuellerMatrix']
#fn = '/Volumes/StepperData/Measurements/Waveguides/Option2/Space Measurements/Off Axis/-2deg/neg2deg_OffAxis_530nm-10-Sep-2021.h5'
#delete_groups(fn, group_list)
