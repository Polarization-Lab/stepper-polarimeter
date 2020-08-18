%this is an example of the functions used for the DMM
exposure =0.5;

prepare_DMM(COMdmm,exposure)
trigger_DMM(COMdmm)
pause(exposure)
read_DMM(COMdmm,exposure)