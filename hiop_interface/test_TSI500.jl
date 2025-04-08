
include("hiop.jl")
include("test_TSACOPF.jl")

PATH_TO_DATA="./data"

model_file = joinpath(PATH_TO_DATA,"model_state_DSPP_500_0.990_0.884_28.9_9.0_4.2_20000_318.0_300_6_0.021_0.00075.pth")
data_record = joinpath(PATH_TO_DATA,"data_record.mat")

ptrTSIGPmodel = load_tsmodel(model_file, data_record)

pf_limit_file = joinpath(PATH_TO_DATA,"pf_new.mat")

case_path=joinpath(PATH_TO_DATA,"ACTIVSg500")

test_TSI(ptrTSIGPmodel, case_path, pf_limit_file)
