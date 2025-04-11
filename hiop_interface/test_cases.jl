
include("hiop.jl")


function test_500bus()

    prob="/p/lustre1/santiago/COSMIN/GITHUB/exajugo/examples/500bus/";
    rawfile=joinpath(prob,"case.raw"); confile=joinpath(prob,"case.con"); ropfile=joinpath(prob,"case.rop");
    inlfile=joinpath(prob,"case.inl")

    return load_ACOPF_dir(prob)

end

function test_9bus()

    prob="/p/lustre1/santiago/COSMIN/GITHUB/exajugo/examples/9bus/";
    rawfile=joinpath(prob,"case.raw"); confile=joinpath(prob,"case.con"); ropfile=joinpath(prob,"case.rop");

    return load_ACOPF(rawfile, ropfile, confile)

end


