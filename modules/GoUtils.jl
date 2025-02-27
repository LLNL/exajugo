#__precompile__()

module GoUtils

## elements to be exported

export indexsets, build_or_split_indexsets, redispatch

## load external modules

using DataFrames, Printf

function computeChecks(L_Nidx, T_Nidx, SSh_Nidx, G_Nidx, Lidxn, Lin, Tidxn, Tin, SShn, Gn)
    sCheck=""
    chk=0.
    for l=1:size(L_Nidx, 1), i=1:2
        chk += L_Nidx[l,i]/16.
    end
    sCheck=@sprintf("%s L_Nidx=%.5e", sCheck, chk)
    chk=0.;
    for t=1:size(T_Nidx, 1), i=1:2
        chk += T_Nidx[t,i]/16.
    end
    sCheck=@sprintf("%s T_Nidx=%.5e", sCheck, chk)

    chk = sum(SSh_Nidx/16)
    sCheck=@sprintf("%s SSh_Nidx=%.5e", sCheck, chk)
    chk = sum(G_Nidx/16)
    sCheck=@sprintf("%s G_Nidx=%.5e", sCheck, chk)

    chk=0.
    for n = 1:size(Lidxn, 1)
        chk += sum(Lidxn[n])/16.
    end
    sCheck=@sprintf("%s Lidxn=%.5e", sCheck, chk)
    chk=0.
    for n = 1:size(Lin, 1)
        chk += sum(Lin[n])/16.
    end
    sCheck=@sprintf("%s Lin=%.5e", sCheck, chk)

    chk=0.
    for n = 1:size(Tidxn, 1)
        chk += sum(Tidxn[n])/16.
    end
    sCheck=@sprintf("%s Tidxn=%.5e", sCheck, chk)
    chk=0.
    for n = 1:size(Tin, 1)
        chk += sum(Tin[n])/16.
    end
    sCheck=@sprintf("%s Tin=%.5e", sCheck, chk)

    chk=0.
    for n = 1:size(SShn, 1)
        chk += sum(SShn[n])/16.
    end
    sCheck=@sprintf("%s SShn=%.5e", sCheck, chk)

    chk=0.
    for n = 1:size(Gn, 1)
        chk += sum(Gn[n])/16.
    end
    sCheck=@sprintf("%s Gn=%.5e", sCheck, chk)

    return sCheck
end

## function for computing indices for efficient formulation
function indexsets(N::DataFrame, L::AbstractDataFrame, T::AbstractDataFrame, SSh::DataFrame,
	           G::AbstractDataFrame, K::Union{DataFrame, Tuple{Symbol, Int}})

	L_Nidx = hcat(convert(Vector{Int}, indexin(L[!, :From], N[!, :Bus])),
		convert(Vector{Int}, indexin(L[!, :To], N[!, :Bus])))
	T_Nidx = hcat(convert(Vector{Int}, indexin(T[!, :From], N[!, :Bus])),
		convert(Vector{Int}, indexin(T[!, :To], N[!, :Bus])))
	SSh_Nidx = convert(Vector{Int}, indexin(SSh[!, :Bus], N[!, :Bus]))
	G_Nidx = convert(Vector{Int}, indexin(G[!, :Bus], N[!, :Bus]))
	Lidxn = Vector{Vector{Int}}(undef, size(N, 1))
	Lin = Vector{Vector{Int}}(undef, size(N, 1))
	Tidxn = Vector{Vector{Int}}(undef, size(N, 1))
	Tin = Vector{Vector{Int}}(undef, size(N, 1))
	SShn = Vector{Vector{Int}}(undef, size(N, 1))
	Gn = Vector{Vector{Int}}(undef, size(N, 1))
	for n = 1:size(N, 1)
		Lidxn[n] = Int[]
		Lin[n] = Int[]
		Tidxn[n] = Int[]
		Tin[n] = Int[]
		SShn[n] = Int[]
		Gn[n] = Int[]
	end
	for l = 1:size(L, 1), i=1:2
		push!(Lidxn[L_Nidx[l,i]], l)
		push!(Lin[L_Nidx[l,i]], i)
	end
	for t = 1:size(T, 1), i=1:2
		push!(Tidxn[T_Nidx[t,i]], t)
		push!(Tin[T_Nidx[t,i]], i)
	end
	for s = 1:size(SSh, 1)
		push!(SShn[SSh_Nidx[s]], s)
	end
	for g = 1:size(G, 1)
		push!(Gn[G_Nidx[g]], g)
	end
	if typeof(K) == DataFrame
		K_outidx = Vector{Int}(undef, size(K, 1))
		Kgen = findall(K[!, :ConType] .== :Generator)
		Klin = findall(K[!, :ConType] .== :Line)
		Ktra = findall(K[!, :ConType] .== :Transformer)
		K_outidx[Kgen] .= convert(Vector{Int}, indexin(K[!, :IDout][Kgen], G[!, :Generator]))
		K_outidx[Klin] .= convert(Vector{Int}, indexin(K[!, :IDout][Klin], L[!, :Line]))
		K_outidx[Ktra] .= convert(Vector{Int}, indexin(K[!, :IDout][Ktra], T[!, :Transformer]))
	else
		if K[1] == :Generator
			K_outidx = Int(findfirst(G[!, :Generator] .== K[2]))
		elseif K[1] == :Line
			K_outidx = Int(findfirst(L[!, :Line] .== K[2]))
		elseif K[1] == :Transformer
			K_outidx = Int(findfirst(T[!, :Transformer] .== K[2]))
		else
			error("unrecognized contingency type ", K[1])
		end
	end

        # @show computeChecks(L_Nidx, T_Nidx, SSh_Nidx, G_Nidx, Lidxn, Lin, Tidxn, Tin, SShn, Gn)
	return L_Nidx, T_Nidx, SSh_Nidx, G_Nidx, Lidxn, Lin, Tidxn, Tin, SShn, Gn, K_outidx
end

function build_or_split_indexsets(IndexSets::Union{Tuple,Nothing},
                                  N::DataFrame, L::AbstractDataFrame, T::AbstractDataFrame, 
	                          SSh::DataFrame, G::AbstractDataFrame, K::Union{DataFrame, Tuple{Symbol,Int}})
    if IndexSets == nothing
        return indexsets(N, L, T, SSh, G, K)
    else
        K_outidx = IndexSets[11]
        if typeof(K) == DataFrame
            if size(K_outidx,1) != size(K,1)
                return indexsets(N, L, T, SSh, G, K)
            end
        else
            @assert typeof(K) == Tuple{Symbol, Int}
            return indexsets(N, L, T, SSh, G, K)
        end        
        return IndexSets[1], IndexSets[2], IndexSets[3], IndexSets[4], IndexSets[5], IndexSets[6], IndexSets[7], IndexSets[8], IndexSets[9], IndexSets[10], IndexSets[11]
    end
end

## re-dispatch function based on delta

function redispatch(delta::Float64, p0::T, alpha::T,
	Plb::T, Pub::T)::Float64 where {T <: AbstractVector{Float64}}
	if length(p0) != length(alpha) || length(alpha) != length(Plb) || length(Plb) != length(Pub)
		DimensionMismatch()
	end
	deltap = 0.0
	for g = 1:length(p0)
		p = p0[g] + delta*alpha[g]
		if p < Plb[g]
			deltap += Plb[g]
		elseif p > Pub[g]
			deltap += Pub[g]
		else
			deltap += p
		end
	end
	return deltap - sum(p0)
end

end
