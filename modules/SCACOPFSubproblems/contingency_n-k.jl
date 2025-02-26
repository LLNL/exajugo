# function to generate a m-element random subset of 1:n
function randsubset(m::Int, n::Int, rng=Random.GLOBAL_RNG)
    m >= 0 || error("subset size must be nonnegative")
    n >= m || error("subset size larger than original set")
    m != 0 || return Int[]
    return sort(Random.randperm(rng, n)[1:m])
end

# function to generate a random contingency of n_elems elements
function random_n_elem_contingency(psd::SCACOPFdata, n_elems::Int, rng=Random.GLOBAL_RNG)
    n_elems >= 0 || error("number of elements to fail must be greater than zero.")
    n_elems != 0 || return GenericContingency()
    subset = randsubset(n_elems, nrow(psd.G) + nrow(psd.L) + nrow(psd.T), rng)
    g_out = Int[]
    l_out = Int[]
    i = 1
    while i <= n_elems && subset[i] <= nrow(psd.G)
        push!(g_out, subset[i])
        i += 1
    end
    while i <= n_elems && subset[i] <= nrow(psd.G) + nrow(psd.L)
        push!(l_out, subset[i] - nrow(psd.G))
        i += 1
    end
    if i <= n_elems && subset[i] <= nrow(psd.G) + nrow(psd.L) + nrow(psd.T)
        t_out = subset[i:end] .- (nrow(psd.G) + nrow(psd.L))
    else
        t_out = Int[]
    end
    return GenericContingency(g_out, l_out, t_out)
end

# function to sample from a binomial conditioned on a maximum number of successes
function n_successes_cond_binomial(n_trials::Int, n_max_successes::Int,
                                   p_success::Float64, rng=Random.GLOBAL_RNG)
    n_trials >= 0 || error("number of trials should be >= 0")
    n_max_successes >=  0 || error("maximum number of successes should be >= 0")
    d = Binomial(n_trials, p_success)
    n_max_successes < n_trials || return rand(d, 1)[1]
    break_points = Float64[]
    for i = 0:n_max_successes
        push!(break_points, pdf(d, i))
    end
    cumsum!(break_points, break_points)
    selector = break_points[end] * rand(rng)
    return findfirst(x -> selector < break_points[x], 1:(n_max_successes+1)) - 1
end

# function to generate a random contingency of up to n_elems elements
function random_up_to_n_elem_contingency(psd::SCACOPFdata, n_elems::Int,
                                         p_failure=0.05, rng=Random.GLOBAL_RNG)
    n_elems >= 0 || error("number of elements to fail must be greater than zero.")
    n_elems != 0 || return GenericContingency()
    n_fails = n_successes_cond_binomial(nrow(psd.G) + nrow(psd.L) + nrow(psd.T),
                                        n_elems, p_failure, rng)
    return random_n_elem_contingency(psd, n_fails, rng)
end
