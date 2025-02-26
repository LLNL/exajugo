# method to generate a penalty function with limited slope, similar to the Huber loss but with
# non-zero slope at zero and smooth gradient (using a cubic function in between quadratic and linear
# parts).
#
# Inputs: - penalty function f(x) = a * x^2 + b * x; specified by a and b
#         - maximum allowed slope: max_slope
#         - interval size to make gradient smooth: delta

struct HuberLikePenalty
    a1::Float64
    b1::Float64
    a2::Float64
    b2::Float64
    c2::Float64
    d2::Float64
    a3::Float64
    b3::Float64
    x_break1::Float64
    x_break2::Float64
    function HuberLikePenalty(a::Float64, b::Float64,
                              max_slope::Float64=2.0 * b,
                              delta::Float64=0.1 * (max_slope - b) / (2.0 * a))
        a > 0.0 || error("penalty function should be strongly convex.")
        b >= 0.0 || error("slope should be nonnegative at 0.")
        max_slope > b || error("max slope smaller or equal than slope at 0.")
        delta > 0.0 || error("need positive delta to ensure hessian is continuous.")
        x_break1_ = (max_slope - b) / (2.0 * a)
        x_break2_ = x_break1_ + delta
        a1_ = a
        b1_ = b
        a2_ = -a / (3.0 * delta)
        b2_ = a * (x_break1_ / delta + 1.0)
        c2_ = (b - (a * x_break1_^2)/delta)
        d2_ = a / (3.0 * delta) * x_break1_^3
        a3_ = (-a / delta) * x_break2_^2 + (2.0 * a * (x_break1_ / delta + 1.0)) * x_break2_ + 
             (b - a/delta * x_break1_^2)
        b3_ = (a2_ * x_break2_^3 + b2_ * x_break2_^2 + c2_ * x_break2_ + d2_) - (a3_ * x_break2_)
        return new(a1_, b1_, a2_, b2_, c2_, d2_, a3_, b3_, x_break1_, x_break2_)
    end
end

function (h::HuberLikePenalty)(x::Float64)::Float64
    if x <= h.x_break1
        return h.a1 * x^2 + h.b1 * x
    elseif x <= h.x_break2
        return h.a2 * x^3 + h.b2 * x^2 + h.c2 * x + h.d2
    else
        return h.a3 * x + h.b3
    end
end

struct HuberLikePenaltyPrime
    a1::Float64
    b1::Float64
    a2::Float64
    b2::Float64
    c2::Float64
    a3::Float64
    x_break1::Float64
    x_break2::Float64
    function HuberLikePenaltyPrime(h::HuberLikePenalty)
        return new(h.a1, h.b1, h.a2, h.b2, h.c2, h.a3, h.x_break1, h.x_break2)
    end
end

function (h::HuberLikePenaltyPrime)(x::Float64)::Float64
    if x <= h.x_break1
        return 2.0 * h.a1 * x + h.b1
    elseif x <= h.x_break2
        return 3.0 * h.a2 * x^2 + 2.0 * h.b2 * x + h.c2
    else
        return h.a3
    end
end

struct HuberLikePenaltyPrimePrime
    a1::Float64
    a2::Float64
    b2::Float64
    x_break1::Float64
    x_break2::Float64
    function HuberLikePenaltyPrimePrime(h::HuberLikePenalty)
        return new(h.a1, h.a2, h.b2, h.x_break1, h.x_break2)
    end
end

function (h::HuberLikePenaltyPrimePrime)(x::Float64)::Float64
    if x <= h.x_break1
        return 2 * h.a1
    elseif x <= h.x_break2
        return 6 * h.a2 * x + 2 * h.b2
    else
        return 0.0
    end
end
