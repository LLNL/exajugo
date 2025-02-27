#__precompile__()

module SmoothApproximations

## elements to be exported

export smoothstep, map2step, smoothclamp

## constants

const DEFAULTMU = 0.0001	# pu

## create smoothstep function and derivatives

function step(x::T)::T where {T <: Real}
	if x < 0
		return 0
	else
		return 1
	end
end

function smoothstep(x::T, mu::Real=DEFAULTMU)::T where {T <: Real}
	if x <= -mu
		return 0
	elseif x >= mu
		return 1
	else
		return 6/(32*mu^5)*x^5-5/(8*mu^3)*x^3+15/(16*mu)*x+1/2
	end
end

function smoothstepprime(x::T, mu::Real=DEFAULTMU)::T where {T <: Real}
	if x <= -mu
		return 0
	elseif x >= mu
		return 0
	else
		return 15/(16*mu^5)*x^4-15/(8*mu^3)*x^2+15/(16*mu)
	end
end

function smoothstepprimeprime(x::T, mu::Real=DEFAULTMU)::T where {T <: Real}
	if x <= -mu
		return 0
	elseif x >= mu
		return 0
	else
		return 15/(4*mu^5)*x^3-15/(4*mu^3)*x
	end
end

## define function mapping pairs (x, y) of the smoothstep curve
## onto the exact step curve -- note we only need x, y is implied
## NOTE: the mapping is performed following the normal vector to
## the smoothstep function

function map2step(x::T, mu::Real=DEFAULTMU)::Tuple{T, T} where {T <: Real}
	
	# compute y coordinate and derivative
	y = smoothstep(x, mu)
	fprime = smoothstepprime(x, mu)
	
	# try each piece of the step function
	gamma = -x - y*fprime
	if gamma >= 0
		return -gamma, 0
	end
	gamma = x + (y-1)*fprime
	if gamma >= 0
		return gamma, 1
	end
	gamma = x/fprime + y
	if 0 <= gamma && gamma <= 1
		return 0, gamma
	end
	
	# if we arrive here there is a problem ...
	error("unable to determine projection of (", x, ", ", y, ") onto unitary step.")
	
end

## define smoothclamp function and derivatives

function clamp(x::T)::T where {T <: Real}
	if x < 0
		return 0
	elseif x >= 0 && x < 1
		return x
	else
		return 1
	end
end

function smoothclamp(x::T, mu::Real=DEFAULTMU)::T where {T <: Real}
	if x <= -mu
		return 0
	elseif -mu < x < mu
		return -1/(16*mu^3)*x^4+3/(8*mu)*x^2+1/2*x+3*mu/16
	elseif mu <= x <= 1-mu
		return x
	elseif 1-mu < x < 1+mu
		return 1/(16*mu^3)*(x-1)^4-3/(8*mu)*(x-1)^2+1/2*(x-1)-3*mu/16+1
	else
		return 1
	end
end

function smoothclampprime(x::T, mu::Real=DEFAULTMU)::T where {T <: Real}
	if x <= -mu
		return 0
	elseif -mu < x < mu
		return -1/(4*mu^3)*x^3+3/(4*mu)*x+1/2
	elseif mu <= x <= 1-mu
		return 1
	elseif 1-mu < x < 1+mu
		return 1/(4*mu^3)*(x-1)^3-3/(4*mu)*(x-1)+1/2
	else
		return 0
	end
end

function smoothclampprimeprime(x::T, mu::Real=DEFAULTMU)::T where {T <: Real}
	if x <= -mu
		return 0
	elseif -mu < x < mu
		return -3/(4*mu^3)*x^2+3/(4*mu)
	elseif mu <= x <= 1-mu
		return 0
	elseif 1-mu < x < 1+mu
		return 3/(4*mu^3)*(x-1)^2-3/(4*mu)
	else
		return 0
	end
end

end
