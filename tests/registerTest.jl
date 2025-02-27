using JuMP, LinearAlgebra
import Ipopt
import JuMP: add_to_expression!

function ft(n,z)
    zz = collect(z)
    zz = reshape(zz,(4,6))
    # println("type zz ", typeof(zz))
    xx = zz[:,1:2]
    yy = zz[:,3:6]
    # println("ft z ", z)
    # println(xx[1:2],"  y ", yy[3:4])
    fret = 0.5*sum(xx.^2)+0.5*sum(yy.^2);
    grad = 1.0.*z
    return fret,grad
end
function f(grad, z...)
    # println("f n= ",n)
    # println("z ", z)
    f, gra = ft(n,z)
    grad .= gra
    return f
end

function ∇f(grad, g, z...)
    # x = z[:,1:2]
    # y = z[:,3:6]
    nele = length(z)
    # println("grad ", grad)
    g .= grad
    # for i=1:(nele/4)
    #     idx = Int(4*(i-1)+1)
    #     g[idx:(idx+3)] .= z[idx:(idx+3)]
    # end
    # g .= z
    # println("g ", g)
    return g
end


function ∇²f(H, z...)
    nele = length(z)
    for i=1:nele
        for j=1:nele
            if i == j
                H[i,j] = 1.0
            end
        end
    end
    return H
end

n=5;
grad = ones(24)
f2(z...) = f(grad, z...)
∇f2(g, z...) = ∇f(grad, g, z...)
m = Model(Ipopt.Optimizer)

znorm = zero(JuMP.AffExpr)
# @variable(m, x[ts=1:4, v=1:2] >= 0)
# @variable(m, y[ts=1:4,v=1:4] >= 1)
@variable(m, z[ts=1:4, v=1:6] >= 0, start = 0.0)
println(start_value.(z))

@constraint(m, con[ts=1:4,v=1:3], z[ts,v] + z[ts,v+3] >= n*v)
# zn = @expression(m, norm(z,1))
# add_to_expression!(znorm, zn)
# register(m, :fobj, 24, f2; autodiff = true)
register(m, :fobj, 24, f2, ∇f2)
@NLobjective(m, Min, fobj(z...))

optimize!(m)

primal_z = JuMP.value.(z)
primal_x = primal_z[:,1:2]
primal_y = primal_z[:,3:6]

println("Primal Solution:")
println("x = ", primal_x)
println("y = ", primal_y)

normz = norm(primal_z, 1)
println("l1 norm z ", normz)
# map = relax_with_penalty!(m)
# Ipopt.GetIpoptCurrentViolations()
