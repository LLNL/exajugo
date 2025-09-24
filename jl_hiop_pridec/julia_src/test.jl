

include("./test_cases.jl")
ptr=test_9bus()
prev_sol=solve_base_case(ptr)

w = solve_contingency_pridec(ptr, 1, prev_sol)

n = length(w[].cont_grad)
H = ones(n)

Pgrad = Ref(w[].cont_grad)
PH= Ref(H)

ptr_rderivaties=get_recourse_derivatives(w[].cont_grad, H, n)


solve_base_case_recourse(ptr, prev_sol, ptr_rderivaties)


