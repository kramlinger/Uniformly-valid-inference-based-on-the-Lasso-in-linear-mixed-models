
using Random, LinearAlgebra, RDatasets, MixedModels, CategoricalArrays, Optim, CSV, Distributions



δβ = -4:0.1:4
Δ = length(δβ)
nsims = 10000

m = 10
n = 100
nᵢ = 10
p = 2

σᵥ₀ = 4
σₑ₀ = 4

Random.seed!(123)
X = randn(n, p)
λ = sqrt(n) / 2

Q = function(β, y, X, V⁻¹, λ)
    (y - X * β)' * V⁻¹ * (y - X * β) + 2 * λ * sum(abs.(β))
end

Coverage = zeros(4, Δ, Δ)
for i in 1:Δ
    for j in 1:Δ

        β₀ = [δβ[i], δβ[j]]
        fixed = X * β₀

        for k in 1:nsims
            # generate data
            v = σᵥ₀ * randn(m)
            y = fixed + repeat(v, inner = nᵢ) + σₑ₀ * randn(n)
            data = hcat(DataFrame(X, :auto), DataFrame(v = repeat(categorical(v), inner = nᵢ), y = y))

            # estimate model
            model = fit(MixedModel, @formula(y ~ -1 + x1 + x2 + (1 | v)), data)
            modelx1 = fit(MixedModel, @formula(y ~ -1 + x1 + (1 | v)), data)
            modelx2 = fit(MixedModel, @formula(y ~ -1 + x2 + (1 | v)), data)
            βₗₛ = coef(model)
            σᵥ = maximum([VarCorr(model).σρ.v.σ[1], 0.001])
            σₑ = maximum([VarCorr(model).s, 0.001])
            V⁻¹ = σₑ^-2 * (I - σᵥ^2 / (σᵥ^2 * nᵢ + σₑ^2) * kron(Diagonal(ones(m)), ones(nᵢ, nᵢ)))
            β = optimize(β -> Q(β, y, X, V⁻¹, λ), βₗₛ).minimizer
            C = X' * V⁻¹ * X / n
            C⁻¹ = inv(C)

            κ = λ^2 / n * maximum([sum(C⁻¹), [-1, 1]' * C⁻¹ * [-1, 1]])

            # test
            uₗₛ = sqrt(n) * (β₀-βₗₛ)
            u = sqrt(n) * (β₀-β)
            LASSO = u' * C * u
            WLS = uₗₛ' * C * uₗₛ #'

            WLSAIC = uₗₛ' * C * uₗₛ #'
            p′ = 2 #number of params

            model_choice = argmin(aic.([model, modelx1, modelx2]))
            if model_choice == 2 # full model
                u′ = sqrt(n) * (β₀[1] - coef(modelx1)[1])
                WLSAIC = u′^2 * C[1, 1]
                p′ = 1
            elseif  model_choice == 3
                u′ = sqrt(n) * (β₀[2] - coef(modelx2)[1])
                WLSAIC = u′^2 * C[2, 2]
                p′ = 1
            end

            Coverage[1, i, j] += (LASSO < quantile(NoncentralChisq(p, κ), 0.95)) / nsims
            Coverage[2, i, j] += (LASSO < quantile(Chisq(p), 0.95)) / nsims
            Coverage[3, i, j] += (WLS < quantile(Chisq(p), 0.95)) / nsims
            Coverage[4, i, j] += (WLSAIC < quantile(Chisq(p′), 0.95)) / nsims
        end

    end
end



# heatmap(δβ, δβ, Coverage[1,:,:])
# heatmap(δβ, δβ, Coverage[2,:,:])
# heatmap(δβ, δβ, Coverage[3,:,:])
# heatmap(δβ, δβ, Coverage[4,:,:])



CSV.write("LASSO.csv",  DataFrame(Coverage[1,:,:], :auto), writeheader = false)
CSV.write("LASSOWLS.csv",  DataFrame(Coverage[2,:,:], :auto), writeheader = false)
CSV.write("WLS.csv",  DataFrame(Coverage[3,:,:], :auto), writeheader = false)
CSV.write("WLSAIC.csv",  DataFrame(Coverage[4,:,:], :auto), writeheader = false)
