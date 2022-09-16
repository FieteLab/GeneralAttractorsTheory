using Plots
using MultivariateStats

# d = 2

# u = rand() 
# # d=np.sum(u**2) **(0.5)
# # (x,y,z) = (u,v,w)/d

d = 4
u = rand(-1:.1:1, (d, 1000))
u = mapslices(
    x -> x/norm(x), u; dims=1
)


pca = fit(PCA, u; maxoutdim=3)
ū = predict(pca, u)

scatter(eachrow(ū)...)

