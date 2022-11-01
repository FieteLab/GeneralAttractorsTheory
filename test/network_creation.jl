import GeneralAttractors.Networks: ShiftOperator

@testset begin
    m = 20
    v = rand(m)
    S = ShiftOperator(m, 1)
    Sᵀ = ShiftOperator(m, 0)

    right = S * v
    left = Sᵀ * v

    for i = 5:10
        @test right[i] == v[i-1]
        @test left[i] == v[i+1]
    end
end
