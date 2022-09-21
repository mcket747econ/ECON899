function gini(dist)
    sorted = sort(dist)
    n = size(dist)
    coef = 2/n
    constant = (n+1)/n
    for i in enumerate(sorted)
        for yi in enumerate(sorted)
            weighted_sum = sum(i+1)*yi
        end
    end
    denom = sum(sorted)
    return coef*weighted_sum/(denom) - constant
end

array = collect(1:1.5:200)

print(sum(array))


gini(array)


