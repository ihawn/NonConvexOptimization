for b = 1:1000
    a = sqrt((b+14)^2+b^2)/2
    c = sqrt(2*a^2)

    if abs(a-floor(a)) <= 1e-10 && (b+14) + b > c
        println((b+14)*b/2)
    end
end
