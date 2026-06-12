# Génère toutes les combinaisons N_INTERIOR et écrit combinations.txt
max_int = 3   # max nœuds intérieurs par cycle
combos = vec(collect(Iterators.product(fill(0:max_int, 5)...)))
open("combinations_3.txt", "w") do io
    for c in combos
        println(io, join(c, " "))
    end
end
println("$(length(combos)) combinaisons générées")
