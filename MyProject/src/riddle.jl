# this code helps me solving a riddle proposed here: https://www.youtube.com/watch?v=6hVPNONm7xw

S = ["a", "b", "c"]
println("______1.___2.___3._")
for i in S
    for j in S
        for k in S 
            println(i * j * k * " | " * j * k * " | " * i * k * " | " * i * j * " | ")
        end
    end
end