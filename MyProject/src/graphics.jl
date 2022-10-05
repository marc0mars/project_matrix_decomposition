using Luxor

function switch_numbers(x)
    return parse(Int64, reverse(string(x)))
end 

function guidelines()
    setline(0.5)
    setdash("dot")
    setcolor("red")
    arrow(Point(-250, 0), Point(250, 0))
    arrow(Point(0, -250), Point(0, 250))
    circle(Point(0,0), 200, :stroke)
    box(Point(0,0), 400, 400, :stroke)
    setdash("solid")
    setline(2)
    setcolor("black")
end

function main(filename)
    Drawing(600, 600, filename)
    origin()
    background("white")
    # guidelines()

    points = ngon(0, 0, 200, 100, -pi/2; vertices = true)
    setcolor("darkblue")
    setline(1)
    circle(0, 0, 200; action = :stroke)
    # c = circle.(points, 2; action = :fill)
    for i in 10:length(points)-1
        line(points[i], points[switch_numbers(i)], :stroke)
    end

    finish()
    preview()
end

main("palindrom-circle.svg")


