using DataStructures

function test()
    t = [1; 3.1; 4]
    n = length(t)
    k = 1
    while (k <= n) && (t[k] > 0)
        k += 1
    end
    print("ok")
end

function test2()
    t = [(1, 1.2); (3.1, 1); (4, 3.2)]
    t2 = [(0, 5.1); (6.1, 2)]
    append!(t, t2)
    print(t)
end

function test3()
    t = [[1; 1.2]; [3.1; 1]; [4; 3.2]]
    t2 = [[0; 5.1]; [6.1; 2]]
    append!(t, t2)
    print(t)
end

function test4()
    t = [(1, 1.2); (3.1, 1); (4, 3.2)]
    sorted_t = sort(t, by=x->x[2])
    print(sorted_t)
end

function test5()
    print(abs(max(-1, -3)))
end

function test6()
    q = Queue{Int}(2)
    enqueue!(q, 1)
    enqueue!(q, 3)
    print(q)
end

#test()
#test2()
#test3()
#test4()
#test5()
test6()
