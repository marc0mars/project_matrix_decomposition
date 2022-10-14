# CHEAT-SHEET STRUCTS

# simple struct
"""
struct Test 
    name::String
end
"""

# mutable struct
"""
mutable struct Test
    name::String
end
"""

# accessing fields
"""
mutable struct Test
    name::String
end

t = Test("Marco")
display(t.name)
# or
display(getfield(t, :name))
"""

# getting field names
"""
fieldnames(Test)
"""

# inner constructors
"""
mutable struct Test
    name::String

    function Test(name)
        new(name)
    end
    function Test()
        Test("Max Mustermann")
    end
end

display(Test("Peter"))
"""

# dynamic field types
"""
mutable struct Student{T}
    age::Int64
    label::T
    function Student(age::Int64, label::Any)
        new{typeof(label)}(age,label)
    end
    function Student()
        new{Int64}(0,0)
    end
end

display(Student(3,"test"))
display(Student())
"""

# Subtyping

abstract type SchoolMember end
mutable struct Student{T} <: SchoolMember # <---- Subtyping
    name::String
    age::Int64
    label::T

    function Student(name::String, age::Int64, label::Any)
        new{typeof(label)}(name,age,label)
    end
end