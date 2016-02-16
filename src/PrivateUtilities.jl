function dictionary_to_variable_names(dict)
    s = ""
    if typeof(collect(keys(eval(dict)))[1]) <: AbstractString
        for i in keys(eval(dict))
            s = string(s ,"global const ",i, " = ", string(dict),"[\"",i,"\"]", " ; ")
        end
    else
        for i in keys(eval(dict))
            s = string(s ,i, " = ", string(dict),"[:",i,"]", " ; ")
        end
    end
    println(s)
    eval(parse(s)) # magic to turn keys of dict to variables 
    nothing
end



