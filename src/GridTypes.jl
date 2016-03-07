 module GridTypes
export GridPatch, GridData

type GridPatch{T<:Real,N}
    LE::Array{T,N} # LeftEdge
    RE::Array{T,N} # RightEdge
    id::Integer    # identifier
    level::Integer # For AMR someday ... 
    dim::Tuple     # dimensions of grid patch
end 

""" Type to hold Grid data, 
"""
type GridData
    d::Dict{}
end

end
