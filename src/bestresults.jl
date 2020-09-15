"""
# Example
resdir = joinpath(outdir,"Results")
DIVAndNN.bestresults(resdir)
"""
function bestresults(resdir)
    files =  sort(glob("*/DIVAndNN.json",resdir));
    minval,ibest = findmin(getindex.(JSON.parsefile.(files),"validation"))
    return files[ibest],minval
end
