
function summary(outdir)
    data = [JSON.parse(String(read(fname))) for fname in sort(glob("*json",outdir))];
    df = DataFrame(name = getindex.(data,"name"), validation = getindex.(data,"validation"));
    CSV.write(joinpath(outdir,"summary.csv"),df)
    score = mean(df.validation)
    @show score
    return mean(df.validation)
end
