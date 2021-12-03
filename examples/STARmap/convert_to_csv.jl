#! /usr/bin/env julia

using CSV
using DataFrames
using MAT
using Statistics

@info "Read $(ARGS[1])"

file = matopen(ARGS[1]); # "goodPoints.mat"

good_points = read(file, "goodPoints");
good_bases = vec(read(file, "goodBases"));

close(file)

@info "Read $(ARGS[2])"

gene_per_base = CSV.read(ARGS[2], DataFrame, header=0); # "cell_barcode_names.csv"
gene_per_base = Dict("$(row[2])" => row[3] for row in eachrow(gene_per_base));

df_spatial = DataFrame(good_points, [:x, :y, :z]);
df_spatial[!, :gene] = [gene_per_base[b] for b in good_bases];

CSV.write("molecules.csv", df_spatial)

if length(ARGS) > 2
    import NPZ
    import Images

    @info "Read $(ARGS[3])"
    Images.save("segmentation.tiff", UInt16.(NPZ.npzread(ARGS[3])["labels"]))
end
