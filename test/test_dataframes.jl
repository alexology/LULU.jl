@testset "DataFrames extension" begin
    df = DataFrame(OTU = ["D", "P", "Z"], a = [1, 10, 0], b = [2, 20, 0], c = [0, 5, 3])
    t = OTUTable(df)
    @test t.ids == ["D", "P", "Z"] && t.samples == ["a", "b", "c"]
    @test t.counts == [1 2 0; 10 20 5; 0 0 3]
    @test OTUTable(df[:, [:a, :OTU, :b, :c]]; id_col = :OTU).samples == ["a", "b", "c"]
    @test_throws ArgumentError OTUTable(DataFrame(id = ["x"], s = ["text"]))
    mdf = DataFrame(q = ["D", "Z"], h = ["P", "P"], m = [90, 90])
    @test MatchList(mdf) == MatchList([("D", "P", 90.0), ("Z", "P", 90.0)])
    @test_throws ArgumentError MatchList(DataFrame(q = ["D"], h = ["P"]))
    r = lulu(df, mdf)
    @test r.curated_otus == ["P"]
    @test lulu(df, MatchList(mdf)).otu_map.parent_id == r.otu_map.parent_id
    @test lulu(t, mdf).otu_map.parent_id == r.otu_map.parent_id
    out = DataFrame(r.curated_table)
    @test names(out) == ["OTUid", "a", "b", "c"] && out.a == [11]
    m = DataFrame(r.otu_map)
    @test names(m) == ["OTUid", "total", "spread", "parent_id", "curated", "rank"]
    @test m.curated == ["parent", "merged", "merged"]
    @test DataFrame(MatchList(mdf)).match == [90.0, 90.0]
    sd = DataFrame(sparse(t))
    @test sd == DataFrame(t) && sd.a isa Vector{Int}
    # Colliding names must fail before an ID or sample column is overwritten.
    colliding = OTUTable(["a", "b"], ["OTUid"], reshape([1, 2], 2, 1))
    @test_throws ArgumentError DataFrame(colliding)
    renamed = DataFrame(colliding; id_col = :feature)
    @test OTUTable(renamed; id_col = :feature) == colliding
    @test_throws ArgumentError DataFrame(colliding; id_col = "OTUid")
    duplicate = OTUTable(["a", "b"], ["s", "s"], [1 3; 2 4])
    @test_throws ArgumentError DataFrame(duplicate)
    @test_throws ArgumentError DataFrame(sparse(duplicate))
    @test colliding.ids == ["a", "b"]
    @test duplicate.counts == [1 3; 2 4]
end
