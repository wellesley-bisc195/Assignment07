using Assignment07
using Test

@testset "Assignment07" begin

@testset "Using Strings"
    
    @testset "normalizeDNA (assignment 4)"
        @test normalizeDNA("aatgn") == "AATGN"
        @test_throws ErrorException normalizeDNA("ZCA")
        @test_throws ErrorException normalizeDNA(42)
        c = normalizeDNA('C') 
        @test c == "C"
        @test typeof(c) == String
    end

    @testset "composition (assignment 4)"
        seq = rand(['A','T','G','C','N'], 20) |> join
        bc = composition(seq)
        @test bc isa Tuple{Int,Int,Int,Int}
        a,c,g,t = bc

        @test bc["A"] == count(x-> x == 'A', seq)
        @test bc["C"] == count(x-> x == 'C', seq)
        @test bc["G"] == count(x-> x == 'G', seq)
        @test bc["T"] == count(x-> x == 'T', seq)
        @test bc["N"] == count(x-> x == 'N', seq)

        a,c,g,t = composition(lowercase(seq))

        @test bc["A"] == count(x-> x == 'A', seq)
        @test bc["C"] == count(x-> x == 'C', seq)
        @test bc["G"] == count(x-> x == 'G', seq)
        @test bc["T"] == count(x-> x == 'T', seq)
        @test bc["N"] == count(x-> x == 'N', seq)

        @test_throws ErrorException composition(a)
    end





end

@testset "Using BioSequences"
    @testset "normalizeDNA (assignment 4)"
        @test normalizeDNA("aatgn") == dna"AATGN"
        @test_throws ErrorException normalizeDNA("ZCA")
        @test_throws ErrorException normalizeDNA(42)
        c = normalizeDNA('C') 
        @test c == dna"c"
        @test c isa LongSequence
    end

end

end # Assignment07
