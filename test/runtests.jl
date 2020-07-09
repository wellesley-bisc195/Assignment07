using Assignment07
using Test

@testset "Assignment07" begin

@testset "Using Strings"
    
    @testset "normalizeDNA"
        @test normalizeDNA("aatgn") == "AATGN"
        @test_throws ErrorException normalizeDNA("ZCA")
        @test_throws ErrorException normalizeDNA(42)
        c = normalizeDNA('C') 
        @test c == "C"
        @test typeof(c) == String
    end # normalizeDNA

    @testset "composition"
        seq = rand(['A','T','G','C','N'], 20) |> join
        bc = composition(seq)
        @test bc isa Dict

        @test bc["A"] == count(x-> x == 'A', seq)
        @test bc["C"] == count(x-> x == 'C', seq)
        @test bc["G"] == count(x-> x == 'G', seq)
        @test bc["T"] == count(x-> x == 'T', seq)
        @test bc["N"] == count(x-> x == 'N', seq)

        bc = composition(lowercase(seq))

        @test bc["A"] == count(x-> x == 'A', seq)
        @test bc["C"] == count(x-> x == 'C', seq)
        @test bc["G"] == count(x-> x == 'G', seq)
        @test bc["T"] == count(x-> x == 'T', seq)
        @test bc["N"] == count(x-> x == 'N', seq)
    end # composition

    @tesetset "gc_content"
        @test gc_content("ANTG") == 0.25
        @test gc_content("cccggg") * 100 == 100.0
        @test gc_content("ATta") == 0.0
        @test_throws ErrorException gc_content("ATty")
    end # gc_content

    @testset "complement"
        @test complement("ATTAN") == "TAATN"
        @test complement("gcta") == "CGAT"
        @test complement("nnnnnnn") == "NNNNNNN"
        @test_throws ErrorException complement("ABC")
    end # complement

    @testset "reverse_complement"
        @test reverse_complement("ATTAN") == "NTAAT"
        @test reverse_complement("gcta") == "TAGC"
        @test reverse_complement("nnnnnnn") == "NNNNNNN"
        @test_throws ErrorException reverse_complement("ABC")
    end # reverse_complement

    @testset "parse_fasta"
        testpath = normpath(joinpath(@__DIR__, "..", "data"))
        genomes = joinpath(testpath, "cov2_genomes.fasta")
        ex1_path = joinpath(testpath, "ex1.fasta")
        ex2_path = joinpath(testpath, "ex2.fasta")

        ex1 = parse_fasta(ex1_path)
        @test ex1 isa Tuple
        @test all(x-> x isa Tuple, ex1[1])
        @test all(x-> x isa String, ex1[2])

        @test ex1[1] == ["ex1.1 | easy"), ("ex1.2 | multiline")]
        @test ex1[2] == ["AATTATAGC", "CGCCCCCCAGTCGGATT"]

        @test_throws parse_fasta(ex2_path)

        cov2 = parse_fasta(genomes)
        @test length(cov2[1]) == 8
        @test length(cov2[2]) == 8
    end #parse_fasta

end # strings

# @testset "Using BioSequences"
    
#     @testset "normalizeDNA"
#         @test normalizeDNA("aatgn") == dna"AATGN"
#         @test_throws ErrorException normalizeDNA("ZCA")
#         @test_throws ErrorException normalizeDNA(42)
#         c = normalizeDNA('C') 
#         @test c == dna"c"
#         @test c isa LongSequence
#     end #  normalizeDNA

#     @tesetset "gc_content"
#         @test gc_content(dna"ANTG") == 0.25
#         @test gc_content(dna"cccggg") * 100 == 100.0
#         @test gc_content(dna"ATta") == 0.0
#         @test_throws ErrorException gc_content(dna"ATty")
#     end #  composition

#     @testset "composition"
#         seq = rand(['A','T','G','C','N'], 20) |> join |> LongDNASeq
#         bc = composition(seq)

#         @test bc[DNA_A] == count(==, dna"A", seq)
#         @test bc[DNA_C] == count(==, dna"C", seq)
#         @test bc[DNA_G] == count(==, dna"G", seq)
#         @test bc[DNA_T] == count(==, dna"T", seq)
#         @test bc[DNA_N] == count(==, dna"N", seq)

#         @test bc[DNA_A] == count(==, dna"A", seq)
#         @test bc[DNA_C] == count(==, dna"C", seq)
#         @test bc[DNA_G] == count(==, dna"G", seq)
#         @test bc[DNA_T] == count(==, dna"T", seq)
#         @test bc[DNA_N] == count(==, dna"N", seq)
#     end #  gc_content

#     @testset "complement"
#         @test complement(dna"ATTAN") == dna"TAATN"
#         @test complement(dna"gcta") == dna"CGAT"
#         @test complement(dna"nnnnnnn") == dna"NNNNNNN"
#         @test_throws ErrorException complement(dna"ABC")
#     end #  complement

#     @testset "reverse_complement"
#         @test reverse_complement(dna"ATTAN") == dna"NTAAT"
#         @test reverse_complement(dna"gcta") == dna"TAGC"
#         @test reverse_complement(dna"nnnnnnn") == dna"NNNNNNN"
#         @test_throws ErrorException reverse_complement(dna"ABC")
#     end #  reverse_complement

#     @testset "parse_fasta"
#         testpath = normpath(joinpath(@__DIR__, "..", "data"))
#         genomes = joinpath(testpath, "cov2_genomes.fasta")
#         ex1_path = joinpath(testpath, "ex1.fasta")
#         ex2_path = joinpath(testpath, "ex2.fasta")

#         ex1 = parse_fasta(ex1_path)
#         @test ex1 isa Tuple
#         @test all(x-> x isa Tuple, ex1[1])
#         @test all(x-> x isa String, ex1[2])

#         @test ex1[1] == ["ex1.1 | easy"), ("ex1.2 | multiline")]
#         @test ex1[2] == [dna"AATTATAGC", dna"CGCCCCCCAGTCGGATT"]

#         @test_throws parse_fasta(ex2_path)

#         cov2 = parse_fasta(genomes)
#         @test length(cov2[1]) == 8
#         @test length(cov2[2]) == 8
#     end # parse_fasta

# end # BioSequences

end # Assignment07
