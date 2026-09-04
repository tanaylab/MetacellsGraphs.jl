nested_test("bar_graphs") do
    # Four base blocks, and a gene per case worth telling apart. The base correlation is what says which blocks a gene
    # is counted in: a zero there means it says nothing about that block.
    #
    # gene | correlated in | improved in | declined in | of the blocks it is counted in
    # A    | 4             | 3           | 0           | 75% improved
    # B    | 4             | 1           | 2           | 25% improved, 50% declined
    # C    | 2             | 0           | 2           | 100% declined
    # D    | 4             | 0           | 0           | neither, and a lateral gene
    # E    | 0             | -           | -           | no fraction at all, so never shown
    base_correlation = Float32[
        0.5 0.5 0.5 0.5
        0.5 0.5 0.5 0.5
        0.5 0.5 0.0 0.0
        0.5 0.5 0.5 0.5
        0.0 0.0 0.0 0.0
    ]
    correlation = Float32[
        0.6 0.6 0.6 0.5
        0.6 0.4 0.4 0.5
        0.4 0.4 0.6 0.6
        0.5 0.5 0.5 0.5
        0.6 0.6 0.6 0.6
    ]

    base_daf = MemoryDaf(; name = "base!")
    add_axis!(base_daf, "gene", ["A", "B", "C", "D", "E"])
    add_axis!(base_daf, "block", ["B1", "B2", "B3", "B4"])
    add_axis!(base_daf, "base_block", ["B1", "B2", "B3", "B4"])
    set_vector!(base_daf, "gene", "is_lateral", [false, false, false, true, false])
    set_vector!(base_daf, "gene", "is_regulator", [false, true, false, false, false])
    set_matrix!(
        base_daf,
        "gene",
        "base_block",
        "correlation_between_base_neighborhood_cells_and_punctuated_metacells",
        base_correlation,
    )

    daf = MemoryDaf(; name = "test!")
    add_axis!(daf, "gene", ["A", "B", "C", "D", "E"])
    add_axis!(daf, "base_block", ["B1", "B2", "B3", "B4"])
    set_matrix!(
        daf,
        "gene",
        "base_block",
        "correlation_between_base_neighborhood_cells_and_punctuated_metacells",
        correlation,
    )

    nested_test("improved") do
        graph = improved_genes_graph(; daf, base_daf)

        # The gene the graph is named for is at the top, which for horizontal bars is the end of the vector.
        @test graph.data.bars_names == ["C", "D", "B", "A"]
        @test graph.data.series_bars_values[1] == Float32[1.0, 0.0, 0.5, 0.0]
        @test graph.data.series_bars_values[2] == Float32[0.0, 0.0, 0.25, 0.75]
        @test graph.data.series_names == ["declined", "improved"]
        @test graph.configuration.mirrored
        @test graph.configuration.values_orientation == HorizontalValues
        @test graph.configuration.value_axis.percent
        return nothing
    end

    nested_test("declined") do
        graph = declined_genes_graph(; daf, base_daf)
        @test graph.data.bars_names == ["A", "D", "B", "C"]
        @test graph.data.series_bars_values[1] == Float32[0.0, 0.0, 0.5, 1.0]
        @test graph.data.series_bars_values[2] == Float32[0.75, 0.0, 0.25, 0.0]
        return nothing
    end

    nested_test("genes_count") do
        graph = improved_genes_graph(; daf, base_daf, genes_count = 2)
        @test graph.data.bars_names == ["B", "A"]
        return nothing
    end

    nested_test("annotations") do
        graph = improved_genes_graph(; daf, base_daf)
        @test length(graph.data.bars_annotations) == 2
        @test graph.data.bars_annotations[1].title == "is lateral"
        @test graph.data.bars_annotations[1].values == ["no", "yes", "no", "no"]
        @test graph.data.bars_annotations[2].title == "is regulator"
        @test graph.data.bars_annotations[2].values == ["no", "no", "yes", "no"]
        return nothing
    end
end
