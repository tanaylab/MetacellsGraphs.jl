nested_test("bar_graphs") do
    # Four base blocks, and a gene per case worth telling apart. The base correlation is what says which blocks a gene
    # is counted in: a zero there means it says nothing about that block. `E` is not a marker, so it is not in the
    # report and is never shown.
    #
    # gene | correlated in | improved in | degraded in | of the blocks it is counted in
    # A    | 4             | 3           | 0           | 75% improved
    # B    | 4             | 1           | 2           | 25% improved, 50% degraded
    # C    | 2             | 0           | 2           | 100% degraded
    # D    | 4             | 0           | 0           | neither, and a lateral gene
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

    # The report names the most correlated marker of each kind, which this graph does not draw; any correlations will
    # do, as long as there are some.
    correlation_between_markers = Float32[
        1.0 0.9 0.8 0.7 0.6
        0.9 1.0 0.7 0.6 0.5
        0.8 0.7 1.0 0.5 0.4
        0.7 0.6 0.5 1.0 0.3
        0.6 0.5 0.4 0.3 1.0
    ]

    # What the hovers say. Only B is a regulator, so only its row is non-zero, and each side is asked about its own
    # base blocks: A shares a module with B where it improved, and C shares one with B where it degraded.
    improved_shared = Float32[
        0.0 0.0 0.0 0.0 0.0
        0.4 0.0 0.0 0.1 0.0
        0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0
    ]
    degraded_shared = Float32[
        0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.6 0.0 0.0
        0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0
    ]

    base_daf = MemoryDaf(; name = "base!")
    add_axis!(base_daf, "gene", ["A", "B", "C", "D", "E"])
    add_axis!(base_daf, "block", ["B1", "B2", "B3", "B4"])
    add_axis!(base_daf, "base_block", ["B1", "B2", "B3", "B4"])
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
    set_vector!(daf, "gene", "is_marker", [true, true, true, true, false])
    set_vector!(daf, "gene", "marker_rank", UInt32[1, 2, 3, 4, 5])
    set_vector!(daf, "gene", "is_lateral", [false, false, false, true, false])
    set_vector!(daf, "gene", "is_transcription_factor", [false, true, false, false, false])
    set_vector!(daf, "gene", "is_regulator", [false, true, false, false, false])
    set_vector!(daf, "gene", "is_skeleton", [true, false, false, false, false])
    set_matrix!(daf, "gene", "gene", "correlation_between_markers", correlation_between_markers)
    set_matrix!(
        daf,
        "gene",
        "base_block",
        "correlation_between_base_neighborhood_cells_and_punctuated_metacells",
        correlation,
    )
    set_vector!(
        daf,
        "gene",
        "mean_no_module_fraction_in_base_neighborhood_cells_at_improved_base_blocks",
        Float32[0.25, 0.5, 1.0, 0.3, 1.0],
    )
    set_vector!(
        daf,
        "gene",
        "mean_no_module_fraction_in_base_neighborhood_cells_at_degraded_base_blocks",
        Float32[1.0, 0.4, 0.2, 1.0, 1.0],
    )
    set_matrix!(
        daf,
        "gene",
        "gene",
        "mean_shared_module_fraction_in_base_neighborhood_cells_at_improved_base_blocks",
        improved_shared,
    )
    set_matrix!(
        daf,
        "gene",
        "gene",
        "mean_shared_module_fraction_in_base_neighborhood_cells_at_degraded_base_blocks",
        degraded_shared,
    )

    nested_test("improved") do
        graph = improved_genes_graph(; daf, base_daf)

        # The gene the graph is named for is at the top, which for horizontal bars is the end of the vector. `C` and `D`
        # both improved in nothing, so they are in marker rank order.
        @test graph.data.bars_names == ["C", "D", "B", "A"]
        @test graph.data.series_bars_values[1] == [1.0, 0.0, 0.5, 0.0]
        @test graph.data.series_bars_values[2] == [0.0, 0.0, 0.25, 0.75]
        @test graph.data.series_names == ["degraded", "improved"]
        @test graph.configuration.mirrored
        @test graph.configuration.values_orientation == HorizontalValues
        @test graph.configuration.value_axis.percent
        return nothing
    end

    nested_test("degraded") do
        graph = degraded_genes_graph(; daf, base_daf)
        @test graph.data.bars_names == ["A", "D", "B", "C"]
        @test graph.data.series_bars_values[1] == [0.0, 0.0, 0.5, 1.0]
        @test graph.data.series_bars_values[2] == [0.75, 0.0, 0.25, 0.0]
        return nothing
    end

    nested_test("genes_count") do
        graph = improved_genes_graph(; daf, base_daf, genes_count = 2)
        @test graph.data.bars_names == ["B", "A"]
        return nothing
    end

    nested_test("hovers") do
        graph = improved_genes_graph(; daf, base_daf)

        # The bars are ["C", "D", "B", "A"], and each wing says what its own side's base blocks say.
        @test graph.data.series_bars_hovers[1][1] == "C<br>degraded: in no module in 20.0% of the cells<br>- B: 60.0%"
        @test graph.data.series_bars_hovers[1][4] == "A<br>degraded: in no module in 100.0% of the cells"
        @test graph.data.series_bars_hovers[2][1] == "C<br>improved: in no module in 100.0% of the cells"
        @test graph.data.series_bars_hovers[2][4] == "A<br>improved: in no module in 25.0% of the cells<br>- B: 40.0%"
        return nothing
    end

    nested_test("regulators_count") do
        graph = improved_genes_graph(; daf, base_daf, regulators_count = 0)
        @test graph.data.series_bars_hovers[2][4] == "A<br>improved: in no module in 25.0% of the cells"
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
