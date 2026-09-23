nested_test("scatter_graphs") do
    daf = MemoryDaf(; name = "test!")

    add_axis!(daf, "gene", ["A", "B", "C"])
    add_axis!(daf, "metacell", ["M1", "M2", "M3", "M4"])
    add_axis!(daf, "block", ["B1", "B2"])

    # A zero fraction is what the regularization exists for, so one of each is included.
    set_matrix!(
        daf,
        "gene",
        "metacell",
        "linear_fraction",
        Float32[
            0.1 0.2 0.3 0.4
            0.4 0.3 0.2 0.1
            0.0 0.1 0.0 0.1
        ],
    )
    set_matrix!(
        daf,
        "gene",
        "block",
        "linear_fraction",
        Float32[
            0.15 0.35
            0.35 0.15
            0.05 0.05
        ],
    )

    set_vector!(daf, "metacell", "umap_x", Float32[0.0, 1.0, 2.0, 3.0])
    set_vector!(daf, "metacell", "umap_y", Float32[3.0, 2.0, 1.0, 0.0])
    set_vector!(daf, "block", "umap_x", Float32[0.5, 2.5])
    set_vector!(daf, "block", "umap_y", Float32[2.5, 0.5])

    nested_test("gene_gene") do
        nested_test("metacells") do
            graph = gene_gene_graph(daf; x_gene = "A", y_gene = "B")
            @test graph.data.x.vector == Float32[0.1, 0.2, 0.3, 0.4]
            @test graph.data.y.vector == Float32[0.4, 0.3, 0.2, 0.1]
            @test graph.configuration.x_axis.title == "A fraction"
            @test graph.configuration.y_axis.title == "B fraction"
            @test graph.configuration.x_axis.scale.log_base == Log2Base
            @test graph.configuration.y_axis.scale.log_regularization == 1e-5
            @test graph.data.points.entities.names == ["M1", "M2", "M3", "M4"]
            @test graph.data.points.entities.hovers[1] == "A fraction: 0.1<br>B fraction: 0.4"

            # The graph says nothing about the colors of the points; the caller does.
            @test graph.data.points.colors.vector === nothing
            return nothing
        end

        nested_test("blocks") do
            graph = gene_gene_graph(daf; axis = "block", x_gene = "A", y_gene = "C")
            @test graph.data.x.vector == Float32[0.15, 0.35]
            @test graph.data.y.vector == Float32[0.05, 0.05]
            @test graph.data.points.entities.names == ["B1", "B2"]
            @test graph.data.points.entities.hovers[2] == "A fraction: 0.35<br>C fraction: 0.05"
            return nothing
        end

        nested_test("entries") do
            graph = gene_gene_graph(daf; x_gene = "A", y_gene = "B", entries = ["M1", "M3"])
            @test graph.data.x.vector == Float32[0.1, 0.3]
            @test graph.data.y.vector == Float32[0.4, 0.2]
            @test graph.data.points.entities.names == ["M1", "M3"]
            return nothing
        end

        nested_test("regularization") do
            graph = gene_gene_graph(daf; x_gene = "A", y_gene = "B", gene_fraction_regularization = 1e-3)
            @test graph.configuration.x_axis.scale.log_regularization == 1e-3
            @test graph.configuration.y_axis.scale.log_regularization == 1e-3
            return nothing
        end
    end

    nested_test("umap") do
        nested_test("metacells") do
            graph = umap_graph(daf)
            @test graph.data.x.vector == Float32[0.0, 1.0, 2.0, 3.0]
            @test graph.data.y.vector == Float32[3.0, 2.0, 1.0, 0.0]
            @test graph.configuration.x_axis.title == "UMAP X"
            @test graph.configuration.y_axis.title == "UMAP Y"
            @test !graph.configuration.x_axis.show_ticks
            @test !graph.configuration.y_axis.show_grid
            @test graph.data.points.entities.names == ["M1", "M2", "M3", "M4"]

            # A UMAP coordinate is the arbitrary output of the projection, so it says nothing in a hover.
            @test graph.data.points.entities.hovers === nothing
            return nothing
        end

        nested_test("blocks") do
            graph = umap_graph(daf; axis = "block")
            @test graph.data.x.vector == Float32[0.5, 2.5]
            @test graph.data.y.vector == Float32[2.5, 0.5]
            @test graph.data.points.entities.names == ["B1", "B2"]
            return nothing
        end
    end

    # Coloring the points is a second call, so any per-entry property can be the color.
    nested_test("colors") do
        add_axis!(daf, "type", ["X", "Y"])
        set_vector!(daf, "type", "color", ["red", "blue"])
        set_vector!(daf, "metacell", "type", ["X", "X", "Y", "Y"])

        nested_test("type") do
            graph = umap_graph(daf)
            fill_type!(points_colors_vector_fields(graph), daf)
            @test graph.data.points.colors.vector == ["X", "X", "Y", "Y"]
            @test graph.configuration.points.colors.title == "type"
            @test graph.configuration.points.colors.palette == Dict("X" => "red", "Y" => "blue", "" => EMPTY_TYPE_COLOR)
            @test graph.configuration.points.colors.show_legend
            @test graph.data.points.entities.hovers == ["type: X", "type: X", "type: Y", "type: Y"]
            return nothing
        end

        nested_test("total_UMIs") do
            set_vector!(daf, "metacell", "total_UMIs", UInt32[100, 200, 300, 400])
            graph = umap_graph(daf)
            fill_total_UMIs!(points_colors_vector_fields(graph), daf)
            @test graph.data.points.colors.vector == UInt32[100, 200, 300, 400]
            @test graph.configuration.points.colors.scale.log_base == Log2Base
            return nothing
        end
    end

    nested_test("gene_base_delta_correlations") do
        # The base is the repository the blocks and the expression come from, and the other repository holds the same
        # genes correlated over the same base blocks. Gene `A` improves in `B1` and degrades in `B2`; gene `B` is
        # correlated in no base block at all, so it has nothing to show.
        add_axis!(daf, "base_block", ["B1", "B2"])
        set_matrix!(
            daf,
            "gene",
            "base_block",
            "correlation_between_base_neighborhood_cells_and_punctuated_metacells",
            Float32[
                0.4 0.6
                0.0 0.0
                0.5 0.5
            ],
        )

        base_daf = MemoryDaf(; name = "base!")
        add_axis!(base_daf, "gene", ["A", "B", "C"])
        add_axis!(base_daf, "block", ["B1", "B2"])
        add_axis!(base_daf, "base_block", ["B1", "B2"])
        set_matrix!(
            base_daf,
            "gene",
            "block",
            "linear_fraction",
            Float32[
                0.15 0.35
                0.35 0.15
                0.05 0.05
            ],
        )
        set_matrix!(
            base_daf,
            "gene",
            "base_block",
            "correlation_between_base_neighborhood_cells_and_punctuated_metacells",
            Float32[
                0.3 0.7
                0.0 0.0
                0.5 0.0
            ],
        )

        nested_test("()") do
            graph = gene_base_delta_correlations_graph(; daf, base_daf, gene = "A")
            @test graph.data.x.vector ≈ Float32[0.1, -0.1]
            @test graph.data.y.vector == Float32[0.15, 0.35]
            @test graph.configuration.x_axis.title == "correlation change"
            @test graph.configuration.y_axis.title == "A fraction"
            @test graph.configuration.y_axis.scale.log_base == Log2Base
            @test !graph.configuration.x_axis.scale.include_hidden
            @test !graph.configuration.y_axis.scale.include_hidden
            @test graph.data.points.entities.names == ["B1", "B2"]
            @test graph.data.points.entities.mask == [true, true]
            @test graph.data.points.colors.vector === nothing
            return nothing
        end

        # A gene correlated in only some of the base blocks is shown in those, since a zero base correlation is the
        # gene saying nothing about the block rather than a correlation of zero. The rest are hidden, not dropped.
        nested_test("uncorrelated") do
            graph = gene_base_delta_correlations_graph(; daf, base_daf, gene = "C")
            @test graph.data.x.vector == Float32[0.0, 0.5]
            @test graph.data.y.vector == Float32[0.05, 0.05]
            @test graph.data.points.entities.mask == [true, false]
            return nothing
        end

        nested_test("!correlated") do
            @test_throws "no base block correlates the gene: B\nof the base daf data: base!" gene_base_delta_correlations_graph(;
                daf,
                base_daf,
                gene = "B",
            )
            return nothing
        end
    end
end
