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

    nested_test("untyped") do
        nested_test("metacells_gene_gene") do
            graph = metacells_gene_gene_graph(daf; x_gene = "A", y_gene = "B")
            @test graph.data.points_xs == Float32[0.1, 0.2, 0.3, 0.4]
            @test graph.data.points_ys == Float32[0.4, 0.3, 0.2, 0.1]
            @test graph.data.points_colors === nothing
            @test graph.data.points_hovers == ["M1", "M2", "M3", "M4"]
            @test graph.configuration.x_axis.log_scale == Log2Scale
            @test graph.configuration.y_axis.log_regularization == 1e-5
            return nothing
        end

        nested_test("blocks_gene_gene") do
            graph = blocks_gene_gene_graph(daf; x_gene = "A", y_gene = "C")
            @test graph.data.points_xs == Float32[0.15, 0.35]
            @test graph.data.points_ys == Float32[0.05, 0.05]
            @test graph.data.points_colors === nothing
            @test graph.data.points_hovers == ["B1", "B2"]
            return nothing
        end

        nested_test("metacells_umap") do
            graph = metacells_umap_graph(daf)
            @test graph.data.points_xs == Float32[0.0, 1.0, 2.0, 3.0]
            @test graph.data.points_ys == Float32[3.0, 2.0, 1.0, 0.0]
            @test graph.data.points_colors === nothing
            @test graph.data.points_hovers == ["M1", "M2", "M3", "M4"]
            @test !graph.configuration.x_axis.show_ticks
            @test !graph.configuration.y_axis.show_grid
            return nothing
        end

        nested_test("blocks_umap") do
            graph = blocks_umap_graph(daf)
            @test graph.data.points_xs == Float32[0.5, 2.5]
            @test graph.data.points_ys == Float32[2.5, 0.5]
            @test graph.data.points_colors === nothing
            @test graph.data.points_hovers == ["B1", "B2"]
            return nothing
        end
    end

    nested_test("typed") do
        add_axis!(daf, "type", ["X", "Y"])
        set_vector!(daf, "type", "color", ["red", "blue"])
        set_vector!(daf, "metacell", "type", ["X", "X", "Y", "Y"])
        set_vector!(daf, "block", "type", ["X", "Y"])

        nested_test("metacells_gene_gene") do
            graph = metacells_gene_gene_graph(daf; x_gene = "A", y_gene = "B")
            @test graph.data.points_colors == ["X", "X", "Y", "Y"]
            @test graph.configuration.points.colors.show_legend
            @test graph.configuration.points.colors.palette == ["red", "blue"]
            return nothing
        end

        nested_test("blocks_gene_gene") do
            graph = blocks_gene_gene_graph(daf; x_gene = "A", y_gene = "B")
            @test graph.data.points_colors == ["X", "Y"]
            @test graph.configuration.points.colors.show_legend
            return nothing
        end

        nested_test("metacells_umap") do
            graph = metacells_umap_graph(daf)
            @test graph.data.points_colors == ["X", "X", "Y", "Y"]
            @test graph.configuration.points.colors.show_legend
            @test graph.configuration.points.colors.palette == ["red", "blue"]
            return nothing
        end

        nested_test("blocks_umap") do
            graph = blocks_umap_graph(daf)
            @test graph.data.points_colors == ["X", "Y"]
            @test graph.configuration.points.colors.show_legend
            return nothing
        end
    end

    nested_test("regularization") do
        graph = metacells_gene_gene_graph(daf; x_gene = "A", y_gene = "B", gene_fraction_regularization = 1e-3)
        @test graph.configuration.x_axis.log_regularization == 1e-3
        @test graph.configuration.y_axis.log_regularization == 1e-3
        return nothing
    end
end
