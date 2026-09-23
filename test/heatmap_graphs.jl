nested_test("heatmap_graphs") do
    daf = MemoryDaf(; name = "test!")

    add_axis!(daf, "gene", ["A", "B", "C", "D", "E"])
    add_axis!(daf, "metacell", ["M1", "M2", "M3", "M4"])
    add_axis!(daf, "block", ["B1", "B2"])

    set_vector!(daf, "metacell", "block", ["B1", "B1", "B2", "B2"])

    # The best marker is `A`, and `E` is no marker at all, so it is ranked above any count that can be asked for.
    set_vector!(daf, "gene", "marker_rank", UInt32[1, 2, 3, 4, typemax(UInt32)])
    set_vector!(daf, "gene", "is_skeleton", [true, true, false, false, false])

    set_matrix!(
        daf,
        "gene",
        "metacell",
        "log_linear_fraction",
        Float32[
            1.0 2.0 2.0 3.0
            3.0 2.0 2.0 1.0
            1.0 1.0 2.0 2.0
            2.0 2.0 2.0 2.0
            0.0 0.0 0.0 0.0
        ],
    )
    set_matrix!(
        daf,
        "gene",
        "block",
        "log_linear_fraction",
        Float32[
            1.5 2.5
            2.5 1.5
            1.0 2.0
            2.0 2.0
            0.0 0.0
        ],
    )

    nested_test("genes") do
        nested_test("markers") do
            genes = get_top_marker_gene_indices(daf; markers_count = 100)
            graph = genes_heatmap_graph(daf; genes, rows_axis_title = "Marker genes")
            @test graph.data.rows.entities.names == ["A", "B", "C", "D"]
            @test graph.data.columns.entities.names == ["M1", "M2", "M3", "M4"]
            @test graph.data.entries.matrix[1, :] == Float32[-1.0, 0.0, 0.0, 1.0]
            @test graph.configuration.rows.title == "Marker genes"
            @test graph.configuration.columns.title == "Metacells"
            @test graph.configuration.rows.reorder == OptimalHclust
            @test graph.configuration.columns.reorder == OptimalHclust
            @test graph.configuration.entries.colors.scale.minimum == -3
            @test graph.configuration.entries.colors.scale.maximum == 3
            @test graph.configuration.entries.colors.palette == "BuWtRd"

            # There are too many columns to label, so they are named in the hovers only.
            @test !graph.configuration.columns.show_ticks
            @test graph.configuration.rows.show_ticks

            # The graph groups nothing and annotates nothing; the caller does.
            @test isempty(graph.data.columns.annotations)
            @test graph.data.columns.groups.vector === nothing
            @test graph.data.columns.subgroups.vector === nothing
            return nothing
        end

        nested_test("markers_count") do
            genes = get_top_marker_gene_indices(daf; markers_count = 2)
            graph = genes_heatmap_graph(daf; genes)
            @test graph.data.rows.entities.names == ["A", "B"]
            @test graph.configuration.rows.title == "Genes"
            return nothing
        end

        nested_test("skeletons") do
            genes = get_skeleton_gene_indices(daf)
            graph = genes_heatmap_graph(daf; genes, rows_axis_title = "Skeleton genes")
            @test graph.data.rows.entities.names == ["A", "B"]
            @test graph.data.entries.matrix[2, :] == Float32[1.0, 0.0, 0.0, -1.0]
            return nothing
        end

        nested_test("blocks") do
            genes = get_top_marker_gene_indices(daf; markers_count = 100)
            graph = genes_heatmap_graph(daf; axis = "block", genes)
            @test graph.data.rows.entities.names == ["A", "B", "C", "D"]
            @test graph.data.columns.entities.names == ["B1", "B2"]
            @test graph.data.entries.matrix[1, :] == Float32[-0.5, 0.5]
            @test graph.configuration.columns.title == "Blocks"
            return nothing
        end

        # A single gene has nothing to be ordered against, so there is nothing to cluster.
        nested_test("one") do
            graph = genes_heatmap_graph(daf; genes = ["A"])
            @test graph.data.rows.entities.names == ["A"]
            @test graph.configuration.rows.reorder === nothing
            @test graph.configuration.columns.reorder == OptimalHclust
            return nothing
        end
    end

    # Grouping and annotating the columns is a second call, so any per-entry property can do either.
    nested_test("columns") do
        add_axis!(daf, "type", ["X", "Y"])
        set_vector!(daf, "type", "color", ["red", "blue"])
        set_vector!(daf, "metacell", "type", ["X", "X", "Y", "Y"])

        # The flow order runs the types in the reverse of the order the type axis holds them in.
        set_vector!(daf, "type", "global_flow_order", UInt32[2, 1])

        genes = get_skeleton_gene_indices(daf)

        nested_test("annotation") do
            graph = genes_heatmap_graph(daf; genes)
            index = add_columns_annotation!(graph)
            fill_type!(columns_annotations_colors_vector_fields(graph, index), daf; show_legend = false)
            @test graph.data.columns.annotations[1].values.vector == ["X", "X", "Y", "Y"]
            @test graph.data.columns.annotations[1].colors.title == "type"
            @test graph.data.columns.annotations[1].colors.palette ==
                  Dict("X" => "red", "Y" => "blue", "" => EMPTY_TYPE_COLOR)
            @test !graph.data.columns.annotations[1].colors.show_legend
            @test graph.data.columns.entities.hovers == ["type: X", "type: X", "type: Y", "type: Y"]
            return nothing
        end

        nested_test("group_by_type") do
            graph = genes_heatmap_graph(daf; genes)
            fill_type!(columns_groups_vector_data_fields(graph), daf)
            @test graph.data.columns.groups.vector == ["X", "X", "Y", "Y"]
            return nothing
        end

        nested_test("group_by_block") do
            graph = genes_heatmap_graph(daf; genes)
            fill_block!(columns_groups_vector_data_fields(graph), daf)
            @test graph.data.columns.groups.vector == ["B1", "B1", "B2", "B2"]
            @test graph.data.columns.entities.hovers == ["block: B1", "block: B1", "block: B2", "block: B2"]
            return nothing
        end

        # Naming both levels nests them: the blocks are grouped inside the types.
        nested_test("group_by_both") do
            graph = genes_heatmap_graph(daf; genes)
            fill_global_flow_order!(columns_groups_vector_data_fields(graph), daf)
            fill_block!(columns_subgroups_vector_data_fields(graph), daf)
            @test graph.data.columns.groups.vector == UInt32[2, 2, 1, 1]
            @test graph.data.columns.subgroups.vector == ["B1", "B1", "B2", "B2"]
            return nothing
        end
    end
end
