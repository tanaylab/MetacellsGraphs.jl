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

    nested_test("untyped") do
        nested_test("markers_metacells") do
            graph = markers_metacells_heatmap_graph(daf)
            @test graph.data.rows_names == ["A", "B", "C", "D"]
            @test graph.data.entries_values[1, :] == Float32[-1.0, 0.0, 0.0, 1.0]
            @test graph.data.columns_hovers[1] == "metacell: M1<br>block: B1"
            @test isempty(graph.data.columns_annotations)
            @test graph.data.columns_groups === nothing
            @test graph.configuration.entries_colors.axis.maximum == 3
            return nothing
        end

        nested_test("markers_count") do
            graph = markers_metacells_heatmap_graph(daf; markers_count = 2)
            @test graph.data.rows_names == ["A", "B"]
            return nothing
        end

        nested_test("skeletons_metacells") do
            graph = skeletons_metacells_heatmap_graph(daf)
            @test graph.data.rows_names == ["A", "B"]
            @test graph.data.entries_values[2, :] == Float32[1.0, 0.0, 0.0, -1.0]
            return nothing
        end

        nested_test("markers_blocks") do
            graph = markers_blocks_heatmap_graph(daf)
            @test graph.data.rows_names == ["A", "B", "C", "D"]
            @test graph.data.entries_values[1, :] == Float32[-0.5, 0.5]
            @test graph.data.columns_hovers == ["block: B1", "block: B2"]
            return nothing
        end

        nested_test("skeletons_blocks") do
            graph = skeletons_blocks_heatmap_graph(daf)
            @test graph.data.rows_names == ["A", "B"]
            @test graph.data.columns_hovers == ["block: B1", "block: B2"]
            return nothing
        end

        nested_test("group_by_block") do
            graph = markers_metacells_heatmap_graph(daf; group_by_block = true)
            @test graph.data.columns_groups == ["B1", "B1", "B2", "B2"]
            @test graph.data.columns_subgroups === nothing
            return nothing
        end

        nested_test("!group_by_type") do
            @test_throws "group_by_type without a type per metacell" markers_metacells_heatmap_graph(
                daf;
                group_by_type = true,
            )
            return nothing
        end
    end

    nested_test("colorless") do
        add_axis!(daf, "type", ["X", "Y"])
        set_vector!(daf, "metacell", "type", ["X", "X", "Y", "Y"])

        # A type with no color to draw it in is not shown, but it is still there to be grouped by.
        graph = markers_metacells_heatmap_graph(daf; group_by_type = true)
        @test isempty(graph.data.columns_annotations)
        @test graph.data.columns_groups == ["X", "X", "Y", "Y"]
        @test graph.data.columns_hovers[1] == "metacell: M1<br>block: B1<br>type: X"
        return nothing
    end

    nested_test("typed") do
        add_axis!(daf, "type", ["X", "Y"])
        set_vector!(daf, "type", "color", ["red", "blue"])
        set_vector!(daf, "metacell", "type", ["X", "X", "Y", "Y"])
        set_vector!(daf, "block", "type", ["X", "Y"])

        # The flow order runs the types in the reverse of the order the type axis holds them in.
        set_vector!(daf, "type", "global_flow_order", UInt32[2, 1])

        nested_test("markers_metacells") do
            graph = markers_metacells_heatmap_graph(daf)
            @test graph.data.columns_hovers[1] == "metacell: M1<br>block: B1<br>type: X"
            @test length(graph.data.columns_annotations) == 1
            @test graph.data.columns_annotations[1].title == "type"
            @test graph.data.columns_annotations[1].values == ["X", "X", "Y", "Y"]
            @test graph.data.columns_annotations[1].colors.palette == ["red", "blue"]
            return nothing
        end

        nested_test("group_by_type") do
            graph = markers_metacells_heatmap_graph(daf; group_by_type = true)
            @test graph.data.columns_groups == ["X", "X", "Y", "Y"]
            @test graph.data.columns_subgroups === nothing
            return nothing
        end

        nested_test("group_by_both") do
            graph = markers_metacells_heatmap_graph(daf; group_by_type = true, group_by_block = true)
            @test graph.data.columns_groups == ["X", "X", "Y", "Y"]
            @test graph.data.columns_subgroups == ["B1", "B1", "B2", "B2"]
            return nothing
        end

        nested_test("use_global_flow_order") do
            graph = markers_metacells_heatmap_graph(daf; group_by_type = true, use_global_flow_order = true)
            @test graph.data.columns_groups == UInt32[2, 2, 1, 1]
            return nothing
        end

        nested_test("!use_global_flow_order") do
            @test_throws "use_global_flow_order without group_by_type" markers_metacells_heatmap_graph(
                daf;
                use_global_flow_order = true,
            )
            return nothing
        end

        nested_test("markers_blocks") do
            graph = markers_blocks_heatmap_graph(daf; group_by_type = true)
            @test graph.data.columns_hovers == ["block: B1<br>type: X", "block: B2<br>type: Y"]
            @test graph.data.columns_groups == ["X", "Y"]
            return nothing
        end

        nested_test("skeletons_blocks") do
            graph = skeletons_blocks_heatmap_graph(daf; group_by_type = true, use_global_flow_order = true)
            @test graph.data.columns_groups == UInt32[2, 1]
            return nothing
        end
    end
end
