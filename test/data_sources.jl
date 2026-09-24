nested_test("data_sources") do
    daf = MemoryDaf(; name = "test!")

    add_axis!(daf, "gene", ["A", "B", "C"])
    add_axis!(daf, "metacell", ["M1", "M2", "M3", "M4"])
    add_axis!(daf, "block", ["B1", "B2"])

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
        "metacell",
        "log_linear_fraction",
        Float32[
            1.0 2.0 2.0 3.0
            3.0 2.0 2.0 1.0
            0.0 0.0 0.0 0.0
        ],
    )
    set_vector!(daf, "metacell", "umap_x", Float32[0.0, 1.0, 2.0, 3.0])
    set_vector!(daf, "metacell", "umap_y", Float32[3.0, 2.0, 1.0, 0.0])
    set_vector!(daf, "metacell", "block", ["B1", "B1", "B2", "B2"])

    nested_test("axis_names") do
        graph = points_graph()

        nested_test("()") do
            fill_axis_names_data!(points_entities(graph), daf; axis = "metacell")
            @test graph.data.points.entities.names == ["M1", "M2", "M3", "M4"]
            @test graph.data.points.entities.hovers === nothing
            return nothing
        end

        nested_test("entries") do
            @test get_axis_entries_vector(daf; axis = "metacell", entries = [1, 3]) == ["M1", "M3"]
            @test get_axis_entries_vector(daf; axis = "metacell", entries = ["M2"]) == ["M2"]
            return nothing
        end

        # Naming the same entities through several sinks names them once.
        nested_test("shared") do
            put_vector_names_data!(
                [x_axis_vector_fields(graph), points_colors_vector_fields(graph)],
                ["M1", "M2", "M3", "M4"],
            )
            @test graph.data.points.entities.names == ["M1", "M2", "M3", "M4"]
            return nothing
        end
    end

    nested_test("axis_vector") do
        graph = points_graph()

        nested_test("()") do
            fill_axis_vector_data!(
                x_axis_vector_fields(graph),
                daf;
                axis = "metacell",
                query_suffix = ": umap_x",
                title = "X",
            )
            @test graph.data.x.vector == Float32[0.0, 1.0, 2.0, 3.0]
            @test graph.data.points.entities.names == ["M1", "M2", "M3", "M4"]
            @test graph.data.points.entities.hovers[1] == "X: 0.0"

            # Nothing is configured, because how to show a value depends on which property it is.
            @test graph.configuration.x_axis.title === nothing
            return nothing
        end

        nested_test("entries") do
            graph = points_graph()
            fill_axis_vector_data!(
                x_axis_vector_fields(graph),
                daf;
                axis = "metacell",
                query_suffix = ": umap_x",
                entries = ["M1", "M3"],
                title = "X",
            )
            @test graph.data.x.vector == Float32[0.0, 2.0]
            @test graph.data.points.entities.names == ["M1", "M3"]
            return nothing
        end

        # A property of what an entry belongs to is reached by naming the way there.
        nested_test("via") do
            set_vector!(daf, "block", "region", ["north", "south"])

            nested_test("()") do
                @test get_axis_vector(daf; axis = "metacell", query_suffix = ": region", via = ["block"]) ==
                      ["north", "north", "south", "south"]
                return nothing
            end

            # An entry which belongs to nothing breaks the chain, and is given the `empty_value` instead.
            nested_test("empty_value") do
                add_axis!(daf, "batch", ["T1", "T2"])
                set_vector!(daf, "batch", "region", ["north", "south"])
                set_vector!(daf, "metacell", "batch", ["T1", "", "T2", ""])
                @test get_axis_vector(
                    daf;
                    axis = "metacell",
                    query_suffix = ": region",
                    via = ["batch"],
                    empty_value = "none",
                ) == ["north", "none", "south", "none"]
                return nothing
            end
        end

        # A value per entry goes to the values of a role and to a hover line, and only the values half is a vector.
        nested_test("halves") do
            put_vector_data!(graph.data.x, Float32[1.0, 2.0, 3.0, 4.0]; title = "ignored")
            @test graph.data.x.vector == Float32[1.0, 2.0, 3.0, 4.0]
            @test graph.data.points.entities.hovers === nothing

            put_vector_data!(points_entities(graph), Float32[1.0, 2.0, 3.0, 4.0]; title = "X")
            @test graph.data.points.entities.hovers[1] == "X: 1.0"
            return nothing
        end
    end

    nested_test("mask") do
        graph = points_graph()
        put_vector_mask_data!(x_axis_vector_fields(graph), [true, false, true, false])
        @test graph.data.points.entities.mask == [true, false, true, false]
        return nothing
    end

    nested_test("gene_expression") do
        nested_test("vector") do
            graph = points_graph()
            fill_gene_expression!(x_axis_vector_fields(graph), daf; gene = "A")
            @test graph.data.x.vector == Float32[0.1, 0.2, 0.3, 0.4]
            @test graph.configuration.x_axis.title == "A fraction"
            @test graph.configuration.x_axis.scale.log_base == Log2Base
            @test graph.configuration.x_axis.scale.log_regularization == GENE_FRACTION_REGULARIZATION_FOR_GRAPHS
            @test graph.data.points.entities.hovers[1] == "A fraction: 0.1"
            return nothing
        end

        nested_test("colors") do
            graph = points_graph()
            fill_gene_expression!(points_colors_vector_fields(graph), daf; gene = "C")
            @test graph.data.points.colors.vector == Float32[0.0, 0.1, 0.0, 0.1]
            @test graph.configuration.points.colors.title == "C fraction"
            @test graph.configuration.points.colors.scale.log_base == Log2Base
            @test graph.configuration.points.colors.show_legend
            return nothing
        end

        nested_test("!regularization") do
            graph = points_graph()
            @test_throws AssertionError fill_gene_expression!(
                x_axis_vector_fields(graph),
                daf;
                gene = "A",
                gene_fraction_regularization = -1,
            )
            return nothing
        end

        nested_test("matrix") do
            nested_test("genes") do
                graph = heatmap_graph()
                fill_genes_expression_matrix!(entries_matrix_fields(graph), daf; genes = [1, 2])
                @test graph.data.entries.matrix == Float32[0.1 0.2 0.3 0.4; 0.4 0.3 0.2 0.1]
                @test graph.data.rows.entities.names == ["A", "B"]
                @test graph.data.columns.entities.names == ["M1", "M2", "M3", "M4"]
                @test graph.data.cells.hovers[1, 1] == "fraction: 0.1"
                @test graph.configuration.entries.colors.scale.log_base == Log2Base
                return nothing
            end

            nested_test("entries") do
                graph = heatmap_graph()
                fill_genes_expression_matrix!(entries_matrix_fields(graph), daf; entries = ["M1", "M3"])
                @test graph.data.entries.matrix == Float32[0.1 0.3; 0.4 0.2; 0.0 0.0]
                @test graph.data.rows.entities.names == ["A", "B", "C"]
                @test graph.data.columns.entities.names == ["M1", "M3"]
                return nothing
            end

            # Asking for everything takes no slice at all, so the matrix is used as the repository holds it.
            nested_test("all") do
                graph = heatmap_graph()
                fill_genes_expression_matrix!(entries_matrix_fields(graph), daf)
                @test graph.data.entries.matrix == Float32[0.1 0.2 0.3 0.4; 0.4 0.3 0.2 0.1; 0.0 0.1 0.0 0.1]
                @test graph.data.rows.entities.names == ["A", "B", "C"]
                return nothing
            end
        end
    end

    nested_test("axes_matrix") do
        graph = heatmap_graph()
        fill_axes_matrix_data!(
            entries_matrix_fields(graph),
            daf;
            rows_axis = "gene",
            columns_axis = "metacell",
            query_suffix = ":: linear_fraction",
            title = "fraction",
        )
        @test graph.data.entries.matrix[1, :] == Float32[0.1, 0.2, 0.3, 0.4]
        @test graph.data.rows.entities.names == ["A", "B", "C"]
        @test graph.data.columns.entities.names == ["M1", "M2", "M3", "M4"]
        @test graph.data.cells.hovers[1, 1] == "fraction: 0.1"

        # Nothing is configured, because how to show a value depends on which property it is.
        @test graph.configuration.entries.colors.palette === nothing
        return nothing
    end

    nested_test("genes_fold") do
        graph = heatmap_graph()

        nested_test("()") do
            fill_genes_fold_matrix!(entries_matrix_fields(graph), daf; genes = [1, 2])
            @test graph.data.entries.matrix[1, :] == Float32[-1.0, 0.0, 0.0, 1.0]
            @test graph.data.entries.matrix[2, :] == Float32[1.0, 0.0, 0.0, -1.0]
            @test graph.data.rows.entities.names == ["A", "B"]
            @test graph.configuration.entries.colors.palette == "BuWtRd"
            @test graph.configuration.entries.colors.show_legend
            @test graph.configuration.entries.colors.scale.minimum == -MAX_FOLD_FOR_GRAPHS
            @test graph.configuration.entries.colors.scale.maximum == MAX_FOLD_FOR_GRAPHS
            return nothing
        end

        nested_test("max_fold") do
            fill_genes_fold_matrix!(entries_matrix_fields(graph), daf; genes = [1, 2], max_fold = 2)
            @test graph.configuration.entries.colors.scale.minimum == -2
            return nothing
        end

        nested_test("!max_fold") do
            @test_throws AssertionError fill_genes_fold_matrix!(
                entries_matrix_fields(graph),
                daf;
                genes = [1, 2],
                max_fold = 0,
            )
            return nothing
        end
    end

    nested_test("umap") do
        graph = points_graph()
        fill_umap!(x_axis_vector_fields(graph), daf; coordinate = "x")
        @test graph.data.x.vector == Float32[0.0, 1.0, 2.0, 3.0]
        @test graph.configuration.x_axis.title == "UMAP X"
        @test !graph.configuration.x_axis.show_ticks
        @test !graph.configuration.x_axis.show_grid
        @test graph.data.points.entities.names == ["M1", "M2", "M3", "M4"]

        # A UMAP coordinate is the arbitrary output of the projection, so it says nothing in a hover.
        @test graph.data.points.entities.hovers === nothing

        fill_umap!(y_axis_vector_fields(graph), daf; coordinate = "y")
        @test graph.data.y.vector == Float32[3.0, 2.0, 1.0, 0.0]
        return nothing
    end

    nested_test("type") do
        add_axis!(daf, "type", ["X", "Y"])
        set_vector!(daf, "metacell", "type", ["X", "X", "Y", "Y"])
        set_vector!(daf, "type", "color", ["red", "blue"])
        set_vector!(daf, "block", "type", ["X", "Y"])

        nested_test("()") do
            graph = points_graph()
            fill_type!(points_colors_vector_fields(graph), daf)
            @test graph.data.points.colors.vector == ["X", "X", "Y", "Y"]
            @test graph.configuration.points.colors.title == "type"
            @test graph.configuration.points.colors.show_legend
            @test graph.data.points.entities.hovers == ["type: X", "type: X", "type: Y", "type: Y"]
            return nothing
        end

        nested_test("!show_legend") do
            graph = points_graph()
            fill_type!(points_colors_vector_fields(graph), daf; show_legend = false)
            @test !graph.configuration.points.colors.show_legend
            return nothing
        end

        # An entry with no type is given a color of its own, so that it is still drawn.
        nested_test("colors") do
            @test get_type_colors(daf) == Dict("X" => "red", "Y" => "blue", "" => EMPTY_TYPE_COLOR)
            @test get_type_colors(daf; empty_type_color = nothing) == Dict("X" => "red", "Y" => "blue")
            return nothing
        end

        # The type of the block of the metacell is reached by naming the way there.
        nested_test("via") do
            @test get_type_vector(daf; axis = "metacell", via = ["block"]) == ["X", "X", "Y", "Y"]
            return nothing
        end

        # The order is a property of the type, so all the entries of a type get the same number. The numbers are used
        # as they are, since the groups are laid out in the order of their numbers.
        nested_test("global_flow_order") do
            set_vector!(daf, "type", "global_flow_order", UInt32[9, 5])
            graph = heatmap_graph()
            fill_global_flow_order!(columns_groups_vector_data_fields(graph), daf)
            @test graph.data.columns.groups.vector == UInt32[9, 9, 5, 5]
            @test graph.data.columns.entities.names == ["M1", "M2", "M3", "M4"]
            return nothing
        end
    end

    nested_test("block") do
        graph = heatmap_graph()
        fill_block!(columns_groups_vector_data_fields(graph), daf)
        @test graph.data.columns.groups.vector == ["B1", "B1", "B2", "B2"]
        @test graph.data.columns.entities.hovers == ["block: B1", "block: B1", "block: B2", "block: B2"]
        return nothing
    end

    nested_test("boolean_annotation") do
        set_vector!(daf, "gene", "is_skeleton", [true, false, true])
        graph = heatmap_graph()
        index = add_rows_annotation!(graph)
        fill_boolean_annotation!(
            rows_annotations_colors_vector_fields(graph, index),
            daf;
            axis = "gene",
            property = "is_skeleton",
        )
        @test graph.data.rows.annotations[1].values.vector == ["true", "false", "true"]
        @test graph.data.rows.annotations[1].colors.title == "is_skeleton"
        @test graph.data.rows.annotations[1].colors.palette == Dict("true" => "black", "false" => "lightgrey")
        @test !graph.data.rows.annotations[1].colors.show_legend
        return nothing
    end

    nested_test("counts") do
        set_vector!(daf, "metacell", "total_UMIs", UInt32[100, 200, 300, 400])
        set_vector!(daf, "metacell", "n_cells", UInt32[10, 20, 30, 40])
        set_vector!(daf, "block", "total_UMIs", UInt32[300, 700])
        set_vector!(daf, "block", "n_cells", UInt32[30, 70])
        set_vector!(daf, "block", "n_metacells", UInt32[2, 2])

        nested_test("total_UMIs") do
            graph = points_graph()
            fill_total_UMIs!(points_colors_vector_fields(graph), daf)
            @test graph.data.points.colors.vector == UInt32[100, 200, 300, 400]
            @test graph.configuration.points.colors.title == "total UMIs"
            @test graph.configuration.points.colors.palette == "YlOrRd"
            @test graph.configuration.points.colors.show_legend
            @test graph.configuration.points.colors.scale.log_base == Log2Base
            @test graph.configuration.points.colors.scale.log_regularization == 0
            return nothing
        end

        nested_test("n_cells") do
            graph = points_graph()
            fill_n_cells!(points_colors_vector_fields(graph), daf)
            @test graph.data.points.colors.vector == UInt32[10, 20, 30, 40]
            return nothing
        end

        nested_test("n_metacells") do
            graph = points_graph()
            fill_n_metacells!(points_colors_vector_fields(graph), daf; axis = "block")
            @test graph.data.points.colors.vector == UInt32[2, 2]
            return nothing
        end

        nested_test("means") do
            nested_test("blocks") do
                graph = points_graph()
                fill_mean_cells_per_metacell!(x_axis_vector_fields(graph), daf; axis = "block")
                fill_mean_total_UMIs_per_metacell!(y_axis_vector_fields(graph), daf; axis = "block")
                @test graph.data.x.vector == Float32[15.0, 35.0]
                @test graph.data.y.vector == Float32[150.0, 350.0]
                @test graph.configuration.x_axis.title == "mean cells per metacell"
                @test graph.configuration.y_axis.title == "mean UMIs per metacell"
                return nothing
            end

            nested_test("metacells") do
                graph = points_graph()
                fill_mean_total_UMIs_per_cell!(points_colors_vector_fields(graph), daf)
                @test graph.data.points.colors.vector == Float32[10.0, 10.0, 10.0, 10.0]
                @test graph.configuration.points.colors.title == "mean UMIs per cell"
                return nothing
            end
        end
    end

    nested_test("gene_selection") do
        set_vector!(daf, "gene", "marker_rank", UInt32[1, 2, typemax(UInt32)])
        set_vector!(daf, "gene", "is_skeleton", [true, false, true])
        @test get_top_marker_gene_indices(daf; markers_count = 100) == [1, 2]
        @test get_top_marker_gene_indices(daf; markers_count = 1) == [1]
        @test get_skeleton_gene_indices(daf) == [1, 3]
        return nothing
    end

    nested_test("gene_correlation") do
        # The base is the repository the blocks come from, and the other repository holds the same genes correlated
        # over the same base blocks. Gene `A` improves in `B1` and degrades in `B2`.
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
        add_axis!(base_daf, "base_block", ["B1", "B2"])
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
            graph = points_graph()
            fill_gene_correlation!(points_entities(graph), base_daf; gene = "A", title = "base correlation")
            fill_gene_correlation!(points_entities(graph), daf; gene = "A")
            @test graph.data.points.entities.names == ["B1", "B2"]
            @test graph.data.points.entities.hovers[1] == "base correlation: 0.3<br>correlation: 0.4"
            return nothing
        end

        nested_test("change") do
            nested_test("axis") do
                graph = points_graph()
                fill_gene_correlation_change!(x_axis_vector_fields(graph), daf, base_daf; gene = "A")
                @test graph.data.x.vector ≈ Float32[0.1, -0.1]
                @test graph.configuration.x_axis.title == "correlation change"
                return nothing
            end

            nested_test("colors") do
                graph = points_graph()
                fill_gene_correlation_change!(points_colors_vector_fields(graph), daf, base_daf; gene = "A")
                @test graph.data.points.colors.vector ≈ Float32[0.1, -0.1]
                @test graph.configuration.points.colors.title == "correlation change"
                @test graph.configuration.points.colors.show_legend
                return nothing
            end
        end
    end

    nested_test("frame") do
        frame = DataFrame(; gene = ["A", "B"], imp_f = [0.75, 0.25], lat = [false, true])

        nested_test("vector") do
            graph = bars_graph()
            fill_column_vector_data!(values_axis_vector_fields(graph), frame; column = "imp_f", title = "improved")
            @test graph.data.values.vector == [0.75, 0.25]
            @test graph.data.bars.hovers == ["improved: 0.75", "improved: 0.25"]

            # A frame has no one column which identifies its rows, so say which one does.
            @test graph.data.bars.names === nothing
            fill_column_names_data!(bars_entities(graph), frame; column = "gene")
            @test graph.data.bars.names == ["A", "B"]
            return nothing
        end

        nested_test("annotation") do
            graph = bars_graph()
            index = add_annotation!(graph)
            fill_column_boolean_annotation!(
                annotations_colors_vector_fields(graph, index),
                frame;
                column = "lat",
                title = "is lateral",
            )
            @test graph.data.annotations[1].values.vector == ["false", "true"]
            @test graph.data.annotations[1].colors.title == "is lateral"
            return nothing
        end
    end

    nested_test("module_regulators") do
        frame = DataFrame(;
            imp_no_mod_f = [0.25, 1.0],
            imp_reg1 = ["B", ""],
            imp_reg1_f = [0.4, 0.0],
            imp_reg2 = ["", ""],
            imp_reg2_f = [0.0, 0.0],
        )

        # A side which gives a gene fewer regulators than asked for pads with empty names, which are left out.
        nested_test("()") do
            hovers = get_module_regulators_hovers(frame; prefix = "imp", side_name = "improved", regulators_count = 2)
            @test hovers[1] == "improved: in no module in 25.0% of the cells<br>- B: 40.0%"
            @test hovers[2] == "improved: in no module in 100.0% of the cells"
            return nothing
        end

        # The lines carry their own labels, so nothing is prefixed to them.
        nested_test("fill") do
            entities = VectorEntitiesData()
            fill_module_regulators_hovers!(
                entities,
                frame;
                prefix = "imp",
                side_name = "improved",
                regulators_count = 2,
            )
            @test entities.hovers[1] == "improved: in no module in 25.0% of the cells<br>- B: 40.0%"
            return nothing
        end
    end

    # A put has a method per leaf it writes and a no-op for the leaves of its kind it ignores. A view whose
    # configuration it ignores is left alone, a leaf of the other kind is walked and nothing is found, and a leaf it
    # says nothing about (the other shape) is an error.
    nested_test("leaves") do
        points = points_graph()
        heatmap = heatmap_graph()
        heatmap.data.entries.matrix = [1.0 2.0; 3.0 4.0]

        nested_test("ignored") do
            sizes = points_sizes_vector_fields(points)
            colors = points_colors_vector_fields(points)
            axis = x_axis_vector_fields(points)
            before = string(points.configuration)

            put_boolean_annotation_configuration!(axis)
            put_count_configuration!(sizes)
            put_gene_correlation_change_configuration!(sizes)
            put_genes_expression_configuration!(sizes)
            put_genes_fold_configuration!(axis)
            put_type_configuration!(axis, Dict("X" => "red"))
            put_umap_configuration!(colors)

            @test string(points.configuration) == before
            return nothing
        end

        nested_test("other") do
            put_count_configuration!(points.data.x)
            put_vector_data!(points.configuration.x_axis, [1.0, 2.0])
            put_matrix_names_data!(heatmap.configuration.entries.colors, ["a", "b"], ["c", "d"])
            @test points.data.x.vector === nothing
            @test heatmap.data.rows.entities.names === nothing
            return nothing
        end

        nested_test("mismatched") do
            @test_throws MethodError put_vector_data!(entries_matrix_fields(heatmap), [1.0, 2.0])
            @test_throws MethodError put_vector_names_data!(entries_matrix_fields(heatmap), ["a", "b"])
            @test_throws MethodError put_matrix_data!(x_axis_vector_fields(points), [1.0 2.0; 3.0 4.0])
            @test_throws "can't name the rows and columns of a vector sink" put_matrix_names_data!(
                x_axis_vector_fields(points),
                ["a", "b"],
                ["c", "d"],
            )
            return nothing
        end

        nested_test("matrix_names") do
            put_matrix_names_data!(entries_matrix_fields(heatmap).data, ["r1", "r2"], ["c1", "c2"])
            @test heatmap.data.rows.entities.names == ["r1", "r2"]
            @test heatmap.data.columns.entities.names == ["c1", "c2"]

            other = heatmap_graph()
            put_matrix_names_data!(
                (entries_matrix_fields(other), heatmap.configuration.entries.colors),
                ["r3", "r4"],
                ["c3", "c4"],
            )
            @test other.data.rows.entities.names == ["r3", "r4"]
            @test other.data.columns.entities.names == ["c3", "c4"]
            return nothing
        end
    end
end
