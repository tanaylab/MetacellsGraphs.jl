# Data Sources

```@docs
MetacellsGraphs.DataSources
```

## Generic

These fill the data half of any role from a `daf` query. A specific data source uses them under the hood and adds only
how to show the result. They hand the data to the generic puts of
[`SomeGraphs.Sources`](@extref SomeGraphs SomeGraphs.Sources).

### Axis Names

```@docs
MetacellsGraphs.DataSources.fill_axis_names_data!
MetacellsGraphs.DataSources.get_axis_entries_vector
```

### Axis Vectors

```@docs
MetacellsGraphs.DataSources.fill_axis_vector_data!
MetacellsGraphs.DataSources.get_axis_vector
```

### Frame Columns

Some graphs show a report computed by `Metacells`, which is a `DataFrame` rather than a property of a `daf` repository.
These fetch from a frame instead, and hand the result to the same `put_` functions.

```@docs
MetacellsGraphs.DataSources.fill_column_vector_data!
MetacellsGraphs.DataSources.fill_column_names_data!
MetacellsGraphs.DataSources.fill_column_boolean_annotation!
MetacellsGraphs.DataSources.get_column_vector
```

## Module Regulators

```@docs
MetacellsGraphs.DataSources.fill_module_regulators_hovers!
MetacellsGraphs.DataSources.get_module_regulators_hovers
```

### Axes Names

```@docs
MetacellsGraphs.DataSources.fill_axes_names_data!
```

### Axes Matrices

```@docs
MetacellsGraphs.DataSources.fill_axes_matrix_data!
MetacellsGraphs.DataSources.get_axes_matrix
```

## Display Buckets

Properties fall into a few buckets by how they are shown. Each bucket has one configuration function, and the data
source of each property uses the function of its bucket.

Most buckets also have generic data sources, taking the name of a property: a get, which fetches it, a fill, which
fetches and shows it, and a put, which shows values already fetched. The data sources of the stored properties of each
bucket (e.g. [`fill_n_cells!`](@ref)) are generated from a table. Each fetches with its own get, and shows the values
with the put of its bucket.

```@docs
MetacellsGraphs.DataSources.put_size_configuration!
MetacellsGraphs.DataSources.get_size_vector
MetacellsGraphs.DataSources.fill_size!
MetacellsGraphs.DataSources.put_size_values!
MetacellsGraphs.DataSources.put_expression_level_configuration!
MetacellsGraphs.DataSources.GENE_FRACTION_REGULARIZATION_FOR_GRAPHS
MetacellsGraphs.DataSources.put_fold_configuration!
MetacellsGraphs.DataSources.MAX_FOLD_FOR_GRAPHS
MetacellsGraphs.DataSources.put_correlation_configuration!
MetacellsGraphs.DataSources.get_correlation_vector
MetacellsGraphs.DataSources.fill_correlation!
MetacellsGraphs.DataSources.put_correlation_values!
MetacellsGraphs.DataSources.put_boolean_configuration!
MetacellsGraphs.DataSources.get_boolean_vector
MetacellsGraphs.DataSources.fill_boolean!
MetacellsGraphs.DataSources.put_boolean_values!
MetacellsGraphs.DataSources.put_category_configuration!
MetacellsGraphs.DataSources.get_category_vector
MetacellsGraphs.DataSources.fill_category!
MetacellsGraphs.DataSources.put_category_values!
MetacellsGraphs.DataSources.get_category_colors
MetacellsGraphs.DataSources.EMPTY_CATEGORY_COLOR
MetacellsGraphs.DataSources.put_coordinates_configuration!
```

## Gene Expression

### Gene Expression Vectors

```@docs
MetacellsGraphs.DataSources.fill_gene_expression!
MetacellsGraphs.DataSources.get_gene_expression_vector
```

### Gene Expression Matrices

```@docs
MetacellsGraphs.DataSources.fill_genes_expression_matrix!
MetacellsGraphs.DataSources.get_genes_expression_matrix
```

### Gene Fold Matrices

```@docs
MetacellsGraphs.DataSources.fill_genes_fold_matrix!
MetacellsGraphs.DataSources.get_genes_fold_matrix
```

## Gene Correlations

```@docs
MetacellsGraphs.DataSources.fill_gene_correlation!
MetacellsGraphs.DataSources.get_gene_correlation_vector
MetacellsGraphs.DataSources.get_gene_is_correlated_vector
MetacellsGraphs.DataSources.fill_correlation_between_cells_and_punctuated_metacells!
MetacellsGraphs.DataSources.get_correlation_between_cells_and_punctuated_metacells_vector
MetacellsGraphs.DataSources.fill_correlation_between_cells_and_projected_punctuated_metacells!
MetacellsGraphs.DataSources.get_correlation_between_cells_and_projected_punctuated_metacells_vector
MetacellsGraphs.DataSources.fill_correlation_between_cells_and_projected_metacells!
MetacellsGraphs.DataSources.get_correlation_between_cells_and_projected_metacells_vector
```

### Gene Correlation Changes

```@docs
MetacellsGraphs.DataSources.fill_gene_correlation_change!
MetacellsGraphs.DataSources.get_gene_correlation_change_vector
```

## UMAP Coordinates

```@docs
MetacellsGraphs.DataSources.fill_umap!
MetacellsGraphs.DataSources.put_umap_data!
MetacellsGraphs.DataSources.get_umap_vector
```

## Gene Selection

These pick which genes a graph shows, and return indices to pass as the `genes` of another data source.

```@docs
MetacellsGraphs.DataSources.get_top_marker_gene_indices
MetacellsGraphs.DataSources.get_skeleton_gene_indices
```

## Entry Masks

These pick which entries of an axis a graph shows, and return a Boolean mask over the whole axis, to pass as the
`entries` of another data source. Masks combine by `.&`.

```@docs
MetacellsGraphs.DataSources.get_is_in_categories_vector
MetacellsGraphs.DataSources.get_is_in_neighborhood_vector
MetacellsGraphs.DataSources.get_is_in_environment_vector
```

## Counts

### Total UMIs

```@docs
MetacellsGraphs.DataSources.fill_total_UMIs!
MetacellsGraphs.DataSources.get_total_UMIs_vector
```

### Cells and Metacells

```@docs
MetacellsGraphs.DataSources.fill_n_cells!
MetacellsGraphs.DataSources.get_n_cells_vector
MetacellsGraphs.DataSources.fill_n_metacells!
MetacellsGraphs.DataSources.get_n_metacells_vector
```

### Means

```@docs
MetacellsGraphs.DataSources.fill_mean_cells_per_metacell!
MetacellsGraphs.DataSources.get_mean_cells_per_metacell_vector
MetacellsGraphs.DataSources.fill_mean_total_UMIs_per_metacell!
MetacellsGraphs.DataSources.get_mean_total_UMIs_per_metacell_vector
MetacellsGraphs.DataSources.fill_mean_total_UMIs_per_cell!
MetacellsGraphs.DataSources.get_mean_total_UMIs_per_cell_vector
```

## Blocks

```@docs
MetacellsGraphs.DataSources.fill_block!
MetacellsGraphs.DataSources.get_block_vector
MetacellsGraphs.DataSources.get_block_colors
```

## Global Flow Order

```@docs
MetacellsGraphs.DataSources.fill_global_flow_order!
MetacellsGraphs.DataSources.get_global_flow_order_vector
```

## Types

```@docs
MetacellsGraphs.DataSources.fill_type!
MetacellsGraphs.DataSources.get_type_vector
MetacellsGraphs.DataSources.get_type_colors
```

## Gene Flags

```@docs
MetacellsGraphs.DataSources.fill_is_lateral!
MetacellsGraphs.DataSources.get_is_lateral_vector
MetacellsGraphs.DataSources.fill_is_marker!
MetacellsGraphs.DataSources.get_is_marker_vector
MetacellsGraphs.DataSources.fill_is_transcription_factor!
MetacellsGraphs.DataSources.get_is_transcription_factor_vector
MetacellsGraphs.DataSources.fill_is_regulator!
MetacellsGraphs.DataSources.get_is_regulator_vector
MetacellsGraphs.DataSources.fill_is_forbidden!
MetacellsGraphs.DataSources.get_is_forbidden_vector
MetacellsGraphs.DataSources.fill_is_skeleton!
MetacellsGraphs.DataSources.get_is_skeleton_vector
MetacellsGraphs.DataSources.fill_is_correlated_with_skeleton!
MetacellsGraphs.DataSources.get_is_correlated_with_skeleton_vector
MetacellsGraphs.DataSources.fill_is_mitochondrial!
MetacellsGraphs.DataSources.get_is_mitochondrial_vector
MetacellsGraphs.DataSources.fill_is_ribosomal!
MetacellsGraphs.DataSources.get_is_ribosomal_vector
MetacellsGraphs.DataSources.fill_is_excluded!
MetacellsGraphs.DataSources.get_is_excluded_vector
```

## Low Level Functions

```@docs
MetacellsGraphs.DataSources.get_vector_query
MetacellsGraphs.DataSources.get_matrix_query
```

## Index

```@index
Pages = ["data_sources.md"]
```
