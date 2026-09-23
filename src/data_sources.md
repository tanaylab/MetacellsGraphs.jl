# Data Sources

```@docs
MetacellsGraphs.DataSources
```

## Generic

These fill the data half of any role from a `daf` query. A specific data source uses them under the hood and adds only
how to show the result.

### Axis Names

```@docs
MetacellsGraphs.DataSources.fill_axis_names_data!
MetacellsGraphs.DataSources.put_vector_names_data!
MetacellsGraphs.DataSources.get_axis_entries_vector
```

### Axis Vectors

```@docs
MetacellsGraphs.DataSources.fill_axis_vector_data!
MetacellsGraphs.DataSources.put_vector_data!
MetacellsGraphs.DataSources.get_axis_vector
```

### Masks

```@docs
MetacellsGraphs.DataSources.put_vector_mask_data!
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
MetacellsGraphs.DataSources.put_matrix_names_data!
```

### Axes Matrices

```@docs
MetacellsGraphs.DataSources.fill_axes_matrix_data!
MetacellsGraphs.DataSources.put_matrix_data!
MetacellsGraphs.DataSources.get_axes_matrix
```

## Gene Expression

```@docs
MetacellsGraphs.DataSources.put_genes_expression_configuration!
MetacellsGraphs.DataSources.GENE_FRACTION_REGULARIZATION_FOR_GRAPHS
```

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
MetacellsGraphs.DataSources.put_genes_fold_configuration!
MetacellsGraphs.DataSources.get_genes_fold_matrix
MetacellsGraphs.DataSources.MAX_FOLD_FOR_GRAPHS
```

## Gene Correlations

```@docs
MetacellsGraphs.DataSources.fill_gene_correlation!
MetacellsGraphs.DataSources.get_gene_correlation_vector
```

### Gene Correlation Changes

```@docs
MetacellsGraphs.DataSources.fill_gene_correlation_change!
MetacellsGraphs.DataSources.put_gene_correlation_change_configuration!
MetacellsGraphs.DataSources.get_gene_correlation_change_vector
```

## UMAP Coordinates

```@docs
MetacellsGraphs.DataSources.fill_umap!
MetacellsGraphs.DataSources.put_umap_data!
MetacellsGraphs.DataSources.put_umap_configuration!
MetacellsGraphs.DataSources.get_umap_vector
```

## Gene Selection

These pick which genes a graph shows, and return indices to pass as the `genes` of another data source.

```@docs
MetacellsGraphs.DataSources.get_top_marker_gene_indices
MetacellsGraphs.DataSources.get_skeleton_gene_indices
```

## Counts

Every count-like quantity is shown the same way, so they share one configuration.

```@docs
MetacellsGraphs.DataSources.put_count_configuration!
```

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
```

## Global Flow Order

```@docs
MetacellsGraphs.DataSources.fill_global_flow_order!
MetacellsGraphs.DataSources.get_global_flow_order_vector
```

## Types

```@docs
MetacellsGraphs.DataSources.fill_type!
MetacellsGraphs.DataSources.put_type_configuration!
MetacellsGraphs.DataSources.get_type_vector
MetacellsGraphs.DataSources.get_type_colors
MetacellsGraphs.DataSources.EMPTY_TYPE_COLOR
```

## Boolean Annotations

```@docs
MetacellsGraphs.DataSources.fill_boolean_annotation!
MetacellsGraphs.DataSources.put_boolean_annotation_configuration!
MetacellsGraphs.DataSources.get_boolean_annotation_vector
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
