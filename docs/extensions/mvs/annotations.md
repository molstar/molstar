# MVS annotations

Annotations are used to define substructures (components) and apply colors, labels, or tooltips to them. In contrast to [selectors](./selectors.md), annotations are defined in a separate file, which can then be referenced in the main MVS file.


## MVS annotation files

MVS annotations can be encoded in multiple different formats, but their logic is always the same and in fact very similar to that of selectors.

### JSON format

The simplest example of an annotation in JSON format is just a JSON-encoded [union component expression](./selectors.md) selector. Here is a simple annotation containing 4 **annotation rows**:

```json
[
    { "label_asym_id": "A" },
    { "label_asym_id": "B" },
    { "label_asym_id": "B", "beg_label_seq_id": 100, "end_label_seq_id": 200 },
    { "label_asym_id": "B", "beg_label_seq_id": 150, "end_label_seq_id": 160 },
]
```

However, in a typical annotation, there is at least one extra field that provides the value of the dependent variable (such as color or label) mapped to each annotation row:

```json
[
    { "label_asym_id": "A", "color": "#00ff00" },
    { "label_asym_id": "B", "color": "blue" },
    { "label_asym_id": "B", "beg_label_seq_id": 100, "end_label_seq_id": 200, "color": "skyblue" }
    { "label_asym_id": "B", "beg_label_seq_id": 150, "end_label_seq_id": 160, "color": "lightblue" }
]
```

This particular annotation (when applied via `color_from_uri` node) will apply green color (#00ff00) to the whole chain A and three shades of blue to the chain B. Later annotation rows override earlier rows, therefore residues 1–99 will be blue, 100–149 skyblue, 150–160 lightblue, 161–200 skyblue, and 201–end blue. (Tip: to color all the rest of the structure in one color, add an annotation row with no selector fields (e.g. `{ "color": "yellow" }`) to the beginning of the annotation.)

Real-life annotation files can include huge numbers of annotation rows. To avoid repeating the same field keys in every row, we can convert the array-of-objects into object-of-arrays. This will result in an equivalent annotation but smaller file size:

```json
{
    "label_asym_id": ["A", "B", "B", "B"],
    "beg_label_seq_id": [null, null, 100, 150],
    "end_label_seq_id": [null, null, 200, 160],
    "color": ["#00ff00", "blue", "skyblue", "lightblue"]
}
```

A more complex example of JSON annotation is provided in [/examples/mvs/1h9t_domains.json](/examples/mvs/1h9t_domains.json).

### CIF format

Annotations can also be encoded using CIF format, a table-based format which is commonly used in structure biology to store structures or any kind of tabular data.

The example from above, encoded as CIF, would look like this:

```cif
data_annotation
loop_
_coloring.label_asym_id
_coloring.beg_label_seq_id
_coloring.end_label_seq_id
_coloring.color
A   .   . '#00ff00'
B   .   . 'blue'
B 100 200 'skyblue'
B 150 160 'lightblue'
```

An advantage of the CIF format is that it can include multiple annotation tables in the same file, organized into blocks and categories. Then the MVS file can reference individual tables using `block_header` (or `block_index`) and `category_name` parameters. The column containing the dependent variable can be specified using `field_name` parameter. In this case, we could use `"block_header": "annotation", "category_name": "coloring", "field_name": "color"`.

### BCIF format

This has exactly the same structure as the CIF format, but encoded using [BinaryCIF](https://github.com/molstar/BinaryCIF).


## Referencing MVS annotations in MVS tree

### From URI

MVS annotations can be referenced in `color_from_uri`, `label_from_uri`, `tooltip_from_uri`, and `component_from_uri` nodes in MVS tree.

For example this part of a MVS tree:

```txt
- representation {type: "cartoon"}
  - color {selector: {label_asym_id: "A"}, color: "#00ff00"}
  - color {selector: {label_asym_id: "B"}, color: "blue"}
  - color {selector: {label_asym_id: "B", beg_label_seq_id: 100, end_label_seq_id: 200}, color: "skyblue"}
  - color {selector: {label_asym_id: "B", beg_label_seq_id: 150, end_label_seq_id: 160}, color: "lightblue"}
```

can be replaced by:

```txt
- representation {type: "cartoon"}
  - color_from_uri {uri: "https://example.org/annotations.json", format: "json", schema: "residue_range"}
```

assuming that the JSON annotation file shown in the previous section is available at `https://example.org/annotations.json`. 

#### Relative URIs

The `uri` parameter can also hold a URI reference (relative URI). In such cases, this URI reference is relative to the URI of the MVS file itself (e.g. if the MVS file is available from `https://example.org/spanish/inquisition/expectations.mvsj`, then the relative URI `./annotations.json` is equivalent to `https://example.org/spanish/inquisition/annotations.json`). This is however not applicable in all cases (e.g. the MVS tree can be constructed ad-hoc within a web application, therefore it has no URI; or the MVS file is loaded from a local disk using drag&drop, therefore the relative location is not accessible by the browser).

### From source

The MVS annotations can in fact be stored within the same mmCIF file from which the structure coordinates are loaded. To reference these annotations, we can use `color_from_source`, `label_from_source`, `tooltip_from_source`, and `component_from_source` nodes. Example:

```txt
- representation {type: "cartoon"}
  - color_from_source {schema: "residue_range", block_header: "annotation", category_name: "coloring"}
```


## Annotation schemas

The `schema` parameter of all `*_from_uri` and `*_from_source` nodes specifies the MVS annotation schema, i.e. a set of fields used to select a substructure. In the example above we are using `residue_range` schema, which uses columns `label_entity_id`, `label_asym_id`, `beg_label_seq_id`, and `end_label_seq_id`. (We didn't provide values for `label_entity_id`, so it is not taken into account even though the schema supports it).


Table of selector field names supported by individual MVS annotation schemas:

|Field \ Schema|whole_structure|entity|chain|residue|residue_range|atom|auth_chain|auth_residue|auth_residue_range|auth_atom|all_atomic|
|:------------------|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| label_entity_id   |   | X | X | X | X | X |   |   |   |   | X |
| label_asym_id     |   |   | X | X | X | X |   |   |   |   | X |
| label_seq_id      |   |   |   | X |   | X |   |   |   |   | X |
| beg_label_seq_id  |   |   |   |   | X |   |   |   |   |   | X |
| end_label_seq_id  |   |   |   |   | X |   |   |   |   |   | X |
| label_atom_id     |   |   |   |   |   | X |   |   |   |   | X |
| auth_asym_id      |   |   |   |   |   |   | X | X | X | X | X |
| auth_seq_id       |   |   |   |   |   |   |   | X |   | X | X |
| pdbx_PDB_ins_code |   |   |   |   |   |   |   | X |   | X | X |
| beg_auth_seq_id   |   |   |   |   |   |   |   |   | X |   | X |
| end_auth_seq_id   |   |   |   |   |   |   |   |   | X |   | X |
| auth_atom_id      |   |   |   |   |   |   |   |   |   | X | X |
| type_symbol       |   |   |   |   |   | X |   |   |   | X | X |
| atom_id           |   |   |   |   |   | X |   |   |   | X | X |
| atom_index        |   |   |   |   |   | X |   |   |   | X | X |

To include all selector field names that are present in the annotation, one can use `"schema": "all_atomic"` (we could use it in the example above and the result would be the same). In future versions of MVS, non-atomic schemas might be added, to select parts of structures that are not composed of atoms, e.g. coarse models or geometric primitives.


## `group_id` field

The `group_id` field is a special field supported by all MVS annotation schemas. It does not change the sets of atoms selected by individual rows but instead groups annotation rows together to create more complex selections. This is useful when adding labels to our visualization.

The following example (when applied via `label_from_uri` node) will create 7 separate labels, each bound to a single residue:

```cif
data_annotation
loop_
_labels.label_asym_id
_labels.label_seq_id
_labels.color
_labels.label
A 100 pink 'Substrate binding site'
A 150 pink 'Substrate binding site'
A 170 pink 'Substrate binding site'
A 200 blue 'Inhibitor binding site'
A 220 blue 'Inhibitor binding site'
A 300 lime 'Glycosylation site'
A 330 lime 'Glycosylation site'
```

On the other hand, the next example will only create 4 labels ("Substrate binding site" label bound to residues 100, 150, and 170; "Inhibitor binding site" label bound to residues 200 and 220; "Glycosylation site" label bound to residue 300; and "Glycosylation site" label bound to residue 330):

```cif
data_annotation
loop_
_labels.group_id
_labels.label_asym_id
_labels.label_seq_id
_labels.color
_labels.label
1 A 100 pink 'Substrate binding site'
1 A 150 pink 'Substrate binding site'
1 A 170 pink 'Substrate binding site'
2 A 200 blue 'Inhibitor binding site'
2 A 220 blue 'Inhibitor binding site'
. A 300 lime 'Glycosylation site'
. A 330 lime 'Glycosylation site'
```

Note: Annotation rows with empty `group_id` field (`.` in CIF, ommitted field or `null` in JSON) are always treated as separate groups.

Note 2: `group_id` field has no effect on colors, tooltips, components. It only makes any difference for labels.
