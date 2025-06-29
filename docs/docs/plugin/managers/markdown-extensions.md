# Markdown Extension Manager

The `markdownExtensions` manager in `PluginContext.manager` allows customizing
the `Markdown` React component to enable executing commands and rendering custom content.

The main use case of this is enriching [MolViewSpec](`https://molstar.org/mol-view-spec`) support.

## API

- `PluginContext.manager.markdownExtensions.register*` functions can be used to register extensions and state/data resolvers to make the the manager work with plugin extension
- `PluginContext.manager.markdownExtensions.remove*` can be used to dynamically remove the above

## Commands

Extends Markdown Hyperlink syntax to support expressions of the form `[title](!c1=v1&c2=v2&...)` into an executable command. The command can be executed either on click, mouse enter, or mouse leave.

### Built-in Commands

- `center-camera` - Centers the camera
- `apply-snapshot=key` - Loads snapshots with the provided key
- `focus-refs=ref1,ref2,...` - On click, focuses nodes with the provided refs
- `highlight-refs=ref1,ref2,...` - On mouse over, highlights the provided refs

## Custom Content

Extends Markdown Image syntax to support expressions of the form `![alt](!c1=v1&c2=v2&...)` to render custom elements instead.

### Built-in Custom Content

- `color-swatch=color` - Renders a box with the provided color
- `color-palette=name` - Renders a gradient with the provivided named color palette (see `mol-util/color/lists.ts` for supported color schemes)
  - `color-palette-width=CCS-value` - Specifies the width of the element, defaults to `150px`
  - `color-palette-height=CCS-value` - Specified the height of the element, defaults to `0.5em`
  - `color-palette-discrete` - Renders discrete color list instead of interpolating


## Example

```markdown
### Highlight/Focus:
- ![blue](!color-swatch=blue) [polymer](!highlight-refs=polymer&focus-refs=polymer)
- ![blue](!color-swatch=red) [ligand](!highlight-refs=ligand&focus-refs=ligand)
- [both](!highlight-refs=polymer,ligand&focus-refs=polymer,ligand)

### Table
|name|visual|
|---:|---|
|viridis|![viridis](!color-palette=viridis)|
|rainbow (discrete)|![simple-rainbow](!color-palette=simple-rainbow&color-palette-discrete)|

### Camera controls
- [center](!center-camera)

### Image embedded in MVSX file
![mvsx image](logo.png)
```

This works with the MolViewSpec state built by:

```py
import molviewspec as mvs

builder = mvs.create_builder()

assets = {
    "1cbs.cif": "https://files.wwpdb.org/download/1cbs.cif",
    "logo.png": "https://molstar.org/img/molstar-logo.png",
}

model = (
    builder.download(url="1cbs.cif")
        .parse(format="mmcif")
        .model_structure()
)
(
    model.component(selector="polymer")
    .representation(ref="polymer")
    .color(color="blue")
)
(
    model.component(selector="ligand")
    .representation(ref="ligand")
    .color(color="red")
)

mvsx = mvs.MVSX(
    data=builder.get_state(
        description="""..."""  # inline the code above
    ),
    assets=assets
)
```