# Markdown Extension Manager

The `markdownExtensions` manager in `PluginContext.manager` allows customizing
the `Markdown` React component to enable executing commands and rendering custom content.

The main use case of this is enriching [MolViewSpec](`https://molstar.org/mol-view-spec`) support.

## API

- `PluginContext.manager.markdownExtensions.register*` functions can be used to register extensions and state/data resolvers to make the the manager work with plugin extension
- `PluginContext.manager.markdownExtensions.remove*` can be used to dynamically remove the above

## Commands

Extends Markdown Hyperlink syntax to support expressions of the form `[title](!c1=v1&c2=v2&...)` into an executable command. The command can be executed either on click, mouse enter, or mouse leave.

Generally, the command should be URL encoded, e.g., `a b` => `a%20b` (in JS, `encodeURIComponent`, in Python `urllib.parse.quote_plus/urlencode`).

### Built-in Commands

- `center-camera` - Centers the camera
- `apply-snapshot=key` - Loads snapshots with the provided key
- `focus-refs=ref1,ref2,...` - On click, focuses nodes with the provided refs
- `highlight-refs=ref1,ref2,...` - On mouse over, highlights the provided refs
- `query=...&lang=...&action=highlight,focus&focus-radius=...`
  - `query` is an expression (e.g., `resn HEM` when using PyMol syntax)
  - (optional) `lang` is one of `mol-script` (default), `pymol`, `vmd`, `jmol`
  - (optional) `action` is an array of `highlight` (default), `focus` (multiple actions can be specified)
  - (optional) `focus-radius` is extra distance applied when focusing the selection (default is `3`)
  - Example: `[HEM](!query%3Dresn%20HEM%26lang%3Dpymol%26action%3Dhighlight%2Cfocus)` highlights or focuses the HEM residue (the command must be URL encoded because it contains spaces and possibly other special characters)
- `play-audio=src`, `toggle-audio[=str]`, `stop-audio`, `pause-audio` - Audio playback support

## Custom Content

Extends Markdown Image syntax to support expressions of the form `![alt](!c1=v1&c2=v2&...)` to render custom elements instead.

### Built-in Custom Content
- `color-swatch=color` - Renders a box with the provided color
-  Color palettes:
  - `color-palette-name=name` - Renders a gradient with the provided named color palette (see `mol-util/color/lists.ts` for supported color schemes)
  - `color-palette-colors=color1,color2` - Renders a gradient with the provided colors
  - `color-palette-width=CCS-value` - Specifies the width of the element, defaults to `150px`
  - `color-palette-height=CCS-value` - Specified the height of the element, defaults to `0.5em`
  - `color-palette-discrete` - Renders discrete color list instead of interpolating


## Example

```markdown
### Highlight/Focus:
- ![blue](!color-swatch=blue) [polymer](!highlight-refs=polymer&focus-refs=polymer)
- ![blue](!color-swatch=red) [ligand](!highlight-refs=ligand&focus-refs=ligand)
- [both](!highlight-refs=polymer,ligand&focus-refs=polymer,ligand)

### Color Palettes
|name|visual|
|---:|---|
|viridis|![viridis](!color-palette-name=viridis)|
|rainbow (discrete)|![simple-rainbow](!color-palette-name=simple-rainbow&color-palette-discrete)|
|custom|![custom](!color-palette-colors=red,#00ff00,rgb(0,0,255))|

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