# Exporting components

Export components data can be useful to reproduce the same view in a different visualization software.
To do that, one would need to loop over all components, extract its selection (for example by using atom indices) and its representations (type, coloring and sizing).

### Getting assets / molecular files

```js
for (const { asset, file } of plugin.managers.asset.assets) {
  const isFile = asset.asset.kind === 'url'
  console.log(asset.asset.id)
  console.log(isFile)
  const data = await file.arrayBuffer()
}
```

### Getting components per structure

```js
import { PluginStateObject as PSO } from 'molstar/lib/mol-plugin-state/objects';
//...

const componentManager = plugin.managers.structure.component;
for (const structure of componentManager.currentStructures) {
  if (!structure.properties) {
      continue;
  }
  const cell = plugin.state.data.select(structure.properties.cell.transform.ref)[0];
  if (!cell || !cell.obj) {
    continue;
  }
  const structureData = (cell.obj as PSO.Molecule.Structure).data;
  for (const component of structure.components) {
    if (!component.cell.obj) {
      continue;
    }
    // For each component in each structure, display the content of the selection
    Structure.eachAtomicHierarchyElement(component.cell.obj.data, {
      atom: location => console.log(location.element)
    });
    for (const rep of component.representations) {
      // For each representation of the component, display its type
      console.log(rep.cell?.transform?.params?.type?.name)

      // Also display the color for each atom
      const colorThemeName = rep.cell.transform.params?.colorTheme.name;
      const colorThemeParams = rep.cell.transform.params?.colorTheme.params;
      const theme = plugin.representation.structure.themes.colorThemeRegistry.create(
        colorThemeName || '',
        { structure: structureData },
        colorThemeParams
      ) as ColorTheme<typeof colorThemeParams>;
      Structure.eachAtomicHierarchyElement(component.cell.obj.data, {
        atom: loc => console.log(theme.color(loc, false))
      });
    }
  }
}
```
