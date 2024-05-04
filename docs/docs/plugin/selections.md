# Selections


Assuming you have a model already loaded into the plugin (see [Creating Plugin Instance](./instance.md)), these are some of the following method you can select structural data.

### Selecting directly from the `hierarchy` manager

One can select a subcomponent's data directly from the plugin manager.

```typescript 
import { Structure } from '../mol-model/structure';

const ligandData = plugin.managers.structure.hierarchy.selection.structures[0]?.components[0]?.cell.obj?.data;
const ligandLoci = Structure.toStructureElementLoci(ligandData as any);

plugin.managers.camera.focusLoci(ligandLoci);
plugin.managers.interactivity.lociSelects.select({ loci: ligandLoci });
```

## Selection callbacks
If you want to subscribe to selection events (e.g. to change external state in your application based on a user selection), you can use: `plugin.behaviors.interaction.click.subscribe`

Here's an example of passing in a React "set" function to update selected residue positions.
```typescript
import {
  Structure,
  StructureProperties,
} from "molstar/lib/mol-model/structure"
// setSelected is assumed to be a "set" function returned by useState
// (selected: any[]) => void
plugin.behaviors.interaction.click.subscribe(
  (event: InteractivityManager.ClickEvent) => {
    const selections = Array.from(
      plugin.managers.structure.selection.entries.values()
    );
    // This bit can be customized to record any piece information you want
    const localSelected: any[] = [];
    for (const { structure } of selections) {
      if (!structure) continue;
      Structure.eachAtomicHierarchyElement(structure, {
        residue: (loc) => {
          const position = StructureProperties.residue.label_seq_id(loc);
          localSelected.push({ position });
        },
      });
    }
    setSelected(localSelected);
  }
)
```

### `Molscript` language

Molscript is a language for addressing crystallographic structures and is a part of the Mol* library found at `https://github.com/molstar/molstar/tree/master/src/mol-script`. It can be used against the Molstar plugin as a query language and transpiled against multiple external molecular visualization libraries(see [here](https://github.com/molstar/molstar/tree/master/src/mol-script/transpilers)).

### Querying a structure for a specific chain and residue range (select residues with 12<res_id<200 of chain with auth_asym_id==A) :

```typescript
import { compileIdListSelection } from 'molstar/lib/mol-script/util/id-list'

const query = compileIdListSelection('A 12-200', 'auth');
window.molstar?.managers.structure.selection.fromCompiledQuery('add',query);
```

## Selection Queries

Another way to create a selection is via a `SelectionQuery` object. This is a more programmatic way to create a selection. The following example shows how to select a chain and a residue range using a `SelectionQuery` object.
This relies on the concept of `Expression` which is basically a intermediate representation between a Molscript statement and a selection query. 
 
### Select residues 10-15 of chains A and F in a structure using a `SelectionQuery` object:

```typescript

import { MolScriptBuilder as MS, MolScriptBuilder } from 'molstar/lib/mol-script/language/builder';
import { Expression } from 'molstar/lib/mol-script/language/expression';
import {  StructureSelectionQuery } from 'molstar/lib/mol-plugin-state/helpers/structure-selection-query'


export function select_multiple() {

 const args = [['A', 10, 15], ['F', 10, 15]]
 const groups: Expression[] = [];
 for (var chain of args) {
   groups.push(MS.struct.generator.atomGroups({
     "chain-test": MS.core.rel.eq([MolScriptBuilder.struct.atomProperty.macromolecular.auth_asym_id(), chain[0]]),
     "residue-test": MS.core.rel.inRange([MolScriptBuilder.struct.atomProperty.macromolecular.label_seq_id(), chain[1], chain[2]])
   }));
 }
 var sq = StructureSelectionQuery('residue_range_10_15_in_A_and_F', MS.struct.combinator.merge(groups))
 mstar.managers.structure.selection.fromSelectionQuery('set', sq)
}
```

Complex queries can be constructed by combining primitive queries at the level of [`chain-test`, `residue-test`, `entity-test`, etc] (https://github.com/molstar/molstar/blob/6edbae80db340134341631f669eec86543a0f1a8/src/mol-script/language/symbol-table/structure-query.ts#L88C4-L94C112) by combining them via logical connectives provided in the `MolscriptBuilder.core.rel` as above.

Inspect these examples to get a better feeling for this syntax: `https://github.com/molstar/molstar/blob/6edbae80db340134341631f669eec86543a0f1a8/src/mol-plugin-state/helpers/structure-selection-query.ts#L88-L580`


Furthermore, a query made this way can be converted to a `Loci` object which is important in many parts of the libary:
```typescript

// Select residue 124 of chain A and convert to Loci
const Q = MolScriptBuilder;
var sel = Script.getStructureSelection(Q => Q.struct.generator.atomGroups({
                'chain-test'  : Q.core.rel.eq([Q.struct.atomProperty.macromolecular.auth_asym_id(), A]),
                "residue-test": Q.core.rel.eq([Q.struct.atomProperty.macromolecular.label_seq_id(), 124]),
              }), objdata)

let loci = StructureSelection.toLociWithSourceUnits(sel);
```
