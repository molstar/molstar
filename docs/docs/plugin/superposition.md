# Structure Superposition

Mol* provides utilities for superposing protein structures, including both sequence-independent (RMSD-based) and structure-based (TM-align) methods.

## RMSD-based Superposition

The basic superposition method uses the Kabsch algorithm to minimize RMSD between corresponding atoms:

```typescript
import { superpose } from 'molstar/lib/mol-model/structure/structure/util/superposition';
import { StructureSelection, QueryContext } from 'molstar/lib/mol-model/structure';
import { compile } from 'molstar/lib/mol-script/runtime/query/compiler';
import { MolScriptBuilder as MS } from 'molstar/lib/mol-script/language/builder';

// Create a query for C-alpha atoms
const caQuery = compile<StructureSelection>(MS.struct.generator.atomGroups({
    'atom-test': MS.core.rel.eq([MS.struct.atomProperty.macromolecular.label_atom_id(), 'CA'])
}));

// Get selections from two structures
const sel1 = StructureSelection.toLociWithCurrentUnits(caQuery(new QueryContext(structure1)));
const sel2 = StructureSelection.toLociWithCurrentUnits(caQuery(new QueryContext(structure2)));

// Compute superposition (returns transformation matrices)
const transforms = superpose([sel1, sel2]);

// transforms[0].bTransform contains the Mat4 to superpose structure2 onto structure1
```

## TM-align Superposition

TM-align is a structure-based alignment algorithm that produces the TM-score, a length-independent metric for comparing protein structures. Unlike RMSD, TM-score is normalized to [0, 1] and is more robust for comparing proteins of different sizes.

### Basic Usage

```typescript
import { tmAlign } from 'molstar/lib/mol-model/structure/structure/util/tm-align';
import { StructureElement } from 'molstar/lib/mol-model/structure';

// Get C-alpha Loci from two structures (see selection examples above)
const loci1: StructureElement.Loci = /* ... */;
const loci2: StructureElement.Loci = /* ... */;

// Run TM-align
const result = tmAlign(loci1, loci2);

console.log('TM-score (normalized by structure 1):', result.tmScoreA);
console.log('TM-score (normalized by structure 2):', result.tmScoreB);
console.log('RMSD:', result.rmsd);
console.log('Aligned residues:', result.alignedLength);

// result.bTransform is a Mat4 to transform structure2 onto structure1
```

### TM-align Result

The `tmAlign` function returns a `TMAlignResult` object with the following properties:

| Property | Type | Description |
|----------|------|-------------|
| `bTransform` | `Mat4` | Transformation matrix to superpose structure B onto A |
| `tmScoreA` | `number` | TM-score normalized by length of structure A |
| `tmScoreB` | `number` | TM-score normalized by length of structure B |
| `rmsd` | `number` | RMSD of aligned residue pairs (in Angstroms) |
| `alignedLength` | `number` | Number of aligned residue pairs |
| `sequenceIdentity` | `number` | Sequence identity of aligned residues (0-1) |
| `alignmentA` | `number[]` | Indices of aligned residues in structure A |
| `alignmentB` | `number[]` | Indices of aligned residues in structure B |

### Understanding TM-score

The TM-score is calculated as:

$$\text{TM-score} = \frac{1}{L} \sum_{i=1}^{L_{ali}} \frac{1}{1 + (d_i/d_0)^2}$$

Where:
- $L$ is the length of the reference protein
- $L_{ali}$ is the number of aligned residues
- $d_i$ is the distance between the $i$-th pair of aligned residues after superposition
- $d_0 = 1.24 \sqrt[3]{L - 15} - 1.8$ is a length-dependent normalization factor

**TM-score interpretation:**
- TM-score > 0.5: Generally indicates proteins with the same fold
- TM-score > 0.17: Generally indicates proteins with random structural similarity

### Low-level API

For direct coordinate-based alignment without structures, use the `TMAlign` namespace:

```typescript
import { TMAlign } from 'molstar/lib/mol-math/linear-algebra/3d/tm-align';

// Create position arrays
const posA = TMAlign.Positions.empty(lengthA);
const posB = TMAlign.Positions.empty(lengthB);

// Fill in coordinates
for (let i = 0; i < lengthA; i++) {
    posA.x[i] = /* x coordinate */;
    posA.y[i] = /* y coordinate */;
    posA.z[i] = /* z coordinate */;
}
// ... similarly for posB

// Compute alignment
const result = TMAlign.compute({ a: posA, b: posB });
```

### Complete Example: Aligning Two PDB Structures

```typescript
import { PluginContext } from 'molstar/lib/mol-plugin/context';
import { MolScriptBuilder as MS } from 'molstar/lib/mol-script/language/builder';
import { compile } from 'molstar/lib/mol-script/runtime/query/compiler';
import { StructureSelection, QueryContext, StructureElement } from 'molstar/lib/mol-model/structure';
import { tmAlign } from 'molstar/lib/mol-model/structure/structure/util/tm-align';
import { StateTransforms } from 'molstar/lib/mol-plugin-state/transforms';
import { Mat4 } from 'molstar/lib/mol-math/linear-algebra';

async function alignStructures(plugin: PluginContext, structure1: any, structure2: any) {
    // Query for C-alpha atoms in chain A
    const caQuery = compile<StructureSelection>(MS.struct.generator.atomGroups({
        'chain-test': MS.core.rel.eq([MS.struct.atomProperty.macromolecular.auth_asym_id(), 'A']),
        'atom-test': MS.core.rel.eq([MS.struct.atomProperty.macromolecular.label_atom_id(), 'CA'])
    }));

    // Get structure data
    const data1 = structure1.cell?.obj?.data;
    const data2 = structure2.cell?.obj?.data;

    // Create selections
    const sel1 = StructureSelection.toLociWithCurrentUnits(caQuery(new QueryContext(data1)));
    const sel2 = StructureSelection.toLociWithCurrentUnits(caQuery(new QueryContext(data2)));

    // Run TM-align
    const result = tmAlign(sel1, sel2);

    // Apply transformation to structure2
    const b = plugin.state.data.build().to(structure2)
        .insert(StateTransforms.Model.TransformStructureConformation, {
            transform: { name: 'matrix', params: { data: result.bTransform, transpose: false } }
        });
    await plugin.runTask(plugin.state.data.updateTree(b));

    return result;
}
```

## References

- Zhang Y, Skolnick J. "TM-align: a protein structure alignment algorithm based on the TM-score." *Nucleic Acids Research* 33, 2302-2309 (2005). DOI: [10.1093/nar/gki524](https://doi.org/10.1093/nar/gki524)
- Kabsch W. "A solution for the best rotation to relate two sets of vectors." *Acta Crystallographica* A32, 922-923 (1976).
