# Bond orders

`@molstar/bond-orders-extension`: double, triple and aromatic bonds perceived from the 3D coordinates where the file
does not state them, as a bond provider ([bond providers plan](../../../.v6/plans/bond-providers.md)).

```ts
import { BondOrders } from '@molstar/bond-orders-extension';

// add the entry to a plugin spec's registry, e.g. the default Prism viewer's
const spec = DefaultPrismViewerSpec();
const withBondOrders = { ...spec, plugin: { ...spec.plugin, registry: [...(spec.plugin?.registry ?? []), BondOrders] } };
```

A structure then offers "Bond order perception" in its `bonds` param (Prism: the structure's Bonds field, or the
`structure.bonds` operation), with two modes:

- `auto` (default): the single bonds whose order the file does not give, in residues other than the canonical ones (a
  file with a bond table, such as SDF or MOL2, is left as it is);
- `forceCompute`: every covalent bond of every residue, stated orders included.

A running viewer: [`examples/bond-providers`](../../../examples/bond-providers/README.md), beside other providers.
Design and decisions: [\_spec](_spec/README.md); open items: [TODO](TODO.md).
