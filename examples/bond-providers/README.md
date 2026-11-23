# Bond providers

A [Prism](../../packages/prism/AGENTS.md) viewer that shows where a structure's bonds can come from
([bond-providers plan](../../.v6/plans/bond-providers.md)): the built-in bonds, the generic providers of the default
registry, three naive providers of the example's own, registered as any registry entry is, and bond-order perception
from the [bond-orders extension](../../packages/extensions/bond-orders/README.md).

```sh
pnpm dev:examples -- bond-providers   # http://localhost:1338/examples/bond-providers/build/, watching the sources
pnpm -C examples/bond-providers build # esbuild → build/
pnpm -C examples/bond-providers typecheck
```

Built like the other examples (`scripts/esbuild/app.mjs`), as an ES module with lazy chunks (`"format": "esm"`, as the
Prism viewer app), from the workspace sources, so nothing needs building first.

## The providers

| Name                    | From                                     | What it does                                                                                                                                                           |
| ----------------------- | ---------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `default`               | always there                             | the built-in bonds: the file's, then computed from distances and the component dictionary                                                                              |
| `none`                  | `DefaultBonds`                           | the file's bonds only (bond tables, CONECT, struct_conn, the dictionary), none from distances                                                                          |
| `per-unit`              | `DefaultBonds`                           | another provider per unit, by entity type, chain or entity; here naive orders for the ligands                                                                          |
| `naive-distance`        | this example (`src/providers.ts`)        | replaces the bonds: every pair of atoms closer than a cutoff; between units, waters and ions too                                                                       |
| `naive-orders`          | this example                             | refines the built-in bonds: short C–C single bonds become double; intra-unit only                                                                                      |
| `random`                | this example                             | random bonds between atoms within a max distance, with random orders; the same for the same seed; none between units                                                   |
| `bond-order-perception` | the bond-orders extension (`BondOrders`) | perceives double, triple and aromatic bonds from the coordinates: `auto` the orders the file does not state, outside the canonical residues; `forceCompute` every bond |

The example's providers use the public interface only (`BondProvider`, the unit, the model, a graph builder): one
replaces both halves (`intra`, and `inter` through the built-in computation with other props), the other refines the
built-in intra-unit bonds (`ctx.builtIn()`), and `random` is deterministic (a hash of the seed and the two atoms), so
the same params give the same bonds, as a provider must. None is meant to be right.

## The page

One page, two viewers side by side, their cameras linked (`syncCameras`): the left one shows a file as served, the right
one the same file with the bonds of the provider under test. Built for the bond-order benchmark's failures (its
[plan](../../packages/extensions/bond-orders/plans/benchmark.md) §7), it shows any provider.

- **Header**: a PDB ID or a CCD code (4 characters starting with a digit read as a PDB ID; the PDB / CCD switch
  overrides the guess), and for a CCD code its ideal or model coordinates. An entry is `models.rcsb.org`'s BinaryCIF, a
  component the CCD's mmCIF file; downloaded once, the bytes given to both viewers.
- **Left**: the built-in bonds, the file's orders. No controls.
- **Right**: the bond provider (bond-order perception by default, the providers above in the list), its mode for
  perception (`forceCompute` by default; `auto` disabled when the model states its orders, as these files do), and two
  removals made on the parsed block, never on the file ([`src/removals.ts`](src/removals.ts), the example's copy of the
  benchmark's `benchmark/src/input.ts`): Remove H (off) and Remove charges (on). Model 1 and a residue's first altloc
  always. The provider and the mode switch the bonds in place (`structure.bonds`, as the Bonds field of the structure's
  inspector); a removal re-parses the right side only, keeping the camera.
- **URL**: follows the state, so a view is a link; the benchmark report's links fill it:

  ```
  ?pdb=1HSG                                   an entry
  ?ccd=197&coords=model                       a CCD component, model coordinates (ideal by default)
  ?pdb=1O4I&asym=B&seq=300                    an entry, focused on a copy (label_asym_id, auth_seq_id)
  &bonds=naive-orders  &mode=auto  &removeH=1  &removeCharges=0
  ```

The right side is the current code: an edit to the algorithm shows on reload, while the report's 2D drawings are the
run's. Ball-and-stick draws Kekulé orders rather than dashed aromatic bonds: the page patches that param's default
before the viewers start, as no plugin setting does it.
