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

[`src/index.ts`](src/index.ts) loads a PDB entry, polymers as cartoon and ligands as ball-and-stick: 1HSG (HIV protease
with an inhibitor), 1NCO (the neocarzinostatin chromophore CHR: double, triple and aromatic bonds, a non-aromatic
five-membered ring) or 5KRE (6WG, a covalent inhibitor bonded to SER A 124). It reads the PDB format: the connectivity
the file states (CONECT, LINK) without orders, so `auto` has orders to perceive; the mmCIF files state the ligands'
orders (`chem_comp_bond`), which `auto` keeps. Ball-and-stick draws Kekulé orders rather than dashed aromatic bonds: the
page patches that param's default before the viewer starts, as no plugin setting does it. The header switches the
structure's bonds through the `structure.bonds` operation, as the Bonds field of the structure's inspector does
(Structure panel, pick the structure). A button unregisters and registers `naive-orders` at run time: the field and the
operation follow the registry, and a structure keeps the bonds it has until it is rebuilt.
