# Bond orders: perception from 3D coordinates

`@molstar/bond-orders-extension`: double, triple and aromatic bonds perceived from the coordinates where the file does
not state them (PDB entries, CONECT records, ligands without `chem_comp_bond`), offered as a bond provider (bond
providers plan, `.v6/plans/bond-providers.md`). This spec is for working on the code; code comments cite it
(`_spec §1.5`). Open items: [../TODO.md](../TODO.md). The algorithm (Sayle 2001, Labute 2005, Zhang 2012) is described
at the top of `src/perception/perceive.ts` and in each step's module.

## 1. Decisions

1. **A bond provider, opt-in**: `BondOrders` is a `PluginRegistryEntry` (`structure.bondProviders`) an app adds to its
   spec; the choice is the structure's `bonds` param (`{ name: 'bond-order-perception', params: { mode } }`). Why: the
   provider interface is how bonds are supplied (plan §1); a plugin without the entry carries none of the code.
2. **Perceived aromaticity is `BondType.Flag.Aromatic`**, the flag files use; no flag of its own. Why: one meaning for
   consumers (rendering, queries, export); an aromatic bond the file states is not perceivable (§1.9), so it is never
   confused with a perceived one during perception. Perception writes orders and `Aromatic` only: the other flags stay
   the input's (`Computed`: found by distance; `OrderUnknown`: order not given), `forceCompute` adding `OrderUnknown`
   (§1.9). Why: an extension does not change what a flag means to its consumers (excluding `computed` bonds in a
   representation hides the distance bonds, never a CONECT bond perception gave an order).
3. **A function of the model, the elements and the mode** (unit.md §4, what makes `cacheKey` sound): coordinates are the
   model's (`conformation.invariantPosition`), never through the unit's operator; a unit's neighbours in other units
   come from `IntraContext.modelStructure()` (`perception/inter-unit-neighbours.ts`, once per model), never from a
   structure. Why: a unit is shared by every structure that holds it and its bonds are cached by its elements.
   Consequence: bonds that exist only between symmetry copies are not seen.
4. **Canonical residues are known** (`perception/known-components.ts`): the 20 standard amino acids and the standard
   nucleotides, whose orders the built-in bonds always give (core order table, backbone C=O and P=O rules). `auto` skips
   them. Why not `AminoAcidNames`: it holds modified residues (PTR, SEP, PCA, the D-forms) whose side-chain orders the
   core table does not give.
5. **A model-wide pattern cache** (`perception/chemcomp-cache.ts`): a residue's perceived multiple and aromatic bonds,
   by atom-name pair, reused for the residues of the same key: the mode, the component id, its sorted atom names (an
   incomplete residue gets its own key) and its links, atoms bonded outside the residue with their partner
   (`C12-SER.OG`). Why: a ligand repeated across chains is perceived once. Known limitation (a deviation from unit.md
   §4): two complete copies with the same key but different geometry get the pattern of whichever unit is perceived
   first.
6. **No closures in the neighbour loops**: `State` collects per-residue adjacency arrays once (intra-unit from the graph
   being perceived, inter-unit neighbours from §1.3); the steps loop over arrays and over the graph's `offset`/`b`. Why:
   `unit.bonds` and `unit.rings` are not available inside `intra` (rings come from `computeRings(unit, bonds)` and
   `createIndex(rings, [])`, without aromaticity), and the core iterators read them.
7. **Output props**: `cacheable: true` (§1.3), and `canRemap: false` once a residue was perceived. Why: perceived orders
   depend on the geometry, so a frame with `dynamicBonds` is perceived again even where the file's bonds would carry
   over; a unit with only canonical residues keeps the built-in props.
8. **`childBonds: 'root'`**: a child unit gets the root unit's perceived bonds, restricted. Why: perception spans rings
   and conjugated systems; a residue alone could differ from its chain.
9. **Modes**: `auto` perceives the single covalent bonds whose order nobody stated (`Computed`: distance bonds;
   `OrderUnknown`: struct_conn/CONECT without an order), outside canonical residues; `forceCompute` resets every
   covalent bond to single, unstated, non-aromatic, and perceives every residue. One instance per mode. In `auto`, a
   residue without a perceivable bond inside it is skipped (`hasPerceivableBond`; its bonds to other residues do not
   count: a modified residue in a chain, its own orders stated and its peptide bonds found by distance, stays as it is),
   and a model with a bond table (`IndexPairBonds`: MOL/SDF, MOL2, CIF core, LAMMPS, topologies, `molstar_bond_site`)
   returns the built-in bonds untouched: the table is then the whole built-in graph, and none of its bonds is flagged
   unknown. Why: a residue whose own orders are stated has nothing to perceive, and perception would only cost time and
   `canRemap`; the inter-unit neighbours (§1.3) are built only once a residue needs perceiving. Consequence: tables
   without orders (LAMMPS, topologies, CIF core without valence, mol2 `un`) are not perceived in `auto` (TODO).
10. **No benchmark here yet**: validation over a large set (the Labute 2005 set, or the whole PDB, possibly compared
    with RDKit) is a project of its own (TODO). The tests are hard cases with their fixtures (§4).

## 2. Layout

| Path                                               | What                                                                        |
| -------------------------------------------------- | --------------------------------------------------------------------------- |
| `src/index.ts`                                     | `BondOrders` (the registry entry), the provider, the types                  |
| `src/provider.ts`                                  | `BondOrderPerceptionProvider`: params, one instance per mode, `intra`       |
| `src/perception/perceive.ts`                       | `perceiveBondOrders`: the copy of the built-in bonds, rings, residue loop   |
| `src/perception/common.ts`                         | `State` (per residue), bond and position helpers, `isPerceivable`           |
| `src/perception/inter-unit-neighbours.ts`          | neighbours in other units, per model (§1.3)                                 |
| `src/perception/hybridization.ts`                  | geometry: PCA dimensionality, angles, ring planarity, ambiguities           |
| `src/perception/functional-groups-connectivity.ts` | functional-group patterns (Sayle step 7)                                    |
| `src/perception/huckel.ts`                         | aromatic ring systems, Hückel 4n+2 (Sayle step 8)                           |
| `src/perception/kekule.ts`                         | double and triple bond assignment, weighted matching (Sayle step 9, Labute) |
| `src/perception/thresholds.ts`                     | bond-length and planarity thresholds                                        |
| `src/perception/chemcomp-cache.ts`                 | the pattern cache (§1.5)                                                    |
| `src/perception/known-components.ts`               | canonical residues (§1.4)                                                   |

## 3. Use

```ts
import { BondOrders } from '@molstar/bond-orders-extension';
const spec = DefaultPrismViewerSpec();
// then a structure's Bonds field offers "Bond order perception"
const withBondOrders = { ...spec, plugin: { ...spec.plugin, registry: [...(spec.plugin?.registry ?? []), BondOrders] } };
```

Or in a state tree: a `StructureFromModel` with `bonds: { name: 'bond-order-perception', params: { mode: 'auto' } }`. A
running example: `examples/bond-orders`.

## 4. Tests

Quick lane, no network (`src/_test/`). Fixtures in `src/_test/fixtures/`, named entry id + the residues kept:

- `1nco-chr.pdb` (CHR extract) against `ccd-chr.cif` (its CCD definition): expected multiple and aromatic bonds read
  from the CCD at test time; missed (forward) and spurious (reverse) bonds by atom-name pair.
- `5kre-ser124-6wg.pdb` (SER A 124, 6WG and their CONECT records): the link `C12-SER.OG` is seen from the ligand unit
  and is in its pattern key.
- Basic cases: benzene, CONECT benzene, exocyclic amino, carboxylate, glycine, altlocs.
- The residue gate: a residue whose orders are all stated is left as is, without building the inter-unit neighbours.
- The provider: modes (a mol2's bond table untouched in `auto`), cache shared by two structures of a model, one instance
  per mode, the registry entry.

## 5. Index names

Four kinds of index meet in the perception code; each has its own names, so a mix-up reads wrong at a glance:

| Names          | Index                                                      | Reads                                                         |
| -------------- | ---------------------------------------------------------- | ------------------------------------------------------------- |
| `i`, `j` (`k`) | local: a residue atom in `State` (`0` to `state.n`)        | `State` arrays (`geometry[i]`, `heavyNeighbours[i]`)          |
| `u`, `v`       | `UnitIndex`: an atom of the unit (`unit.elements[u]`)      | the bond graph (`offset[u]`, `b[t]`), `getOrder(bonds, u, v)` |
| `eI`           | `ElementIndex`: an atom of the model                       | model columns, coordinates, `InterUnitNeighbours`             |
| `t`            | an edge of the bond graph (`offset[u]` to `offset[u + 1]`) | `order[t]`, `flags[t]`, `assignedBonds[t]`                    |

Conversions: `u = state.unitIndices[i]`, `eI = unit.elements[u]`. `e` is an edge object (`interUnitBonds.edges`), not an
index. Why: `State` holds one residue, its steps read local arrays and write the graph, and every function takes one
kind or the other; Mol\* has no shared convention (`rings/compute.ts` uses `i` and `a`/`b` for local and unit indices
alike).
