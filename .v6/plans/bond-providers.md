# Bond providers: design and plan

Status: implemented on `v6-bond-providers` (2026-10-20), steps 1–6; decisions in §7. Origin: molstar/molstar#1939 (Paul Pillot's
POC) and the review there: a provider is a static parameter of a `Unit`, never looked up at `unit.bonds` time; the
registry is on `PluginContext`; a transform resolves a provider by name and passes the instance to the `Structure`,
which passes it to its units. This plan covers the interface only. Bond-order perception itself stays in an extension.

Paths: `model/` is `packages/model/src/model/structure/structure/`, `plugin/` is `packages/plugin/core/src/`.

## 1. Goals and non-goals

- Anyone can supply intra-unit bonds, inter-unit bonds, or both, and keep the built-in computation as a baseline to
  refine. Bond-order perception is the first use; bonds read from a sidecar file, or no bonds at all, work the same way.
- The built-in computation is the default and the fallback: no provider, an unknown name, a provider that does not
  apply to the model, or an empty registry all mean the built-in bonds, with a warning where a snapshot names one that
  is not registered.
- **By default nothing changes.** Every structure gets exactly today's bonds unless something names another provider:
  - the `bonds` param is `PD.Optional` and absent by default;
  - absent, `'default'` and unknown names all resolve to `BondProvider.DefaultInstance`, which runs today's
    computation through today's cache;
  - registering `DefaultBonds` only adds options;
  - snapshots and MVS files without the field load as before, and a structure built outside the state tree
    (`Structure.ofModel(model)`) gets the default instance.

  The existing bond tests stay untouched and pass as the proof.
- **Static.** Two reads of `unit.bonds` return the same object. The provider is fixed when a unit is constructed and
  never changes; switching provider makes a new `Structure`.
- **Serializable selection.** The state tree holds a provider name and its JSON params. Provider objects are never
  serialized.
- Per-unit choice (§5), serializable, without a second mechanism in the model layer.
- Not in scope: perception algorithms, MVS syntax (§8), changing bonds at run time.

## 2. Model layer: `@molstar/model`

New `model/unit/bonds/provider.ts`:

```ts
/** A registered kind of bonds: named, with params; makes instances. */
export interface BondProvider<P extends PD.Params = PD.Params> {
  readonly name: string; // registry key; also what the state tree stores
  readonly label: string;
  readonly description?: string;
  /** Params for a model; called with no model when a param definition is built without data (must not throw). */
  getParams(model: Model | undefined): P;
  /** Not applicable: not offered for the model, and a stored choice falls back to the built-in bonds. Default: true. */
  isApplicable?(model: Model): boolean;
  create(props: PD.Values<P>): BondProvider.Instance;
}

export namespace BondProvider {
  /** What a Unit and a Structure hold. Model-independent: survives remapModel. */
  export interface Instance {
    readonly provider: BondProvider;
    /** Equal keys ⇒ equal bonds for equal (model, elements). Keys the intra cache; undefined disables caching. */
    readonly cacheKey: string | undefined;
    /** How a child unit gets its bonds when its parent has none yet (§2.3). Default 'self'. */
    readonly childBonds?: 'self' | 'root';
    intra?(unit: Unit.Atomic, ctx: IntraContext): IntraUnitBonds;
    inter?(structure: Structure, ctx: InterContext): InterUnitBonds;
  }
  export interface IntraContext {
    /** The built-in computation for this unit (computeIntraUnitBonds), computed on demand; `noDistance` §2.4. */
    builtIn(props?: { noDistance?: boolean }): IntraUnitBonds;
    /** The unit's model as a structure with the built-in bonds, once per model (§2.5). */
    modelStructure(): Structure;
  }
  export interface InterContext {
    builtIn(props?: Partial<InterBondComputationProps>): InterUnitBonds;
  }

  /** The built-in bonds: the default and the fallback. Name 'default'. */
  export const Default: BondProvider<{}>;
  export const DefaultInstance: Instance; // shared; what every unit gets unless told otherwise
}
```

- Both hooks are in from the start (Q1). A missing `intra` or `inter` means the built-in computation for that half. A
  provider that only refines orders implements `intra` as `refine(ctx.builtIn())`; one that perceives across
  struct_conn links (glycans, covalent ligands) implements `inter` too.
- `Instance` is not bound to a model or a structure (the POC's `context.structure` assignment after building is what
  this replaces). A provider that needs neighbours outside the unit reads them from the unit's model
  (`ctx.modelStructure()`, §2.5); one that needs the structure implements `inter`.
- `getParams(model)` may depend on the model (Q5): the transform has the model when it builds its params (§4).

### 2.1 Unit

- `Unit.Atomic` gets a constructor parameter `bondProvider: BondProvider.Instance` (a readonly field, not in `_props`,
  which symmetry copies share by object but which `remapModel` rebuilds). `Unit.create(..., bondProvider =
  DefaultInstance)`. Coarse units: no bonds getter today; leave them out.
- Propagation: `getChild`, `getCopy`, `applyOperator` and `remapModel` pass `this.bondProvider` on. Every other path
  goes through `Unit.create`, so the callers below pass it.
- `bonds` getter (unit.ts:441), in order:
  1. the parent's bonds restricted (unit.md §4) when they exist and the parent has the same instance;
  2. with `childBonds: 'root'` and a parent: the parent's `bonds` getter (computing the root unit's bonds), then
     restricted;
  3. `ElementSetIntraBondCache.get(model, instance.cacheKey)`, skipped when `cacheKey` is undefined;
  4. `instance.intra ? instance.intra(this, ctx) : computeIntraUnitBonds(this)`; cached when `props.cacheable`, as
     today.
- `ElementSetIntraBondCache` is keyed by model and elements only (in `model._dynamicPropertyData`). It becomes one
  cache per `cacheKey`, with the built-in under `'default'`. Without this, two structures of one model with different
  providers would share bonds.
- `tryRemapBonds` keeps working on `IntraUnitBondProps.canRemap`; a provider sets it on what it returns.

### 2.2 Structure

- `Structure.Props.bondProvider?: BondProvider.Instance`, held in `state`. `Structure.create` inherits it from
  `parent` like `dynamicBonds` (structure.ts:747); `asParent` copies state; `remapModel` copies it explicitly.
- `StructureBuilder.addUnit` passes `props.bondProvider` to `Unit.create`. That covers `ofModel` and the symmetry
  builders, provided the props carry it.
- Paths that build a new root without a parent and must pass it explicitly: `StructureSymmetry.buildAssembly`,
  `buildSymmetryAssembly`, `assembleOperators`, `findMatesRadius` (symmetry.ts:33, 66, 256, 317), `mergeUnitsByOperator`
  (structure.ts:945) and `mergeUnitsWithSameOperator` / `partitionUntransformedUnits` (util/unit-merging.ts:59, 75),
  and `ofTrajectory` (structure.ts:773, which takes props for the first time).
- `interUnitBonds` getter: restrict the parent's bonds only when the instance is the parent's too (structure.ts:240,
  next to the `dynamicBonds` check); else `instance.inter ? instance.inter(this, ctx) : computeInterUnitBonds(...)`.
- A structure's instance decides inter-unit bonds; each unit's instance decides its intra-unit bonds. A root structure
  built by the builder has one instance everywhere. A `Structure.create(units, ...)` over units of different
  instances is allowed and is what §5 builds on.

### 2.3 Contract for a provider (goes into unit.md §4 and structure.md §5)

- **Child bonds are the provider's choice (Q2).** `'self'` (the default, today's rule): a child restricts the parent's
  bonds when they exist and otherwise computes its own; right for providers whose decisions are per atom pair, and
  cheap. `'root'`: a child computes the root unit's bonds and restricts them; right for perception over rings or
  conjugated systems, where a residue alone can differ from its chain. Cost: a one-residue child pays for the whole
  chain, once (the root unit keeps them).
- **Determinism.** `intra(unit)` depends on the model, the elements and the props only, never on call order: that is
  what makes `cacheKey` sound.

### 2.4 The built-in bonds: `'default'` and `'none'`

The built-in computation stays as it is: format bonds and distance bonds are interleaved in one pass, and splitting it
into separate providers would fight that. One flag is added instead: `noDistance` in `BondComputationProps`
(`unit/bonds/common.ts:13`). The rest of the computation runs as today; only the distance step over the unit's grid is
skipped:

- in `findBonds` (intra-compute.ts:176): the grid pass after the dictionary and struct_conn bonds;
- in `findPairBonds` (inter-compute.ts:90): the distance branch.

Index pairs need nothing, because they never use distance (`forceCompute` aside). What still loads under `noDistance`
is what the format states:

- index pairs (`IndexPairBonds`: SDF/mol/mol2 bond tables, CONECT, `molstar_bond_site`);
- `struct_conn` records;
- `chem_comp_bond` (`ComponentBond`), matched per residue within `maxRadius` as today.

Two built-in providers, neither with params:

| Name        | Label                     | Bonds                                 | Where                                    |
| ----------- | ------------------------- | ------------------------------------- | ---------------------------------------- |
| `'default'` | Default                   | today's computation, unchanged        | the fallback; always there, not an entry |
| `'none'`    | None (format bonds only)  | the same computation with `noDistance` | `DefaultBonds` (§3)                      |

- **No provider never means no bonds.** `'none'` removes the computed bonds, not the format's: an SDF, a PDB with
  CONECT records or an mmCIF's struct_conn links keep their bonds. A scene of cartoons or surfaces pays for no distance
  search. There is no fully empty option; a provider that wants one returns `IntraUnitBonds.Empty`.
- With `'none'`, a PDB entry without `chem_comp_bond` has no peptide bonds within a chain, since those come from
  distance today.
- Each has its own `cacheKey`, so cached bonds of one never serve the other.
- `ctx.builtIn({ noDistance: true })` gives a perceiving provider the format's bonds alone, e.g. to assign orders
  without the distance guesses.

### 2.5 What a perceiving provider reads (molstar/molstar#2016 review)

Perceiving bond orders (#1939) needs, besides the unit's bonds:

- **The unit's links to other units** (a covalent ligand's to its protein). `intra` has no structure: a unit is shared
  by every structure that holds it and its bonds are cached by its elements, so a structure passed in would make the
  bonds depend on which structure asked first. `ctx.modelStructure()` is the unit's model as a structure
  (`Structure.ofModel`, the built-in bonds, default instances), built once per model and kept in its dynamic data: its
  `interUnitBonds` are the links, and a function of the model keeps `intra` one of the model, the elements and the
  props. Links made by symmetry operators are not in it. `inter` already receives the structure.
- **Which orders are unknown.** `BondType.Flag.OrderUnknown`: a covalent struct_conn bond whose order the file does not
  state and the inter-residue table does not know (order 1 then). Not one of `BondType.Names`: a property of a bond, not
  a type to select.
- **Rings over the graph being built.** `computeRings(unit, bonds)` and `UnitRings.create(unit, bonds)` (and
  `UnitRing.isAromatic(unit, ring, bonds)`): reading `unit.bonds` from inside `intra` would compute them again.

## 3. Plugin layer: registry and spec

- `BondProviderRegistry` (`plugin/state/bonds/registry.ts`), counted and name-keyed like `DataFormatRegistry`
  (`plugin/state/formats/registry.ts:27`): `add`, `remove`, `findConflict`, `has`, `get(name)` returning
  `BondProvider.Default` for `'default'` and `undefined` for an unknown name, `getApplicable(model?)` and `options(model?)`
  for the select. `'default'` is always present and cannot be registered or removed.
- **On `PluginContext` as `plugin.bondProviders`**, a top-level field next to `dataFormats` and the
  `customModelProperties` / `customStructureProperties` registries. Why there: those are the registries that decide how
  data is built (not how it is drawn, which is `plugin.representation.*`, or picked, which is `plugin.query.*`), and they
  are flat. A `plugin.structure.*` namespace does not exist and would be one field. The spec side is
  `PluginRegistryEntry.structure.bondProviders`, because spec entries are grouped by data kind.
- `PluginRegistryEntry.structure.bondProviders?: readonly BondProvider[]` (`plugin/spec.ts`), wired into
  `registry-entry.ts` `targets()` and `entryRegistrations()` before presets (plugin-composition §4.1). Pure data, like
  every other entry.
- **Dynamic add and remove**, as for representation providers (`RepresentationRegistry.add/remove`,
  graphics/repr/representation.ts:133): `plugin.bondProviders.add(provider)` / `remove(provider)` at any time, directly
  or through `plugin.register` (which returns the undo), from a behavior or an app. Same object again counts up, a
  different object under a taken name throws, removal counts down and drops at zero, removing an unknown provider is a
  no-op (plugin-composition §4.2).
  - Adding makes the provider an option the next time a `StructureFromModel` param definition is built (the params are
    a function of the plugin, built on demand), so a UI opened after the add offers it.
  - Removing never touches existing structures: their units hold the instance, which keeps working (§1, static). The
    next apply or update of a node that names it normalizes to `'default'` and recreates the structure with the
    built-in bonds, as a removed representation type falls back to the registry default. Nothing walks the state tree
    on removal.
  - Unlike the representation registries, the registry has no default to lose: `'default'` is not an entry, so removing
    every provider leaves a working select.
- **`DefaultRegistry` gets a `DefaultBonds` entry (Q4)** with the generic providers: `'none'` (§2.4; a scene of
  cartoons, surfaces or spacefill needs no computed bonds, and mesoscale scenes pay for them today) and `'per-unit'`
  (§5). Its own named entry, so an app that needs none of them leaves it out (plugin-composition §5.3 lists the
  entries; this adds one).
  `'default'` is not in it: it is the fallback and needs no registration.

## 4. Transforms and params

- `CommonStructureParams` (`plugin/state/helpers/root-structure.ts:20`) gets a third optional field, next to
  `dynamicBonds` (Q3):

  ```ts
  bonds: PD.Optional(
    PD.Mapped('default', registry.options(model), (name) => registry.get(name)?.getParams(model) ?? {}),
  );
  ```

  Stored as `{ name, params }`. The options are the registry's providers applicable to the model, `'default'` first;
  each case's params are that provider's `getParams(model)` (Q5). Placing it in the common params keeps it in every
  mode (`model`, `assembly`, the symmetry modes), in the hierarchy presets that take `RootStructureDefinition` params,
  and in Prism's `structure.assembly` and `symmetry` ops, which carry the common params over already
  (`prism/viewer/src/api/structure.ts:484`).
- `RootStructureDefinition.getParams(model, defaultValue, plugin?)`: without a plugin (`cli/state-docs`, data-less
  helpers) the select offers `'default'` only; without a model, every registered provider with `getParams(undefined)`.
  `StructureFromModel.params(a, plugin)` passes both.
- **Fallback is normalization.** `PD.normalizeParam` maps a `mapped` name that is not an option to the default
  (param-definition.ts:839), and the state re-runs it against the current definition on apply and on restore
  (core/state/state.ts:1371). So an unregistered provider, one not applicable to this model, and a stale params shape
  all land on `'default'` with no code of ours. `reportUnregisteredNames` (`plugin/state.ts:181`) gets a branch for
  `structure-from-model`'s `type.params.bonds.name`, so a snapshot's fallback is reported, not silent. Resolution in
  `create` still checks `get(name)` and `isApplicable` and falls back, for callers that bypass normalization.
- `RootStructureDefinition.create` resolves `props.bonds` once to an instance and passes it in the props of
  `Structure.ofModel` and the symmetry builders. One instance per apply.
- `StructureFromModel.update`: a change of `bonds` is a param change, so it recreates already. `remapModel` keeps the
  instance (it is model-independent). A frame whose model the provider does not apply to keeps the instance: the
  choice is per structure, made at apply.
- `StructureFromTrajectory` gets the same optional `bonds` param (its first) and passes the instance to `ofTrajectory`.
- Out of the state tree: the raw `Structure.ofModel(model)` callers (particle targets, json-cif, mesoscale) keep the
  built-in bonds. Prism's sequence source is the exception (§4.1).
- `atomic-detail.ts:63` decides "bonds given" from the format providers. With a non-default instance on the structure
  it counts bonds as given too (step 2).

### 4.1 Prism

- **An API operation, `structure.bonds`** (`{ structure, provider, params? }`), next to `structure.assembly` and
  `structure.symmetry` (`prism/viewer/src/api/structure.ts`; `_spec/operations.md`). It updates the
  `StructureFromModel` node's `type.params.bonds` and keeps the rest; `commonTypeParams` (structure.ts:484) already
  carries it over through `structure.assembly` and `symmetry`. An unknown or non-applicable provider is refused with a
  message, not normalized silently, as the other operations refuse. Being an operation, it is in the history, the
  journal and the agent's tools for free.
- **A "Bonds" field in the source form** (`prism/viewer/src/sources/scene.ts`, `sourceForm` ~:1840): a select over
  `plugin.bondProviders.options(model)`, hidden when only `'default'` is there (so a viewer without `DefaultBonds` or
  extra providers looks as today). The chosen provider's params show under it from its `getParams(model)`. Edits go
  through `sourceChanges` to `structure.bonds`, with the slider debounce the radius uses.
- **The sequence source** (`prism/viewer/src/sequence/molstar/source.ts`) needs nothing: it builds structures from the
  files for Seq*'s kernels (sequences, pLDDT, secondary structure), never shown and never reading bonds.
- `DefaultPrismViewerSpec` gets `DefaultBonds` through `DefaultRegistry` (`sources/base.ts`), so nothing to add.
- Spec: `_spec/operations.md` (the operation), `parity.md` (the row for structure source).

### 4.2 MolViewSpec

Only through a Mol\*-specific custom param for now, no MVS schema change:

```json
{ "kind": "structure", "params": { "type": "assembly" }, "custom": { "molstar_bond_provider": { "name": "none" } } }
```

- `structureProps` (`mvs/runtime/src/load-helpers.ts:432`) merges `node.custom.molstar_bond_provider` (`{ name,
  params? }`) into `type.params.bonds`, as `representationProps` merges `molstar_representation_params` (:615).
- An unregistered name falls back to `'default'` by normalization, and the load logs a warning naming it.
- The MVS builder gets nothing: `custom` is already free-form there.
- Spec: the runtime's `_spec` section on `custom` keys.

## 5. Per-unit providers

The rule: per-unit choice is a provider, not a model-layer feature. The model layer already supports units with
different instances (§2.2); only the selection needs to be serializable.

**A composite provider** `'per-unit'` (label "By unit"), in `DefaultBonds`, with params:

```ts
{
  rules: PD.ObjectList({
    match: PD.MappedStatic('entity-type', {
      'entity-type': PD.MultiSelect(['non-polymer'], [['polymer', ...], ['non-polymer', ...], ['branched', ...], ['water', ...]]),
      'label-asym-id': PD.Text(''), // comma-separated
      'entity-id': PD.Text(''),
    }),
    provider: <the registry's Mapped select for the model, without 'per-unit'>,
  }, ...),
  inter: <the registry's Mapped select>, // who computes inter-unit bonds; 'default' by default
  // units no rule matches get the built-in bonds
}
```

- `create(props)` makes the child instances once, and its `intra(unit)` hands the unit to the first rule that matches.
  `cacheKey` combines the children's keys; `childBonds` is per child (the composite asks the matched child).
- Rules match the unit's **first element** on properties that are the same for every element of a chain and every
  operator: entity type, entity id, label_asym_id. A function of the model and the elements alone, as the intra-unit
  cache requires (equal elements, equal bonds; matching on the root broke that for children of multi-chain units). A
  symmetry copy and a remapped frame have the elements and so the match; a child of a unit with several chains (merged
  waters) matches by its own first chain.
- Why not per-element selections (MolScript, MVS selectors): bonds are a per-unit computation, so a selection that
  splits a unit has no meaning, and expression rules would need a structure to evaluate against before the structure
  exists.
- Rejected: a `(unit) => Instance` function on `Structure.Props`. Serializable only if the rule list is plugin data,
  which is this provider again, with a function on the `Structure` state that every derivation must carry.

## 6. Steps

1. **Model:** `BondProvider`, `Default`, the `Unit` field and getter (including `childBonds`), the per-key
   `ElementSetIntraBondCache`, `Structure.Props.bondProvider` and its propagation through every path in §2.2. Quick
   tests in `model/.../_test/`: a counting provider is called once per unit; children, symmetry copies and
   `remapModel` keep the instance; two structures of one model with different providers do not share cached bonds;
   `childBonds: 'root'` restricts the root's bonds; `inter` is called and restricted only under the same instance; the
   default path is identical to today (existing bond tests unchanged).
2. **Plugin:** the registry on `plugin.bondProviders`, the spec field, `registry-entry.ts` wiring, the `bonds` common
   param with model-dependent options and params, resolution in `RootStructureDefinition.create` and
   `StructureFromTrajectory`, the `reportUnregisteredNames` branch. Tests: an unknown name and a non-applicable
   provider fall back (with a warning on restore); a snapshot round-trips the name and params; registering twice
   counts; params build with no plugin and no model; with `DefaultBonds` registered, a structure built with no `bonds`
   param, and a snapshot saved before this change, both get `DefaultInstance`; a provider added after init appears in the options, and one
   removed leaves existing structures' bonds untouched and the next update falls back to `'default'`.
3. **`DefaultBonds`:** `noDistance` in the intra and inter computations; `'none'` and `'per-unit'` (§2.4, §5). Tests:
   `'none'` on an SDF and on CONECT records equals `'default'`; on an mmCIF it keeps struct_conn and dictionary bonds
   and drops the distance ones, intra and inter; `'default'` is unchanged; each per-unit match kind, and a child of a
   matched
   unit.
4. **Prism example, `examples/bond-providers`.** A Vite app modelled on `examples/minimal-prism-viewer`, added to
   `scripts/workspace/inventory.json` and the root `tsconfig.json` references.
   - **Composition:** the minimal spec plus `DefaultBonds` and the example's own providers, registered through
     `registry.structure.bondProviders` in the spec.
   - **Naive providers**, written as a third party would, against the public interface only:
     - `'naive-distance'` (param `cutoff`, default 1.9 Å): every atom pair under the cutoff is bonded, ignoring the
       dictionary and struct_conn. It shows a provider that replaces both halves, `intra` and `inter`.
     - `'naive-orders'`: takes `ctx.builtIn()` and marks C–C bonds shorter than 1.42 Å as double. It shows a provider
       that refines the built-in bonds, `intra` only.
   - **The page:**
     - a header with the entry field (default a ligand-rich entry, e.g. 1HSG);
     - the provider is picked where a Prism user would pick it: the structure's "Bonds" field in the inspector (§4.1),
       so the spec keeps the Structure feature's inspector;
     - header buttons that run `structure.bonds` for the presets worth comparing: default, none, naive distance, naive
       orders, and "By unit" (`'per-unit'`: naive orders for non-polymer entities, the default for the rest).
   - **Run-time registration:** a button adds and removes `'naive-orders'`, to show that the field follows the
     registry and that removing a provider falls back to `'default'` on the next update.
   - **Not shown:** Paul's perception code. It stays his PR, and can move onto this interface there.
5. **Prism and MVS** (§4.1, §4.2): the `structure.bonds` operation, the "Bonds" field, the sequence source, `DefaultBonds`
   in `DefaultPrismViewerSpec`; `molstar_bond_provider` in `structureProps`. Tests: the operation updates only `bonds`
   and survives `structure.assembly`; the field is hidden with only `'default'`; an MVS file with the custom key builds
   with that provider, and one without it builds as before.
6. **Docs**: unit.md §1/§4, structure.md §2/§5, the plugin-composition design §3 and §5.3 (the new entry and
   `DefaultBonds`), CHANGELOG `[Unreleased]`, an entry file in `.v6/changelogs/`, and `breaking-v6-changes.md` if the
   `Unit` constructor change is counted as public.

## 7. Decisions (2026-10-20)

- **Q1. Both hooks.** `intra` and `inter` from the start.
- **Q2. Child bonds per provider.** `childBonds: 'self' | 'root'` on the instance, default `'self'`.
- **Q3. In `CommonStructureParams`**, next to `dynamicBonds`.
- **Q4. In `DefaultRegistry`** as its own `DefaultBonds` entry, with `'none'` and `'per-unit'`.
- **Q5. Model-dependent params.** `getParams(model | undefined)` and `isApplicable(model)`; the select is a `PD.Mapped`
  derived from the registry for the model, and normalization is the fallback to `'default'`.
- **Q6. The built-in computation is not split.** `'none'` is the built-in with `noDistance`: format bonds (index
  pairs, struct_conn, chem_comp_bond) always load; only distance bonds are dropped (§2.4).
- **Q7. Run-time add and remove**, as for representation providers (§3).

- **Q8. Prism now**: an operation and a field in the source form (§4.1).
- **Q9. MVS through `custom.molstar_bond_provider` only** (§4.2); a schema field when the interface has settled.

## 8. Deferred

Each deferred item is an entry in the TODO of the area that needs it:

- `packages/prism/TODO.md`: the "Bonds" field refreshes its options when the registry changes. None of the registries
  has a change event today.
- `.v6/plans/checklist.md`, "Remaining architecture and release work": a `bonds` field in the MVS schema in place of
  the custom key. The MVS packages have no TODO file, and the checklist holds their open items.

Out by design, so not a TODO: switching bonds without a new structure (§1).
