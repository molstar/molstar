## Bond providers

- Pluggable bonds: a `BondProvider` (`@molstar/model/.../unit/bonds/provider`) supplies a structure's intra-unit bonds,
  inter-unit bonds or both, starting from the built-in computation (`ctx.builtIn()`) or replacing it. Its instance is a
  constructor argument of a unit and a prop of a structure (`Structure.Props.bondProvider`), fixed for their life and
  passed on to children, symmetry copies, assemblies, merged units and remapped frames; the intra-unit bond cache is
  one per provider. Without one, the bonds are exactly the built-in ones, as before (#2016)
- `PluginContext.bondProviders`: the plugin's providers, registered by a spec (`registry.structure.bondProviders`) or
  added and removed at run time. `StructureFromModel` (in every structure type) and `StructureFromTrajectory` take an
  optional `bonds` param (`{ name, params }`); a name the plugin lacks, or a provider not applicable to the model, falls
  back to the built-in bonds, with a warning when a snapshot names it (#2016)
- `DefaultRegistry` gets a `DefaultBonds` entry: `'none'` (only the bonds the file states: bond tables, CONECT,
  struct_conn, the component dictionary; none from distances) and `'per-unit'` (other providers chosen per unit by
  entity type, chain or entity). `BondComputationProps.noDistance` is the flag behind `'none'` (#2016)
- Prism: the `structure.bonds` operation and a Bonds field in a structure's inspector (hidden when the built-in bonds
  are the only choice). MolViewSpec: a structure node's `custom.molstar_bond_provider` (#2016)
- `examples/bond-providers`: a Prism viewer with three naive providers of its own (distance, bond orders, random) (#2016)
- For providers that perceive bond orders: `ctx.modelStructure()` in `intra` (the unit's model as a structure with the
  built-in bonds, once per model, for a unit's links to other units), `BondType.Flag.OrderUnknown` on covalent
  struct_conn bonds whose order the file does not state, and `computeRings(unit, bonds)` / `UnitRings.create(unit,
  bonds)` over a graph being built. Bond hashes pinned in tests change with the new flag (#2016)
