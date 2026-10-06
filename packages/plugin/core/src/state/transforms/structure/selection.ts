/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Adam Midlik <midlik@gmail.com>
 * @author Ludovic Autin <ludovic.autin@gmail.com>
 */

import { PluginStateTransform, PluginStateObject as SO } from '@molstar/plugin/state/objects';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { MolScriptBuilder } from '@molstar/model/script/language/builder';
import type { Expression } from '@molstar/model/script/language/expression';
import { StructureQueryHelper } from '@molstar/plugin/state/helpers/structure-query';
import {
  StructureSelection as Sel,
  Structure,
  StructureElement,
  type StructureQuery,
  Queries,
  QueryContext,
} from '@molstar/model/model/structure';
import { StateObject, StateTransformer } from '@molstar/core/state';
import { Script } from '@molstar/model/script/script';
import {
  polymer,
  protein,
  nucleic,
  branchedPlusConnected,
  ligandPlusConnected,
  coarse,
} from '@molstar/plugin/state/queries/structure/type';
import { nonStandardPolymer } from '@molstar/plugin/state/queries/structure/residue';
import { assertUnreachable } from '@molstar/core/util/type-helpers';
import {
  StructureComponentParams,
  createStructureComponent,
  updateStructureComponent,
} from '@molstar/plugin/state/helpers/structure-component';

export { StructureSelectionFromExpression };
type StructureSelectionFromExpression = typeof StructureSelectionFromExpression;
const StructureSelectionFromExpression = PluginStateTransform.BuiltIn({
  name: 'structure-selection-from-expression',
  display: { name: 'Selection', description: 'Create a molecular structure from the specified expression.' },
  from: SO.Molecule.Structure,
  to: SO.Molecule.Structure,
  params: () => ({
    expression: PD.Value<Expression>(MolScriptBuilder.struct.generator.all, { isHidden: true }),
    label: PD.Optional(PD.Text('', { isHidden: true })),
  }),
})({
  apply({ a, params, cache }) {
    const { selection, entry } = StructureQueryHelper.createAndRun(a.data, params.expression);
    (cache as any).entry = entry;

    if (Sel.isEmpty(selection)) return StateObject.Null;
    const s = Sel.unionStructure(selection);
    const props = { label: `${params.label || 'Selection'}`, description: Structure.elementDescription(s) };
    return new SO.Molecule.Structure(s, props);
  },
  update: ({ a, b, oldParams, newParams, cache }) => {
    if (oldParams.expression !== newParams.expression) return StateTransformer.UpdateResult.Recreate;

    const entry = (cache as { entry: StructureQueryHelper.CacheEntry }).entry;

    if (entry.currentStructure === a.data) {
      return StateTransformer.UpdateResult.Unchanged;
    }

    const selection = StructureQueryHelper.updateStructure(entry, a.data);
    if (Sel.isEmpty(selection)) return StateTransformer.UpdateResult.Null;

    StructureQueryHelper.updateStructureObject(b, selection, newParams.label);
    return StateTransformer.UpdateResult.Updated;
  },
  dispose({ b }) {
    b?.data.customPropertyDescriptors.dispose();
  },
});

export { MultiStructureSelectionFromExpression };
type MultiStructureSelectionFromExpression = typeof MultiStructureSelectionFromExpression;
const MultiStructureSelectionFromExpression = PluginStateTransform.BuiltIn({
  name: 'structure-multi-selection-from-expression',
  display: {
    name: 'Multi-structure Measurement Selection',
    description: 'Create selection object from multiple structures.',
  },
  from: SO.Root,
  to: SO.Molecule.Structure.Selections,
  params: () => ({
    selections: PD.ObjectList(
      {
        key: PD.Text(void 0, { description: 'A unique key.' }),
        ref: PD.Text(),
        groupId: PD.Optional(PD.Text()),
        expression: PD.Value<Expression>(MolScriptBuilder.struct.generator.empty),
      },
      (e) => e.ref,
      { isHidden: true },
    ),
    isTransitive: PD.Optional(
      PD.Boolean(false, {
        isHidden: true,
        description: 'Remap the selections from the original structure if structurally equivalent.',
      }),
    ),
    label: PD.Optional(PD.Text('', { isHidden: true })),
  }),
})({
  apply({ params, cache, dependencies }) {
    const entries = new Map<string, StructureQueryHelper.CacheEntry>();

    const selections: SO.Molecule.Structure.SelectionEntry[] = [];
    let totalSize = 0;

    for (const sel of params.selections) {
      const { selection, entry } = StructureQueryHelper.createAndRun(
        dependencies![sel.ref].data as Structure,
        sel.expression,
      );
      entries.set(sel.key, entry);
      const loci = Sel.toLociWithSourceUnits(selection);
      selections.push({ key: sel.key, structureRef: sel.ref, loci, groupId: sel.groupId });
      totalSize += StructureElement.Loci.size(loci);
    }

    (cache as object as any).entries = entries;

    const props = {
      label: `${params.label || 'Multi-selection'}`,
      description: `${params.selections.length} source(s), ${totalSize} element(s) total`,
    };
    return new SO.Molecule.Structure.Selections(selections, props);
  },
  update: ({ b, oldParams, newParams, cache, dependencies }) => {
    if (!!oldParams.isTransitive !== !!newParams.isTransitive) return StateTransformer.UpdateResult.Recreate;

    const cacheEntries = (cache as any).entries as Map<string, StructureQueryHelper.CacheEntry>;
    const entries = new Map<string, StructureQueryHelper.CacheEntry>();

    const current = new Map<string, SO.Molecule.Structure.SelectionEntry>();
    for (const e of b.data) current.set(e.key, e);

    let changed = false;
    let totalSize = 0;

    const selections: SO.Molecule.Structure.SelectionEntry[] = [];
    for (const sel of newParams.selections) {
      const structure = dependencies![sel.ref].data as Structure;

      let recreate = false;

      if (cacheEntries.has(sel.key)) {
        const entry = cacheEntries.get(sel.key)!;
        if (StructureQueryHelper.isUnchanged(entry, sel.expression, structure) && current.has(sel.key)) {
          const loci = current.get(sel.key)!;
          if (loci.groupId !== sel.groupId) {
            loci.groupId = sel.groupId;
            changed = true;
          }
          entries.set(sel.key, entry);
          selections.push(loci);
          totalSize += StructureElement.Loci.size(loci.loci);

          continue;
        }
        if (entry.expression !== sel.expression) {
          recreate = true;
        } else {
          // TODO: properly support "transitive" queries. For that Structure.areUnitAndIndicesEqual needs to be fixed;
          let update = false;

          if (!!newParams.isTransitive) {
            if (Structure.areUnitIdsAndIndicesEqual(entry.originalStructure, structure)) {
              const selection = StructureQueryHelper.run(entry, entry.originalStructure);
              entry.currentStructure = structure;
              entries.set(sel.key, entry);
              const loci = StructureElement.Loci.remap(Sel.toLociWithSourceUnits(selection), structure);
              selections.push({ key: sel.key, structureRef: sel.ref, loci, groupId: sel.groupId });
              totalSize += StructureElement.Loci.size(loci);
              changed = true;
            } else {
              update = true;
            }
          } else {
            update = true;
          }

          if (update) {
            changed = true;
            const selection = StructureQueryHelper.updateStructure(entry, structure);
            entries.set(sel.key, entry);
            const loci = Sel.toLociWithSourceUnits(selection);
            selections.push({ key: sel.key, structureRef: sel.ref, loci, groupId: sel.groupId });
            totalSize += StructureElement.Loci.size(loci);
          }
        }
      } else {
        recreate = true;
      }

      if (recreate) {
        changed = true;

        // create new selection
        const { selection, entry } = StructureQueryHelper.createAndRun(structure, sel.expression);
        entries.set(sel.key, entry);
        const loci = Sel.toLociWithSourceUnits(selection);
        selections.push({ key: sel.key, structureRef: sel.ref, loci, groupId: sel.groupId });
        totalSize += StructureElement.Loci.size(loci);
      }
    }

    if (!changed) return StateTransformer.UpdateResult.Unchanged;

    (cache as object as any).entries = entries;
    b.data = selections;
    b.label = `${newParams.label || 'Multi-selection'}`;
    b.description = `${selections.length} source(s), ${totalSize} element(s) total`;

    return StateTransformer.UpdateResult.Updated;
  },
});

export { MultiStructureSelectionFromBundle };
type MultiStructureSelectionFromBundle = typeof MultiStructureSelectionFromBundle;
const MultiStructureSelectionFromBundle = PluginStateTransform.BuiltIn({
  name: 'structure-multi-selection-from-bundle',
  display: {
    name: 'Multi-structure Measurement Selection',
    description: 'Create selection object from multiple structures.',
  },
  from: SO.Root,
  to: SO.Molecule.Structure.Selections,
  params: () => ({
    selections: PD.ObjectList(
      {
        key: PD.Text(void 0, { description: 'A unique key.' }),
        ref: PD.Text(),
        groupId: PD.Optional(PD.Text()),
        bundle: PD.Value<StructureElement.Bundle>(StructureElement.Bundle.Empty),
      },
      (e) => e.ref,
      { isHidden: true },
    ),
    isTransitive: PD.Optional(
      PD.Boolean(false, {
        isHidden: true,
        description: 'Remap the selections from the original structure if structurally equivalent.',
      }),
    ),
    label: PD.Optional(PD.Text('', { isHidden: true })),
  }),
})({
  apply({ params, cache, dependencies }) {
    const entries = new Map<string, { source: Structure }>();

    const selections: SO.Molecule.Structure.SelectionEntry[] = [];
    let totalSize = 0;

    for (const sel of params.selections) {
      const source = dependencies![sel.ref].data as Structure;
      const loci = StructureElement.Bundle.toLoci(sel.bundle, source);
      selections.push({ key: sel.key, structureRef: sel.ref, loci, groupId: sel.groupId });
      totalSize += StructureElement.Loci.size(loci);
      entries.set(sel.key, { source });
    }

    (cache as object as any).entries = entries;

    const props = {
      label: `${params.label || 'Multi-selection'}`,
      description: `${params.selections.length} source(s), ${totalSize} element(s) total`,
    };
    return new SO.Molecule.Structure.Selections(selections, props);
  },
  update: ({ b, oldParams, newParams, cache, dependencies }) => {
    if (!!oldParams.isTransitive !== !!newParams.isTransitive) return StateTransformer.UpdateResult.Recreate;

    const cacheEntries = (cache as any).entries as Map<string, { source: Structure }>;
    const entries = new Map<string, { source: Structure }>();

    const prevBundles = new Map<string, StructureElement.Bundle>();
    for (const sel of oldParams.selections) {
      prevBundles.set(sel.key, sel.bundle);
    }

    const current = new Map<string, SO.Molecule.Structure.SelectionEntry>();
    for (const e of b.data) current.set(e.key, e);

    let changed = false;
    let totalSize = 0;

    const selections: SO.Molecule.Structure.SelectionEntry[] = [];
    for (const sel of newParams.selections) {
      let recreate = false;

      if (cacheEntries.has(sel.key)) {
        const source = dependencies![sel.ref].data as Structure;
        const entry = cacheEntries.get(sel.key)!;
        const prev = prevBundles.get(sel.key);
        if (
          prev &&
          source === entry.source &&
          sel.bundle.hash === entry.source.hashCode &&
          StructureElement.Bundle.areEqual(sel.bundle, prev)
        ) {
          const loci = current.get(sel.key)!;
          if (loci.groupId !== sel.groupId) {
            loci.groupId = sel.groupId;
            changed = true;
          }
          entries.set(sel.key, entry);
          selections.push(loci);
          totalSize += StructureElement.Loci.size(loci.loci);
          continue;
        }
        recreate = true;
      } else {
        recreate = true;
      }

      if (recreate) {
        changed = true;

        // create new selection
        const source = dependencies![sel.ref].data as Structure;
        const loci = StructureElement.Bundle.toLoci(sel.bundle, source);
        selections.push({ key: sel.key, structureRef: sel.ref, loci, groupId: sel.groupId });
        totalSize += StructureElement.Loci.size(loci);
        entries.set(sel.key, { source });
      }
    }

    if (!changed) return StateTransformer.UpdateResult.Unchanged;

    (cache as object as any).entries = entries;
    b.data = selections;
    b.label = `${newParams.label || 'Multi-selection'}`;
    b.description = `${selections.length} source(s), ${totalSize} element(s) total`;

    return StateTransformer.UpdateResult.Updated;
  },
});

export { StructureSelectionFromScript };
type StructureSelectionFromScript = typeof StructureSelectionFromScript;
const StructureSelectionFromScript = PluginStateTransform.BuiltIn({
  name: 'structure-selection-from-script',
  display: { name: 'Selection', description: 'Create a molecular structure from the specified script.' },
  from: SO.Molecule.Structure,
  to: SO.Molecule.Structure,
  params: () => ({
    script: PD.Script({
      language: 'mol-script',
      expression: '(sel.atom.atom-groups :residue-test (= atom.resname ALA))',
    }),
    label: PD.Optional(PD.Text('')),
  }),
})({
  apply({ a, params, cache }) {
    const { selection, entry } = StructureQueryHelper.createAndRun(a.data, params.script);
    (cache as any).entry = entry;

    const s = Sel.unionStructure(selection);
    const props = { label: `${params.label || 'Selection'}`, description: Structure.elementDescription(s) };
    return new SO.Molecule.Structure(s, props);
  },
  update: ({ a, b, oldParams, newParams, cache }) => {
    if (!Script.areEqual(oldParams.script, newParams.script)) {
      return StateTransformer.UpdateResult.Recreate;
    }

    const entry = (cache as { entry: StructureQueryHelper.CacheEntry }).entry;

    if (entry.currentStructure === a.data) {
      return StateTransformer.UpdateResult.Unchanged;
    }

    const selection = StructureQueryHelper.updateStructure(entry, a.data);
    StructureQueryHelper.updateStructureObject(b, selection, newParams.label);
    return StateTransformer.UpdateResult.Updated;
  },
  dispose({ b }) {
    b?.data.customPropertyDescriptors.dispose();
  },
});

export { StructureSelectionFromBundle };
type StructureSelectionFromBundle = typeof StructureSelectionFromBundle;
const StructureSelectionFromBundle = PluginStateTransform.BuiltIn({
  name: 'structure-selection-from-bundle',
  display: {
    name: 'Selection',
    description: 'Create a molecular structure from the specified structure-element bundle.',
  },
  from: SO.Molecule.Structure,
  to: SO.Molecule.Structure,
  params: () => ({
    bundle: PD.Value<StructureElement.Bundle>(StructureElement.Bundle.Empty, { isHidden: true }),
    label: PD.Optional(PD.Text('', { isHidden: true })),
  }),
})({
  apply({ a, params, cache }) {
    if (params.bundle.hash !== a.data.hashCode) {
      return StateObject.Null;
    }

    (cache as { source: Structure }).source = a.data;

    const s = StructureElement.Bundle.toStructure(params.bundle, a.data);
    if (s.elementCount === 0) return StateObject.Null;

    const props = { label: `${params.label || 'Selection'}`, description: Structure.elementDescription(s) };
    return new SO.Molecule.Structure(s, props);
  },
  update: ({ a, b, oldParams, newParams, cache }) => {
    if (!StructureElement.Bundle.areEqual(oldParams.bundle, newParams.bundle)) {
      return StateTransformer.UpdateResult.Recreate;
    }

    if (newParams.bundle.hash !== a.data.hashCode) {
      return StateTransformer.UpdateResult.Null;
    }

    if ((cache as { source: Structure }).source === a.data) {
      return StateTransformer.UpdateResult.Unchanged;
    }
    (cache as { source: Structure }).source = a.data;

    const s = StructureElement.Bundle.toStructure(newParams.bundle, a.data);
    if (s.elementCount === 0) return StateTransformer.UpdateResult.Null;

    b.label = `${newParams.label || 'Selection'}`;
    b.description = Structure.elementDescription(s);
    b.data = s;
    return StateTransformer.UpdateResult.Updated;
  },
  dispose({ b }) {
    b?.data.customPropertyDescriptors.dispose();
  },
});

export const StructureComplexElementTypes = {
  polymer: 'polymer',

  protein: 'protein',
  nucleic: 'nucleic',
  water: 'water',

  branched: 'branched', // = carbs
  ligand: 'ligand',
  'non-standard': 'non-standard',

  coarse: 'coarse',

  // Legacy
  'atomic-sequence': 'atomic-sequence',
  'atomic-het': 'atomic-het',
  spheres: 'spheres',
} as const;

export type StructureComplexElementTypes = keyof typeof StructureComplexElementTypes;

const StructureComplexElementTypeTuples = PD.objectToOptions(StructureComplexElementTypes);

export { StructureComplexElement };
type StructureComplexElement = typeof StructureComplexElement;
const StructureComplexElement = PluginStateTransform.BuiltIn({
  name: 'structure-complex-element',
  display: { name: 'Complex Element', description: 'Create a molecular structure from the specified model.' },
  from: SO.Molecule.Structure,
  to: SO.Molecule.Structure,
  params: {
    type: PD.Select<StructureComplexElementTypes>('atomic-sequence', StructureComplexElementTypeTuples, {
      isHidden: true,
    }),
  },
})({
  apply({ a, params }) {
    // TODO: update function.

    let query: StructureQuery, label: string;
    switch (params.type) {
      case 'polymer':
        query = polymer.query;
        label = 'Polymer';
        break;

      case 'protein':
        query = protein.query;
        label = 'Protein';
        break;
      case 'nucleic':
        query = nucleic.query;
        label = 'Nucleic';
        break;
      case 'water':
        query = Queries.internal.water();
        label = 'Water';
        break;

      case 'branched':
        query = branchedPlusConnected.query;
        label = 'Branched';
        break;
      case 'ligand':
        query = ligandPlusConnected.query;
        label = 'Ligand';
        break;

      case 'non-standard':
        query = nonStandardPolymer.query;
        label = 'Non-standard';
        break;

      case 'coarse':
        query = coarse.query;
        label = 'Coarse';
        break;

      case 'atomic-sequence':
        query = Queries.internal.atomicSequence();
        label = 'Sequence';
        break;
      case 'atomic-het':
        query = Queries.internal.atomicHet();
        label = 'HET Groups/Ligands';
        break;
      case 'spheres':
        query = Queries.internal.spheres();
        label = 'Coarse Spheres';
        break;

      default:
        assertUnreachable(params.type);
    }

    const result = query(new QueryContext(a.data));
    const s = Sel.unionStructure(result);

    if (s.elementCount === 0) return StateObject.Null;
    return new SO.Molecule.Structure(s, { label, description: Structure.elementDescription(s) });
  },
  dispose({ b }) {
    b?.data.customPropertyDescriptors.dispose();
  },
});

export { StructureComponent };
type StructureComponent = typeof StructureComponent;
const StructureComponent = PluginStateTransform.BuiltIn({
  name: 'structure-component',
  display: { name: 'Component', description: 'A molecular structure component.' },
  from: SO.Molecule.Structure,
  to: SO.Molecule.Structure,
  params: StructureComponentParams,
})({
  apply({ a, params, cache }) {
    return createStructureComponent(a.data, params, cache as any);
  },
  update: ({ a, b, oldParams, newParams, cache }) => {
    return updateStructureComponent(a.data, b, oldParams, newParams, cache as any);
  },
  dispose({ b }) {
    b?.data.customPropertyDescriptors.dispose();
  },
});
