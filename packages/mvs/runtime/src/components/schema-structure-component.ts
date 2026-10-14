/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Structure } from '@molstar/model/model/structure';
import { PluginStateObject as SO } from '@molstar/plugin/state/objects';
import { StateObject, StateTransformer } from '@molstar/core/state';
import { deepEqual } from '@molstar/core/util';
import { omitObjectKeys } from '@molstar/core/util/object';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { substructureFromRows } from '@molstar/model/model/structure/structure/element/schema-selection/row-structure';
import { describeRows } from '@molstar/mvs/helpers/row-label';
import type { MVSAnnotationRow } from '@molstar/mvs/helpers/schemas';
import { MVSTransform } from './annotation-structure-component.js';

/** Parameter definition for `MVSSchemaStructureComponent` transformer */
export const MVSSchemaStructureComponentParams = {
  schema: PD.Value<MVSAnnotationRow[]>([], {
    isHidden: true,
    description: 'Selection criteria; an element belongs to the component if it satisfies at least one of them',
  }),
  nullIfEmpty: PD.Optional(PD.Boolean(true, { isHidden: false })),
  label: PD.Text('', { isHidden: false }),
};

/** Parameter values for `MVSSchemaStructureComponent` transformer */
export type MVSSchemaStructureComponentProps = PD.ValuesFor<typeof MVSSchemaStructureComponentParams>;

interface MVSSchemaComponentCache {
  source?: Structure;
  /** The component is the source structure itself (everything is selected), so it must not be disposed with the component */
  isSource?: boolean;
}

/** Transformer for creating a structure component from selection criteria ("schema" rows, e.g. `{ label_asym_id: 'A', beg_label_seq_id: 10 }`),
 * resolved through the indices of the model (not by running a query on every atom).
 * Unlike the standard `StructureComponent` it does not take a MolScript expression, which cannot be evaluated faster than by testing each atom. */
export type MVSSchemaStructureComponent = typeof MVSSchemaStructureComponent;
export const MVSSchemaStructureComponent = MVSTransform({
  name: 'mvs-structure-component-from-schema',
  display: {
    name: 'MVS Schema Component',
    description: 'A molecular structure component defined by selection criteria (rows of label_*, auth_*, ... fields).',
  },
  from: SO.Molecule.Structure,
  to: SO.Molecule.Structure,
  params: MVSSchemaStructureComponentParams,
})({
  apply({ a, params, cache }) {
    return createMVSSchemaStructureComponent(a.data, params, cache as MVSSchemaComponentCache);
  },
  update: ({ a, b, oldParams, newParams, cache }) => {
    return updateMVSSchemaStructureComponent(a.data, b, oldParams, newParams, cache as MVSSchemaComponentCache);
  },
  dispose({ b, cache }) {
    // (when everything is selected the component is the source structure, which has other users)
    if (!(cache as MVSSchemaComponentCache | undefined)?.isSource) b?.data.customPropertyDescriptors.dispose();
  },
});

/** Create a substructure PSO based on `MVSSchemaStructureComponentProps` */
export function createMVSSchemaStructureComponent(
  structure: Structure,
  params: MVSSchemaStructureComponentProps,
  cache: MVSSchemaComponentCache,
) {
  cache.source = structure;
  const component = substructureFromRows(structure, params.schema);
  cache.isSource = component === structure;
  if (params.nullIfEmpty && component.elementCount === 0) return StateObject.Null;
  const label = params.label || describeRows(params.schema);
  return new SO.Molecule.Structure(component, { label, description: Structure.elementDescription(component) });
}

/** Update a substructure PSO based on `MVSSchemaStructureComponentProps` */
export function updateMVSSchemaStructureComponent(
  a: Structure,
  b: SO.Molecule.Structure,
  oldParams: MVSSchemaStructureComponentProps,
  newParams: MVSSchemaStructureComponentProps,
  cache: MVSSchemaComponentCache,
) {
  const structureChanged = !cache.source || !Structure.areEquivalent(a, cache.source);
  cache.source = a;
  if (structureChanged) return StateTransformer.UpdateResult.Recreate;

  const coreParamsChanged = !deepEqual(omitObjectKeys(newParams, ['label']), omitObjectKeys(oldParams, ['label']));
  if (coreParamsChanged) return StateTransformer.UpdateResult.Recreate;

  if (newParams.label !== oldParams.label) {
    b.label = newParams.label || b.label;
    return StateTransformer.UpdateResult.Updated;
  }
  return StateTransformer.UpdateResult.Unchanged;
}
