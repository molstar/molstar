/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { deepEqual } from '@molstar/core/util';
import { PluginStateTransform, PluginStateObject as SO } from '@molstar/plugin/state/objects';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Script } from '@molstar/model/script/script';
import { Clipping } from '@molstar/graphics/theme/clipping';
import { ObjectKeys } from '@molstar/core/util/type-helpers';
import { StructureElement } from '@molstar/model/model/structure';
import { StateTransformer } from '@molstar/core/state';
import { registerTransformLayer } from '@molstar/plugin/state/precompile';

export { ClippingStructureRepresentation3DFromScript };
type ClippingStructureRepresentation3DFromScript = typeof ClippingStructureRepresentation3DFromScript;
const ClippingStructureRepresentation3DFromScript = PluginStateTransform.BuiltIn({
  name: 'clipping-structure-representation-3d-from-script',
  display: 'Clipping 3D Representation',
  from: SO.Molecule.Structure.Representation3D,
  to: SO.Molecule.Structure.Representation3DState,
  params: () => ({
    layers: PD.ObjectList(
      {
        script: PD.Script(Script('(sel.atom.all)', 'mol-script')),
        groups: PD.Converted(
          (g: Clipping.Groups) => Clipping.Groups.toNames(g),
          (n) => Clipping.Groups.fromNames(n),
          PD.MultiSelect(ObjectKeys(Clipping.Groups.Names), PD.objectToOptions(Clipping.Groups.Names)),
        ),
      },
      (e) => `${Clipping.Groups.toNames(e.groups).length} group(s)`,
      {
        defaultValue: [
          {
            script: Script('(sel.atom.all)', 'mol-script'),
            groups: Clipping.Groups.Flag.None,
          },
        ],
      },
    ),
  }),
})({
  canAutoUpdate() {
    return true;
  },
  apply({ a, params }) {
    const structure = a.data.sourceData;
    const clipping = Clipping.ofScript(params.layers, structure);

    return new SO.Molecule.Structure.Representation3DState(
      {
        state: { clipping },
        initialState: { clipping: Clipping.Empty },
        info: structure,
        repr: a.data.repr,
      },
      { label: `Clipping (${clipping.layers.length} Layers)` },
    );
  },
  update({ a, b, newParams, oldParams }) {
    // a trajectory frame: the layers again on the new structure, in place (T3)
    const structure = a.data.sourceData;
    const structureChanged = structure !== b.data.info;
    if (a.data.repr !== b.data.repr) return StateTransformer.UpdateResult.Recreate;

    const oldClipping = b.data.state.clipping!;
    const newClipping = Clipping.ofScript(newParams.layers, structure);
    if (!structureChanged && Clipping.areEqual(oldClipping, newClipping))
      return StateTransformer.UpdateResult.Unchanged;

    b.data.state.clipping = newClipping;
    b.data.info = structure;
    b.data.repr = a.data.repr;
    b.label = `Clipping (${newClipping.layers.length} Layers)`;
    return StateTransformer.UpdateResult.Updated;
  },
});

export { ClippingStructureRepresentation3DFromBundle };
type ClippingStructureRepresentation3DFromBundle = typeof ClippingStructureRepresentation3DFromBundle;
const ClippingStructureRepresentation3DFromBundle = PluginStateTransform.BuiltIn({
  name: 'clipping-structure-representation-3d-from-bundle',
  display: 'Clipping 3D Representation',
  from: SO.Molecule.Structure.Representation3D,
  to: SO.Molecule.Structure.Representation3DState,
  params: () => ({
    layers: PD.ObjectList(
      {
        bundle: PD.Value<StructureElement.Bundle>(StructureElement.Bundle.Empty),
        groups: PD.Converted(
          (g: Clipping.Groups) => Clipping.Groups.toNames(g),
          (n) => Clipping.Groups.fromNames(n),
          PD.MultiSelect(ObjectKeys(Clipping.Groups.Names), PD.objectToOptions(Clipping.Groups.Names)),
        ),
      },
      (e) => `${Clipping.Groups.toNames(e.groups).length} group(s)`,
      {
        defaultValue: [
          {
            bundle: StructureElement.Bundle.Empty,
            groups: Clipping.Groups.Flag.None,
          },
        ],
        isHidden: true,
      },
    ),
  }),
})({
  canAutoUpdate() {
    return true;
  },
  apply({ a, params }) {
    const structure = a.data.sourceData;
    const clipping = Clipping.ofBundle(params.layers, structure);

    return new SO.Molecule.Structure.Representation3DState(
      {
        state: { clipping },
        initialState: { clipping: Clipping.Empty },
        info: structure,
        repr: a.data.repr,
      },
      { label: `Clipping (${clipping.layers.length} Layers)` },
    );
  },
  update({ a, b, newParams, oldParams }) {
    // a trajectory frame: the layers again on the new structure, in place (T3)
    const structure = a.data.sourceData;
    const structureChanged = structure !== b.data.info;
    if (a.data.repr !== b.data.repr) return StateTransformer.UpdateResult.Recreate;

    const oldClipping = b.data.state.clipping!;
    // bundles and the structure it last ran on: the same result (a bundle reads topology only; a script may read a
    // custom property, which an ancestor's update computes again on the same structure), without the work (S7)
    if (!structureChanged && deepEqual(oldParams.layers, newParams.layers))
      return StateTransformer.UpdateResult.Unchanged;
    const newClipping = Clipping.ofBundle(newParams.layers, structure);
    if (!structureChanged && Clipping.areEqual(oldClipping, newClipping))
      return StateTransformer.UpdateResult.Unchanged;

    b.data.state.clipping = newClipping;
    b.data.info = structure;
    b.data.repr = a.data.repr;
    b.label = `Clipping (${newClipping.layers.length} Layers)`;
    return StateTransformer.UpdateResult.Updated;
  },
});

registerTransformLayer(ClippingStructureRepresentation3DFromScript, 'clipping');
registerTransformLayer(ClippingStructureRepresentation3DFromBundle, 'clipping');
