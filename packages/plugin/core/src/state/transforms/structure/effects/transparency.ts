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
import { Transparency } from '@molstar/graphics/theme/transparency';
import { type Structure, StructureElement } from '@molstar/model/model/structure';
import { StateTransformer } from '@molstar/core/state';
import { hasColorSmoothingProp } from '@molstar/graphics/geo/geometry/base';
import { registerTransformLayer } from '@molstar/plugin/state/precompile';

export { TransparencyStructureRepresentation3DFromScript };
type TransparencyStructureRepresentation3DFromScript = typeof TransparencyStructureRepresentation3DFromScript;
const TransparencyStructureRepresentation3DFromScript = PluginStateTransform.BuiltIn({
  name: 'transparency-structure-representation-3d-from-script',
  display: 'Transparency 3D Representation',
  from: SO.Molecule.Structure.Representation3D,
  to: SO.Molecule.Structure.Representation3DState,
  params: () => ({
    layers: PD.ObjectList(
      {
        script: PD.Script(Script('(sel.atom.all)', 'mol-script')),
        value: PD.Numeric(0.5, { min: 0, max: 1, step: 0.01 }, { label: 'Transparency' }),
      },
      (e) => `Transparency (${e.value})`,
      {
        defaultValue: [
          {
            script: Script('(sel.atom.all)', 'mol-script'),
            value: 0.5,
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
    const geometryVersion = a.data.repr.geometryVersion;
    const transparency = Transparency.ofScript(params.layers, structure);

    return new SO.Molecule.Structure.Representation3DState(
      {
        state: { transparency },
        initialState: { transparency: Transparency.Empty },
        info: { structure, geometryVersion },
        repr: a.data.repr,
      },
      { label: `Transparency (${transparency.layers.length} Layers)` },
    );
  },
  update({ a, b, newParams, oldParams }) {
    const info = b.data.info as { structure: Structure; geometryVersion: number };
    const newStructure = a.data.sourceData;
    // a trajectory frame: the layers again on the new structure, in place
    const structureChanged = newStructure !== info.structure;
    if (a.data.repr !== b.data.repr) return StateTransformer.UpdateResult.Recreate;

    const newGeometryVersion = a.data.repr.geometryVersion;
    // smoothing needs to be re-calculated when geometry changes
    if (newGeometryVersion !== info.geometryVersion && hasColorSmoothingProp(a.data.repr.props))
      return StateTransformer.UpdateResult.Recreate;

    const oldTransparency = b.data.state.transparency!;
    const newTransparency = Transparency.ofScript(newParams.layers, newStructure);
    if (!structureChanged && Transparency.areEqual(oldTransparency, newTransparency))
      return StateTransformer.UpdateResult.Unchanged;

    info.geometryVersion = newGeometryVersion;
    info.structure = newStructure;
    b.data.state.transparency = newTransparency;
    b.data.repr = a.data.repr;
    b.label = `Transparency (${newTransparency.layers.length} Layers)`;
    return StateTransformer.UpdateResult.Updated;
  },
});

export { TransparencyStructureRepresentation3DFromBundle };
type TransparencyStructureRepresentation3DFromBundle = typeof TransparencyStructureRepresentation3DFromBundle;
const TransparencyStructureRepresentation3DFromBundle = PluginStateTransform.BuiltIn({
  name: 'transparency-structure-representation-3d-from-bundle',
  display: 'Transparency 3D Representation',
  from: SO.Molecule.Structure.Representation3D,
  to: SO.Molecule.Structure.Representation3DState,
  params: () => ({
    layers: PD.ObjectList(
      {
        bundle: PD.Value<StructureElement.Bundle>(StructureElement.Bundle.Empty),
        value: PD.Numeric(0.5, { min: 0, max: 1, step: 0.01 }, { label: 'Transparency' }),
      },
      (e) => `Transparency (${e.value})`,
      {
        defaultValue: [
          {
            bundle: StructureElement.Bundle.Empty,
            value: 0.5,
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
    const geometryVersion = a.data.repr.geometryVersion;
    const transparency = Transparency.ofBundle(params.layers, structure);

    return new SO.Molecule.Structure.Representation3DState(
      {
        state: { transparency },
        initialState: { transparency: Transparency.Empty },
        info: { structure, geometryVersion },
        repr: a.data.repr,
      },
      { label: `Transparency (${transparency.layers.length} Layers)` },
    );
  },
  update({ a, b, newParams, oldParams }) {
    const info = b.data.info as { structure: Structure; geometryVersion: number };
    const newStructure = a.data.sourceData;
    // a trajectory frame: the layers again on the new structure, in place
    const structureChanged = newStructure !== info.structure;
    if (a.data.repr !== b.data.repr) return StateTransformer.UpdateResult.Recreate;

    const newGeometryVersion = a.data.repr.geometryVersion;
    // smoothing needs to be re-calculated when geometry changes
    if (newGeometryVersion !== info.geometryVersion && hasColorSmoothingProp(a.data.repr.props))
      return StateTransformer.UpdateResult.Recreate;

    const oldTransparency = b.data.state.transparency!;
    // bundles and the structure it last ran on: the same result (a bundle reads topology only; a script may read a
    // custom property, which an ancestor's update computes again on the same structure), without the work
    if (!structureChanged && deepEqual(oldParams.layers, newParams.layers))
      return StateTransformer.UpdateResult.Unchanged;
    const newTransparency = Transparency.ofBundle(newParams.layers, newStructure);
    if (!structureChanged && Transparency.areEqual(oldTransparency, newTransparency))
      return StateTransformer.UpdateResult.Unchanged;

    info.geometryVersion = newGeometryVersion;
    info.structure = newStructure;
    b.data.state.transparency = newTransparency;
    b.data.repr = a.data.repr;
    b.label = `Transparency (${newTransparency.layers.length} Layers)`;
    return StateTransformer.UpdateResult.Updated;
  },
});

registerTransformLayer(TransparencyStructureRepresentation3DFromScript, 'transparency');
registerTransformLayer(TransparencyStructureRepresentation3DFromBundle, 'transparency');
