/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginStateTransform, PluginStateObject as SO } from '@molstar/plugin/state/objects';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Script } from '@molstar/model/script/script';
import { Emissive } from '@molstar/graphics/theme/emissive';
import { type Structure, StructureElement } from '@molstar/model/model/structure';
import { StateTransformer } from '@molstar/core/state';
import { hasColorSmoothingProp } from '@molstar/graphics/geo/geometry/base';

export { EmissiveStructureRepresentation3DFromScript };
type EmissiveStructureRepresentation3DFromScript = typeof EmissiveStructureRepresentation3DFromScript;
const EmissiveStructureRepresentation3DFromScript = PluginStateTransform.BuiltIn({
  name: 'emissive-structure-representation-3d-from-script',
  display: 'Emissive 3D Representation',
  from: SO.Molecule.Structure.Representation3D,
  to: SO.Molecule.Structure.Representation3DState,
  params: () => ({
    layers: PD.ObjectList(
      {
        script: PD.Script(Script('(sel.atom.all)', 'mol-script')),
        value: PD.Numeric(0.5, { min: 0, max: 1, step: 0.01 }, { label: 'Emissive' }),
      },
      (e) => `Emissive (${e.value})`,
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
    const emissive = Emissive.ofScript(params.layers, structure);

    return new SO.Molecule.Structure.Representation3DState(
      {
        state: { emissive },
        initialState: { emissive: Emissive.Empty },
        info: { structure, geometryVersion },
        repr: a.data.repr,
      },
      { label: `Emissive (${emissive.layers.length} Layers)` },
    );
  },
  update({ a, b, newParams, oldParams }) {
    const info = b.data.info as { structure: Structure; geometryVersion: number };
    const newStructure = a.data.sourceData;
    if (newStructure !== info.structure) return StateTransformer.UpdateResult.Recreate;
    if (a.data.repr !== b.data.repr) return StateTransformer.UpdateResult.Recreate;

    const newGeometryVersion = a.data.repr.geometryVersion;
    // smoothing needs to be re-calculated when geometry changes
    if (newGeometryVersion !== info.geometryVersion && hasColorSmoothingProp(a.data.repr.props))
      return StateTransformer.UpdateResult.Recreate;

    const oldEmissive = b.data.state.emissive!;
    const newEmissive = Emissive.ofScript(newParams.layers, newStructure);
    if (Emissive.areEqual(oldEmissive, newEmissive)) return StateTransformer.UpdateResult.Unchanged;

    info.geometryVersion = newGeometryVersion;
    b.data.state.emissive = newEmissive;
    b.data.repr = a.data.repr;
    b.label = `Emissive (${newEmissive.layers.length} Layers)`;
    return StateTransformer.UpdateResult.Updated;
  },
});

export { EmissiveStructureRepresentation3DFromBundle };
type EmissiveStructureRepresentation3DFromBundle = typeof EmissiveStructureRepresentation3DFromBundle;
const EmissiveStructureRepresentation3DFromBundle = PluginStateTransform.BuiltIn({
  name: 'emissive-structure-representation-3d-from-bundle',
  display: 'Emissive 3D Representation',
  from: SO.Molecule.Structure.Representation3D,
  to: SO.Molecule.Structure.Representation3DState,
  params: () => ({
    layers: PD.ObjectList(
      {
        bundle: PD.Value<StructureElement.Bundle>(StructureElement.Bundle.Empty),
        value: PD.Numeric(0.5, { min: 0, max: 1, step: 0.01 }, { label: 'Emissive' }),
      },
      (e) => `Emissive (${e.value})`,
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
    const emissive = Emissive.ofBundle(params.layers, structure);

    return new SO.Molecule.Structure.Representation3DState(
      {
        state: { emissive },
        initialState: { emissive: Emissive.Empty },
        info: { structure, geometryVersion },
        repr: a.data.repr,
      },
      { label: `Emissive (${emissive.layers.length} Layers)` },
    );
  },
  update({ a, b, newParams, oldParams }) {
    const info = b.data.info as { structure: Structure; geometryVersion: number };
    const newStructure = a.data.sourceData;
    if (newStructure !== info.structure) return StateTransformer.UpdateResult.Recreate;
    if (a.data.repr !== b.data.repr) return StateTransformer.UpdateResult.Recreate;

    const newGeometryVersion = a.data.repr.geometryVersion;
    // smoothing needs to be re-calculated when geometry changes
    if (newGeometryVersion !== info.geometryVersion && hasColorSmoothingProp(a.data.repr.props))
      return StateTransformer.UpdateResult.Recreate;

    const oldEmissive = b.data.state.emissive!;
    const newEmissive = Emissive.ofBundle(newParams.layers, newStructure);
    if (Emissive.areEqual(oldEmissive, newEmissive)) return StateTransformer.UpdateResult.Unchanged;

    info.geometryVersion = newGeometryVersion;
    b.data.state.emissive = newEmissive;
    b.data.repr = a.data.repr;
    b.label = `Emissive (${newEmissive.layers.length} Layers)`;
    return StateTransformer.UpdateResult.Updated;
  },
});
