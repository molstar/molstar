/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginStateTransform, PluginStateObject as SO } from '@molstar/plugin/state/objects';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Script } from '@molstar/model/script/script';
import { Wiggle } from '@molstar/graphics/theme/wiggle';
import { type Structure, StructureElement } from '@molstar/model/model/structure';
import { StateTransformer } from '@molstar/core/state';

export { WiggleStructureRepresentation3DFromScript };
type WiggleStructureRepresentation3DFromScript = typeof WiggleStructureRepresentation3DFromScript;
const WiggleStructureRepresentation3DFromScript = PluginStateTransform.BuiltIn({
  name: 'wiggle-structure-representation-3d-from-script',
  display: 'Wiggle 3D Representation',
  from: SO.Molecule.Structure.Representation3D,
  to: SO.Molecule.Structure.Representation3DState,
  params: () => ({
    layers: PD.ObjectList(
      {
        script: PD.Script(Script('(sel.atom.all)', 'mol-script')),
        value: PD.Numeric(5, { min: 0, max: 1, step: 0.01 }, { label: 'Wiggle' }),
      },
      (e) => `Wiggle (${e.value})`,
      {
        defaultValue: [
          {
            script: Script('(sel.atom.all)', 'mol-script'),
            value: 1,
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
    const wiggle = Wiggle.ofScript(params.layers, structure);

    return new SO.Molecule.Structure.Representation3DState(
      {
        state: { wiggle },
        initialState: { wiggle: Wiggle.Empty },
        info: { structure, geometryVersion },
        repr: a.data.repr,
      },
      { label: `Wiggle (${wiggle.layers.length} Layers)` },
    );
  },
  update({ a, b, newParams, oldParams }) {
    const info = b.data.info as { structure: Structure; geometryVersion: number };
    const newStructure = a.data.sourceData;
    if (newStructure !== info.structure) return StateTransformer.UpdateResult.Recreate;
    if (a.data.repr !== b.data.repr) return StateTransformer.UpdateResult.Recreate;

    const oldWiggle = b.data.state.wiggle!;
    const newWiggle = Wiggle.ofScript(newParams.layers, newStructure);
    if (Wiggle.areEqual(oldWiggle, newWiggle)) return StateTransformer.UpdateResult.Unchanged;

    info.geometryVersion = a.data.repr.geometryVersion;
    b.data.state.wiggle = newWiggle;
    b.data.repr = a.data.repr;
    b.label = `Wiggle (${newWiggle.layers.length} Layers)`;
    return StateTransformer.UpdateResult.Updated;
  },
});

export { WiggleStructureRepresentation3DFromBundle };
type WiggleStructureRepresentation3DFromBundle = typeof WiggleStructureRepresentation3DFromBundle;
const WiggleStructureRepresentation3DFromBundle = PluginStateTransform.BuiltIn({
  name: 'wiggle-structure-representation-3d-from-bundle',
  display: 'Wiggle 3D Representation',
  from: SO.Molecule.Structure.Representation3D,
  to: SO.Molecule.Structure.Representation3DState,
  params: () => ({
    layers: PD.ObjectList(
      {
        bundle: PD.Value<StructureElement.Bundle>(StructureElement.Bundle.Empty),
        value: PD.Numeric(5, { min: 0, max: 1, step: 0.01 }, { label: 'Wiggle' }),
      },
      (e) => `Wiggle (${e.value})`,
      {
        defaultValue: [
          {
            bundle: StructureElement.Bundle.Empty,
            value: 1,
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
    const wiggle = Wiggle.ofBundle(params.layers, structure);

    return new SO.Molecule.Structure.Representation3DState(
      {
        state: { wiggle },
        initialState: { wiggle: Wiggle.Empty },
        info: { structure, geometryVersion },
        repr: a.data.repr,
      },
      { label: `Wiggle (${wiggle.layers.length} Layers)` },
    );
  },
  update({ a, b, newParams, oldParams }) {
    const info = b.data.info as { structure: Structure; geometryVersion: number };
    const newStructure = a.data.sourceData;
    if (newStructure !== info.structure) return StateTransformer.UpdateResult.Recreate;
    if (a.data.repr !== b.data.repr) return StateTransformer.UpdateResult.Recreate;

    const oldWiggle = b.data.state.wiggle!;
    const newWiggle = Wiggle.ofBundle(newParams.layers, newStructure);
    if (Wiggle.areEqual(oldWiggle, newWiggle)) return StateTransformer.UpdateResult.Unchanged;

    info.geometryVersion = a.data.repr.geometryVersion;
    b.data.state.wiggle = newWiggle;
    b.data.repr = a.data.repr;
    b.label = `Wiggle (${newWiggle.layers.length} Layers)`;
    return StateTransformer.UpdateResult.Updated;
  },
});
