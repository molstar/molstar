/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginStateTransform, PluginStateObject as SO } from '@molstar/plugin/state/objects';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Script } from '@molstar/model/script/script';
import { Material } from '@molstar/core/util/material';
import { Substance } from '@molstar/graphics/theme/substance';
import { type Structure, StructureElement } from '@molstar/model/model/structure';
import { StateTransformer } from '@molstar/core/state';
import { hasColorSmoothingProp } from '@molstar/graphics/geo/geometry/base';

export { SubstanceStructureRepresentation3DFromScript };
type SubstanceStructureRepresentation3DFromScript = typeof SubstanceStructureRepresentation3DFromScript;
const SubstanceStructureRepresentation3DFromScript = PluginStateTransform.BuiltIn({
  name: 'substance-structure-representation-3d-from-script',
  display: 'Substance 3D Representation',
  from: SO.Molecule.Structure.Representation3D,
  to: SO.Molecule.Structure.Representation3DState,
  params: () => ({
    layers: PD.ObjectList(
      {
        script: PD.Script(Script('(sel.atom.all)', 'mol-script')),
        material: Material.getParam(),
        clear: PD.Boolean(false),
      },
      (e) => `${e.clear ? 'Clear' : Material.toString(e.material)}`,
      {
        defaultValue: [
          {
            script: Script('(sel.atom.all)', 'mol-script'),
            material: Material({ roughness: 1 }),
            clear: false,
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
    const substance = Substance.ofScript(params.layers, structure);

    return new SO.Molecule.Structure.Representation3DState(
      {
        state: { substance },
        initialState: { substance: Substance.Empty },
        info: { structure, geometryVersion },
        repr: a.data.repr,
      },
      { label: `Substance (${substance.layers.length} Layers)` },
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

    const oldSubstance = b.data.state.substance!;
    const newSubstance = Substance.ofScript(newParams.layers, newStructure);
    if (Substance.areEqual(oldSubstance, newSubstance)) return StateTransformer.UpdateResult.Unchanged;

    info.geometryVersion = newGeometryVersion;
    b.data.state.substance = newSubstance;
    b.data.repr = a.data.repr;
    b.label = `Substance (${newSubstance.layers.length} Layers)`;
    return StateTransformer.UpdateResult.Updated;
  },
});

export { SubstanceStructureRepresentation3DFromBundle };
type SubstanceStructureRepresentation3DFromBundle = typeof SubstanceStructureRepresentation3DFromBundle;
const SubstanceStructureRepresentation3DFromBundle = PluginStateTransform.BuiltIn({
  name: 'substance-structure-representation-3d-from-bundle',
  display: 'Substance 3D Representation',
  from: SO.Molecule.Structure.Representation3D,
  to: SO.Molecule.Structure.Representation3DState,
  params: () => ({
    layers: PD.ObjectList(
      {
        bundle: PD.Value<StructureElement.Bundle>(StructureElement.Bundle.Empty),
        material: Material.getParam(),
        clear: PD.Boolean(false),
      },
      (e) => `${e.clear ? 'Clear' : Material.toString(e.material)}`,
      {
        defaultValue: [
          {
            bundle: StructureElement.Bundle.Empty,
            material: Material({ roughness: 1 }),
            clear: false,
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
    const substance = Substance.ofBundle(params.layers, structure);

    return new SO.Molecule.Structure.Representation3DState(
      {
        state: { substance },
        initialState: { substance: Substance.Empty },
        info: { structure, geometryVersion },
        repr: a.data.repr,
      },
      { label: `Substance (${substance.layers.length} Layers)` },
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

    const oldSubstance = b.data.state.substance!;
    const newSubstance = Substance.ofBundle(newParams.layers, newStructure);
    if (Substance.areEqual(oldSubstance, newSubstance)) return StateTransformer.UpdateResult.Unchanged;

    info.geometryVersion = newGeometryVersion;
    b.data.state.substance = newSubstance;
    b.data.repr = a.data.repr;
    b.label = `Substance (${newSubstance.layers.length} Layers)`;
    return StateTransformer.UpdateResult.Updated;
  },
});
