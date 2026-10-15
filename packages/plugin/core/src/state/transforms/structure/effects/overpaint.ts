/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginStateTransform, PluginStateObject as SO } from '@molstar/plugin/state/objects';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Script } from '@molstar/model/script/script';
import { ColorNames } from '@molstar/core/util/color/names';
import { Color } from '@molstar/core/util/color';
import { Overpaint } from '@molstar/graphics/theme/overpaint';
import { type Structure, StructureElement } from '@molstar/model/model/structure';
import { StateTransformer } from '@molstar/core/state';
import { hasColorSmoothingProp } from '@molstar/graphics/geo/geometry/base';

export { OverpaintStructureRepresentation3DFromScript };
type OverpaintStructureRepresentation3DFromScript = typeof OverpaintStructureRepresentation3DFromScript;
const OverpaintStructureRepresentation3DFromScript = PluginStateTransform.BuiltIn({
  name: 'overpaint-structure-representation-3d-from-script',
  display: 'Overpaint 3D Representation',
  from: SO.Molecule.Structure.Representation3D,
  to: SO.Molecule.Structure.Representation3DState,
  params: () => ({
    layers: PD.ObjectList(
      {
        script: PD.Script(Script('(sel.atom.all)', 'mol-script')),
        color: PD.Color(ColorNames.blueviolet),
        clear: PD.Boolean(false),
      },
      (e) => `${e.clear ? 'Clear' : Color.toRgbString(e.color)}`,
      {
        defaultValue: [
          {
            script: Script('(sel.atom.all)', 'mol-script'),
            color: ColorNames.blueviolet,
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
    const overpaint = Overpaint.ofScript(params.layers, structure);

    return new SO.Molecule.Structure.Representation3DState(
      {
        state: { overpaint },
        initialState: { overpaint: Overpaint.Empty },
        info: { structure, geometryVersion },
        repr: a.data.repr,
      },
      { label: `Overpaint (${overpaint.layers.length} Layers)` },
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

    const oldOverpaint = b.data.state.overpaint!;
    const newOverpaint = Overpaint.ofScript(newParams.layers, newStructure);
    if (Overpaint.areEqual(oldOverpaint, newOverpaint)) return StateTransformer.UpdateResult.Unchanged;

    info.geometryVersion = newGeometryVersion;
    b.data.state.overpaint = newOverpaint;
    b.data.repr = a.data.repr;
    b.label = `Overpaint (${newOverpaint.layers.length} Layers)`;
    return StateTransformer.UpdateResult.Updated;
  },
});

export { OverpaintStructureRepresentation3DFromBundle };
type OverpaintStructureRepresentation3DFromBundle = typeof OverpaintStructureRepresentation3DFromBundle;
const OverpaintStructureRepresentation3DFromBundle = PluginStateTransform.BuiltIn({
  name: 'overpaint-structure-representation-3d-from-bundle',
  display: 'Overpaint 3D Representation',
  from: SO.Molecule.Structure.Representation3D,
  to: SO.Molecule.Structure.Representation3DState,
  params: () => ({
    layers: PD.ObjectList(
      {
        bundle: PD.Value<StructureElement.Bundle>(StructureElement.Bundle.Empty),
        color: PD.Color(ColorNames.blueviolet),
        clear: PD.Boolean(false),
      },
      (e) => `${e.clear ? 'Clear' : Color.toRgbString(e.color)}`,
      {
        defaultValue: [
          {
            bundle: StructureElement.Bundle.Empty,
            color: ColorNames.blueviolet,
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
    const overpaint = Overpaint.ofBundle(params.layers, structure);

    return new SO.Molecule.Structure.Representation3DState(
      {
        state: { overpaint },
        initialState: { overpaint: Overpaint.Empty },
        info: { structure, geometryVersion },
        repr: a.data.repr,
      },
      { label: `Overpaint (${overpaint.layers.length} Layers)` },
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

    const oldOverpaint = b.data.state.overpaint!;
    const newOverpaint = Overpaint.ofBundle(newParams.layers, newStructure);
    if (Overpaint.areEqual(oldOverpaint, newOverpaint)) return StateTransformer.UpdateResult.Unchanged;

    info.geometryVersion = newGeometryVersion;
    b.data.state.overpaint = newOverpaint;
    b.data.repr = a.data.repr;
    b.label = `Overpaint (${newOverpaint.layers.length} Layers)`;
    return StateTransformer.UpdateResult.Updated;
  },
});
