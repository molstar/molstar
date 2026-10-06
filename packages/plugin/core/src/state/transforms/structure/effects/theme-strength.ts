/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginStateTransform, PluginStateObject as SO } from '@molstar/plugin/state/objects';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { StateTransformer } from '@molstar/core/state';
import { lerp } from '@molstar/core/math/interpolate';

export { ThemeStrengthRepresentation3D };
type ThemeStrengthRepresentation3D = typeof ThemeStrengthRepresentation3D;
const ThemeStrengthRepresentation3D = PluginStateTransform.BuiltIn({
  name: 'theme-strength-representation-3d',
  display: 'Theme Strength 3D Representation',
  from: SO.Molecule.Structure.Representation3D,
  to: SO.Molecule.Structure.Representation3DState,
  params: () => ({
    overpaintStrength: PD.Numeric(1, { min: 0, max: 1, step: 0.01 }),
    transparencyStrength: PD.Numeric(1, { min: 0, max: 1, step: 0.01 }),
    emissiveStrength: PD.Numeric(1, { min: 0, max: 1, step: 0.01 }),
    substanceStrength: PD.Numeric(1, { min: 0, max: 1, step: 0.01 }),
    wiggleStrength: PD.Numeric(1, { min: 0, max: 1, step: 0.01 }),
  }),
})({
  canAutoUpdate() {
    return true;
  },
  apply({ a, params }) {
    return new SO.Molecule.Structure.Representation3DState(
      {
        state: {
          themeStrength: {
            overpaint: params.overpaintStrength,
            transparency: params.transparencyStrength,
            emissive: params.emissiveStrength,
            substance: params.substanceStrength,
            wiggle: params.wiggleStrength,
          },
        },
        initialState: {
          themeStrength: { overpaint: 1, transparency: 1, emissive: 1, substance: 1, wiggle: 1 },
        },
        info: {},
        repr: a.data.repr,
      },
      {
        label: 'Theme Strength',
        description: `${params.overpaintStrength.toFixed(2)}, ${params.transparencyStrength.toFixed(2)}, ${params.emissiveStrength.toFixed(2)}, ${params.substanceStrength.toFixed(2)}, ${params.wiggleStrength.toFixed(2)}`,
      },
    );
  },
  update({ a, b, newParams, oldParams }) {
    if (
      newParams.overpaintStrength === b.data.state.themeStrength?.overpaint &&
      newParams.transparencyStrength === b.data.state.themeStrength?.transparency &&
      newParams.emissiveStrength === b.data.state.themeStrength?.emissive &&
      newParams.substanceStrength === b.data.state.themeStrength?.substance &&
      newParams.wiggleStrength === b.data.state.themeStrength?.wiggle
    )
      return StateTransformer.UpdateResult.Unchanged;

    b.data.state.themeStrength = {
      overpaint: newParams.overpaintStrength,
      transparency: newParams.transparencyStrength,
      emissive: newParams.emissiveStrength,
      substance: newParams.substanceStrength,
      wiggle: newParams.wiggleStrength,
    };
    b.data.repr = a.data.repr;
    b.label = 'Theme Strength';
    b.description = `${newParams.overpaintStrength.toFixed(2)}, ${newParams.transparencyStrength.toFixed(2)}, ${newParams.emissiveStrength.toFixed(2)}, ${newParams.substanceStrength.toFixed(2)}, ${newParams.wiggleStrength.toFixed(2)}`;
    return StateTransformer.UpdateResult.Updated;
  },
  interpolate(src, tar, t) {
    return {
      overpaintStrength: lerp(src.overpaintStrength, tar.overpaintStrength, t),
      transparencyStrength: lerp(src.transparencyStrength, tar.transparencyStrength, t),
      emissiveStrength: lerp(src.emissiveStrength, tar.emissiveStrength, t),
      substanceStrength: lerp(src.substanceStrength, tar.substanceStrength, t),
      wiggleStrength: lerp(src.wiggleStrength, tar.wiggleStrength, t),
    };
  },
});
