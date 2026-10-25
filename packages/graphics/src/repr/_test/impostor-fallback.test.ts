/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import type { GraphicsContext } from '../../gl/graphics-context.js';
import { BridgeVisual } from '../../props/computed/representations/interactions-bridge-cylinder.js';
import { InteractionsInterUnitVisual } from '../../props/computed/representations/interactions-inter-unit-cylinder.js';
import { InteractionsIntraUnitVisual } from '../../props/computed/representations/interactions-intra-unit-cylinder.js';
import { SpacefillParticlesVisual } from '../particles/representation/spacefill.js';
import { InterUnitBondCylinderVisual } from '../structure/visual/bond-inter-unit-cylinder.js';
import {
  IntraUnitBondCylinderVisual,
  StructureIntraUnitBondCylinderVisual,
} from '../structure/visual/bond-intra-unit-cylinder.js';
import { ElementSphereVisual, StructureElementSphereVisual } from '../structure/visual/element-sphere.js';
import { NucleotideAtomicBondVisual } from '../structure/visual/nucleotide-atomic-bond.js';
import { NucleotideAtomicElementVisual } from '../structure/visual/nucleotide-atomic-element.js';
import { PolymerBackboneCylinderVisual } from '../structure/visual/polymer-backbone-cylinder.js';
import { PolymerBackboneSphereVisual } from '../structure/visual/polymer-backbone-sphere.js';
import type { SelectingVisualFactory } from '../visual-factory.js';
import { VolumeSphereVisual } from '../volume/dot.js';

const factories: [string, SelectingVisualFactory<any, any>][] = [
  ['element sphere', ElementSphereVisual],
  ['structure element sphere', StructureElementSphereVisual],
  ['intra-unit bond', IntraUnitBondCylinderVisual],
  ['structure intra-unit bond', StructureIntraUnitBondCylinderVisual],
  ['inter-unit bond', InterUnitBondCylinderVisual],
  ['nucleotide atomic element', NucleotideAtomicElementVisual],
  ['nucleotide atomic bond', NucleotideAtomicBondVisual],
  ['polymer backbone sphere', PolymerBackboneSphereVisual],
  ['polymer backbone cylinder', PolymerBackboneCylinderVisual],
  ['interactions intra-unit', InteractionsIntraUnitVisual],
  ['interactions inter-unit', InteractionsInterUnitVisual],
  ['interactions bridge', BridgeVisual],
  ['particles spacefill', SpacefillParticlesVisual],
  ['volume dot', VolumeSphereVisual],
];

const graphics = (fragDepth: boolean) => ({ capabilities: { fragDepth } }) as unknown as GraphicsContext;

// A context without `fragDepth` gets the mesh variant of an impostor-capable visual: it must keep it, not recreate it on
// every update, and the impostor variant must give way to it (perf-audit-2 09/B5). The variant is the factory's choice
// (`chooseVisual`), and so is when it is recreated.
describe('impostor visuals without fragDepth', () => {
  for (const [name, factory] of factories) {
    it(name, () => {
      const props = { tryUseImpostor: true };
      const mesh = factory(0, {}, props, graphics(false));
      const impostor = factory(0, {}, props, graphics(true));
      expect(mesh.mustRecreate!({}, props, graphics(false))).toBe(false);
      expect(mesh.mustRecreate!({}, props, graphics(true))).toBe(true);
      expect(mesh.mustRecreate!({}, { tryUseImpostor: false }, graphics(true))).toBe(false);
      expect(impostor.mustRecreate!({}, props, graphics(false))).toBe(true);
      expect(impostor.mustRecreate!({}, props, graphics(true))).toBe(false);
      expect(impostor.mustRecreate!({}, props, undefined)).toBe(true);
    });
  }
});
