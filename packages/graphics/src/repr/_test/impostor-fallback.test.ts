/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import type { GraphicsContext } from '../../gl/graphics-context.js';
import {
  BridgeImpostorVisual,
  BridgeMeshVisual,
} from '../../props/computed/representations/interactions-bridge-cylinder.js';
import {
  InteractionsInterUnitImpostorVisual,
  InteractionsInterUnitMeshVisual,
} from '../../props/computed/representations/interactions-inter-unit-cylinder.js';
import {
  InteractionsIntraUnitImpostorVisual,
  InteractionsIntraUnitMeshVisual,
} from '../../props/computed/representations/interactions-intra-unit-cylinder.js';
import {
  SpacefillParticlesImpostorVisual,
  SpacefillParticlesMeshVisual,
} from '../particles/representation/spacefill.js';
import {
  InterUnitBondCylinderImpostorVisual,
  InterUnitBondCylinderMeshVisual,
} from '../structure/visual/bond-inter-unit-cylinder.js';
import {
  IntraUnitBondCylinderImpostorVisual,
  IntraUnitBondCylinderMeshVisual,
  StructureIntraUnitBondCylinderImpostorVisual,
  StructureIntraUnitBondCylinderMeshVisual,
} from '../structure/visual/bond-intra-unit-cylinder.js';
import {
  ElementSphereImpostorVisual,
  ElementSphereMeshVisual,
  StructureElementSphereImpostorVisual,
  StructureElementSphereMeshVisual,
} from '../structure/visual/element-sphere.js';
import {
  NucleotideAtomicBondImpostorVisual,
  NucleotideAtomicBondMeshVisual,
} from '../structure/visual/nucleotide-atomic-bond.js';
import {
  NucleotideAtomicElementImpostorVisual,
  NucleotideAtomicElementMeshVisual,
} from '../structure/visual/nucleotide-atomic-element.js';
import {
  PolymerBackboneCylinderImpostorVisual,
  PolymerBackboneCylinderMeshVisual,
} from '../structure/visual/polymer-backbone-cylinder.js';
import {
  PolymerBackboneSphereImpostorVisual,
  PolymerBackboneSphereMeshVisual,
} from '../structure/visual/polymer-backbone-sphere.js';
import { VolumeSphereImpostorVisual, VolumeSphereMeshVisual } from '../volume/dot.js';

type Recreatable = { mustRecreate?: (key: any, props: any, graphics?: GraphicsContext) => boolean };
const pairs: [string, (materialId: number) => Recreatable, (materialId: number) => Recreatable][] = [
  ['element sphere', ElementSphereImpostorVisual, ElementSphereMeshVisual],
  ['structure element sphere', StructureElementSphereImpostorVisual, StructureElementSphereMeshVisual],
  ['intra-unit bond', IntraUnitBondCylinderImpostorVisual, IntraUnitBondCylinderMeshVisual],
  ['structure intra-unit bond', StructureIntraUnitBondCylinderImpostorVisual, StructureIntraUnitBondCylinderMeshVisual],
  ['inter-unit bond', InterUnitBondCylinderImpostorVisual, InterUnitBondCylinderMeshVisual],
  ['nucleotide atomic element', NucleotideAtomicElementImpostorVisual, NucleotideAtomicElementMeshVisual],
  ['nucleotide atomic bond', NucleotideAtomicBondImpostorVisual, NucleotideAtomicBondMeshVisual],
  ['polymer backbone sphere', PolymerBackboneSphereImpostorVisual, PolymerBackboneSphereMeshVisual],
  ['polymer backbone cylinder', PolymerBackboneCylinderImpostorVisual, PolymerBackboneCylinderMeshVisual],
  ['interactions intra-unit', InteractionsIntraUnitImpostorVisual, InteractionsIntraUnitMeshVisual],
  ['interactions inter-unit', InteractionsInterUnitImpostorVisual, InteractionsInterUnitMeshVisual],
  ['interactions bridge', BridgeImpostorVisual, BridgeMeshVisual],
  ['particles spacefill', SpacefillParticlesImpostorVisual, SpacefillParticlesMeshVisual],
  ['volume dot', VolumeSphereImpostorVisual, VolumeSphereMeshVisual],
];

const graphics = (fragDepth: boolean) => ({ capabilities: { fragDepth } }) as unknown as GraphicsContext;

// A context without `fragDepth` gets the mesh variant of an impostor-capable visual: it must keep it, not recreate it on
// every update, and the impostor variant must give way to it (perf-audit-2 09/B5).
describe('impostor visuals without fragDepth', () => {
  for (const [name, impostor, mesh] of pairs) {
    it(name, () => {
      const props = { tryUseImpostor: true };
      expect(mesh(0).mustRecreate!({}, props, graphics(false))).toBe(false);
      expect(mesh(0).mustRecreate!({}, props, graphics(true))).toBe(true);
      expect(mesh(0).mustRecreate!({}, { tryUseImpostor: false }, graphics(true))).toBe(false);
      expect(impostor(0).mustRecreate!({}, props, graphics(false))).toBe(true);
      expect(impostor(0).mustRecreate!({}, props, graphics(true))).toBe(false);
      expect(impostor(0).mustRecreate!({}, props, undefined)).toBe(true);
    });
  }
});
