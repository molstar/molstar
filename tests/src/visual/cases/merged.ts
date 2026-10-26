/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { BackboneRepresentationProvider } from '@molstar/graphics/repr/structure/representation/backbone';
import { GaussianSurfaceRepresentationProvider } from '@molstar/graphics/repr/structure/representation/gaussian-surface';
import { LineRepresentationProvider } from '@molstar/graphics/repr/structure/representation/line';
import { PointRepresentationProvider } from '@molstar/graphics/repr/structure/representation/point';
import { SpacefillRepresentationProvider } from '@molstar/graphics/repr/structure/representation/spacefill';
import type { RepresentationProvider } from '@molstar/graphics/repr/representation';
import type { Structure } from '@molstar/model/model/structure';
import { type VisualTestCase } from '../framework/test-case.js';
import { atomStructure, structureRepr } from '../scenes/structure.js';
import { Impostor, ObliqueView } from './geometry.js';

/**
 * `count` short chains (a unit each) on a ring, one carbon atom per residue on a zigzag; 6 to 9 atoms per chain, so
 * that the group-instance data of the segments are offset (`aSegment.y` is zero for equally sized segments of one
 * instance each)
 */
async function chainRing(count: number): Promise<Structure> {
  const lengths = Array.from({ length: count }, (_, c) => 6 + (c % 4));
  const total = lengths.reduce((a, b) => a + b, 0);
  const x = new Float32Array(total);
  const y = new Float32Array(total);
  const z = new Float32Array(total);
  const chainOf = new Int32Array(total);
  let k = 0;
  for (let c = 0; c < count; ++c) {
    const a = (c / count) * 2 * Math.PI;
    for (let i = 0; i < lengths[c]; ++i, ++k) {
      const r = 10 + i * 1.5;
      x[k] = r * Math.cos(a);
      y[k] = r * Math.sin(a);
      z[k] = (i % 2) * 1.2;
      chainOf[k] = c;
    }
  }
  const chains = Array.from({ length: count }, (_, c) => String.fromCharCode(65 + c));
  return atomStructure('chain-ring', x, y, z, chains, chainOf);
}

const Units = 12;

/**
 * The units of a representation merged into render items per geometry type (`renderMerged`), colored per atom: by
 * `chain-id` (group data, offset per segment), or by `element-index` (group-instance data); compared with WebGL
 */
function mergedCase(
  name: string,
  provider: RepresentationProvider<Structure, any, any>,
  props: { [k: string]: unknown },
  options: { requires?: VisualTestCase['requires']; colorTheme?: 'chain-id' | 'element-index' } = {},
): VisualTestCase {
  const colorTheme = options.colorTheme ?? 'chain-id';
  return {
    id: `geometry/merged/${name}`,
    description: `${Units} units as ${name}, merged into render items (renderMerged), colored by ${colorTheme}.`,
    requires: options.requires,
    view: ObliqueView,
    scene: async ({ graphics }) => {
      const structure = await chainRing(Units);
      const repr = await structureRepr(graphics, provider, structure, { ...props, renderMerged: true } as any, {
        name: colorTheme,
      });
      // the case is about merged items (a representation can have several visuals per geometry type)
      if (repr.renderObjects.length >= Units) {
        throw new Error(`${repr.renderObjects.length} render objects of ${Units} units, the units were not merged`);
      }
      return [repr];
    },
  };
}

export const MergedCases: VisualTestCase[] = [
  mergedCase('spacefill', SpacefillRepresentationProvider, {}, { requires: Impostor }),
  mergedCase(
    'spacefill-group-instance',
    SpacefillRepresentationProvider,
    {},
    {
      requires: Impostor,
      colorTheme: 'element-index',
    },
  ),
  mergedCase('spacefill-mesh', SpacefillRepresentationProvider, { tryUseImpostor: false }),
  mergedCase('backbone', BackboneRepresentationProvider, {}, { requires: Impostor }),
  mergedCase('point', PointRepresentationProvider, {}),
  mergedCase('line', LineRepresentationProvider, {}),
  {
    // near cuts the units in front (a camera radius below the scene's), which get their caps (`solidInterior`, on for
    // Gaussian surfaces): decided per unit, not by the merged item's sphere
    ...mergedCase('gaussian-surface-near-cut', GaussianSurfaceRepresentationProvider, { tryUseGpu: false }),
    view: { direction: [0, 1, -0.4] },
    camera: { radius: 14 },
  },
];
