/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { OrderedSet } from '@molstar/core/data/int';
import { SpacefillRepresentationProvider } from '@molstar/graphics/repr/structure/representation/spacefill';
import { IsosurfaceRepresentationProvider } from '@molstar/graphics/repr/volume/isosurface';
import { StructureElement, type Structure } from '@molstar/model/model/structure';
import { Volume } from '@molstar/model/model/volume';
import { type VisualTestCase } from '../framework/test-case.js';
import { getThemeRegistries, helixStructure, structureRepr } from '../scenes/structure.js';
import { Tensor } from '@molstar/core/math/linear-algebra';
import { Theme } from '@molstar/graphics/theme/theme';
import { blobVolume, createVolume } from '../scenes/volume.js';

import { Impostor, ObliqueView } from './geometry.js';

/** A linear gradient along x on the grid of `blobVolume(n)`, from 0 to 1 */
function gradientVolume(n: number) {
  const space = Tensor.Space([n, n, n], [0, 1, 2], Float32Array);
  const data = new Float32Array(n * n * n);
  for (let x = 0; x < n; ++x) {
    for (let y = 0; y < n; ++y) {
      for (let z = 0; z < n; ++z) data[space.dataOffset(x, y, z)] = x / (n - 1);
    }
  }
  return createVolume('gradient', space, data);
}

/**
 * An isosurface of the blobs colored by the value of a gradient volume at its vertices (`volume-value`, a palette
 * lookup, `dUsePalette`)
 */
function paletteCase(name: string, tryUseGpu: boolean): VisualTestCase {
  return {
    id: `geometry/palette/${name}`,
    description: `An isosurface (${tryUseGpu ? 'GPU texture mesh' : 'CPU mesh'}) colored by a gradient volume through a palette.`,
    requires: tryUseGpu ? ['floatRenderTarget', 'multipleRenderTargets'] : undefined,
    view: ObliqueView,
    scene: async ({ graphics }) => {
      const volume = blobVolume(32);
      const registries = getThemeRegistries();
      const provider = IsosurfaceRepresentationProvider;
      const repr = provider.factory({ graphics, ...registries }, provider.getParams);
      repr.setTheme(
        Theme.create(
          registries,
          { volume: gradientVolume(32) },
          {
            colorTheme: { name: 'volume-value', params: { domain: { name: 'custom', params: [0.2, 0.8] } } },
            sizeTheme: { name: 'uniform', params: {} },
          },
        ),
      );
      await repr.createOrUpdate({ isoValue: Volume.IsoValue.absolute(0.5), tryUseGpu }, volume).run();
      return [repr];
    },
  };
}

/** The atoms of chain A */
function firstUnitLoci(structure: Structure) {
  const unit = structure.units[0];
  const indices = OrderedSet.ofBounds(
    0 as StructureElement.UnitIndex,
    unit.elements.length as StructureElement.UnitIndex,
  );
  return StructureElement.Loci(structure, [{ unit, indices }]);
}

/**
 * Wiggle data (`tWiggle`, per group instance) on chain A with animation enabled: the image is rendered at time 0, where
 * the wiggle noise displaces the atoms by a fixed amount (chain B, without wiggle, stays in place)
 */
function wiggleCase(
  name: string,
  props: { [k: string]: unknown },
  requires?: VisualTestCase['requires'],
): VisualTestCase {
  return {
    id: `renderer/wiggle/${name}`,
    description: `Spacefill (${name}) with wiggle data on chain A, animation enabled, rendered at time 0.`,
    requires,
    view: ObliqueView,
    canvasProps: { renderer: { enableAnimation: true } },
    scene: async ({ graphics }) => {
      const structure = await helixStructure(24);
      const repr = await structureRepr(graphics, SpacefillRepresentationProvider, structure, {
        ...props,
        animation: { wiggleMode: 'position', wiggleSpeed: 1, wiggleAmplitude: 0, wiggleFrequency: 0.3 },
      } as any);
      repr.setState({ wiggle: { kind: 'element-loci', layers: [{ loci: firstUnitLoci(structure), value: 1 }] } });
      return [repr];
    },
  };
}

export const ThemeDataCases: VisualTestCase[] = [
  paletteCase('isosurface', false),
  paletteCase('isosurface-gpu', true),
  wiggleCase('impostor', {}, Impostor),
  wiggleCase('mesh', { tryUseImpostor: false }),
];
