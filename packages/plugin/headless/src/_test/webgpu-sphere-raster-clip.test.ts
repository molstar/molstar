/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { WebGPU } from './webgpu-module.js';

const Pdb = `ATOM      1  N   ALA A   1      -1.458   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C
ATOM      3  C   ALA A   1       0.551   1.420   0.000  1.00  0.00           C
ATOM      4  O   ALA A   1      -0.200   2.390   0.000  1.00  0.00           O
ATOM      5  CB  ALA A   1       0.540  -0.760  -1.210  1.00  0.00           C
ATOM      6  N   GLY A   2       1.870   1.550   0.000  1.00  0.00           N
ATOM      7  CA  GLY A   2       2.520   2.860   0.000  1.00  0.00           C
ATOM      8  C   GLY A   2       4.030   2.720   0.000  1.00  0.00           C
ATOM      9  O   GLY A   2       4.570   1.610   0.000  1.00  0.00           O
END
`;

const Size = 96;

/** A 4×4×4 grid of carbon atoms 2.6 Å apart, centered: atoms behind atoms along every view */
function gridPdb() {
  const f = (v: number) => v.toFixed(3).padStart(8);
  const lines: string[] = [];
  let n = 0;
  for (let i = 0; i < 4; ++i) {
    for (let j = 0; j < 4; ++j) {
      for (let k = 0; k < 4; ++k) {
        n += 1;
        const [x, y, z] = [i, j, k].map((v) => (v - 1.5) * 2.6);
        lines.push(
          `ATOM  ${String(n).padStart(5)}  CA  ALA A${String(n).padStart(4)}    ${f(x)}${f(y)}${f(z)}  1.00  0.00           C`,
        );
      }
    }
  }
  return `${lines.join('\n')}\nEND\n`;
}

const object = (type: 'plane' | 'sphere', invert: boolean) => ({
  type,
  invert,
  position: Vec3.create(1.5, 1.4, 0),
  rotation: { axis: Vec3.create(1, 0.3, 0), angle: type === 'plane' ? 70 : 0 },
  scale: Vec3.create(4, 4, 4),
  transform: Mat4.identity(),
});

/** The structure as spacefill clipped by pixel, the spheres by the raster (`only`) or as impostors */
async function render(sphereRaster: 'off' | 'only', typeParams: { [k: string]: unknown }, pdb = Pdb) {
  const plugin = await HeadlessPluginContext.create(
    { webgpu: WebGPU.module! },
    DefaultPluginSpec(),
    { width: Size, height: Size },
    { backend: 'webgpu' },
  );
  try {
    await plugin.init();
    const { graphics } = plugin.canvas3d!;
    if (!graphics.capabilities.atomicMinMax64) return undefined;
    plugin.canvas3d!.setProps({
      renderer: { sphereRaster, sphereRasterRadius: 32 },
      camera: { mode: 'perspective' } as never,
      postprocessing: { occlusion: { name: 'off', params: {} }, outline: { name: 'off', params: {} } },
    });
    const data = await plugin.builders.data.rawData({ data: pdb, label: 'test' });
    const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
    const model = await plugin.builders.structure.createModel(trajectory);
    const structure = await plugin.builders.structure.createStructure(model);
    await plugin.builders.structure.representation.addRepresentation(structure, {
      type: 'spacefill',
      typeParams: typeParams as never,
      color: 'element-index',
    });
    // the grid's front face toward the camera (its own view would be along a principal axis of either sign)
    if (pdb !== Pdb) {
      const { camera } = plugin.canvas3d!;
      camera.setState(camera.getFocus(Vec3(), 9, Vec3.create(0, 1, 0), Vec3.create(0, 0, -1)), 0);
    }
    const image = () => plugin.getImageRaw({ width: Size, height: Size });
    // pipelines and kernels compile in the background: frames drawn before they are ready draw impostors
    for (let i = 0; i < 3; ++i) {
      await image();
      await graphics.resources.pending.whenReady();
    }
    const resolves = graphics.stats.sphereRasterResolves;
    const pixels = Array.from((await image()).data);
    const rasterDrawn = graphics.stats.sphereRasterResolves > resolves;
    expect(graphics.takeErrors?.() ?? []).toEqual([]);
    return { pixels, rasterDrawn };
  } finally {
    plugin.dispose();
  }
}

/** pixels whose color differs by more than 8 levels */
function countDiffer(a: number[], b: number[]) {
  let differ = 0;
  for (let i = 0; i < a.length; i += 4) {
    const d = Math.max(Math.abs(a[i] - b[i]), Math.abs(a[i + 1] - b[i + 1]), Math.abs(a[i + 2] - b[i + 2]));
    if (d > 8) differ += 1;
  }
  return differ;
}

/** The grid's front face (toward the camera) holed by a clip sphere removing its inside */
const hole = {
  type: 'sphere' as const,
  invert: false,
  position: Vec3.create(0, 0, 3.9),
  rotation: { axis: Vec3.create(1, 0, 0), angle: 0 },
  scale: Vec3.create(8, 8, 8),
  transform: Mat4.identity(),
};

/**
 * The clip modes; a cap of a convex clip region (a plane, an inverted sphere) is the nearest point no clip object removes
 * along its ray, so a raster writing the clipped fronts would draw the same (the resolve finds the cap): those cases
 * check that the raster's clipping agrees, the others (backs, nothing, a hole, whole spheres) differ without it
 */
const Cases: { name: string; params: { [k: string]: unknown }; pdb?: string }[] = [
  {
    name: 'a plane, the caps of the solid interiors',
    params: { clip: { variant: 'pixel', objects: [object('plane', false)] } },
  },
  {
    name: 'a plane, the backs (double sided, no solid interior)',
    params: { clip: { variant: 'pixel', objects: [object('plane', false)] }, solidInterior: false, doubleSided: true },
  },
  {
    name: 'a plane, nothing behind the clipped fronts (no solid interior, single sided)',
    params: { clip: { variant: 'pixel', objects: [object('plane', false)] }, solidInterior: false },
  },
  {
    name: 'an inverted sphere, the caps of the solid interiors',
    params: { clip: { variant: 'pixel', objects: [object('sphere', true)] } },
  },
  {
    // atoms within the clip sphere clipped through along their rays: those behind them show (with the fronts written,
    // the resolve would find no hit there and leave the pixels empty)
    name: 'a sphere holing a grid of atoms, spheres clipped through and those behind',
    params: { clip: { variant: 'pixel', objects: [hole] } },
    pdb: gridPdb(),
  },
  {
    name: 'a plane, spheres whole by their centers (clipPrimitive)',
    params: { clip: { variant: 'pixel', objects: [object('plane', false)] }, clipPrimitive: true },
  },
];

// pixel clipping in the raster as the impostors clip: per pixel the point a sphere shows (its front, a cap, its back,
// none), or spheres whole by their centers
describe.skipIf(!WebGPU.module)('sphere raster pixel clipping on WebGPU', () => {
  for (const { name, params, pdb } of Cases) {
    it(`draws the image of the impostors: ${name}`, async () => {
      const impostors = await render('off', params, pdb);
      if (!impostors) return;
      const unclipped = (await render('off', {}, pdb))!.pixels;
      // the clip objects cut what is drawn
      expect(countDiffer(impostors.pixels, unclipped)).toBeGreaterThan(50);
      const raster = (await render('only', params, pdb))!;
      expect(raster.rasterDrawn).toBe(true);
      // but where a pixel's center lies on an edge (the raster's arithmetic is not the fragment stage's)
      expect(countDiffer(impostors.pixels, raster.pixels)).toBeLessThanOrEqual(Size * Size * 0.01);
    }, 60_000);
  }
});
