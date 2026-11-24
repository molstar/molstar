/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// a smoothed layer on a merged surface: drawn as the per-group visuals' (their grids), not as the plain layer
import { describe, expect, it } from 'vitest';
import { Color } from '@molstar/core/util/color';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import { Interval } from '@molstar/core/data/int';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { type Structure, StructureElement } from '@molstar/model/model/structure';
import { Renderer } from '@molstar/graphics/gl/renderer';
import { createFrame } from '@molstar/graphics/gl/renderable';
import { isMergedRenderObject } from '@molstar/graphics/gl/render-object';
import { Scene } from '@molstar/graphics/gl/scene';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import { Camera } from '@molstar/graphics/canvas3d/camera';
import { Viewport } from '@molstar/graphics/canvas3d/camera/util';
import { Theme } from '@molstar/graphics/theme/theme';
import { UnitsRepresentation } from '@molstar/graphics/repr/structure/units-representation';
import {
  GaussianSurfaceMeshParams,
  GaussianSurfaceMeshVisual,
} from '@molstar/graphics/repr/structure/visual/gaussian-surface-mesh';
import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { WebGPU } from './webgpu-module.js';

const Size = 128;
const TwoChainsPdb = `${Array.from({ length: 48 }, (_, i) => {
  const f = (v: number) => v.toFixed(3).padStart(8);
  const chain = i < 24 ? 'A' : 'B';
  const k = i % 24;
  const [x, y, z] = [(i < 24 ? 0 : 6) + Math.cos(k * 0.6) * 2, Math.sin(k * 0.6) * 2, k * 0.6];
  return `ATOM  ${String(i + 1).padStart(5)}  CA  ALA ${chain}${String(k + 1).padStart(4)}    ${f(x)}${f(y)}${f(z)}  1.00  0.00           C`;
}).join('\n')}\nEND\n`;

const Smooth = { name: 'on', params: { resolutionFactor: 2, sampleStride: 3 } };

describe.skipIf(!WebGPU.module)('smoothed layers on merged surfaces (WebGPU)', () => {
  it('draws an overpaint on one chain as the per-group visuals do', async () => {
    const plugin = await HeadlessPluginContext.create(
      { webgpu: WebGPU.module! },
      DefaultPluginSpec(),
      { width: Size, height: Size },
      { backend: 'webgpu' },
    );
    try {
      await plugin.init();
      const graphics = plugin.canvas3d!.graphics;
      const data = await plugin.builders.data.rawData({ data: TwoChainsPdb, label: 'test' });
      const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
      const model = await plugin.builders.structure.createModel(trajectory);
      const structure = (await plugin.builders.structure.createStructure(model)).obj!.data;
      const registries = {
        colorThemeRegistry: plugin.representation.structure.themes.colorThemeRegistry,
        sizeThemeRegistry: plugin.representation.structure.themes.sizeThemeRegistry,
      };
      const unit = structure.units[0];
      // the half of chain A toward the camera (along the helices' axis): the smoothing seen at the layer's edge
      const n = unit.elements.length;
      const loci = StructureElement.Loci(structure, [{ unit, indices: Interval.ofBounds(n / 2, n) }]);
      const overpaint = { kind: 'element-loci', layers: [{ loci, color: Color(0xff0000), clear: false }] } as const;

      const merged = await image(graphics, registries, structure, 'merged', overpaint);
      const perGroup = await image(graphics, registries, structure, 'per-group', overpaint);
      const plain = await image(graphics, registries, structure, 'merged', overpaint, { name: 'off', params: {} });
      const toPerGroup = differing(merged.pixels, perGroup.pixels);
      const plainToPerGroup = differing(plain.pixels, perGroup.pixels);
      // a part of its own for the smoothed chain beside the plain one
      expect(merged.objects).toBe(2);
      // the plain layer's hard edge apart from the grids' soft one; merged per vertex near them (measured 92 against 264
      // pixels: the vertices' sampling, not the grid's)
      expect(plainToPerGroup).toBeGreaterThan(150);
      expect(toPerGroup).toBeLessThan(plainToPerGroup / 2);
    } finally {
      plugin.dispose();
    }
  }, 60_000);
});

/** pixels of `a` and `b` that differ by more than 24 summed over rgb */
function differing(a: Uint8Array, b: Uint8Array) {
  let n = 0;
  for (let i = 0; i < a.length; i += 4) {
    if (Math.abs(a[i] - b[i]) + Math.abs(a[i + 1] - b[i + 1]) + Math.abs(a[i + 2] - b[i + 2]) > 24) n++;
  }
  return n;
}

/** a two-chain Gaussian surface with `overpaint`, merged or per group (a builder-less visual: grids), drawn opaque */
async function image(
  graphics: GraphicsContext,
  registries: any,
  structure: Structure,
  mode: 'merged' | 'per-group',
  overpaint: any,
  smoothColors: any = Smooth,
) {
  const visualCtor = (materialId: number) => {
    const visual = GaussianSurfaceMeshVisual(materialId);
    return mode === 'merged'
      ? visual
      : new Proxy(visual, { get: (t, k) => (k === 'builder' ? undefined : Reflect.get(t, k)) });
  };
  const repr = UnitsRepresentation(
    'test',
    { graphics, ...registries },
    () => GaussianSurfaceMeshParams as any,
    visualCtor as any,
  );
  repr.setTheme(
    Theme.create(registries, { structure }, { colorTheme: { name: 'uniform' }, sizeTheme: { name: 'physical' } }),
  );
  const props = { ...PD.getDefaultValues(GaussianSurfaceMeshParams), resolution: 0.5, smoothColors };
  await repr.createOrUpdate(props as any, structure).run();
  repr.setState({ overpaint });
  const scene = Scene.create(graphics);
  for (const o of repr.renderObjects) scene.add(o);
  scene.commit();
  const camera = new Camera();
  Viewport.set(camera.viewport, 0, 0, Size, Size);
  const { center, radius } = scene.boundingSphere;
  // from the side: the helices' axis (z) across the view
  camera.setState(camera.getFocus(center, radius, Vec3.create(0, 0, 1), Vec3.create(0, 1, 0)), 0);
  camera.update();
  const renderer = Renderer.create(graphics);
  const target = graphics.createRenderTarget(Size, Size, 'depth-stencil');
  await graphics.resources.preparePipelines(
    renderer.collect({ name: 'opaque' }, scene, graphics.describeTarget(target)),
  );
  let pixels = new Uint8Array(0);
  for (let i = 0; i < 2; ++i) {
    renderer.update(camera, scene, createFrame());
    graphics.beginRenderPass(target);
    renderer.clear(true);
    renderer.setViewport(0, 0, Size, Size);
    renderer.renderOpaque(scene, camera);
    graphics.flush();
    await graphics.resources.pending.whenReady();
    pixels = new Uint8Array(await graphics.readRenderTarget(target, 0, 0, Size, Size));
  }
  const objects = repr.renderObjects.filter((o) => isMergedRenderObject(o)).length;
  for (const o of repr.renderObjects) scene.remove(o);
  repr.destroy();
  return { pixels, objects };
}
