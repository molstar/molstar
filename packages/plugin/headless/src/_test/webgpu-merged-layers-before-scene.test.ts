/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// a merged surface's layer set before a scene adds its objects: drawn, as one set after
import { describe, expect, it } from 'vitest';
import { Color } from '@molstar/core/util/color';
import { Interval } from '@molstar/core/data/int';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { StructureElement } from '@molstar/model/model/structure';
import { Renderer } from '@molstar/graphics/gl/renderer';
import { createFrame } from '@molstar/graphics/gl/renderable';
import { Scene } from '@molstar/graphics/gl/scene';
import { Camera } from '@molstar/graphics/canvas3d/camera';
import { Viewport } from '@molstar/graphics/canvas3d/camera/util';
import { Theme } from '@molstar/graphics/theme/theme';
import { GaussianSurfaceRepresentationProvider } from '@molstar/graphics/repr/structure/representation/gaussian-surface';
import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { createWebGPUPlugin } from './webgpu-plugin.js';
import { WebGPU } from './webgpu-module.js';

const Size = 96;
const TwoChainsPdb = `${Array.from({ length: 48 }, (_, i) => {
  const f = (v: number) => v.toFixed(3).padStart(8);
  const chain = i < 24 ? 'A' : 'B';
  const k = i % 24;
  const [x, y, z] = [(i < 24 ? 0 : 6) + Math.cos(k * 0.6) * 2, Math.sin(k * 0.6) * 2, k * 0.6];
  return `ATOM  ${String(i + 1).padStart(5)}  CA  ALA ${chain}${String(k + 1).padStart(4)}    ${f(x)}${f(y)}${f(z)}  1.00  0.00           C`;
}).join('\n')}\nEND\n`;

describe.skipIf(!WebGPU.module)('merged layers set before a scene adds the objects (WebGPU)', () => {
  for (const order of ['state then scene', 'scene then state'] as const) {
    it(`draws the overpaint: ${order}`, async () => {
      const plugin = await createWebGPUPlugin(DefaultPluginSpec(), { width: Size, height: Size });
      try {
        await plugin.init();
        const graphics = plugin.canvas3d!.graphics;
        const data = await plugin.builders.data.rawData({ data: TwoChainsPdb, label: 'test' });
        const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
        const model = await plugin.builders.structure.createModel(trajectory);
        const sc = await plugin.builders.structure.createStructure(model);
        const s = sc.obj!.data;
        const registries = {
          colorThemeRegistry: plugin.representation.structure.themes.colorThemeRegistry,
          sizeThemeRegistry: plugin.representation.structure.themes.sizeThemeRegistry,
        };
        const provider = GaussianSurfaceRepresentationProvider;
        const repr = provider.factory({ graphics, ...registries } as any, provider.getParams);
        repr.setTheme(
          Theme.create(
            registries as any,
            { structure: s },
            { colorTheme: { name: 'uniform' }, sizeTheme: { name: 'physical' } },
          ),
        );
        await repr
          .createOrUpdate(
            {
              ...PD.getDefaultValues(provider.getParams(registries as any, s)),
              resolution: 0.5,
              tryUseGpu: false,
              smoothColors: { name: 'on', params: { resolutionFactor: 2, sampleStride: 3 } },
            } as any,
            s,
          )
          .run();
        const unit = s.units[0];
        const loci = StructureElement.Loci(s, [{ unit, indices: Interval.ofLength(unit.elements.length) }]);
        const state: any = {
          overpaint: { kind: 'element-loci', layers: [{ loci, color: Color(0xff0000), clear: false }] },
        };
        const scene = Scene.create(graphics);
        if (order === 'state then scene') repr.setState(state);
        for (const o of repr.renderObjects) scene.add(o);
        scene.commit();
        if (order === 'scene then state') {
          const before = [...repr.renderObjects];
          repr.setState(state);
          // as Canvas3D's `update(repr)`: objects a smoothed layer split off added, those gone removed
          for (const o of repr.renderObjects) if (!before.includes(o)) scene.add(o);
          for (const o of before) if (!repr.renderObjects.includes(o)) scene.remove(o);
          scene.update(repr.renderObjects, false);
          scene.commit();
        }
        const camera = new Camera();
        Viewport.set(camera.viewport, 0, 0, Size, Size);
        const { center, radius } = scene.boundingSphere;
        camera.setState(camera.getFocus(center, radius), 0);
        camera.update();
        const renderer = Renderer.create(graphics);
        const target = graphics.createRenderTarget(Size, Size, 'depth-stencil');
        await graphics.resources.preparePipelines(
          renderer.collect({ name: 'opaque' }, scene, graphics.describeTarget(target)),
        );
        let pixels: number[] = [];
        for (let i = 0; i < 3; ++i) {
          renderer.update(camera, scene, createFrame());
          graphics.beginRenderPass(target);
          renderer.setViewport(0, 0, Size, Size);
          renderer.clear(true);
          renderer.renderOpaque(scene, camera);
          graphics.flush();
          await graphics.resources.pending.whenReady();
          pixels = Array.from(await graphics.readRenderTarget(target, 0, 0, Size, Size));
        }
        let red = 0;
        for (let i = 0; i < pixels.length; i += 4) if (pixels[i] > 180 && pixels[i + 1] < 90) red++;
        expect(red).toBeGreaterThan(50);
      } finally {
        plugin.dispose();
      }
    }, 60_000);
  }
});
