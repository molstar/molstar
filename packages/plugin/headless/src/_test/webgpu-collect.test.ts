/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { Renderer } from '@molstar/graphics/gl/renderer';
import { createFrame } from '@molstar/graphics/gl/renderable';
import type { RenderState } from '@molstar/graphics/gl/render-state';
import { Scene } from '@molstar/graphics/gl/scene';
import { Camera } from '@molstar/graphics/canvas3d/camera';
import { Viewport } from '@molstar/graphics/canvas3d/camera/util';
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

describe.skipIf(!WebGPU.module)('Renderer.collect on WebGPU', () => {
  it('collects the pipelines each renderer pass draws with', async () => {
    const plugin = await HeadlessPluginContext.create(
      { webgpu: WebGPU.module! },
      DefaultPluginSpec(),
      { width: 64, height: 64 },
      { backend: 'webgpu' },
    );
    try {
      await plugin.init();
      const data = await plugin.builders.data.rawData({ data: Pdb, label: 'test' });
      const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
      const model = await plugin.builders.structure.createModel(trajectory);
      const structure = await plugin.builders.structure.createStructure(model);
      const { representation } = plugin.builders.structure;
      // meshes, impostors and points; opaque and transparent, single and double sided
      await representation.addRepresentation(structure, { type: 'cartoon' });
      await representation.addRepresentation(structure, { type: 'ball-and-stick' });
      await representation.addRepresentation(structure, { type: 'point' });
      await representation.addRepresentation(structure, {
        type: 'molecular-surface',
        typeParams: { alpha: 0.5, doubleSided: true, transparentBackfaces: 'on' },
      });
      await representation.addRepresentation(structure, {
        type: 'gaussian-surface',
        typeParams: { alpha: 0.4, doubleSided: false, emissive: 0.5 },
      });
      plugin.canvas3d!.commit(true);

      const graphics = plugin.canvas3d!.graphics;
      const scene = Scene.create(graphics);
      for (const ro of plugin.canvas3d!.getRenderObjects()) scene.add(ro);
      scene.commit();
      expect(scene.renderables.length).toBeGreaterThan(4);
      // the whole scene in view, so that frustum culling leaves every object drawn
      const camera = new Camera();
      Viewport.set(camera.viewport, 0, 0, 64, 64);
      const { center, radius } = scene.boundingSphere;
      camera.setState(camera.getFocus(center, radius), 0);
      camera.update();

      const renderer = Renderer.create(graphics);
      const target = graphics.createRenderTarget(64, 64, 'depth-stencil');
      const depth = graphics.resources.texture('image-depth', 'depth', 'ushort', 'nearest');
      depth.define(64, 64);
      const descriptor = graphics.describeTarget(target);
      const { pending } = graphics.resources;

      const passes: [Renderer.Pass, () => void][] = [
        [{ name: 'pick', variant: 'pick' }, () => renderer.renderPick(scene, camera, 'pick', 0)],
        [{ name: 'pick', variant: 'depth' }, () => renderer.renderPick(scene, camera, 'depth', 0)],
        [{ name: 'depth' }, () => renderer.renderDepth(scene, camera)],
        [{ name: 'depthOpaque' }, () => renderer.renderDepthOpaque(scene, camera)],
        [{ name: 'depthOpaqueBack' }, () => renderer.renderDepthOpaqueBack(scene, camera)],
        [{ name: 'depthTransparent' }, () => renderer.renderDepthTransparent(scene, camera, depth)],
        [{ name: 'markingDepth' }, () => renderer.renderMarkingDepth(scene, camera)],
        [{ name: 'markingMask' }, () => renderer.renderMarkingMask(scene, camera, null)],
        [
          { name: 'emissiveOpaque', occludeWithOpaqueDepth: false },
          () => renderer.renderEmissiveOpaque(scene, camera, depth, false),
        ],
        [
          { name: 'emissiveOpaque', occludeWithOpaqueDepth: true },
          () => renderer.renderEmissiveOpaque(scene, camera, depth, true),
        ],
        [{ name: 'emissiveTransparent' }, () => renderer.renderEmissiveTransparent(scene, camera, depth)],
        [{ name: 'tracing' }, () => renderer.renderTracing(scene, camera)],
        [{ name: 'opaque' }, () => renderer.renderOpaque(scene, camera)],
        [{ name: 'blendedTransparent' }, () => renderer.renderBlendedTransparent(scene, camera)],
        [{ name: 'volume' }, () => renderer.renderVolume(scene, camera, depth)],
        [{ name: 'wboitTransparent' }, () => renderer.renderWboitTransparent(scene, camera, depth)],
      ];

      // the draws the renderer issues: render item, variant and the render state set for the draw
      const key = (item: { id: number }, variant: string, state: RenderState) =>
        `${item.id}|${variant}|${JSON.stringify(state)}`;
      let currentState: RenderState = {};
      const setRenderState = graphics.setRenderState;
      graphics.setRenderState = (state) => {
        currentState = state;
        setRenderState(state);
      };
      const drawn: string[] = [];
      for (const r of scene.renderables) {
        const item = r.renderItem;
        const render = item.render;
        item.render = (variant, ...args) => {
          drawn.push(key(item, variant, currentState));
          render(variant, ...args);
        };
      }

      const culls = new Set<string>();
      for (const transparentBackground of [false, true]) {
        renderer.setTransparentBackground(transparentBackground);
        for (const [pass, render] of passes) {
          const requests = renderer.collect(pass, scene, descriptor);
          await graphics.resources.preparePipelines(requests);

          drawn.length = 0;
          const tracked = pending.trackedCount;
          graphics.beginRenderPass(target);
          renderer.update(camera, scene, createFrame());
          renderer.setViewport(0, 0, 64, 64);
          render();
          graphics.flush();
          // every draw of the pass was collected (objects out of view are collected too), and found its pipeline
          const collected = new Set(requests.map((r) => key(r.item, r.variant, r.state)));
          expect({ pass, missing: drawn.filter((d) => !collected.has(d)) }).toEqual({ pass, missing: [] });
          expect({ pass, started: pending.trackedCount - tracked }).toEqual({ pass, started: 0 });
          for (const r of requests) culls.add(`${pass.name}|${r.state.cull ?? 'none'}`);
        }
      }
      // the scene covers per-object states: back and front faces of the double-sided transparent surface, single-sided
      // objects
      expect(culls).toContain('blendedTransparent|front');
      expect(culls).toContain('blendedTransparent|back');
      expect(culls).toContain('depthTransparent|back');

      // with `all`, also what the objects' state leaves out: the marking mask of unmarked objects
      expect(renderer.collect({ name: 'markingMask' }, scene, descriptor)).toEqual([]);
      expect(renderer.collect({ name: 'markingMask' }, scene, descriptor, { all: true }).length).toBeGreaterThan(0);

      renderer.dispose();
      scene.clear();
    } finally {
      plugin.dispose();
    }
  });
});
