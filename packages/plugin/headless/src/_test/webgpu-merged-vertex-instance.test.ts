/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import { Sphere3D } from '@molstar/core/math/geometry/primitives/sphere3d';
import { Color } from '@molstar/core/util/color';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { NullLocation } from '@molstar/model/model/location';
import { Mesh } from '@molstar/graphics/geo/geometry/mesh/mesh';
import { BaseGeometry } from '@molstar/graphics/geo/geometry/base';
import { createTransform } from '@molstar/graphics/geo/geometry/transform-data';
import { createOverpaint } from '@molstar/graphics/geo/geometry/overpaint-data';
import { LocationIterator } from '@molstar/graphics/geo/util/location-iterator';
import {
  canMergeRenderObjects,
  createMergedRenderObject,
  createRenderObject,
  type GraphicsRenderObject,
} from '@molstar/graphics/gl/render-object';
import { Renderer } from '@molstar/graphics/gl/renderer';
import { createFrame } from '@molstar/graphics/gl/renderable';
import { Scene } from '@molstar/graphics/gl/scene';
import { Camera } from '@molstar/graphics/canvas3d/camera';
import { Viewport } from '@molstar/graphics/canvas3d/camera/util';
import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { WebGPU } from './webgpu-module.js';

const Size = 64;

/**
 * A grid of `quads` quads (two triangles each) in the xy plane, as `copies` instances along y from `offset`, colored per
 * vertex and instance by its placed position (`vertexInstance`: each copy its own colors); `overpaint`: a layer per
 * vertex and instance too, as smoothing makes (half of each copy's vertices painted)
 */
function createMember(quads: number, copies: number, offset: Vec3, overpaint = false) {
  const vertices: number[] = [];
  const indices: number[] = [];
  for (let q = 0; q < quads; ++q) {
    const x = q * 1.2;
    const b = vertices.length / 3;
    vertices.push(x, 0, 0, x + 1, 0, 0, x + 1, 1, 0, x, 1, 0);
    indices.push(b, b + 1, b + 2, b, b + 2, b + 3);
  }
  const vertexCount = vertices.length / 3;
  const normals = new Float32Array(vertexCount * 3);
  for (let i = 0; i < vertexCount; ++i) normals[i * 3 + 2] = 1;
  const mesh = Mesh.create(
    new Float32Array(vertices),
    new Uint32Array(indices),
    normals,
    new Float32Array(vertexCount),
    vertexCount,
    indices.length / 3,
  );
  const transforms = new Float32Array(copies * 16);
  for (let i = 0; i < copies; ++i) {
    Mat4.toArray(
      Mat4.fromTranslation(Mat4(), Vec3.add(Vec3(), offset, Vec3.create(0, i * 1.5, 0))),
      transforms,
      i * 16,
    );
  }
  const props = PD.getDefaultValues(Mesh.Params);
  const sphere = Sphere3D.create(Vec3.create(quads * 0.6, 0.5, 0), quads * 0.6 + 1);
  const transform = createTransform(transforms, copies, sphere, props.cellSize, props.batchSize);
  const simple = BaseGeometry.createSimple(Color(0xffffff), 1, transform);
  const color = {
    ...simple.theme.color,
    granularity: 'vertexInstance' as const,
    color: (l: any) => {
      const p = l.position as Vec3;
      return Color.fromRgb((p[0] * 37) % 256, (p[1] * 53) % 256, 128);
    },
  };
  const locationIt = LocationIterator(1, copies, 1, () => NullLocation);
  const values = Mesh.Utils.createValues(mesh, transform, locationIt, { ...simple.theme, color }, props);
  if (overpaint) {
    createOverpaint(vertexCount * copies, 'vertexInstance', values);
    const { array } = values.tOverpaint.ref.value;
    for (let i = 0; i < vertexCount * copies; ++i) {
      if (i % 4 >= 2) continue;
      array.set([255, (i * 29) % 256, 0, 255], i * 4);
    }
  }
  return createRenderObject('mesh', values, Mesh.Utils.createRenderableState(props), -1);
}

describe.skipIf(!WebGPU.module)('merged objects of data per vertex and instance', () => {
  it.each([
    ['colors', false],
    ['colors and an overpaint layer', true],
  ])(
    'draw as their members alone (%s): each segment its own vertices and instances',
    async (_, overpaint) => {
      const plugin = await HeadlessPluginContext.create(
        { webgpu: WebGPU.module! },
        DefaultPluginSpec(),
        { width: Size, height: Size },
        { backend: 'webgpu' },
      );
      try {
        await plugin.init();
        const graphics = plugin.canvas3d!.graphics;
        // other vertex and instance counts: the second segment's offset not the first's
        const members = () => [
          createMember(3, 2, Vec3.create(-4, -3, 0), overpaint),
          createMember(2, 3, Vec3.create(0, -1, 0), overpaint),
        ];
        const unmerged = members();
        expect(unmerged[0].values.colorType.ref.value).toBe('vertexInstance');
        if (overpaint) expect(unmerged[0].values.overpaintType.ref.value).toBe('vertexInstance');
        const merging = members();
        expect(canMergeRenderObjects(merging)).toBe(true);
        const scene = (objects: GraphicsRenderObject[]) => {
          const s = Scene.create(graphics);
          for (const o of objects) s.add(o);
          s.commit();
          return s;
        };
        const unmergedScene = scene(unmerged);
        const mergedScene = scene([createMergedRenderObject(merging)]);

        const camera = new Camera();
        Viewport.set(camera.viewport, 0, 0, Size, Size);
        const { center, radius } = unmergedScene.boundingSphere;
        camera.setState(camera.getFocus(center, radius), 0);
        camera.update();
        const renderer = Renderer.create(graphics);
        const target = graphics.createRenderTarget(Size, Size, 'depth-stencil');
        for (const s of [unmergedScene, mergedScene]) {
          await graphics.resources.preparePipelines(
            renderer.collect({ name: 'opaque' }, s, graphics.describeTarget(target)),
          );
        }
        const frame = async (s: Scene) => {
          renderer.update(camera, s, createFrame());
          renderer.prepareCull(s);
          graphics.beginRenderPass(target);
          renderer.setViewport(0, 0, Size, Size);
          renderer.clear(true);
          renderer.renderOpaque(s, camera);
          graphics.flush();
          const pixels = Array.from(await graphics.readRenderTarget(target, 0, 0, Size, Size));
          await graphics.resources.pending.whenReady();
          return pixels;
        };
        // the second frame of each (the first clears the background otherwise)
        await frame(unmergedScene);
        const expected = await frame(unmergedScene);
        expect(new Set(expected).size).toBeGreaterThan(8);
        await frame(mergedScene);
        expect(await frame(mergedScene)).toEqual(expected);
        target.destroy();
        expect(graphics.takeErrors?.() ?? []).toEqual([]);
      } finally {
        plugin.dispose();
      }
    },
    60_000,
  );
});
