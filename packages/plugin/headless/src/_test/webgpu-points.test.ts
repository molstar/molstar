/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { Vec2, Vec3 } from '@molstar/core/math/linear-algebra';
import { Color } from '@molstar/core/util/color';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Points } from '@molstar/graphics/geo/geometry/points/points';
import { PointsBuilder } from '@molstar/graphics/geo/geometry/points/points-builder';
import * as Shape from '@molstar/graphics/geo/shape/shape';
import { ShapeRepresentation } from '@molstar/graphics/repr/shape/representation';
import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { WebGPU } from './webgpu-module.js';

const Size = 128;
/** The point's size on screen, CSS pixels across */
const Diameter = 60;

/** A point at the origin, `color`, as a representation of its own */
function pointRepr(color: number) {
  const builder = PointsBuilder.create(1, 1);
  builder.add(0, 0, 0, 0);
  return ShapeRepresentation(
    () =>
      Shape.create(
        'point',
        {},
        builder.getPoints(),
        () => Color(color),
        () => 1,
        () => 'point',
      ),
    Points.Utils,
  );
}

describe.skipIf(!WebGPU.module)('points on WebGPU', () => {
  // a `ring` point (a ring, its width a fraction of the radius) is drawn and picked on its ring only: not in its hole,
  // not outside it (the fragment's discard in every variant), so a wider pick-only ring is a hit band
  it('picks a ring point on its ring, not in its hole nor outside', async () => {
    const plugin = await HeadlessPluginContext.create(
      { webgpu: WebGPU.module! },
      DefaultPluginSpec(),
      { width: Size, height: Size },
      { backend: 'webgpu' },
    );
    try {
      await plugin.init();
      const canvas = plugin.canvas3d!;
      const builder = PointsBuilder.create(1, 1);
      builder.add(0, 0, 0, 0);
      const ring = ShapeRepresentation(
        () =>
          Shape.create(
            'ring',
            {},
            builder.getPoints(),
            () => Color(0x377eb8),
            () => 1,
            () => 'ring',
          ),
        Points.Utils,
      );
      const props = {
        ...PD.getDefaultValues(Points.Params),
        pointStyle: 'ring' as const,
        pointRingWidth: 0.3,
        sizeFactor: Diameter,
        pointSizeAttenuation: false,
      };
      await ring.createOrUpdate(props, {}).run();
      canvas.add(ring);
      canvas.commit(true);
      canvas.camera.setState({ position: Vec3.create(0, 0, 30), target: Vec3(), up: Vec3.unitY }, 0);

      const pixelRatio = canvas.graphics.pixelRatio;
      const c = Size / 2;
      const pickAt = async (dx: number) => {
        const pick = await canvas.asyncIdentify(Vec2.create((c + dx) / pixelRatio, c / pixelRatio))?.wait();
        return pick ? canvas.getLoci(pick.id).repr : undefined;
      };
      // the pick's pipelines compile in the background
      const onRing = (Diameter / 2) * 0.85;
      for (let i = 0; i < 20 && !(await pickAt(onRing)); i++);

      expect(await pickAt(onRing)).toBe(ring);
      expect(await pickAt(-onRing)).toBe(ring);
      // in the hole (inside 70% of the radius) and outside the point
      expect(await pickAt(0)).toBeUndefined();
      expect(await pickAt((Diameter / 2) * 0.5)).toBeUndefined();
      expect(await pickAt((Diameter / 2) * 1.3)).toBeUndefined();

      // a full disk (width 1): its center picked too
      await ring.createOrUpdate({ ...props, pointRingWidth: 1 }, {}).run();
      canvas.update(ring);
      canvas.commit(true);
      for (let i = 0; i < 20 && !(await pickAt(0)); i++);
      expect(await pickAt(0)).toBe(ring);
      expect(canvas.graphics.takeErrors?.() ?? []).toEqual([]);
    } finally {
      plugin.dispose();
    }
  }, 60_000);
  // a point's depth offset (CSS pixels toward the camera, negative away) in depth only: of two points at one place, the
  // one moved toward the camera picked, moved away the other; where it is on screen kept
  it('picks a point moved in depth by its offset over, or under, another at its place', async () => {
    const plugin = await HeadlessPluginContext.create(
      { webgpu: WebGPU.module! },
      DefaultPluginSpec(),
      { width: Size, height: Size },
      { backend: 'webgpu' },
    );
    try {
      await plugin.init();
      const canvas = plugin.canvas3d!;
      const props = {
        ...PD.getDefaultValues(Points.Params),
        pointStyle: 'circle' as const,
        sizeFactor: Diameter,
        pointSizeAttenuation: false,
      };
      const a = pointRepr(0x377eb8);
      const b = pointRepr(0xe41a1c);
      await a.createOrUpdate(props, {}).run();
      await b.createOrUpdate(props, {}).run();
      canvas.add(a);
      canvas.add(b);
      canvas.commit(true);
      canvas.camera.setState({ position: Vec3.create(0, 0, 30), target: Vec3(), up: Vec3.unitY }, 0);

      const pixelRatio = canvas.graphics.pixelRatio;
      const c = Size / 2;
      let time = 0;
      const pickAt = async (dx: number) => {
        // a frame first: the pick renders what was drawn
        canvas.requestDraw();
        canvas.tick((time += 1000) as never, { isSynchronous: true });
        await canvas.graphics.resources.pending.whenReady();
        const pick = await canvas.asyncIdentify(Vec2.create((c + dx) / pixelRatio, c / pixelRatio))?.wait();
        return pick ? canvas.getLoci(pick.id).repr : undefined;
      };
      const offset = async (r: typeof a, pointDepthOffset: number) => {
        await r.createOrUpdate({ ...props, pointDepthOffset }, {}).run();
        canvas.update(r);
        canvas.commit(true);
      };
      for (let i = 0; i < 20 && !(await pickAt(0)); i++);

      for (const [front, back] of [
        [a, b],
        [b, a],
      ]) {
        await offset(front, 20);
        await offset(back, 0);
        expect(await pickAt(0)).toBe(front);
        // its place on screen kept: picked to its edge, not beyond
        expect(await pickAt((Diameter / 2) * 0.8)).toBe(front);
        expect(await pickAt((Diameter / 2) * 1.3)).toBeUndefined();
        // far past the depth range the planes fit to the bounds: at the near plane, still in front
        await offset(front, 200);
        expect(await pickAt(0)).toBe(front);
        // moved away: the other one
        await offset(front, -20);
        expect(await pickAt(0)).toBe(back);
      }
      expect(canvas.graphics.takeErrors?.() ?? []).toEqual([]);
    } finally {
      plugin.dispose();
    }
  }, 60_000);
});
