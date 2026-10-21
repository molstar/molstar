/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Mat4, Quat, Vec3 } from '@molstar/core/math/linear-algebra';
import { afterEach, describe, expect, it } from 'vitest';
import { cbs } from '../../../_test/harness.js';
import { DefaultPrismViewerSpec } from '../../../default-spec.js';
import { createPrismViewerHeadless } from '../../../headless.js';
import {
  applyToPlacement,
  placementPose,
  poseMatrix,
  type StructurePlacement,
} from '../../../sources/structure-move/gizmos.js';
import type { PrismViewer } from '../../../viewer.js';
import { StructureMoveFeature } from '../index.js';

// A structure moved by hand: its placement as the handles' pose and back (a turn about its center, then a move), and
// the feature in a viewer (opt-in: not in the default spec).

const close = (a: ArrayLike<number>, b: ArrayLike<number>) =>
  Array.from(b).forEach((x, i) => expect(a[i]).toBeCloseTo(x, 6));

const placement = (matrix = Mat4.identity()): StructurePlacement => ({
  matrix,
  center: Vec3.create(10, 20, 30),
  radius: 25,
});

let viewer: PrismViewer | undefined;
afterEach(async () => {
  await viewer?.dispose();
  viewer = undefined;
});

describe('structure-move', () => {
  it('a placement is its pose: the center as placed, the rotation', () => {
    const p = placement();
    const pose = placementPose(p);
    close(pose.position, [10, 20, 30]);
    close(poseMatrix(pose, p.center), Mat4.identity());
  });

  it('a turn keeps the center where it is; a move adds a translation', () => {
    const p = placement();
    const turned = { ...placementPose(p), rotation: Quat.setAxisAngle(Quat(), Vec3.create(0, 1, 0), Math.PI / 2) };
    const m = poseMatrix(turned, p.center);
    close(Vec3.transformMat4(Vec3(), p.center, m), [10, 20, 30]);
    // and the pose of that placement is the turn again
    const back = placementPose({ ...p, matrix: m });
    expect(Math.abs(Quat.dot(back.rotation, turned.rotation))).toBeCloseTo(1, 6);

    const moved = { ...placementPose(p), position: Vec3.create(15, 20, 30) };
    close(Vec3.transformMat4(Vec3(), Vec3.create(0, 0, 0), poseMatrix(moved, p.center)), [5, 0, 0]);
  });

  it('a drag keeps what a placement has beyond a rigid move (a scale given as a matrix)', () => {
    const p = placement(Mat4.fromScaling(Mat4(), Vec3.create(2, 2, 2)));
    const frame = { right: Vec3.create(1, 0, 0), up: Vec3.create(0, 1, 0), perPixel: 0.5 };
    const next = applyToPlacement(p, { kind: 'move' }, frame, 10, 0);
    // still twice the size, its center (as placed: twice the loaded one) moved 5 Å right
    expect(Mat4.getMaxScaleOnAxis(next.matrix)).toBeCloseTo(2, 6);
    close(Vec3.transformMat4(Vec3(), p.center, next.matrix), [25, 40, 60]);
  });

  it('is opt-in, and lists the structures it can move', async () => {
    const plain = await createPrismViewerHeadless({ spec: DefaultPrismViewerSpec() });
    expect(plain.api.move).toBeUndefined();
    await plain.dispose();

    const base = DefaultPrismViewerSpec();
    viewer = await createPrismViewerHeadless({
      spec: { ...base, extensions: [...(base.extensions ?? []), StructureMoveFeature] },
    });
    expect((await viewer.api.load.files({ files: [cbs()] })).ok).toBe(true);
    const move = viewer.api.move!;
    expect(move.structures.value.map((s) => s.label)).toEqual(['1CBS']);
    // headless: no canvas, no handles
    expect(move.showGizmos(viewer.api)).toBe(false);
  });
});
