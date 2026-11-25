/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Mat4, Quat, Vec3 } from '@molstar/core/math/linear-algebra';
import { afterEach, describe, expect, it } from 'vitest';
import { cbs } from '../../../_test/harness.js';
import { mapFile, particlesFile } from '../../../api/_test/volume-fixtures.js';
import { MeshFiles } from '../../../testing/meshes.js';
import { DefaultPrismViewerSpec } from '../../../default-spec.js';
import { createPrismViewerHeadless } from '../../../headless.js';
import { applyToPlacement, placementPose, poseMatrix, type Placement } from '../../../sources/move/gizmos.js';
import type { PrismViewer } from '../../../viewer.js';
import { MoveFeature } from '../index.js';

// A model moved by hand: its placement as the handles' pose and back (a turn about its center, then a move), and the
// feature in a viewer (opt-in: not in the default spec).

const close = (a: ArrayLike<number>, b: ArrayLike<number>) =>
  Array.from(b).forEach((x, i) => expect(a[i]).toBeCloseTo(x, 6));

const placement = (matrix = Mat4.identity()): Placement => ({
  matrix,
  center: Vec3.create(10, 20, 30),
});

let viewer: PrismViewer | undefined;
afterEach(async () => {
  await viewer?.dispose();
  viewer = undefined;
});

describe('move', () => {
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

  it('is opt-in, and lists the structures, particle lists, shapes and volumes it can move', async () => {
    const plain = await createPrismViewerHeadless({ spec: DefaultPrismViewerSpec() });
    expect(plain.api.move).toBeUndefined();
    await plain.dispose();

    const base = DefaultPrismViewerSpec();
    viewer = await createPrismViewerHeadless({
      spec: { ...base, extensions: [...(base.extensions ?? []), MoveFeature] },
    });
    expect((await viewer.api.load.files({ files: [cbs()] })).ok).toBe(true);
    const move = viewer.api.move!;
    expect(move.movable.value.map((s) => [s.label, s.kind])).toEqual([['1CBS', 'structure']]);
    // headless: no canvas, no handles
    expect(move.showGizmos(viewer.api)).toBe(false);

    // a particle list too, after the structures; Reset by its kind
    expect((await viewer.api.load.files({ files: [particlesFile()] })).ok).toBe(true);
    const [, list] = move.movable.value;
    expect(list.kind).toBe('particles');
    expect(move.gizmos.kindOf(list.ref)).toBe('particles');
    expect((await viewer.api.run('particles.transform', { particles: list.ref, translation: [1, 2, 3] })).ok).toBe(
      true,
    );
    expect((await move.reset(viewer.api, list.ref)).ok).toBe(true);

    // a shape and a volume, after them; each moved and put back by its kind
    expect((await viewer.api.load.files({ files: [MeshFiles.torus(), mapFile()] })).ok).toBe(true);
    const kinds = move.movable.value.map((s) => s.kind);
    expect(kinds).toEqual(['structure', 'particles', 'shape', 'volume']);
    for (const m of move.movable.value.slice(2)) {
      expect(move.gizmos.kindOf(m.ref)).toBe(m.kind);
      const arg = m.kind === 'shape' ? 'node' : 'volume';
      expect((await viewer.api.run(`${m.kind}.transform`, { [arg]: m.ref, translation: [1, 2, 3] })).ok).toBe(true);
      expect((await move.reset(viewer.api, m.ref)).ok).toBe(true);
    }
    // no new rows for the transforms
    expect(move.movable.value.map((s) => s.kind)).toEqual(kinds);
  });
});
