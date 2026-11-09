/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { Sphere3D } from '@molstar/core/math/geometry';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import { LineStride, type Lines } from '@molstar/graphics/geo/geometry/lines/lines';
import type { Mesh } from '@molstar/graphics/geo/geometry/mesh/mesh';
import { PositionGroupStride } from '@molstar/graphics/geo/geometry/center-group';
import type { Points } from '@molstar/graphics/geo/geometry/points/points';
import { CutShapes, type CutObject } from '../../cuts.js';
import type { GizmoPart } from '@molstar/graphics/util/gizmo/controller';
import { gizmoGeometry, MoveScreenSize } from '../geometry.js';

const cut = (type: CutObject['type']): CutObject => ({
  type,
  invert: false,
  position: Vec3.create(0, 0, 0),
  rotation: { axis: Vec3.create(1, 0, 0), angle: 0 },
  scale: Vec3.create(10, 10, 10),
  transform: undefined as never,
});
/** What the cut cuts: 20 Å of reach */
const Bounds = Sphere3D.create(Vec3(), 20);

function groupsOf(mesh: Mesh) {
  return new Set(Array.from(mesh.groupBuffer.ref.value));
}

/** The groups of a part's geometry */
function partGroupsOf(part: GizmoPart) {
  const g = part.geometry;
  if (g.kind === 'mesh') return groupsOf(g);
  if (g.kind === 'lines') return lineGroupsOf(g);
  return pointGroupsOf(g);
}

function lineGroupsOf(lines: Lines | undefined) {
  const groups = new Set<number>();
  if (!lines) return groups;
  const { array } = lines.lineData.ref.value;
  for (let i = 0; i < lines.lineCount; i++) groups.add(array[i * LineStride + 3]);
  return groups;
}

/** The largest distance of a mesh's vertices from the origin */
function extent(mesh: Mesh) {
  const v = mesh.vertexBuffer.ref.value;
  let max = 0;
  for (let i = 0; i < mesh.vertexCount; i++) max = Math.max(max, Math.hypot(v[i * 3], v[i * 3 + 1], v[i * 3 + 2]));
  return max;
}

/** A point's center and group: `PositionGroupStride` floats a point, xyz and the group */
function pointAt(points: Points, i: number) {
  const a = points.positionGroup.ref.value.array;
  const o = i * PositionGroupStride;
  return { center: [a[o], a[o + 1], a[o + 2]], group: a[o + 3] };
}

function pointGroupsOf(points: Points | undefined) {
  const groups = new Set<number>();
  if (points) for (let i = 0; i < points.pointCount; i++) groups.add(pointAt(points, i).group);
  return groups;
}

describe('the clipping gizmos geometry', () => {
  // the transform handles (Blender's gizmo) by cut: what changes it
  const kinds = (type: CutObject['type'], scale?: Vec3) => {
    const g = gizmoGeometry({ ...cut(type), ...(scale ? { scale } : {}) }, Bounds);
    const move = new Set(g.moveParts.flatMap((p) => [...partGroupsOf(p)]));
    return g.groups
      .filter((_, i) => move.has(i))
      .map(({ handle }) => {
        if (handle.kind !== 'rotate' && handle.kind !== 'slide') return handle.kind;
        return `${handle.kind} ${'xyz'[[0, 1, 2].find((a) => Math.abs(handle.axis[a]) > 0.5)!]}`;
      });
  };

  it("gives each cut the move and turn handles that change it, a plane's center a slide along its normal", () => {
    const all = [
      'move',
      'slide x',
      'slide y',
      'slide z',
      'plane',
      'plane',
      'plane',
      'rotate x',
      'rotate y',
      'rotate z',
    ];
    expect(kinds('sphere')).toEqual([...all, 'roll', 'turn']);
    expect(kinds('cube')).toEqual([...all, 'roll', 'turn']);
    // no squares, one arrow: the plane's normal; its rings about the two axes in it
    expect(kinds('plane')).toEqual(['slide y', 'slide y', 'rotate x', 'rotate z', 'roll', 'turn']);
    // no ring about their own axis (y, a cone's opening z), round about it
    expect(kinds('cylinder')).toEqual([...all.filter((k) => k !== 'rotate y'), 'roll', 'turn']);
    // an elliptic cylinder's radii (x, z) turn with it: its ring too
    expect(kinds('cylinder', Vec3.create(10, 20, 30))).toEqual([...all, 'roll', 'turn']);
    expect(kinds('infiniteCone')).toEqual([...all.filter((k) => k !== 'rotate z'), 'roll', 'turn']);
    // a plane's arrow toward the side it cuts away
    const flipped = gizmoGeometry({ ...cut('plane'), invert: true }, Bounds);
    expect((flipped.groups[1].handle as { axis: Vec3 }).axis[1]).toBe(-1);
  });

  // every handle in the move parts (at a fixed size, picked through a hit part) or the shape (lines and points)
  for (const type of CutShapes) {
    it(`gives a ${type} each handle in a part, the move handles first, each picked`, () => {
      const g = gizmoGeometry(cut(type), Bounds);
      const move = new Set(g.moveParts.flatMap((p) => [...partGroupsOf(p)]));
      const hit = new Set(g.moveParts.filter((p) => p.role === 'hit').flatMap((p) => [...partGroupsOf(p)]));
      const shape = new Set([...lineGroupsOf(g.handleLines), ...pointGroupsOf(g.handlePoints)]);
      expect(new Set([...move, ...shape])).toEqual(new Set(g.groups.map((_, i) => i)));
      expect(hit).toEqual(move);
      for (const i of shape) expect(move.has(i)).toBe(false);
      expect(Math.max(...move)).toBeLessThan(Math.min(...shape, Infinity));
    });
  }

  // the move handles in units of their size on screen, the shape where the clip is: zoomed out on a large structure
  // the first stay as large on screen, the second as large as the clip
  it('sizes the move handles on screen and puts the shape where the clip is', () => {
    const small = gizmoGeometry(cut('cube'), Bounds);
    const large = gizmoGeometry(cut('cube'), Sphere3D.create(Vec3(), 200));
    for (const g of [small, large]) for (const p of g.moveParts) expect(p.screenSize).toBe(MoveScreenSize);
    const extentOf = (g: typeof small) =>
      Math.max(...g.moveParts.map((p) => (p.geometry.kind === 'mesh' ? extent(p.geometry) : 0)));
    expect(extentOf(large)).toBeCloseTo(extentOf(small), 6);
    // the box's 12 edges, its 6 face knobs and 8 corner knobs at its size (10 Å), whatever it cuts
    expect(small.handleLines!.lineCount).toBe(12);
    const points = small.handlePoints!;
    expect(points.pointCount).toBe(14);
    let corners = 0;
    for (let i = 0; i < points.pointCount; i++) {
      if (pointAt(points, i).center.every((x) => Math.abs(Math.abs(x) - 5) < 1e-6)) corners += 1;
    }
    expect(corners).toBe(8);
  });

  // the sphere is an ellipsoid (the shader's): the box's handles on its surface, each axis its own size
  it("puts an ellipsoid's handles on its surface: tips, sections, octant knobs", () => {
    const g = gizmoGeometry({ ...cut('sphere'), scale: Vec3.create(10, 20, 30) }, Bounds);
    const half = [5, 10, 15];
    const points = g.handlePoints!;
    expect(points.pointCount).toBe(6 + 8);
    for (let i = 0; i < points.pointCount; i++) {
      const { center, group } = pointAt(points, i);
      // on the surface
      expect(Math.hypot(center[0] / half[0], center[1] / half[1], center[2] / half[2])).toBeCloseTo(1, 4);
      const handle = g.groups[group].handle;
      if (handle.kind === 'faces') {
        // a tip, moved along its axis
        const { axis, sign } = handle.faces[0];
        expect(center[axis]).toBeCloseTo(sign * half[axis], 4);
      } else {
        // an octant's middle, sizing it whole
        expect(handle).toEqual({ kind: 'scale' });
        for (const [a, x] of center.entries()) expect(Math.abs(x)).toBeCloseTo(half[a] / Math.sqrt(3), 4);
      }
    }
    // three sections, each sizing the two axes it lies across, every line on its ellipse
    const { array } = g.handleLines!.lineData.ref.value;
    const sections = new Map<number, [number, number]>();
    for (const [i, h] of g.groups.entries()) {
      if (h.handle.kind === 'scale' && h.handle.axes) sections.set(i, h.handle.axes as [number, number]);
    }
    expect([...sections.values()]).toEqual([
      [0, 2],
      [0, 1],
      [1, 2],
    ]);
    for (let i = 0; i < g.handleLines!.lineCount; i++) {
      const o = i * LineStride;
      const [b, c] = sections.get(array[o + 3])!;
      const p = [array[o], array[o + 1], array[o + 2]];
      const off = [0, 1, 2].find((a) => a !== b && a !== c)!;
      expect(p[off]).toBeCloseTo(0, 4);
      expect(Math.hypot(p[b] / half[b], p[c] / half[c])).toBeCloseTo(1, 4);
    }
  });

  // the cylinder is elliptic (the shader's): its rims ellipses, a knob on each side moving it across its axis, one on
  // each cap along it
  it("puts an elliptic cylinder's rims and side knobs where its sides are", () => {
    const g = gizmoGeometry({ ...cut('cylinder'), scale: Vec3.create(10, 20, 30) }, Bounds);
    const [a, h, c] = [5, 20, 15];
    const { array } = g.handleLines!.lineData.ref.value;
    for (let i = 0; i < g.handleLines!.lineCount; i++) {
      const o = i * LineStride;
      expect(Math.abs(array[o + 1])).toBeCloseTo(h / 2, 4);
      expect(Math.hypot(array[o] / a, array[o + 2] / c)).toBeCloseTo(1, 4);
    }
    const points = g.handlePoints!;
    const sides: string[] = [];
    for (let i = 0; i < points.pointCount; i++) {
      const { center, group } = pointAt(points, i);
      const handle = g.groups[group].handle;
      if (handle.kind !== 'faces') continue;
      const { axis, sign } = handle.faces[0];
      // a cap's center, or a side halfway up
      expect(center[axis]).toBeCloseTo(sign * [a, h / 2, c][axis], 4);
      sides.push(`${sign > 0 ? '+' : '-'}${'xyz'[axis]}`);
    }
    expect(sides.sort()).toEqual(['+x', '+y', '+z', '-x', '-y', '-z']);
  });
});
