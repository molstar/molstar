/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Vec3 } from '@molstar/core/math/linear-algebra';
import { describe, expect, it } from 'vitest';
import { LineStride, type Lines } from '../../../geo/geometry/lines/lines.js';
import type { Mesh } from '../../../geo/geometry/mesh/mesh.js';
import type { Points } from '../../../geo/geometry/points/points.js';
import type { GizmoPart } from '../controller.js';
import { CenterRadius, RingRadius, RollRadius, transformHandles } from '../transform.js';

// Blender's move and rotate gizmo combined: which part holds which handle, how each is shown, what wins the pick.

/** The groups a part's geometry holds */
function groupsOf(part: GizmoPart | Mesh | Lines | Points) {
  const g = 'role' in part ? part.geometry : part;
  if (g.kind === 'lines') {
    const { array } = g.lineData.ref.value;
    const out = new Set<number>();
    for (let i = 0; i < g.lineCount; i++) out.add(array[i * LineStride + 3]);
    return out;
  }
  if (g.kind === 'mesh') return new Set(Array.from(g.groupBuffer.ref.value));
  const { array } = g.positionGroup.ref.value;
  return new Set(Array.from({ length: g.pointCount }, (_, i) => array[i * 4 + 3]));
}

/** The vertices of a mesh, as points */
function vertices(m: Mesh) {
  const v = m.vertexBuffer.ref.value;
  return Array.from({ length: m.vertexCount }, (_, i) => Vec3.create(v[i * 3], v[i * 3 + 1], v[i * 3 + 2]));
}

// groups: center 0, slides 1-3, planes 4-6, rings 7-9, roll 10, trackball 11
describe('transform handles', () => {
  it('has the center, arrows, squares, rings, roll and trackball, in that order', () => {
    const { groups } = transformHandles(1);
    expect(groups.map((g) => g.handle.kind)).toEqual([
      'move',
      'slide',
      'slide',
      'slide',
      'plane',
      'plane',
      'plane',
      'rotate',
      'rotate',
      'rotate',
      'roll',
      'turn',
    ]);
  });

  it('puts each handle in its part: heads and squares, shafts while moving, front rings, circles facing the camera', () => {
    const { parts, handles } = transformHandles(1);
    const drawn = parts.filter((p) => p.role === 'handle');
    const find = (f: (p: GizmoPart) => boolean) => drawn.find(f)!;
    expect(groupsOf(handles)).toEqual(new Set([1, 2, 3, 4, 5, 6]));
    // the shafts: shown during a move only
    const shafts = find((p) => !!p.when?.drag);
    expect(groupsOf(shafts)).toEqual(new Set([1, 2, 3]));
    expect(shafts.when!.drag).toEqual(['move', 'slide', 'plane']);
    expect(shafts.geometry.kind).toBe('mesh');
    // the rings, their front halves
    const rings = find((p) => !!p.front);
    expect(groupsOf(rings)).toEqual(new Set([7, 8, 9]));
    // the center and the roll circle facing the camera: ring points
    const circles = drawn.filter((p) => p.point?.style === 'ring');
    expect(circles.flatMap((p) => [...groupsOf(p)])).toEqual([0, 10]);
    expect(circles.every((p) => p.geometry.kind === 'points')).toBe(true);
    // the trackball's fill: a disc on hover, faint
    const fill = find((p) => !!p.when?.hover);
    expect(fill.when!.hover).toBe(11);
    expect(fill.alpha).toBeLessThan(0.5);
    expect(fill.point?.style ?? 'circle').toBe('circle');
    // every handle picked through a hit part
    const hit = new Set(parts.filter((p) => p.role === 'hit').flatMap((p) => [...groupsOf(p)]));
    expect(hit).toEqual(new Set([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]));
    // the rings' hit area their front half too
    expect(parts.find((p) => p.role === 'hit' && p.front)).toBeDefined();
  });

  it('draws the arrows as heads past the rings and the circles at their radii', () => {
    const { handles, parts } = transformHandles(2);
    // the heads from 0.85 to 1 of their length
    const x = vertices(handles).filter((v) => Math.abs(v[1]) < 0.3 && Math.abs(v[2]) < 0.3 && v[0] > 0.9);
    expect(Math.min(...x.map((v) => v[0]))).toBeCloseTo(2 * 0.85, 4);
    expect(Math.max(...x.map((v) => v[0]))).toBeCloseTo(2, 4);
    expect(2 * RingRadius).toBeLessThan(2 * 0.85);
    // the circles facing the camera: ring points whose band's middle is the center's and the roll's radius
    const circles = parts.filter((p) => p.role === 'handle' && p.point?.style === 'ring');
    const middles = circles.map((p) => p.point!.radius * (1 - p.point!.ringWidth! / 2));
    expect(middles[0]).toBeCloseTo(2 * CenterRadius, 6);
    expect(middles[1]).toBeCloseTo(2 * RollRadius, 6);
  });

  // the hit discs moved in depth (toward the camera, away): the center's over every part, the trackball's under
  it("picks the center over every part, and every part over the trackball's area", () => {
    const { parts, groups } = transformHandles(2, { scale: true });
    const hitOf = (group: number) =>
      parts.find((p) => p.role === 'hit' && p.point?.style !== 'ring' && groupsOf(p).has(group))!;
    // its disc past every other part (within 1.16, a scale knob, of the origin)
    expect(hitOf(0).point!.depth).toBeGreaterThan(2 * 1.16);
    const trackball = groups.findIndex((h) => h.handle.kind === 'turn');
    expect(hitOf(trackball).point!.depth).toBeLessThan(-2 * 1.16);
  });

  it('takes the axes it is given: the arrows, the rings, a sliding center', () => {
    const up = Vec3.create(0, -1, 0);
    const { groups } = transformHandles(1, { arrows: [up], rings: [0, 2], planes: false, centerAxis: up });
    expect(groups.map((g) => g.handle.kind)).toEqual(['slide', 'slide', 'rotate', 'rotate', 'roll', 'turn']);
    expect(groups[1].handle).toEqual({ kind: 'slide', axis: up });
    expect(groups.slice(2, 4).map((g) => (g.handle as { axis: Vec3 }).axis[0])).toEqual([1, 0]);
    // no rotation: no ring, roll or trackball, nothing facing the camera but the center
    const still = transformHandles(1, { rotate: false });
    expect(still.groups.some((g) => ['rotate', 'roll', 'turn'].includes(g.handle.kind))).toBe(false);
    expect(still.parts.filter((p) => p.front)).toHaveLength(0);
  });
});
