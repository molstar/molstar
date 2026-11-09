/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Ludovic Autin <ludovic.autin@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import { Color } from '@molstar/core/util/color';
import type { Mesh } from '../../geo/geometry/mesh/mesh.js';
import { MeshBuilder } from '../../geo/geometry/mesh/mesh-builder.js';
import { PointsBuilder } from '../../geo/geometry/points/points-builder.js';
import { Box } from '../../geo/primitive/box.js';
import type { GizmoPart } from './controller.js';
import { AxisColors, type Handle, type HandleGroup } from './handle.js';
import { arrowHead, arrowShaftEnd, ring, tube } from './primitives.js';

// The standard transform handles, as Blender's move and rotate gizmo combined, in the item's frame (`size` the arrows'
// length; X red, Y green, Z blue): the center, a small circle facing the camera (move in the screen's plane); an arrow
// head on each axis past the rings (slide along it; its shaft drawn only while it or another move is dragged); a small
// square in each pair of axes' plane (move in it); each axis's ring, its half toward the viewer (turn about it); a grey
// circle facing the camera around them (roll: turn about the view axis); the area inside the rings, not drawn (a faint
// fill on hover: turn as the camera's trackball turns the scene); a knob past each arrow's tip sizes along its axis.
// What faces the camera is points (turned to it by nature, a set size on screen: the circles `ring` points, the discs
// `circle` points). What is picked over what: the center's hit disk in front of every part (a ring seen edge-on crosses
// it, a depth offset), the trackball's behind every part.

export interface TransformHandlesOptions {
  /** The center and the arrows (by default) */
  readonly move?: boolean;
  /** The center slides along this local axis instead of moving in the screen's plane (a plane cut's normal) */
  readonly centerAxis?: Vec3;
  /** The arrows (local directions, each sliding along its own): the three axes by default */
  readonly arrows?: readonly Vec3[];
  /** The plane squares (by default) */
  readonly planes?: boolean;
  /** The rings (by default) */
  readonly rotate?: boolean;
  /** The local axes the rings turn about: the three by default */
  readonly rings?: readonly (0 | 1 | 2)[];
  /** The roll circle (with `rotate` by default) */
  readonly roll?: boolean;
  /** Free rotation, the trackball's area (with `rotate` by default) */
  readonly free?: boolean;
  /** The scale knobs (not by default) */
  readonly scale?: boolean;
}

const Axes = [Vec3.create(1, 0, 0), Vec3.create(0, 1, 0), Vec3.create(0, 0, 1)] as const;
const Names = ['X', 'Y', 'Z'] as const;

/** The hit area's rings and circles, times the drawn ones' thickness */
const HitThickness = 4;

/** Radii, times `size`: the center circle, the rings, the roll circle */
export const CenterRadius = 0.1;
export const RingRadius = 0.65;
export const RollRadius = 0.8;
/** An arrow head from here to 1 (times `size`) */
const HeadFrom = 0.85;

/** The center's hit disk, this far toward the camera (times `size`): in front of every part, a scale knob's included */
const CenterHitDepth = 1.25;
/** The trackball's area, this far behind the origin (times `size`): behind every part, a scale knob's included */
const FreeAreaDepth = 1.25;
/** The trackball's fill on hover */
const FillAlpha = 0.35;

/** The circles' color: grey, not Blender's white (white on a white background or structure is lost) */
const CircleColor = Color(0x8c8c8c);

/** The handles moving the item: the arrows' shafts are drawn while one of them is dragged */
const MoveKinds: Handle['kind'][] = ['move', 'slide', 'plane'];

/** One point at `at` (the frame's origin by default) of `group`: a part facing the camera (`GizmoPart.point`) */
function pointAt(group: number, at: Vec3 = Vec3()) {
  const points = PointsBuilder.create(1, 1);
  points.add(at[0], at[1], at[2], group);
  return points.getPoints();
}

/**
 * The transform handles' parts (`GizmoLayout.parts`; one group a handle) and what each does: `handles` the heads and
 * squares, `hitArea` theirs thick; the parts facing the camera and the knobs points
 */
export function transformHandles(
  size: number,
  options: TransformHandlesOptions = {},
): { handles: Mesh; groups: HandleGroup[]; hitArea: Mesh; parts: GizmoPart[] } {
  const { move = true, planes = true, rotate = true, scale = false, roll = rotate, free = rotate } = options;
  const groups: HandleGroup[] = [];
  const group = (handle: Handle, label: string, color: Color) => {
    groups.push({ handle, label, color });
    return groups.length - 1;
  };
  // the heads and squares; their hit area
  const heads = MeshBuilder.createState(4096, 2048);
  const headsHit = MeshBuilder.createState(4096, 2048);
  // the shafts, shown while moving
  const shaftMesh = MeshBuilder.createState(1024, 512);
  // the rings, their front halves; their hit area
  const ringMesh = MeshBuilder.createState(4096, 2048);
  const ringsHit = MeshBuilder.createState(4096, 2048);
  /** the parts facing the camera and the knobs: points, each its own radius and look */
  const points: GizmoPart[] = [];
  const at = (states: MeshBuilder.State[], g: number) => {
    for (const s of states) s.currentGroup = g;
  };
  const width = size * 0.024;
  const shaft = size * 0.025;
  // a circle facing the camera, of radius `r`: a `ring` point (a band `width` either side of `r`), its hit band wider
  const circle = (r: number, g: number) => {
    const drawn = r + width;
    points.push({
      role: 'handle',
      geometry: pointAt(g),
      point: { radius: drawn, style: 'ring', ringWidth: (2 * width) / drawn },
    });
    const hit = r + width * HitThickness;
    points.push({
      role: 'hit',
      geometry: pointAt(g),
      point: { radius: hit, style: 'ring', ringWidth: (2 * width * HitThickness) / hit },
    });
  };

  if (move) {
    const centerAxis = options.centerAxis;
    const center = group(
      centerAxis ? { kind: 'slide', axis: centerAxis } : { kind: 'move' },
      centerAxis ? 'Drag to slide it' : 'Drag to move it',
      CircleColor,
    );
    circle(size * CenterRadius, center);
    // its hit a disk in front of every part
    points.push({
      role: 'hit',
      geometry: pointAt(center),
      point: { radius: size * CenterRadius * 1.4, depth: size * CenterHitDepth },
    });
    for (const dir of options.arrows ?? Axes) {
      const a = Vec3.normalize(Vec3(), dir);
      const axis = dominantAxis(a);
      const g = group({ kind: 'slide', axis: a }, `Drag to move it along ${Names[axis]}`, AxisColors[axis]);
      at([heads, headsHit, shaftMesh], g);
      // `arrowHead` draws from the shaft's end (`arrowShaftEnd`: 0.8 of the way) to the tip
      const tip = Vec3.scale(Vec3(), a, size);
      const from = Vec3.scale(Vec3(), a, size * (1 - (1 - HeadFrom) / 0.2));
      arrowHead(heads, from, tip, shaft);
      arrowHead(headsHit, from, tip, shaft * 2);
      const start = Vec3.scale(Vec3(), a, size * CenterRadius);
      const end = arrowShaftEnd(from, tip);
      tube(shaftMesh, start, end, width);
    }
  }
  if (planes) {
    // the plane across axis i, between the other two: a small square near the center
    for (let i = 0; i < 3; i++) {
      const [b, c] = [0, 1, 2].filter((x) => x !== i);
      const g = group({ kind: 'plane', normal: Axes[i] }, `Drag to move it in ${Names[b]}${Names[c]}`, AxisColors[i]);
      at([heads, headsHit], g);
      const center = Vec3();
      center[b] = center[c] = size * 0.28;
      const extent = Vec3.create(size * 0.1, size * 0.1, size * 0.1);
      extent[i] = size * 0.01;
      const square = (state: MeshBuilder.State, grow: number, thick: number) => {
        const e = Vec3.scale(Vec3(), extent, grow);
        e[i] = thick;
        MeshBuilder.addPrimitive(
          state,
          Mat4.mul(Mat4(), Mat4.fromTranslation(Mat4(), center), Mat4.fromScaling(Mat4(), e)),
          Box(),
        );
      };
      square(heads, 1, size * 0.01);
      square(headsHit, 1.3, size * 0.06);
    }
  }
  if (rotate) {
    // a ring around axis i: around +Y as built, turned onto X (about Z) and Z (about X)
    const onto = [Mat4.rotZ90 as Mat4, Mat4.identity(), Mat4.rotX90 as Mat4];
    for (const i of options.rings ?? ([0, 1, 2] as const)) {
      const g = group({ kind: 'rotate', axis: Axes[i] }, `Drag to turn it about ${Names[i]}`, AxisColors[i]);
      at([ringMesh, ringsHit], g);
      const t = Mat4.mul(Mat4(), onto[i], Mat4.rotX90 as Mat4);
      ring(ringMesh, size * RingRadius, width, 0, t);
      ring(ringsHit, size * RingRadius, width * HitThickness, 0, t);
    }
  }
  if (scale) {
    for (let i = 0; i < 3; i++) {
      const g = group({ kind: 'scale', axes: [i as 0 | 1 | 2] }, `Drag to size it along ${Names[i]}`, AxisColors[i]);
      // a knob: in front of what meets it, drawn and picked (`front`)
      const p = Vec3.scale(Vec3(), Axes[i], size * 1.12);
      points.push({ role: 'handle', geometry: pointAt(g, p), point: { radius: size * 0.04, front: true } });
      points.push({ role: 'hit', geometry: pointAt(g, p), point: { radius: size * 0.06, front: true } });
    }
  }
  if (roll) {
    const g = group({ kind: 'roll' }, 'Drag to turn it about the view', CircleColor);
    circle(size * RollRadius, g);
  }
  let trackball = -1;
  if (free) {
    trackball = group({ kind: 'turn' }, 'Drag to turn it freely', CircleColor);
    // the area inside the rings, behind every part; its fill just behind the center
    points.push({
      role: 'hit',
      geometry: pointAt(trackball),
      point: { radius: size * RingRadius, depth: -size * FreeAreaDepth },
    });
    points.push({
      role: 'handle',
      geometry: pointAt(trackball),
      point: { radius: size * RingRadius, depth: -size * 0.01 },
      alpha: FillAlpha,
      when: { hover: trackball },
    });
  }

  const handles = MeshBuilder.getMesh(heads);
  const hitArea = MeshBuilder.getMesh(headsHit);
  const parts: GizmoPart[] = [
    { role: 'handle', geometry: handles },
    { role: 'hit', geometry: hitArea },
  ];
  if (move) {
    parts.push({ role: 'handle', geometry: MeshBuilder.getMesh(shaftMesh), when: { drag: MoveKinds } });
  }
  if (rotate) {
    parts.push(
      { role: 'handle', geometry: MeshBuilder.getMesh(ringMesh), front: true },
      { role: 'hit', geometry: MeshBuilder.getMesh(ringsHit), front: true },
    );
  }
  parts.push(...points);
  return { handles, groups, hitArea, parts };
}

/** The axis a unit direction is nearest to */
function dominantAxis(a: Vec3): 0 | 1 | 2 {
  const [x, y, z] = [Math.abs(a[0]), Math.abs(a[1]), Math.abs(a[2])];
  return x >= y && x >= z ? 0 : y >= z ? 1 : 2;
}
