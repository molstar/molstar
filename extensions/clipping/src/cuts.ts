/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Sphere3D } from '@molstar/core/math/geometry';
import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import { type State, type StateObjectCell, StateTransform } from '@molstar/core/state';
import type { Clip } from '@molstar/core/util/clip';
import type { PluginContext } from '@molstar/plugin/context';
import { RepresentationClip } from '@molstar/plugin/state/helpers/representation-clip';

// Cuts (_spec/README.md §1): the clip objects of the representations, read from the state tree (its record of what is
// clipped, so a clip set by code or a file shows too) and grouped: equal objects on several representations are one
// cut. A cut has no id of its own in the state; its number is its place in the list (tree order of first use).

export type CutObject = RepresentationClip.Object;

export const CutShapes = ['plane', 'sphere', 'cube', 'cylinder', 'infiniteCone'] as const;
export type CutShape = (typeof CutShapes)[number];

export const CutShapeLabels: Readonly<Record<CutShape, string>> = {
  plane: 'Plane',
  sphere: 'Sphere',
  cube: 'Box',
  cylinder: 'Cylinder',
  infiniteCone: 'Cone',
};

export interface CutTarget {
  /** The representation cell's ref */
  readonly ref: string;
  /** "Cartoon of Polymer" */
  readonly label: string;
  /** Its representation type (`cartoon`, `isosurface`), `shape` for a shape's */
  readonly type: string;
}

export interface Cut {
  /** Its number, from 1, in the list's order */
  readonly index: number;
  readonly object: CutObject;
  readonly targets: readonly CutTarget[];
}

/** The representations whose geometry can be clipped, in tree order (a walk from the root) */
export function clippableCells(state: State): StateObjectCell[] {
  const out: StateObjectCell[] = [];
  const walk = (ref: string) => {
    const cell = state.cells.get(ref);
    if (cell && RepresentationClip.isClippable(cell)) out.push(cell);
    for (const child of state.tree.children.get(ref).values()) walk(child);
  };
  walk(state.tree.root.ref);
  return out;
}

export function representationType(cell: StateObjectCell): string {
  return (cell.transform.params as { type?: { name?: string } } | undefined)?.type?.name ?? 'shape';
}

/** A representation in words: its label and what it draws ("Cartoon of Polymer") */
export function representationLabel(state: State, cell: StateObjectCell): string {
  const own = cell.obj?.label ?? representationType(cell);
  for (let ref = cell.transform.parent; ref && ref !== StateTransform.RootRef; ) {
    const parent = state.cells.get(ref);
    if (!parent) break;
    if (!parent.transform.transformer.definition.isDecorator && parent.obj?.label)
      return `${own} of ${parent.obj.label}`;
    ref = parent.transform.parent;
  }
  return own;
}

/** The objects a representation's clip holds that are not the component manager's (its "clip everything" objects) */
function ownObjects(plugin: PluginContext, clip: Clip.Props): CutObject[] {
  const global = plugin.managers.structure.component.state.options.clipObjects.objects;
  return clip.objects.filter((o) => !global.some((g) => RepresentationClip.objectsEqual(o, g)));
}

/** The cuts in the scene: equal objects on several representations are one cut */
export function findCuts(plugin: PluginContext): Cut[] {
  const state = plugin.state.data;
  const cuts: { object: CutObject; targets: CutTarget[] }[] = [];
  for (const cell of clippableCells(state)) {
    const clip = RepresentationClip.get(cell.transform.params)!;
    for (const object of ownObjects(plugin, clip)) {
      let cut = cuts.find((c) => RepresentationClip.objectsEqual(c.object, object));
      if (!cut) cuts.push((cut = { object, targets: [] }));
      if (!cut.targets.some((t) => t.ref === cell.transform.ref)) {
        cut.targets.push({
          ref: cell.transform.ref,
          label: representationLabel(state, cell),
          type: representationType(cell),
        });
      }
    }
  }
  return cuts.map((c, i) => ({ index: i + 1, ...c }));
}

/** The bounding sphere of what the representations draw (their render objects, their transforms included) */
export function representationsSphere(cells: readonly StateObjectCell[]): Sphere3D | undefined {
  let out: Sphere3D | undefined;
  for (const cell of cells) {
    const repr = (cell.obj?.data as { repr?: { renderObjects: readonly { values: any }[] } } | undefined)?.repr;
    for (const ro of repr?.renderObjects ?? []) {
      const s = ro.values.boundingSphere?.ref.value as Sphere3D | undefined;
      if (!s || s.radius <= 0) continue;
      out = out ? Sphere3D.expandBySphere(Sphere3D(), out, s) : Sphere3D.clone(s);
    }
  }
  return out;
}

/** The rotation (axis, degrees) that turns `from` onto `to`; half a turn about a perpendicular when they are opposite */
export function rotationTo(from: Vec3, to: Vec3): Clip.Rotation {
  const a = Vec3.normalize(Vec3(), from);
  const b = Vec3.normalize(Vec3(), to);
  const axis = Vec3.cross(Vec3(), a, b);
  const dot = Vec3.dot(a, b);
  if (Vec3.magnitude(axis) < 1e-6) {
    // parallel: none; opposite: 180° about any axis perpendicular to `from` (a plane's normal facing down, H005)
    const perpendicular = Math.abs(a[0]) < 0.9 ? Vec3.create(1, 0, 0) : Vec3.create(0, 0, 1);
    const p = Vec3.normalize(Vec3(), Vec3.cross(Vec3(), a, perpendicular));
    return dot > 0 ? { axis: Vec3.create(1, 0, 0), angle: 0 } : { axis: p, angle: 180 };
  }
  return { axis: Vec3.normalize(axis, axis), angle: (Math.acos(Math.max(-1, Math.min(1, dot))) * 180) / Math.PI };
}

/** A direction turned by a clip object's rotation (its local axis in the world: a plane's normal is local +Y) */
export function rotated(rotation: Clip.Rotation, local: Vec3): Vec3 {
  const m = Mat4.fromRotation(Mat4(), (rotation.angle * Math.PI) / 180, rotation.axis);
  return Vec3.transformDirection(Vec3(), local, m);
}

export const LocalUp = Vec3.create(0, 1, 0);
/** An infinite cone opens along its local −Z */
export const ConeOpening = Vec3.create(0, 0, -1);

/** The half-angle (degrees) of a cone with this scale, and back (shader: `size.xy`, the opening's cosine and sine) */
export const coneAngle = (scale: Vec3) => (Math.atan2(scale[1], scale[0]) * 180) / Math.PI;
export const coneScale = (angle: number) =>
  Vec3.create(2 * Math.cos((angle * Math.PI) / 180), 2 * Math.sin((angle * Math.PI) / 180), 1);

/**
 * A new cut of a shape over what it is for (`sphere`), seen from `view` (toward the viewer): a plane through the center
 * facing the viewer (it removes the near half, so the inside shows); closed shapes on the near side, where they cut
 * something seen (one centered in a hollow capsid cuts nothing); a cylinder drilled toward the viewer; a cone from the
 * center opening toward the viewer.
 */
export function placeCut(shape: CutShape, sphere: Sphere3D, view: Vec3): CutObject {
  const { center, radius: r } = sphere;
  const toward = Vec3.normalize(Vec3(), Vec3.magnitude(view) > 0 ? view : Vec3.create(0, 0, 1));
  const at = (k: number) => Vec3.scaleAndAdd(Vec3(), center, toward, k * r);
  const base = { type: shape, invert: false, transform: Mat4.identity() };
  switch (shape) {
    case 'plane':
      return {
        ...base,
        position: Vec3.clone(center),
        rotation: rotationTo(LocalUp, toward),
        scale: Vec3.create(1, 1, 1),
      };
    case 'sphere':
      return {
        ...base,
        position: at(0.85),
        rotation: { axis: Vec3.create(1, 0, 0), angle: 0 },
        scale: Vec3.create(r, r, r),
      };
    case 'cube':
      return { ...base, position: at(0.85), rotation: rotationTo(LocalUp, toward), scale: Vec3.create(r, r / 2, r) };
    case 'cylinder':
      return {
        ...base,
        position: at(0.5),
        rotation: rotationTo(LocalUp, toward),
        scale: Vec3.create(0.4 * r, r, 0.4 * r),
      };
    case 'infiniteCone':
      return { ...base, position: Vec3.clone(center), rotation: rotationTo(ConeOpening, toward), scale: coneScale(25) };
  }
}

/** Each representation's params with `next` in place of `previous` (or appended when absent; removed when undefined) */
export function withObject(params: unknown, previous: CutObject | undefined, next: CutObject | undefined): unknown {
  const clip = RepresentationClip.get(params);
  if (!clip) return params;
  const at = previous ? clip.objects.findIndex((o) => RepresentationClip.objectsEqual(o, previous)) : -1;
  const objects = [...clip.objects];
  if (at >= 0) {
    if (next) objects[at] = next;
    else objects.splice(at, 1);
  } else if (next) objects.push(next);
  return RepresentationClip.set(params, { ...clip, objects });
}
