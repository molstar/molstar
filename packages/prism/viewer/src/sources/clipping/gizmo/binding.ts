/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Ludovic Autin <ludovic.autin@gmail.com>
 */

import type { Sphere3D } from '@molstar/core/math/geometry';
import { Quat, Vec3 } from '@molstar/core/math/linear-algebra';
import type { Clip } from '@molstar/core/util/clip';
import type { GizmoBinding, GizmoItem, GizmoLayout } from '@molstar/graphics/util/gizmo/controller';
import { type DragFrame, type GizmoPose, poseGesture, sizeRatio } from '@molstar/graphics/util/gizmo/gestures';
import type { Handle } from '@molstar/graphics/util/gizmo/handle';
import type { Representation } from '@molstar/graphics/repr/representation';
import type { PluginContext } from '@molstar/plugin/context';
import { RepresentationClip } from '@molstar/plugin/state/helpers/representation-clip';
import { type Observable, Subject } from 'rxjs';
import { type CutObject, coneAngle, coneScale } from '../cuts.js';
import { ConeAngleHandle, gizmoGeometry } from './geometry.js';

// The clipping gizmos' binding (clipping.md §4) to the gizmo library (`@molstar/graphics/util/gizmo`): an item is a cut,
// its value the cut's clip object. A drag's gesture works on the object's pose (position, rotation, scale; a cone's
// angle is the one gesture of its own); a preview is the object on the live representations (a clip-only update
// touches uniforms only, repr/_test/clip-update.test.ts); a cancel puts them back to what the state holds; a commit is
// the app's (one `clipping.update` through the person's API). What a drag is on changing in the state (an undo, another
// edit) invalidates it.

/** One `clipping.update` of the cut that is `from` (by its object: its number may shift) to `next`; false: not taken */
export type GizmoCommit = (from: CutObject, next: CutObject) => Promise<boolean>;

/** A clip object's pose: its position, its rotation (axis and degrees) as a quaternion, its scale */
export function cutPose(o: CutObject): GizmoPose {
  const axis = Vec3.normalize(Vec3(), o.rotation.axis);
  return {
    position: o.position,
    rotation: Quat.setAxisAngle(Quat(), axis, (o.rotation.angle * Math.PI) / 180),
    scale: o.scale,
  };
}

/** The clip object a pose makes from `o`: only what the gesture changed (an untouched rotation keeps its axis and angle) */
export function poseCut(o: CutObject, base: GizmoPose, pose: GizmoPose): CutObject {
  return {
    ...o,
    position: pose.position === base.position ? o.position : Vec3.clone(pose.position),
    rotation: pose.rotation === base.rotation ? o.rotation : rotationOf(pose.rotation),
    scale: pose.scale === base.scale ? o.scale : Vec3.clone(pose.scale),
  };
}

/** A unit quaternion as Mol*'s clip rotation (axis, degrees) */
export function rotationOf(q: Quat): Clip.Rotation {
  const axis = Vec3();
  const rad = Quat.getAxisAngle(axis, q);
  if (!Number.isFinite(rad) || Vec3.magnitude(axis) < 1e-9) return { axis: Vec3.create(1, 0, 0), angle: 0 };
  return { axis: Vec3.normalize(axis, axis), angle: (rad * 180) / Math.PI };
}

/** The clip object a drag on `handle` makes from `from` */
export function applyToCut(from: CutObject, handle: Handle, frame: DragFrame, dx: number, dy: number): CutObject {
  if (handle.kind === 'custom' && handle.id === ConeAngleHandle) {
    // the rim's radius goes as the angle's tangent
    const tan = Math.tan((coneAngle(from.scale) * Math.PI) / 180) * sizeRatio(frame, dx, dy);
    const angle = Math.max(1, Math.min(89, (Math.atan(tan) * 180) / Math.PI));
    return { ...from, scale: coneScale(angle) };
  }
  if (from.type === 'plane' && handle.kind === 'scale') return from;
  const base = cutPose(from);
  const pose = poseGesture(base, handle, frame, dx, dy);
  return pose ? poseCut(from, base, pose) : from;
}

interface Target {
  readonly ref: string;
  readonly params: unknown;
  readonly repr: Representation.Any;
}

/** A cut as the gizmos see it: its object, what it cuts (the representations' bounds) */
export interface ClipItem extends GizmoItem<CutObject> {
  readonly bounds: Sphere3D;
}

export class ClipBinding implements GizmoBinding<CutObject> {
  readonly #invalidated = new Subject<string>();
  /** The representations a drag found (by its item's key): what a cancel restores, what a change under it is told by */
  #drag: { key: string; targets: readonly Target[] } | undefined;
  #commit: GizmoCommit | undefined;
  readonly #unsubscribe: () => void;

  constructor(
    private readonly plugin: PluginContext,
    readonly items: Observable<readonly ClipItem[]>,
  ) {
    const sub = plugin.state.data.events.changed.subscribe(() => {
      const d = this.#drag;
      if (d?.targets.some((t) => plugin.state.data.cells.get(t.ref)?.transform.params !== t.params)) {
        this.#invalidated.next(d.key);
      }
    });
    this.#unsubscribe = () => sub.unsubscribe();
  }

  get invalidated(): Observable<string> {
    return this.#invalidated;
  }

  /** How a drag's release is committed (the person's API) */
  setCommit(commit: GizmoCommit) {
    this.#commit = commit;
  }

  layout(item: GizmoItem<CutObject>): GizmoLayout {
    const g = gizmoGeometry(item.value, (item as ClipItem).bounds);
    return { body: g.body, handles: g.handles, groups: g.groups, frame: g.transform, facesCamera: g.facesCamera };
  }

  sameGeometry(a: CutObject, b: CutObject): boolean {
    return a.type === b.type && a.invert === b.invert && Vec3.exactEquals(a.scale, b.scale);
  }

  equal(a: CutObject, b: CutObject): boolean {
    return RepresentationClip.objectsEqual(a, b);
  }

  apply(from: CutObject, handle: Handle, frame: DragFrame, dx: number, dy: number): CutObject {
    return applyToCut(from, handle, frame, dx, dy);
  }

  begin(key: string, value: CutObject) {
    const targets: Target[] = [];
    for (const cell of this.plugin.state.data.cells.values()) {
      const clip = RepresentationClip.get(cell.transform.params);
      const repr = (cell.obj?.data as { repr?: Representation.Any } | undefined)?.repr;
      if (!clip || !repr || !clip.objects.some((o) => RepresentationClip.sameCut(o, value))) continue;
      targets.push({ ref: cell.transform.ref, params: cell.transform.params, repr });
    }
    this.#drag = { key, targets };
  }

  async preview(_key: string, from: CutObject, next: CutObject): Promise<void> {
    await Promise.all(
      (this.#drag?.targets ?? []).map((t) => {
        const clip = RepresentationClip.get(t.params);
        if (!clip) return undefined;
        const objects = clip.objects.map((o) => (RepresentationClip.sameCut(o, from) ? next : o));
        return t.repr.createOrUpdate({ ...t.repr.props, clip: { ...clip, objects } }).run();
      }),
    );
  }

  async commit(_key: string, from: CutObject, next: CutObject): Promise<boolean> {
    const done = (await this.#commit?.(from, next)) ?? false;
    if (done) this.#drag = undefined;
    return done;
  }

  cancel(): void {
    // what the state holds now (a change under the drag cancels it: the old params would undo that change on screen)
    for (const t of this.#drag?.targets ?? []) {
      const cell = this.plugin.state.data.cells.get(t.ref);
      const clip = RepresentationClip.get(cell?.transform.params);
      const repr = (cell?.obj?.data as { repr?: Representation.Any } | undefined)?.repr;
      if (clip && repr) void repr.createOrUpdate({ ...repr.props, clip }).run();
    }
    this.#drag = undefined;
  }

  dispose() {
    this.cancel();
    this.#unsubscribe();
    this.#invalidated.complete();
  }
}
