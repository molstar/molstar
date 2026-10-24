/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Ludovic Autin <ludovic.autin@gmail.com>
 */

import { Mat4, Quat, Vec3 } from '@molstar/core/math/linear-algebra';
import type { StateObjectCell } from '@molstar/core/state';
import { type GizmoBinding, GizmoController, type GizmoItem } from '@molstar/graphics/util/gizmo/controller';
import { type DragFrame, type GizmoPose, poseGesture } from '@molstar/graphics/util/gizmo/gestures';
import type { Handle } from '@molstar/graphics/util/gizmo/handle';
import { transformHandles } from '@molstar/graphics/util/gizmo/transform';
import type { Representation } from '@molstar/graphics/repr/representation';
import type { PluginContext } from '@molstar/plugin/context';
import { PluginStateObject } from '@molstar/plugin/state/objects';
import { observableState, type ReadonlyState } from '@molstar/prism-ui/core';
import { BehaviorSubject, type Observable, Subject } from 'rxjs';
import { currentTransform } from '../../api/structure.js';

// A structure moved by hand (structure-move.md): the gizmo library's transform handles on a structure, through its
// conformation transform (Mol*'s `TransformStructureConformation`, as `structure.transform` writes it). A rigid move
// is a turn about the structure's center then a move: the gizmo's frame is that center as placed now, turned as now.
// A drag's preview moves the structure's representations as drawn (their state transforms, on top of what they had):
// recomputing a structure a frame would not keep up; only the structure's own representations (what the transform
// moves: not a volume streamed under it). Its release is one `structure.transform` with the new matrix (`commit`); the
// representations' own transforms are put back once the state has the move, frames held meanwhile (`holdFrames`).

/** A structure's placement: the transform on top of where it was loaded, and its center there */
export interface StructurePlacement {
  readonly matrix: Mat4;
  readonly center: Vec3;
  /** Its reach (the bounds' radius), for the handles' size */
  readonly radius: number;
}

/** One `structure.transform` (the structure by ref, the matrix from where it was loaded); false: not taken */
export type MoveCommit = (structure: string, matrix: Mat4) => Promise<boolean>;

/** A placement's pose: its center as placed, its rotation */
export function placementPose(p: StructurePlacement): GizmoPose {
  const position = Vec3.transformMat4(Vec3(), p.center, p.matrix);
  const rotation = Quat.normalize(Quat(), Quat.fromMat4(Quat(), p.matrix));
  return { position, rotation, scale: Vec3.create(1, 1, 1) };
}

/** The matrix a pose makes of a structure: turned about its center, then moved there (`T(position) R T(−center)`) */
export function poseMatrix(pose: GizmoPose, center: Vec3): Mat4 {
  const out = Mat4.fromTranslation(Mat4(), pose.position);
  Mat4.mul(out, out, Mat4.fromQuat(Mat4(), pose.rotation));
  return Mat4.mul(out, out, Mat4.fromTranslation(Mat4(), Vec3.negate(Vec3(), center)));
}

/** The placement a drag on `handle` makes from `from`: what it has beyond a rigid move (a scale or shear given as a
 * matrix) kept under the new pose */
export function applyToPlacement(
  from: StructurePlacement,
  handle: Handle,
  frame: DragFrame,
  dx: number,
  dy: number,
): StructurePlacement {
  const base = placementPose(from);
  const pose = poseGesture(base, handle, frame, dx, dy);
  if (!pose) return from;
  const residual = Mat4.mul(Mat4(), Mat4.invert(Mat4(), poseMatrix(base, from.center)), from.matrix);
  return { ...from, matrix: Mat4.mul(Mat4(), poseMatrix(pose, from.center), residual) };
}

/** The frame the handles are drawn in */
function frameOf(pose: GizmoPose) {
  return Mat4.mul(Mat4(), Mat4.fromTranslation(Mat4(), pose.position), Mat4.fromQuat(Mat4(), pose.rotation));
}

/** The longest a commit holds the canvas's frames (ms) */
const MaxCommitHoldMs = 3000;

interface Drawn {
  readonly repr: Representation.Any;
  readonly transform: Mat4;
}

class MoveBinding implements GizmoBinding<StructurePlacement> {
  readonly #invalidated = new Subject<string>();
  #drag: { key: string; drawn: readonly Drawn[]; params: readonly unknown[] } | undefined;
  readonly #sub: { unsubscribe(): void };

  constructor(
    private readonly plugin: PluginContext,
    readonly items: Observable<readonly GizmoItem<StructurePlacement>[]>,
    private readonly commitMove: MoveCommit,
  ) {
    // the structure's cells changed under a drag (an undo, another edit): the drag is cancelled
    this.#sub = plugin.state.data.events.changed.subscribe(() => {
      const d = this.#drag;
      if (d && this.#paramsOf(d.key).some((p, i) => p !== d.params[i])) this.#invalidated.next(d.key);
    });
  }

  get invalidated(): Observable<string> {
    return this.#invalidated;
  }

  layout(item: GizmoItem<StructurePlacement>) {
    return { ...transformHandles(item.value.radius * 0.6), frame: frameOf(placementPose(item.value)) };
  }

  sameGeometry(a: StructurePlacement, b: StructurePlacement) {
    return a.radius === b.radius;
  }

  equal(a: StructurePlacement, b: StructurePlacement) {
    return Mat4.areEqual(a.matrix, b.matrix, 1e-6);
  }

  apply(from: StructurePlacement, handle: Handle, frame: DragFrame, dx: number, dy: number): StructurePlacement {
    return applyToPlacement(from, handle, frame, dx, dy);
  }

  begin(key: string) {
    const drawn: Drawn[] = [];
    for (const cell of this.#cellsUnder(key)) {
      if (!PluginStateObject.Molecule.Structure.Representation3D.is(cell.obj)) continue;
      const repr = (cell.obj?.data as { repr?: Representation.Any } | undefined)?.repr;
      if (repr?.state?.transform) drawn.push({ repr, transform: Mat4.clone(repr.state.transform) });
    }
    this.#drag = { key, drawn, params: this.#paramsOf(key) };
  }

  async preview(_key: string, from: StructurePlacement, next: StructurePlacement) {
    // what the drag adds, on top of where the state has it
    const delta = Mat4.mul(Mat4(), next.matrix, Mat4.invert(Mat4(), from.matrix));
    for (const d of this.#drag?.drawn ?? []) this.#draw(d.repr, Mat4.mul(Mat4(), delta, d.transform));
    this.plugin.canvas3d?.requestDraw();
  }

  async commit(key: string, _from: StructurePlacement, next: StructurePlacement) {
    // this drag's representations (the controller starts no other drag until the commit is done)
    const drag = this.#drag;
    this.#drag = undefined;
    const done = (async () => {
      try {
        return await this.commitMove(key, next.matrix);
      } finally {
        // the state has the move (or not): the representations as they were drawn before the drag
        this.#restore(drag);
      }
    })();
    // Frames held until the state has the move and the preview is off: the representations are rebuilt where the move
    // puts them while still carrying the preview's transform (moved twice), or, the preview off first, drawn where
    // they were (a blink either way). At most a few seconds: a commit that hangs does not freeze the canvas.
    this.plugin.canvas3d?.holdFrames(done, MaxCommitHoldMs);
    return done;
  }

  cancel() {
    this.#restore(this.#drag);
    this.#drag = undefined;
  }

  dispose() {
    this.cancel();
    this.#sub.unsubscribe();
    this.#invalidated.complete();
  }

  #restore(drag: { drawn: readonly Drawn[] } | undefined) {
    for (const d of drag?.drawn ?? []) this.#draw(d.repr, d.transform);
    this.plugin.canvas3d?.requestDraw();
  }

  #draw(repr: Representation.Any, transform: Mat4) {
    repr.setState({ transform });
    this.plugin.canvas3d?.update(repr);
  }

  /** The structure's cell and every cell under it (its transform, components, representations) */
  #cellsUnder(ref: string): StateObjectCell[] {
    const state = this.plugin.state.data;
    const out: StateObjectCell[] = [];
    const walk = (r: string) => {
      const cell = state.cells.get(r);
      if (cell) out.push(cell);
      // a ref gone (the structure removed under a drag) has no children
      const children = state.tree.children.get(r);
      if (children) for (const child of children.values()) walk(child);
    };
    walk(ref);
    return out;
  }

  #paramsOf(ref: string): unknown[] {
    return this.#cellsUnder(ref).map((c) => c.transform.params);
  }
}

/** The structures that can be moved, and the handles on one of them (the target) while shown */
export class StructureMoveGizmos {
  readonly #target = observableState<string | undefined>(undefined);
  readonly #shown = observableState(false);
  readonly #items = new BehaviorSubject<readonly GizmoItem<StructurePlacement>[]>([]);
  #controller: GizmoController<StructurePlacement> | undefined;
  #binding: MoveBinding | undefined;
  #release: (() => void) | undefined;
  #sub: { unsubscribe(): void } | undefined;

  constructor(private readonly plugin: PluginContext) {}

  /** The structure the handles are on (its ref) */
  get target(): ReadonlyState<string | undefined> {
    return this.#target.asReadonly;
  }

  get shown(): ReadonlyState<boolean> {
    return this.#shown.asReadonly;
  }

  setTarget(ref: string | undefined) {
    this.#target.set(ref);
    this.#publish();
  }

  /** The handles on the target; a drag's release commits through `commit`. Without a canvas (headless): false */
  show(commit: MoveCommit): boolean {
    if (this.#controller) return true;
    const canvas = this.plugin.canvas3d;
    if (!canvas) return false;
    const binding = new MoveBinding(this.plugin, this.#items, commit);
    const controller = new GizmoController(canvas, binding);
    this.#binding = binding;
    this.#controller = controller;
    this.#release = this.plugin.managers.interactivity.clickOwners.claim((r) => controller.owns(r));
    this.#sub = this.plugin.state.data.events.changed.subscribe(() => this.#publish());
    this.#publish();
    this.#shown.set(true);
    return true;
  }

  hide() {
    if (!this.#controller) return;
    this.#sub?.unsubscribe();
    this.#release?.();
    this.#controller.dispose();
    this.#binding?.dispose();
    this.#controller = undefined;
    this.#binding = undefined;
    this.#shown.set(false);
  }

  /** The target as an item: its placement now */
  #publish() {
    const ref = this.#target.value;
    const s = ref
      ? this.plugin.managers.structure.hierarchy.current.structures.find((x) => x.cell.transform.ref === ref)
      : undefined;
    const data = s?.cell.obj?.data;
    if (!s || !data) {
      this.#items.next([]);
      return;
    }
    const sphere = data.boundary.sphere;
    const value: StructurePlacement = {
      matrix: currentTransform(s, data),
      center: Vec3.clone(sphere.center),
      radius: Math.max(sphere.radius, 1),
    };
    this.#items.next([{ key: s.cell.transform.ref, value, label: s.cell.obj?.label ?? 'Structure' }]);
  }
}
