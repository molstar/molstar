/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Ludovic Autin <ludovic.autin@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
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
import { Particle } from '@molstar/model/model/particles/particle-list';
import { observableState, type ReadonlyState } from '@molstar/prism-ui/core';
import { BehaviorSubject, type Observable, Subject } from 'rxjs';
import { currentTransform } from '../../api/structure.js';
import { volumeNames } from '../../api/volume.js';
import { currentParticleTransform } from '../particle-transform.js';
import { currentShapeMatrix, shapeDrawnCenter } from '../shape-transform.js';
import { shapeEntries } from '../shapes.js';
import { currentVolumeTransform, isStreamedForStructure, volumeCenter } from '../volume-transform.js';

// A model moved by hand (move.md): a structure, a particle list, a shape or a volume, the gizmo library's transform
// handles on it, through its transform decorator (Mol*'s `TransformStructureConformation`, as `structure.transform`
// writes it; `TransformParticleList`, as `particles.transform`; `TransformShapeProvider`, as `shape.transform`;
// `VolumeTransform`, as `volume.transform`). A rigid move is a turn about the object's center then a move:
// the gizmo's frame is that center as placed now, turned as now. A drag's preview moves the object's representations
// as drawn (their transient transform, `previewTransform`: applied on the GPU, no instance touched); recomputing it a
// frame would not keep up; only its own representations (what the transform moves: not a volume streamed under a
// structure, not a list's unit cell). Its release is one transform operation with the new matrix (`commit`); the
// transient transform is put back once the state has the move, frames held meanwhile (`holdFrames`).

/** What can be moved: a structure, a particle list, a shape, a volume */
export type MovableKind = 'structure' | 'particles' | 'shape' | 'volume';

/** The representations a kind's transform moves (not those of another kind under it: a volume streamed for a structure) */
const MovedRepresentation: Record<MovableKind, (o: unknown) => boolean> = {
  structure: (o) => PluginStateObject.Molecule.Structure.Representation3D.is(o as never),
  particles: (o) => PluginStateObject.Particle.Representation3D.is(o as never),
  shape: (o) => PluginStateObject.Shape.Representation3D.is(o as never),
  volume: (o) => PluginStateObject.Volume.Representation3D.is(o as never),
};

/** The volumes loaded on their own (not streamed for a structure): what can be moved, with their names */
export function movableVolumes(plugin: PluginContext) {
  const state = plugin.state.data;
  const list = plugin.managers.volume.hierarchy.current.volumes.filter((v) => !isStreamedForStructure(state, v.cell));
  return { list, names: volumeNames(list) };
}

/** An object's placement: the transform on top of where it was loaded, and its center there */
export interface Placement {
  readonly matrix: Mat4;
  readonly center: Vec3;
}

/** The handles' size on screen at any zoom: CSS pixels an arrow is long (the rings 0.65 of it, the roll circle 0.8) */
const HandleScreenSize = 90;
const Handles = transformHandles(1);

/** One transform operation (the object by ref, the matrix from where it was loaded); false: not taken */
export type MoveCommit = (ref: string, kind: MovableKind, matrix: Mat4) => Promise<boolean>;

/** A placement's pose: its center as placed, its rotation */
export function placementPose(p: Placement): GizmoPose {
  const position = Vec3.transformMat4(Vec3(), p.center, p.matrix);
  const rotation = Quat.normalize(Quat(), Quat.fromMat4(Quat(), p.matrix));
  return { position, rotation, scale: Vec3.create(1, 1, 1) };
}

/** The matrix a pose makes of an object: turned about its center, then moved there (`T(position) R T(−center)`) */
export function poseMatrix(pose: GizmoPose, center: Vec3): Mat4 {
  const out = Mat4.fromTranslation(Mat4(), pose.position);
  Mat4.mul(out, out, Mat4.fromQuat(Mat4(), pose.rotation));
  return Mat4.mul(out, out, Mat4.fromTranslation(Mat4(), Vec3.negate(Vec3(), center)));
}

/** The placement a drag on `handle` makes from `from`: what it has beyond a rigid move (a scale or shear given as a
 * matrix) kept under the new pose */
export function applyToPlacement(from: Placement, handle: Handle, frame: DragFrame, dx: number, dy: number): Placement {
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

/** No transient transform: the representations drawn where the state has them */
const NoPreview = Mat4.identity();

class MoveBinding implements GizmoBinding<Placement> {
  readonly #invalidated = new Subject<string>();
  #drag: { key: string; drawn: readonly Representation.Any[]; params: readonly unknown[] } | undefined;
  readonly #sub: { unsubscribe(): void };

  constructor(
    private readonly plugin: PluginContext,
    readonly items: Observable<readonly GizmoItem<Placement>[]>,
    private readonly commitMove: MoveCommit,
    private readonly kindOf: (ref: string) => MovableKind,
  ) {
    // the object's cells changed under a drag (an undo, another edit): the drag is cancelled
    this.#sub = plugin.state.data.events.changed.subscribe(() => {
      const d = this.#drag;
      if (d && this.#paramsOf(d.key).some((p, i) => p !== d.params[i])) this.#invalidated.next(d.key);
    });
  }

  get invalidated(): Observable<string> {
    return this.#invalidated;
  }

  /** The handles at a fixed size on screen (a large object's handles are not off screen, a small one's not tiny) */
  layout(item: GizmoItem<Placement>) {
    const { parts, groups } = Handles;
    return { parts, groups, frame: frameOf(placementPose(item.value)), screenSize: HandleScreenSize };
  }

  /** The handles are the same for any placement: only their frame moves */
  sameGeometry() {
    return true;
  }

  /** The same frame: the center moves without the matrix (a trajectory's frame, a model's change) */
  equal(a: Placement, b: Placement) {
    return Mat4.areEqual(a.matrix, b.matrix, 1e-6) && Vec3.distance(a.center, b.center) < 1e-6;
  }

  apply(from: Placement, handle: Handle, frame: DragFrame, dx: number, dy: number): Placement {
    return applyToPlacement(from, handle, frame, dx, dy);
  }

  begin(key: string) {
    const drawn: Representation.Any[] = [];
    const isMoved = MovedRepresentation[this.kindOf(key)];
    for (const cell of this.#cellsUnder(key)) {
      if (!isMoved(cell.obj)) continue;
      const repr = (cell.obj?.data as { repr?: Representation.Any } | undefined)?.repr;
      if (repr?.state) drawn.push(repr);
    }
    this.#drag = { key, drawn, params: this.#paramsOf(key) };
  }

  async preview(_key: string, from: Placement, next: Placement) {
    // what the drag adds, on top of where the state has it: a transient transform (a uniform, no instance touched)
    const delta = Mat4.mul(Mat4(), next.matrix, Mat4.invert(Mat4(), from.matrix));
    for (const repr of this.#drag?.drawn ?? []) this.#draw(repr, delta);
    this.plugin.canvas3d?.requestDraw();
  }

  async commit(key: string, _from: Placement, next: Placement) {
    // this drag's representations (the controller starts no other drag until the commit is done)
    const drag = this.#drag;
    this.#drag = undefined;
    const done = (async () => {
      try {
        return await this.commitMove(key, this.kindOf(key), next.matrix);
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

  #restore(drag: { drawn: readonly Representation.Any[] } | undefined) {
    for (const repr of drag?.drawn ?? []) this.#draw(repr, NoPreview);
    this.plugin.canvas3d?.requestDraw();
  }

  #draw(repr: Representation.Any, previewTransform: Mat4) {
    repr.setState({ previewTransform });
    this.plugin.canvas3d?.update(repr);
  }

  /** The object's cell and every cell under it (its transform, components, representations) */
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

/** The models that can be moved, and the handles on one of them (the target) while shown */
export class MoveGizmos {
  readonly #target = observableState<string | undefined>(undefined);
  readonly #shown = observableState(false);
  readonly #items = new BehaviorSubject<readonly GizmoItem<Placement>[]>([]);
  #controller: GizmoController<Placement> | undefined;
  #binding: MoveBinding | undefined;
  #release: (() => void) | undefined;
  #sub: { unsubscribe(): void } | undefined;

  constructor(private readonly plugin: PluginContext) {}

  /** The model the handles are on (its ref) */
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
    const binding = new MoveBinding(this.plugin, this.#items, commit, (ref) => this.kindOf(ref));
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

  /** What a ref is: a particle list's, a volume's, a shape's, else a structure's */
  kindOf(ref: string): MovableKind {
    const { managers, state } = this.plugin;
    if (managers.particles.hierarchy.current.refs.get(ref)?.kind === 'particle-list') return 'particles';
    if (managers.volume.hierarchy.current.volumes.some((v) => v.cell.transform.ref === ref)) return 'volume';
    if (shapeEntries(state.data).some((e) => e.cell.transform.ref === ref)) return 'shape';
    return 'structure';
  }

  /** The target as an item: its placement now */
  #publish() {
    const item = this.#target.value ? this.#itemOf(this.#target.value) : undefined;
    this.#items.next(item ? [item] : []);
  }

  #itemOf(ref: string): GizmoItem<Placement> | undefined {
    const { plugin } = this;
    const list = plugin.managers.particles.hierarchy.current.lists.find((x) => x.cell.transform.ref === ref);
    if (list) {
      const data = list.cell.obj?.data;
      if (!data?.count) return undefined;
      const value: Placement = {
        matrix: currentParticleTransform(plugin.state.data, ref, data),
        center: Vec3.clone(Particle.getBoundary(data).sphere.center),
      };
      return { key: ref, value, label: list.cell.obj?.label ?? 'Particles' };
    }
    const { list: volumes, names } = movableVolumes(plugin);
    const volume = volumes.find((v) => v.cell.transform.ref === ref);
    if (volume) {
      const data = volume.cell.obj?.data;
      if (!data) return undefined;
      const value: Placement = {
        matrix: currentVolumeTransform(plugin.state.data, ref, data),
        center: Vec3.clone(volumeCenter(data)),
      };
      return { key: ref, value, label: names.get(ref) ?? 'Volume' };
    }
    const shape = shapeEntries(plugin.state.data).find((e) => e.cell.transform.ref === ref);
    if (shape) {
      const matrix = currentShapeMatrix(plugin.state.data, ref);
      // known once drawn: the handles come with its first representation
      const center = shapeDrawnCenter(shape, matrix);
      return center ? { key: ref, value: { matrix, center }, label: shape.label } : undefined;
    }
    const s = plugin.managers.structure.hierarchy.current.structures.find((x) => x.cell.transform.ref === ref);
    const data = s?.cell.obj?.data;
    if (!s || !data) return undefined;
    const value: Placement = { matrix: currentTransform(s, data), center: Vec3.clone(data.boundary.sphere.center) };
    return { key: ref, value, label: s.cell.obj?.label ?? 'Structure' };
  }
}
