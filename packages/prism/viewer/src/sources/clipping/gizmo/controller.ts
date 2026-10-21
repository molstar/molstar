/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { OrderedSet } from '@molstar/core/data/int';
import type { Sphere3D } from '@molstar/core/math/geometry';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import { Binding } from '@molstar/core/util/binding';
import { Color } from '@molstar/core/util/color';
import { MarkerActions } from '@molstar/core/util/marker-action';
import type { Canvas3D } from '@molstar/graphics/canvas3d/canvas3d';
import { Mesh } from '@molstar/graphics/geo/geometry/mesh/mesh';
import type { Representation } from '@molstar/graphics/repr/representation';
import { VisualSet } from '@molstar/graphics/util/visual-set';
import { Shape, ShapeGroup } from '@molstar/model/model/shape';
import type { PluginContext } from '@molstar/plugin/context';
import { RepresentationClip } from '@molstar/plugin/state/helpers/representation-clip';
import { observableState, type ReadonlyState } from '@molstar/prism-ui/core';
import { type Cut, type CutObject, representationsSphere } from '../cuts.js';
import { gizmoGeometry, type GizmoGeometry, HandleGroups, HandleLabels, slideAxis } from './geometry.js';
import { type DragFrame, gesture, type Handle } from './gestures.js';

// The gizmos of the cuts (clipping.md §4): one per cut the person handles (`handled`), drawn in the canvas outside the state tree (a `VisualSet`),
// their clicks their own (`PluginContext.clickOwners`). Hovering a handle takes the trackball's drags; a drag moves the
// cut on the live representations (a clip-only update touches uniforms only, repr/_test/clip-update.test.ts), and its
// release is one `clipping.update` through the person's API (`commit`). Esc, or the cut's representations changed
// meanwhile (an undo, another edit), cancels: the live representations go back to what the state holds. Which cuts
// are handled outlives `hide` (the section closed): opened again, the same gizmos come back. A cut is known by its
// number, which shifts when one before it goes: a handled cut is followed by its object, or by its number when its
// object changed (an edit).

/**
 * One `clipping.update` of the cut that is `from` (by its object: its number may have shifted during the drag) to
 * `next`, through the person's API; false when it did not happen (the cut is gone, the call failed)
 */
export type GizmoCommit = (from: CutObject, next: CutObject) => Promise<boolean>;

const Colors: Readonly<Record<Handle, Color>> = {
  move: Color(0xf5a623),
  slide: Color(0x4a90e2),
  turn: Color(0x50b860),
  size: Color(0xbd10e0),
  'size-x': Color(0xbd10e0),
  'size-y': Color(0xbd10e0),
  'size-z': Color(0xbd10e0),
  radius: Color(0xbd10e0),
};
const BodyColor = Color(0x8fb8ec);
/** The handles' glow (emissive: the canvas's bloom picks it up) */
const Glow = 0.35;

interface Gizmo {
  cut: Cut;
  bounds: Sphere3D;
  geometry: GizmoGeometry;
  body: Representation.Any;
  handles: Representation.Any;
}

interface Drag {
  readonly index: number;
  readonly handle: Handle;
  readonly from: CutObject;
  readonly axis: Vec3;
  readonly frame: DragFrame;
  /** The cut's representations as the drag found them: what a cancel restores, what a commit checks is unchanged */
  readonly targets: readonly { readonly ref: string; readonly params: unknown; readonly repr: Representation.Any }[];
  /** Where the drag started and how far it went, in CSS pixels (the input's `x`, `y`; its `dx` is a fraction of the canvas) */
  readonly x0: number;
  readonly y0: number;
  dx: number;
  dy: number;
  next: CutObject;
}

export class ClipGizmos {
  readonly #shown = observableState(false);
  /** The handled cuts by number, and their objects as last seen (to follow them when the numbers shift) */
  readonly #handled = observableState<readonly number[]>([]);
  #handledObjects = new Map<number, CutObject>();
  /** How many cuts there were when the handled ones were last followed: an edit keeps the count, a removal does not */
  #followedCount = 0;
  #set: VisualSet | undefined;
  #canvas: Canvas3D | undefined;
  #commit: GizmoCommit | undefined;
  #subs: { unsubscribe(): void }[] = [];
  #releaseClicks: (() => void) | undefined;
  readonly #gizmos = new Map<number, Gizmo>();
  #hovered: { index: number; handle: Handle } | undefined;
  #drag: Drag | undefined;
  /** The trackball's bindings while a handle holds the drags (restored whole: a partial set keeps only its part) */
  #bindings: Canvas3D['attribs']['trackball']['bindings'] | undefined;
  #syncing = 0;

  constructor(
    private readonly plugin: PluginContext,
    private readonly cuts: ReadonlyState<readonly Cut[]>,
  ) {}

  get shown(): ReadonlyState<boolean> {
    return this.#shown.asReadonly;
  }

  /** The cuts whose gizmos show (while `shown`), by number */
  get handled(): ReadonlyState<readonly number[]> {
    return this.#handled.asReadonly;
  }

  /** A cut's gizmo on or off */
  setHandled(index: number, on: boolean): void {
    const cut = this.cuts.value.find((c) => c.index === index);
    if (!cut || on === this.#handled.value.includes(index)) return;
    if (on) this.#handledObjects.set(index, cut.object);
    else this.#handledObjects.delete(index);
    this.#publishHandled();
  }

  /** Every cut's gizmo on, or none */
  handleAll(on: boolean): void {
    this.#handledObjects = new Map(on ? this.cuts.value.map((c) => [c.index, c.object] as const) : []);
    this.#publishHandled();
  }

  #publishHandled() {
    this.#followedCount = this.cuts.value.length;
    const next = [...this.#handledObjects.keys()].sort((a, b) => a - b);
    const prev = this.#handled.value;
    if (next.length !== prev.length || next.some((x, i) => x !== prev[i])) this.#handled.set(next);
    void this.#sync(this.cuts.value);
  }

  /**
   * The handled cuts in `cuts`: each by its object, else by its number when its object is gone and the count is the
   * same (it was edited; a removal shifts the numbers, and the next cut is not the removed one)
   */
  #followHandled(cuts: readonly Cut[]) {
    const before = this.#handledObjects;
    const next = new Map<number, CutObject>();
    const sameCount = cuts.length === this.#followedCount;
    this.#followedCount = cuts.length;
    for (const cut of cuts) {
      const same = [...before.values()].some((o) => RepresentationClip.objectsEqual(o, cut.object));
      const old = before.get(cut.index);
      const edited = sameCount && !!old && !cuts.some((c) => RepresentationClip.objectsEqual(c.object, old));
      if (same || edited) next.set(cut.index, cut.object);
    }
    this.#handledObjects = next;
    const ids = [...next.keys()];
    const prev = this.#handled.value;
    if (ids.length !== prev.length || ids.some((x, i) => x !== prev[i])) this.#handled.set(ids);
  }

  /** Shows the gizmos; a drag's release commits through `commit`. Without a canvas (headless): false */
  show(commit: GizmoCommit): boolean {
    this.#commit = commit;
    if (this.#shown.value) return true;
    const canvas = this.plugin.canvas3d;
    if (!canvas) return false;
    this.#canvas = canvas;
    const set = new VisualSet(canvas);
    this.#set = set;
    this.#releaseClicks = this.plugin.clickOwners.claim((r) => set.has(r));
    this.#subs.push(
      this.cuts.subscribe((cuts) => void this.#sync(cuts)),
      canvas.interaction.hover.subscribe(({ current }) => this.#hover(current)),
      canvas.input.drag.subscribe((e) => this.#onDrag(e.isStart, e.x, e.y, e.dx, e.dy)),
      canvas.input.interactionEnd.subscribe(() => void this.#end()),
      canvas.input.keyDown.subscribe((e) => {
        if (e.key === 'Escape' && this.#drag) this.cancel();
      }),
      // what the cut is on changed under the drag (an undo, another edit): the drag goes, nothing committed over it
      this.plugin.state.data.events.changed.subscribe(() => {
        const d = this.#drag;
        if (d?.targets.some((t) => this.plugin.state.data.cells.get(t.ref)?.transform.params !== t.params)) {
          this.cancel();
        }
      }),
    );
    this.#shown.set(true);
    return true;
  }

  hide(): void {
    if (!this.#shown.value) return;
    this.cancel();
    this.#restoreBindings();
    for (const s of this.#subs) s.unsubscribe();
    this.#subs = [];
    this.#releaseClicks?.();
    this.#releaseClicks = undefined;
    this.#set?.dispose();
    this.#set = undefined;
    this.#gizmos.clear();
    this.#hovered = undefined;
    this.#shown.set(false);
  }

  /** The drag in progress, undone: the live representations back to what the state holds */
  cancel(): void {
    const d = this.#drag;
    if (!d) return;
    this.#drag = undefined;
    for (const t of d.targets) {
      const clip = RepresentationClip.get(t.params);
      if (clip) void t.repr.createOrUpdate({ ...t.repr.props, clip }).run();
    }
    if (!this.#hovered) this.#restoreBindings();
    // the gizmos as the cuts are now (the drag's skipped its own, and the state may have moved on under it)
    void this.#sync(this.cuts.value);
    this.#canvas?.requestDraw();
  }

  /** The gizmos follow the cuts: made, remade where a cut changed, removed with it */
  async #sync(cuts: readonly Cut[]) {
    if (this.#drag === undefined) this.#followHandled(cuts);
    const set = this.#set;
    if (!set) return;
    const version = ++this.#syncing;
    cuts = cuts.filter((c) => this.#handledObjects.has(c.index));
    for (const [index, g] of [...this.#gizmos]) {
      if (cuts.some((c) => c.index === index)) continue;
      set.remove(g.body);
      set.remove(g.handles);
      this.#gizmos.delete(index);
    }
    for (const cut of cuts) {
      // the cut being dragged is the drag's
      if (this.#drag?.index === cut.index) continue;
      const g = this.#gizmos.get(cut.index);
      if (g && RepresentationClip.objectsEqual(g.cut.object, cut.object) && sameTargets(g.cut, cut)) {
        g.cut = cut;
        continue;
      }
      const bounds = representationsSphere(
        cut.targets.flatMap((t) => {
          const cell = this.plugin.state.data.cells.get(t.ref);
          return cell ? [cell] : [];
        }),
      );
      if (!bounds) continue;
      const geometry = gizmoGeometry(cut.object, bounds);
      const made = await this.#make(cut, geometry);
      if (version !== this.#syncing || this.#set !== set) {
        set.remove(made.body);
        set.remove(made.handles);
        return;
      }
      if (g) {
        set.remove(g.body);
        set.remove(g.handles);
      }
      this.#gizmos.set(cut.index, { cut, bounds, geometry, ...made });
    }
    this.#canvas?.requestDraw();
  }

  async #make(cut: Cut, geometry: GizmoGeometry) {
    const set = this.#set!;
    const name = `Cut ${cut.index}`;
    const bodyShape = Shape.create(
      name,
      {},
      geometry.body,
      () => BodyColor,
      () => 1,
      () => name,
      1,
    );
    const handleShape = Shape.create(
      name,
      {},
      geometry.handles,
      (g) => Colors[HandleGroups[g]] ?? BodyColor,
      () => 1,
      (g) => `${name}: ${HandleLabels[HandleGroups[g]]}`,
      HandleGroups.length,
    );
    const helper = { markerActions: MarkerActions.Highlighting };
    const [body, handles] = await Promise.all([
      set.addShape(
        (_, s: Shape<Mesh>) => s,
        Mesh.Utils,
        bodyShape,
        { alpha: 0.14, ignoreLight: true },
        {
          ...helper,
          pickable: false,
        },
      ),
      set.addShape((_, s: Shape<Mesh>) => s, Mesh.Utils, handleShape, { alpha: 1, emissive: Glow }, helper),
    ]);
    set.setTransform(body, geometry.transform);
    set.setTransform(handles, geometry.transform);
    return { body: body as unknown as Representation.Any, handles: handles as unknown as Representation.Any };
  }

  /** A gizmo where `object` is: moved and turned as a whole; remade when its size changed */
  async #place(g: Gizmo, object: CutObject) {
    const set = this.#set;
    if (!set) return;
    const geometry = gizmoGeometry(object, g.bounds);
    if (sameSize(g.cut.object, object) && g.cut.object.invert === object.invert) {
      g.geometry = geometry;
      set.setTransform(g.body, geometry.transform);
      set.setTransform(g.handles, geometry.transform);
      return;
    }
    const made = await this.#make({ ...g.cut, object }, geometry);
    // remade meanwhile (a sync after the commit) or gone: this one is not drawn
    if (this.#set !== set || ![...this.#gizmos.values()].includes(g)) {
      set.remove(made.body);
      set.remove(made.handles);
      return;
    }
    set.remove(g.body);
    set.remove(g.handles);
    Object.assign(g, { geometry, ...made, cut: { ...g.cut, object } });
  }

  #hover(current: Representation.Loci) {
    if (this.#drag) return;
    const index = [...this.#gizmos].find(([, g]) => g.handles === current.repr)?.[0];
    if (index === undefined || !ShapeGroup.isLoci(current.loci) || !current.loci.groups.length) {
      if (this.#hovered) {
        this.#hovered = undefined;
        this.#restoreBindings();
      }
      return;
    }
    const handle = HandleGroups[OrderedSet.start(current.loci.groups[0].ids)];
    if (!handle) return;
    this.#hovered = { index, handle };
    this.#takeDrags();
  }

  #onDrag(isStart: boolean, x: number, y: number, fx: number, fy: number) {
    if (isStart) {
      if (!this.#hovered || this.#drag) return;
      // the first event comes a few pixels after the press: its `dx`, `dy` (fractions of the canvas) lead back to it
      const canvas = this.#canvas!;
      const ratio = canvas.input.pixelRatio || 1;
      const w = canvas.camera.viewport.width / ratio;
      const h = canvas.camera.viewport.height / ratio;
      this.#begin(this.#hovered.index, this.#hovered.handle, x - fx * w, y - fy * h);
    }
    const d = this.#drag;
    if (!d) return;
    d.dx = x - d.x0;
    d.dy = y - d.y0;
    d.next = gesture(d.handle, d.from, d.axis, d.frame, d.dx, d.dy);
    this.#live();
  }

  #begin(index: number, handle: Handle, x0: number, y0: number) {
    const g = this.#gizmos.get(index);
    const canvas = this.#canvas;
    if (!g || !canvas) return;
    const camera = canvas.camera;
    const dir = Vec3.normalize(Vec3(), Vec3.sub(Vec3(), camera.state.target, camera.state.position));
    const right = Vec3.normalize(Vec3(), Vec3.cross(Vec3(), dir, camera.state.up));
    const up = Vec3.normalize(Vec3(), Vec3.cross(Vec3(), right, dir));
    const at = Vec3.fromArray(Vec3(), g.geometry.transform, 12);
    const targets = g.cut.targets.flatMap((t) => {
      const cell = this.plugin.state.data.cells.get(t.ref);
      const repr = (cell?.obj?.data as { repr?: Representation.Any } | undefined)?.repr;
      return cell && repr ? [{ ref: t.ref, params: cell.transform.params, repr }] : [];
    });
    this.#drag = {
      index,
      handle,
      from: g.cut.object,
      axis: slideAxis(g.cut.object),
      frame: { right, up, perPixel: perPixel(canvas, at) },
      targets,
      x0,
      y0,
      dx: 0,
      dy: 0,
      next: g.cut.object,
    };
  }

  /**
   * The drag's cut on the live representations and the gizmo, at most one update in flight: each round takes the drag
   * as it is then (the latest move; a new drag's first move too, while the last one's round finishes)
   */
  #running: Promise<void> | undefined;
  #again = false;
  #live() {
    if (this.#running) {
      this.#again = true;
      return;
    }
    const run = async () => {
      do {
        this.#again = false;
        const d = this.#drag;
        if (!d) break;
        const next = d.next;
        await Promise.all(
          d.targets.map((t) => {
            const clip = RepresentationClip.get(t.params);
            if (!clip) return undefined;
            const objects = clip.objects.map((o) => (RepresentationClip.objectsEqual(o, d.from) ? next : o));
            return t.repr.createOrUpdate({ ...t.repr.props, clip: { ...clip, objects } }).run();
          }),
        );
        const g = this.#gizmos.get(d.index);
        if (g && this.#drag === d) await this.#place(g, next);
        this.#canvas?.requestDraw();
      } while (this.#again);
    };
    this.#running = run().finally(() => {
      this.#running = undefined;
    });
  }

  async #end() {
    const d = this.#drag;
    if (!d) return;
    this.#drag = undefined;
    // a live update still running would land after the commit's
    await this.#running;
    const unchanged = d.targets.every((t) => this.plugin.state.data.cells.get(t.ref)?.transform.params === t.params);
    if (!unchanged || RepresentationClip.objectsEqual(d.from, d.next)) {
      this.#drag = d;
      this.cancel();
      return;
    }
    if (!this.#hovered) this.#restoreBindings();
    // the state takes the new cut and the gizmos follow it (sync); not taken, the live representations go back
    const done = (await this.#commit?.(d.from, d.next)) ?? false;
    if (!done) {
      this.#drag = d;
      this.cancel();
    }
  }

  #takeDrags() {
    const canvas = this.#canvas;
    if (!canvas || this.#bindings) return;
    const bindings = { ...canvas.attribs.trackball.bindings };
    this.#bindings = bindings;
    canvas.setAttribs({
      trackball: {
        bindings: {
          ...bindings,
          dragRotate: Binding.Empty,
          dragRotateZ: Binding.Empty,
          dragPan: Binding.Empty,
          dragFocus: Binding.Empty,
          dragFocusZoom: Binding.Empty,
        },
      },
    });
  }

  #restoreBindings() {
    const bindings = this.#bindings;
    if (!bindings || this.#drag) return;
    this.#bindings = undefined;
    this.#canvas?.setAttribs({ trackball: { bindings } });
  }
}

function sameTargets(a: Cut, b: Cut) {
  return a.targets.length === b.targets.length && a.targets.every((t, i) => t.ref === b.targets[i].ref);
}

function sameSize(a: CutObject, b: CutObject) {
  return Vec3.exactEquals(a.scale, b.scale) && a.type === b.type;
}

/** Å a CSS pixel is at `at`: the height the camera sees there over the viewport's height */
function perPixel(canvas: Canvas3D, at: Vec3): number {
  const { camera } = canvas;
  const s = camera.state;
  const ratio = canvas.input.pixelRatio || 1;
  const height =
    s.mode === 'orthographic'
      ? 2 * Math.tan(s.fov / 2) * Vec3.distance(s.position, s.target)
      : 2 *
        Math.tan(s.fov / 2) *
        Math.max(
          1e-3,
          Vec3.dot(Vec3.sub(Vec3(), at, s.position), Vec3.normalize(Vec3(), Vec3.sub(Vec3(), s.target, s.position))),
        );
  return height / (camera.viewport.height / ratio);
}
