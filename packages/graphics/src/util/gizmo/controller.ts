/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { OrderedSet } from '@molstar/core/data/int';
import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import { Binding } from '@molstar/core/util/binding';
import { Color } from '@molstar/core/util/color';
import { MarkerActions } from '@molstar/core/util/marker-action';
import { Shape, ShapeGroup } from '@molstar/model/model/shape';
import { BehaviorSubject, type Observable, type Subscription } from 'rxjs';
import type { Canvas3D } from '../../canvas3d/canvas3d.js';
import { Mesh } from '../../geo/geometry/mesh/mesh.js';
import type { Representation } from '../../repr/representation.js';
import { VisualSet } from '../visual-set.js';
import { dragFrame, pressPoint } from './drag-frame.js';
import type { DragFrame } from './gestures.js';
import type { Handle, HandleGroup } from './handle.js';
import { facingCamera } from './primitives.js';

// Gizmos: handles in the 3D view for items an app edits (a clip object, a component's place), drawn outside the state
// tree (a `VisualSet`) and dragged live. What an item is and how it changes is the app's (`GizmoBinding`); the controller
// draws each item's gizmo, takes the trackball's drags while a handle is hovered (its bindings saved whole and restored
// at every end), runs a drag's previews (one in flight, the latest wins), and commits on release: one step. Esc, the
// binding saying what the drag started from changed (`invalidated`), or a commit not taken cancels: the binding puts the
// item back, the gizmos follow what the items are. Only `Canvas3D`'s public API (no gl, no Canvas3D internals).

/** An item a binding shows a gizmo for: its key (stable through edits), its value, its name */
export interface GizmoItem<T> {
  readonly key: string;
  readonly value: T;
  readonly label: string;
}

/** An item's gizmo: its meshes in its own frame, what each handle group does, the frame in the world */
export interface GizmoLayout {
  /** Translucent, not picked (a click needs alpha ≥ 0.5, and a body covers what is behind it) */
  readonly body?: Mesh;
  /** Opaque, one group a handle (`groups[i]` is group i's) */
  readonly handles: Mesh;
  readonly groups: readonly HandleGroup[];
  readonly frame: Mat4;
  /** The frame turned to face the camera, at its place (a sphere's outline: always round) */
  readonly facesCamera?: boolean;
}

export interface GizmoBinding<T> {
  /** The items shown */
  readonly items: Observable<readonly GizmoItem<T>[]>;
  /** An item's gizmo (none: not shown) */
  layout(item: GizmoItem<T>): GizmoLayout | undefined;
  /** The same meshes for both values: only the frame moves (a move, a turn), not rebuilt (a size) */
  sameGeometry(a: T, b: T): boolean;
  /** Equal values: a drag that ends here changed nothing (nothing committed) */
  equal(a: T, b: T): boolean;
  /** The value a drag on `handle` makes from `from` (the value at the press) */
  apply(from: T, handle: Handle, frame: DragFrame, dx: number, dy: number): T;
  /** A drag on `key` starts from `value` (what the binding remembers to tell a change under it) */
  begin?(key: string, value: T): void;
  /** Live, outside the state: the item drawn as `next`; at most one in flight */
  preview(key: string, from: T, next: T): Promise<void>;
  /** On release: `next` taken (one step); false: not taken, the controller cancels */
  commit(key: string, from: T, next: T): Promise<boolean>;
  /** The drag undone: the item back to what the state holds */
  cancel(key: string, from: T): void;
  /** Fires with an item's key when what its drag started from changed (an undo, another edit): the drag is cancelled */
  readonly invalidated?: Observable<string>;
}

export interface GizmoControllerOptions {
  readonly bodyColor?: Color;
  readonly bodyAlpha?: number;
  /** The handles' glow (emissive: the canvas's bloom picks it up) */
  readonly glow?: number;
}

interface Gizmo<T> {
  item: GizmoItem<T>;
  layout: GizmoLayout;
  body?: Representation.Any;
  handles: Representation.Any;
}

interface Drag<T> {
  readonly key: string;
  readonly handle: Handle;
  readonly from: T;
  readonly frame: DragFrame;
  /** Where it was pressed (CSS px) */
  readonly x0: number;
  readonly y0: number;
  next: T;
}

export class GizmoController<T> {
  /** The key of the item being dragged */
  readonly dragging = new BehaviorSubject<string | undefined>(undefined);

  readonly #canvas: Canvas3D;
  readonly #binding: GizmoBinding<T>;
  readonly #options: Required<GizmoControllerOptions>;
  readonly #set: VisualSet;
  readonly #subs: Subscription[] = [];
  readonly #gizmos = new Map<string, Gizmo<T>>();
  #items: readonly GizmoItem<T>[] = [];
  #hovered: { key: string; handle: Handle } | undefined;
  #drag: Drag<T> | undefined;
  /** The trackball's bindings while a handle holds the drags (restored whole: a partial set keeps only its part) */
  #bindings: Canvas3D['attribs']['trackball']['bindings'] | undefined;
  #syncing = 0;
  #disposed = false;

  constructor(canvas: Canvas3D, binding: GizmoBinding<T>, options: GizmoControllerOptions = {}) {
    this.#canvas = canvas;
    this.#binding = binding;
    this.#options = { bodyColor: Color(0x8fb8ec), bodyAlpha: 0.14, glow: 0.35, ...options };
    this.#set = new VisualSet(canvas);
    this.#subs.push(
      binding.items.subscribe((items) => {
        this.#items = items;
        void this.#sync();
      }),
      canvas.interaction.hover.subscribe(({ current }) => this.#hover(current)),
      canvas.input.drag.subscribe((e) => this.#onDrag(e.isStart, e.x, e.y, e.dx, e.dy)),
      canvas.input.interactionEnd.subscribe(() => void this.#end()),
      canvas.input.keyDown.subscribe((e) => {
        if (e.key === 'Escape' && this.#drag) this.cancel();
      }),
      canvas.didDraw.subscribe(() => this.#faceCamera()),
    );
    if (binding.invalidated) {
      this.#subs.push(
        binding.invalidated.subscribe((key) => {
          if (this.#drag?.key === key) this.cancel();
        }),
      );
    }
  }

  /** A representation of a gizmo (its clicks are the controller's: an app claims them) */
  owns(repr: Representation.Any | undefined): boolean {
    return !!repr && this.#set.has(repr);
  }

  /** The drag in progress, undone: the item back to what the state holds */
  cancel(): void {
    const d = this.#drag;
    if (!d) return;
    this.#drag = undefined;
    this.dragging.next(undefined);
    this.#binding.cancel(d.key, d.from);
    if (!this.#hovered) this.#restoreBindings();
    // the gizmos as the items are now (the drag's skipped its own, and they may have moved on under it)
    void this.#sync();
    this.#canvas.requestDraw();
  }

  dispose(): void {
    if (this.#disposed) return;
    this.cancel();
    this.#hovered = undefined;
    this.#restoreBindings();
    for (const s of this.#subs) s.unsubscribe();
    this.#set.dispose();
    this.#gizmos.clear();
    this.dragging.complete();
    this.#disposed = true;
  }

  /** The gizmos follow the items: made, moved or remade where an item changed, removed with it */
  async #sync() {
    if (this.#disposed) return;
    const version = ++this.#syncing;
    const items = this.#items;
    for (const [key, g] of [...this.#gizmos]) {
      if (items.some((i) => i.key === key)) continue;
      this.#remove(g);
      this.#gizmos.delete(key);
    }
    for (const item of items) {
      // the item being dragged is the drag's
      if (this.#drag?.key === item.key) continue;
      const g = this.#gizmos.get(item.key);
      if (g && this.#binding.equal(g.item.value, item.value)) {
        g.item = item;
        continue;
      }
      if (g && this.#binding.sameGeometry(g.item.value, item.value)) {
        const layout = this.#binding.layout(item);
        if (layout) {
          this.#place(g, item, layout);
          continue;
        }
      }
      const layout = this.#binding.layout(item);
      if (!layout) continue;
      const made = await this.#make(item, layout);
      if (version !== this.#syncing || this.#disposed) {
        this.#remove(made);
        return;
      }
      if (g) this.#remove(g);
      this.#gizmos.set(item.key, made);
    }
    this.#canvas.requestDraw();
  }

  async #make(item: GizmoItem<T>, layout: GizmoLayout): Promise<Gizmo<T>> {
    const set = this.#set;
    const { bodyColor, bodyAlpha, glow } = this.#options;
    const helper = { markerActions: MarkerActions.Highlighting };
    const handleShape = Shape.create(
      item.label,
      {},
      layout.handles,
      (g) => layout.groups[g]?.color ?? bodyColor,
      () => 1,
      (g) => `${item.label}: ${layout.groups[g]?.label ?? ''}`,
      Math.max(layout.groups.length, 1),
    );
    const handles = set.addShape(
      (_, s: Shape<Mesh>) => s,
      Mesh.Utils,
      handleShape,
      { alpha: 1, emissive: glow },
      helper,
    );
    const body = layout.body
      ? set.addShape(
          (_, s: Shape<Mesh>) => s,
          Mesh.Utils,
          Shape.create(
            item.label,
            {},
            layout.body,
            () => bodyColor,
            () => 1,
            () => item.label,
            1,
          ),
          { alpha: bodyAlpha, ignoreLight: true },
          { ...helper, pickable: false },
        )
      : undefined;
    const [h, b] = await Promise.all([handles, body]);
    const g: Gizmo<T> = {
      item,
      layout,
      handles: h as unknown as Representation.Any,
      ...(b ? { body: b as unknown as Representation.Any } : {}),
    };
    this.#setFrame(g);
    return g;
  }

  #remove(g: Gizmo<T>) {
    if (g.body) this.#set.remove(g.body);
    this.#set.remove(g.handles);
  }

  /** Where a gizmo is drawn: its frame, or (facing the camera) that place turned to it */
  #frameOf(layout: GizmoLayout): Mat4 {
    const eye = this.#canvas.camera.state.position;
    return layout.facesCamera ? facingCamera(layout.frame, eye) : layout.frame;
  }

  #setFrame(g: Gizmo<T>) {
    const t = this.#frameOf(g.layout);
    if (g.body) this.#set.setTransform(g.body, t);
    this.#set.setTransform(g.handles, t);
  }

  /** A gizmo moved and turned as a whole (its meshes the same) */
  #place(g: Gizmo<T>, item: GizmoItem<T>, layout: GizmoLayout) {
    g.item = item;
    g.layout = layout;
    this.#setFrame(g);
  }

  /** A gizmo drawn as `value` during a drag: moved, or remade when its meshes differ */
  async #show(g: Gizmo<T>, value: T) {
    const item = { ...g.item, value };
    if (this.#binding.sameGeometry(g.item.value, value)) {
      const layout = this.#binding.layout(item);
      if (layout) this.#place(g, item, layout);
      return;
    }
    const layout = this.#binding.layout(item);
    if (!layout) return;
    const made = await this.#make(item, layout);
    // remade meanwhile (a sync after the commit) or gone: this one is not drawn
    if (this.#disposed || this.#gizmos.get(item.key) !== g) {
      this.#remove(made);
      return;
    }
    this.#remove(g);
    this.#gizmos.set(item.key, made);
  }

  /** After a frame: the gizmos facing the camera turned to it again where it moved (one more frame, then still) */
  readonly #lastEye = Vec3.create(Number.NaN, 0, 0);
  #faceCamera() {
    const eye = this.#canvas.camera.state.position;
    if (Vec3.exactEquals(eye, this.#lastEye)) return;
    Vec3.copy(this.#lastEye, eye);
    for (const g of this.#gizmos.values()) if (g.layout.facesCamera) this.#setFrame(g);
  }

  #hover(current: Representation.Loci) {
    if (this.#drag) return;
    const entry = [...this.#gizmos].find(([, g]) => g.handles === current.repr);
    if (!entry || !ShapeGroup.isLoci(current.loci) || !current.loci.groups.length) {
      if (this.#hovered) {
        this.#hovered = undefined;
        this.#restoreBindings();
      }
      return;
    }
    const [key, g] = entry;
    const handle = g.layout.groups[OrderedSet.start(current.loci.groups[0].ids)]?.handle;
    if (!handle) return;
    this.#hovered = { key, handle };
    this.#takeDrags();
  }

  #onDrag(isStart: boolean, x: number, y: number, fx: number, fy: number) {
    if (isStart) {
      if (!this.#hovered || this.#drag) return;
      const [x0, y0] = pressPoint(this.#canvas, x, y, fx, fy);
      this.#begin(this.#hovered.key, this.#hovered.handle, x0, y0);
    }
    const d = this.#drag;
    if (!d) return;
    d.next = this.#binding.apply(d.from, d.handle, d.frame, x - d.x0, y - d.y0);
    this.#live();
  }

  #begin(key: string, handle: Handle, x0: number, y0: number) {
    const g = this.#gizmos.get(key);
    if (!g) return;
    const transform = this.#frameOf(g.layout);
    const at = Vec3.fromArray(Vec3(), transform, 12);
    // a handle with a pivot is grabbed from it (a rim's center), the rest from the item's origin
    const local = handle.kind === 'scale' || handle.kind === 'custom' ? handle.pivot : undefined;
    const pivot = local ? Vec3.transformMat4(Vec3(), local, transform) : at;
    const from = g.item.value;
    this.#binding.begin?.(key, from);
    this.#drag = { key, handle, from, frame: dragFrame(this.#canvas, at, pivot, x0, y0), x0, y0, next: from };
    this.dragging.next(key);
  }

  /**
   * The drag's value previewed and its gizmo drawn there, at most one update in flight: each round takes the drag as it
   * is then (the latest move; a new drag's first move too, while the last one's round finishes)
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
        await this.#binding.preview(d.key, d.from, next);
        const g = this.#gizmos.get(d.key);
        if (g && this.#drag === d) await this.#show(g, next);
        this.#canvas.requestDraw();
      } while (this.#again);
    };
    this.#running = run().finally(() => {
      this.#running = undefined;
    });
  }

  async #end() {
    const d = this.#drag;
    if (!d) return;
    // a preview still running would land after the commit's
    await this.#running;
    if (this.#drag !== d) return;
    if (this.#binding.equal(d.from, d.next)) {
      this.cancel();
      return;
    }
    this.#drag = undefined;
    this.dragging.next(undefined);
    if (!this.#hovered) this.#restoreBindings();
    // taken: the items follow (sync); not: the binding puts it back
    const done = await this.#binding.commit(d.key, d.from, d.next);
    if (!done) {
      this.#drag = d;
      this.cancel();
    }
  }

  #takeDrags() {
    if (this.#bindings) return;
    const bindings = { ...this.#canvas.attribs.trackball.bindings };
    this.#bindings = bindings;
    this.#canvas.setAttribs({
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
    this.#canvas.setAttribs({ trackball: { bindings } });
  }
}
