/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Ludovic Autin <ludovic.autin@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { OrderedSet } from '@molstar/core/data/int';
import { smoothstep } from '@molstar/core/math/interpolate';
import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import { Binding } from '@molstar/core/util/binding';
import { Color } from '@molstar/core/util/color';
import { MarkerActions } from '@molstar/core/util/marker-action';
import { Shape, ShapeGroup } from '@molstar/model/model/shape';
import { BehaviorSubject, type Observable, type Subscription } from 'rxjs';
import type { Canvas3D } from '../../canvas3d/canvas3d.js';
import type { GeometryUtils } from '../../geo/geometry/geometry.js';
import { Lines } from '../../geo/geometry/lines/lines.js';
import { LinesBuilder } from '../../geo/geometry/lines/lines-builder.js';
import { Points } from '../../geo/geometry/points/points.js';
import { Mesh } from '../../geo/geometry/mesh/mesh.js';
import type { GraphicsContext } from '../../gl/context.js';
import type { Representation } from '../../repr/representation.js';
import { ShapeRepresentation } from '../../repr/shape/representation.js';
import { Transparency } from '../../theme/transparency.js';
import { VisualSet } from '../visual-set.js';
import { dragFrame, grabFrom, pressPoint } from './drag-frame.js';
import type { DragFrame } from './gestures.js';
import type { Handle, HandleGroup } from './handle.js';

// Gizmos: handles in the 3D view for items an app edits (a clip object, a component's place), drawn outside the state
// tree (a `VisualSet`) and dragged live. What an item is and how it changes is the app's (`GizmoBinding`); the controller
// draws each item's gizmo, takes the trackball's drags while a handle is hovered (its bindings saved whole, once a
// canvas whatever number of controllers hold them, and restored when the last lets go: `TrackballDrags`), runs a drag's
// previews (one in flight, the latest wins), and commits on release: one step, no new drag until it is done. Esc, the
// binding saying what the drag started from changed (`invalidated`), the item going away, or a commit not taken (or
// failed) cancels: the binding puts the item back (after a preview in flight), the gizmos follow what the items are.
// Gizmos are in the overlay layer (above the scene, out of its bounds): the body x-ray, the handles on top by default;
// a handle turned edge-on fades out (not picked once hidden). During a drag a turn hides every other handle, a move the
// handles that do not move (the others dim), and guides show what the drag does: its axis, a ring's radius at the grab,
// a dashed line to the cursor, a grey copy of a moved handle where it started. Parts may show only in a state (a fill
// on hover, shafts during a move) or only their half toward the viewer (a ring's front).
// Only `Canvas3D`'s public API (no gl, no Canvas3D internals).

type TrackballBindings = Canvas3D['attribs']['trackball']['bindings'];

const tmpAxis = Vec3();
const Up = Vec3.create(0, 1, 0);

/**
 * The trackball's drags held by gizmo controllers, per canvas: the bindings saved by the first holder and put back by
 * the last (two controllers on one canvas, the pointer going straight from one's handle to the other's, would otherwise
 * save each other's emptied bindings)
 */
const TrackballDrags = {
  held: new WeakMap<Canvas3D, { saved: TrackballBindings; holders: Set<object> }>(),
  take(canvas: Canvas3D, holder: object) {
    const h = this.held.get(canvas);
    if (h) {
      h.holders.add(holder);
      return;
    }
    const saved = { ...canvas.attribs.trackball.bindings };
    this.held.set(canvas, { saved, holders: new Set([holder]) });
    canvas.setAttribs({
      trackball: {
        bindings: {
          ...saved,
          dragRotate: Binding.Empty,
          dragRotateZ: Binding.Empty,
          dragPan: Binding.Empty,
          dragFocus: Binding.Empty,
          dragFocusZoom: Binding.Empty,
        },
      },
    });
  },
  release(canvas: Canvas3D, holder: object) {
    const h = this.held.get(canvas);
    if (!h || !h.holders.delete(holder) || h.holders.size) return;
    this.held.delete(canvas);
    canvas.setAttribs({ trackball: { bindings: h.saved } });
  },
};

/** An item a binding shows a gizmo for: its key (stable through edits), its value, its name */
export interface GizmoItem<T> {
  readonly key: string;
  readonly value: T;
  readonly label: string;
}

/**
 * A part of an item's gizmo: a geometry in the item's frame, one group a handle (`GizmoLayout.groups[i]` is group i's).
 * Every handle and hit part is one representation whose loci are the gizmo's handles' (`ShapeBuilder.lociShape`): a
 * pick on any is the handle's (hover, highlight, drags), a fade of a handle fades it in every part.
 */
export interface GizmoPart {
  /**
   * `handle`: drawn, opaque, picked; `hit`: picked, not drawn (larger than what is drawn: the pick is drawn at a fraction
   * of the pixels, a thin shaft falls between them); `body`: translucent, not picked (a click needs alpha ≥ 0.5, and a
   * body covers what is behind it), not a handle
   */
  readonly role: 'handle' | 'hit' | 'body';
  /**
   * A mesh; lines a set width on screen (`GizmoControllerOptions.lineWidth`: shafts, rims); points a set size on
   * screen (`point`, else `GizmoControllerOptions.pointSize`: knobs at their places, round)
   */
  readonly geometry: Mesh | Lines | Points;
  /**
   * Points of their own size and look, facing the camera by nature (circles and discs around the view axis): `radius`
   * and `depth` in the frame's units, times the screen size (`screenSize`; a unit a CSS pixel without one); `depth`
   * toward the camera in depth only, negative away (`Points.Params.pointDepthOffset`: a hit area over or under every
   * other part); `ringWidth` a `ring`'s, a fraction of the radius; `front` twice the radius toward the camera (a knob
   * over the edges ending in it: a small point, not clamped to the depth range as `depth` is)
   */
  readonly point?: {
    readonly radius: number;
    readonly style?: 'circle' | 'ring';
    readonly ringWidth?: number;
    readonly depth?: number;
    readonly front?: boolean;
  };
  /**
   * CSS pixels a unit is on screen at any zoom (`Mesh.Params.screenSize`; a point's radius and depth in units of it);
   * none: the layout's. Not for lines (a set width on screen, their lengths the world's)
   */
  readonly screenSize?: number;
  /**
   * Shown only in a state, hidden otherwise: while the group `hover` is hovered and nothing dragged (a trackball's fill),
   * while a handle of a kind in `drag` is dragged (the arrows' shafts during a move)
   */
  readonly when?: { readonly hover?: number; readonly drag?: readonly Handle['kind'][] };
  /** Its alpha, 1 by default (a faint fill: below 0.5 it is not picked). Not for a body (`bodyAlpha`) */
  readonly alpha?: number;
  /**
   * Its half toward the viewer only (a ring's front): cut by a plane through the frame's origin across the view, set as
   * the camera or the gizmo moves; whole while a `rotate` handle of the gizmo is dragged
   */
  readonly front?: boolean;
}

/**
 * A representation of a gizmo beyond its parts (a label, a guide): built with the canvas's graphics, placed in the
 * gizmo's frame and layer, removed with it. A handle where its loci are `handles`' (`ShapeBuilder.lociShape`);
 * otherwise drawn only (its picks no handle's)
 */
export interface GizmoExtra {
  build(graphics: GraphicsContext, handles: Shape): Representation.Any;
  readonly props?: {};
  readonly data?: unknown;
}

/** An item's gizmo: its parts in its own frame, what each handle group does, the frame in the world */
export interface GizmoLayout {
  readonly parts: readonly GizmoPart[];
  readonly groups: readonly HandleGroup[];
  readonly extras?: readonly GizmoExtra[];
  readonly frame: Mat4;
  /** CSS pixels a unit of the parts is on screen, at any zoom (`Mesh.Params.screenSize`); none: world units */
  readonly screenSize?: number;
  /**
   * How the scene's depth applies to each part (`occludedAlpha`: 1 on top, 0 depth tested, between x-ray); by default
   * the body x-ray (`DefaultOccludedAlpha`) and the handles on top, so a handle inside a structure can be grabbed
   */
  readonly occludedAlpha?: { readonly body?: number; readonly handles?: number };
}

/** The parts' occluded alpha where a layout gives none */
export const DefaultOccludedAlpha = { body: 0.25, handles: 1 } as const;

/**
 * The transparency of the handles faded as they turn edge-on (a slide's axis toward the viewer, a plane across the
 * view; hidden and not picked at 1: `uPickingAlphaThreshold`), and of the others during a drag
 */
const Fade = {
  /** |cos| of a slide's axis and the view: fading from, hidden at */
  slide: [0.9, 0.98],
  /** |cos| of a plane's normal and the view: hidden at, fading until */
  plane: [0.08, 0.25],
  dimmed: 0.75,
  /** Steps of the transparency set (fewer state updates as the camera turns) */
  step: 0.05,
} as const;

/** The handles that turn an item and those that move it: a drag on one hides every other, on the other those not moving */
const Turning: ReadonlySet<Handle['kind']> = new Set(['rotate', 'roll', 'turn']);
const Moving: ReadonlySet<Handle['kind']> = new Set(['move', 'slide', 'plane']);
/** The handles whose drags draw guides (`guideLines`) */
const Guided: ReadonlySet<Handle['kind']> = new Set(['slide', 'plane', 'rotate', 'roll', 'turn']);
/** The handles whose drags go by steps, each from the last value (their path matters) */
const Stepped: ReadonlySet<Handle['kind']> = new Set(['turn', 'roll', 'rotate']);

/** What a drag frame's scale depends on of the camera: where it is, what it looks at, its field of view, its frame */
function viewKey(canvas: Canvas3D, out: number[]) {
  const { camera } = canvas;
  const { position: p, target: t, fov } = camera.state;
  out.length = 0;
  out.push(p[0], p[1], p[2], t[0], t[1], t[2], fov, camera.viewport.width, camera.viewport.height);
  return out;
}

function isSameView(a: readonly number[], b: readonly number[]) {
  for (let i = 0, il = a.length; i < il; ++i) if (a[i] !== b[i]) return false;
  return a.length === b.length;
}

/** A moved handle's copy where the drag started, and the guides */
const GhostColor = Color(0x9a9a9a);
const GhostAlpha = 0.6;
/** Grey, not white: read on a white background or structure as on a dark one */
const GuideColor = Color(0x8c8c8c);
/** CSS pixels a dash of a guide is, and its gap */
const DashPixels = 6;
/** CSS pixels the cursor's double arrow and the trackball's cross reach either way */
const ArrowPixels = 9;
/** How far (Å) behind a front part's origin its plane goes to show it whole (one clip object either way: no other shader) */
const WholeOffset = 1e6;

/** Where a drag's guides go: in the world, `perPixel` Å a CSS pixel at `origin` */
export interface GuideFrame {
  /** The item's frame at the press */
  readonly frame: Mat4;
  readonly origin: Vec3;
  /** How far the axis lines go either way (toward the viewer, at most half of `depth`) */
  readonly reach: number;
  /** The origin's distance from the eye along the view */
  readonly depth: number;
  /** The cursor at the origin's depth */
  readonly cursor: Vec3;
  /** The screen's right and up */
  readonly right: Vec3;
  readonly up: Vec3;
  readonly perPixel: number;
  /** Where a turn about an axis started (the grab on a ring, the press for a roll), and where that is now */
  readonly start?: Vec3;
  readonly current?: Vec3;
}

/**
 * A drag's guides in the world: a slide's axis, a plane's two axes, a ring's axis (across the view, through the
 * origin, `reach` either way but in front of the eye; none end-on; group 0, the handle's color); a turn about an axis
 * (a ring, a roll) its radius at the
 * start (group 1) and now (group 0), a dashed line from the origin to the cursor and a double arrow there, the way the
 * cursor turns it (group 1); a trackball's turn a cross at the origin, what it turns about (group 1). Built into
 * `into` where given (a shape representation keeps its geometry's buffers: the guides drawn are updated in place)
 */
export function guideLines(handle: Handle, f: GuideFrame, into?: Lines): Lines {
  const { frame, origin, reach, depth, cursor, right, up, perPixel, start, current } = f;
  const lines = LinesBuilder.create(64, 64, into);
  const toViewer = Vec3.cross(Vec3(), right, up);
  const along = (local: Vec3, group: number) => {
    const a = Vec3.normalize(Vec3(), Vec3.transformDirection(Vec3(), local, frame));
    const c = Vec3.dot(a, toViewer);
    // end-on: a dot, left out (as a slide is hidden, `Fade.slide`)
    if (Math.abs(c) >= Fade.slide[1]) return;
    // toward the viewer no further than half the origin's depth: past the eye it streaks across the view; on the
    // screen it ends tan(angle to the view) from the origin, short as the axis turns end-on
    const near = Math.abs(c) > 1e-6 ? Math.min(reach, Math.max(0, (0.5 * depth) / Math.abs(c))) : reach;
    lines.addVec(
      Vec3.scaleAndAdd(Vec3(), origin, a, c < 0 ? -near : -reach),
      Vec3.scaleAndAdd(Vec3(), origin, a, c > 0 ? near : reach),
      group,
    );
  };
  switch (handle.kind) {
    case 'slide':
      along(handle.axis, 0);
      break;
    case 'plane': {
      const n = Vec3.normalize(Vec3(), handle.normal);
      for (const e of [Vec3.create(1, 0, 0), Vec3.create(0, 1, 0), Vec3.create(0, 0, 1)]) {
        if (Math.abs(Vec3.dot(e, n)) < 0.5) along(e, 0);
      }
      break;
    }
    case 'rotate':
    case 'roll': {
      if (handle.kind === 'rotate') along(handle.axis, 0);
      if (start) lines.addVec(origin, start, 1);
      if (current) lines.addVec(origin, current, 0);
      const count = Math.max(1, Math.round(Vec3.distance(origin, cursor) / (DashPixels * perPixel)));
      lines.addFixedCountDashes(origin, cursor, count, 1);
      // a double arrow at the cursor, across the dashes on the screen
      const out = Vec3.sub(Vec3(), cursor, origin);
      if (Vec3.magnitude(out) > 1e-9) {
        Vec3.normalize(out, out);
        const way = Vec3.normalize(Vec3(), Vec3.cross(Vec3(), toViewer, out));
        const arm = ArrowPixels * perPixel;
        const a = Vec3.scaleAndAdd(Vec3(), cursor, way, arm);
        const b = Vec3.scaleAndAdd(Vec3(), cursor, way, -arm);
        lines.addVec(a, b, 1);
        for (const [tip, sign] of [
          [a, 1],
          [b, -1],
        ] as const) {
          const back = Vec3.scaleAndAdd(Vec3(), tip, way, -sign * arm * 0.5);
          lines.addVec(tip, Vec3.scaleAndAdd(Vec3(), back, out, arm * 0.4), 1);
          lines.addVec(tip, Vec3.scaleAndAdd(Vec3(), back, out, -arm * 0.4), 1);
        }
      }
      break;
    }
    case 'turn': {
      const arm = ArrowPixels * perPixel;
      for (const dir of [right, up]) {
        lines.addVec(Vec3.scaleAndAdd(Vec3(), origin, dir, -arm), Vec3.scaleAndAdd(Vec3(), origin, dir, arm), 1);
      }
      break;
    }
  }
  return lines.getLines();
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
  /**
   * The value a drag on `handle` makes from `from`: the value at the press and the drag's total, for a `turn`, a
   * `roll` and a `rotate` the last value and the step (their path matters: `gestures.turn`, `roll`, `rotateAbout`;
   * `frame.grab` the last cursor from the origin)
   */
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
  /** The width of lines parts (CSS pixels) */
  readonly lineWidth?: number;
  /** The size of points parts (CSS pixels across) */
  readonly pointSize?: number;
  /**
   * A slide fades out as its axis turns toward the viewer, hidden and not picked at the end (by default). Off: always
   * drawn; seen end-on, a pick on it is the gizmo's move handle where it has one (its head covers the knob)
   */
  readonly fadeSlides?: boolean;
}

interface Gizmo<T> {
  item: GizmoItem<T>;
  layout: GizmoLayout;
  /** Every representation of it: its parts, its extras */
  readonly all: readonly Representation.Any[];
  /** Those picked as its handles and faded with them: the handle and hit parts, the extras (`#hover` tells by loci) */
  readonly handles: readonly Representation.Any[];
  /** The handles' shape: the loci of every handle and hit part (`ShapeBuilder.lociShape`), its geometry none */
  readonly shape: Shape<Mesh>;
  /** The representations of its parts, as `layout.parts` */
  readonly parts: readonly Representation.Any[];
  /** The handles' transparency per group as last set (`#fade`) */
  faded: readonly number[];
  /** Whether each part shown in a state (`GizmoPart.when`) is, as last set */
  readonly shown: boolean[];
  /** The front parts' clip as last set (`#front`) */
  front: string;
}

interface Drag<T> {
  readonly key: string;
  readonly handle: Handle;
  /** The handle's group (a layout made again during the drag has handles of its own) */
  readonly group: number;
  readonly from: T;
  /** The frame the moves are measured in: at the press, or at the last rebase (`base`) */
  frame: DragFrame;
  /** Where it was pressed (CSS px) */
  readonly x0: number;
  readonly y0: number;
  /**
   * What a move (not a step) applies to, and from where (CSS px): the press, or where the camera last moved during the
   * drag (a zoom): from then on measured in a frame made there (`#rebase`), the item kept under the cursor
   */
  base: T;
  bx: number;
  by: number;
  /** The camera the frame was made with (`viewKey`) */
  readonly view: number[];
  /** Where the last move was (CSS px): a `turn` goes by steps */
  x: number;
  y: number;
  next: T;
  /** The frame at the press: where the guides go */
  readonly pressFrame: Mat4;
  /** Its ghosts and guides (`#adorn`), once made; done: the drag ended (whatever is made then is removed) */
  adorned: { readonly reprs: Representation.Any[]; guides?: Representation.Any; done: boolean };
}

export class GizmoController<T> {
  /** The key of the item being dragged */
  readonly dragging = new BehaviorSubject<string | undefined>(undefined);

  readonly #canvas: Canvas3D;
  readonly #binding: GizmoBinding<T>;
  #options: Required<GizmoControllerOptions>;
  readonly #set: VisualSet;
  readonly #subs: Subscription[] = [];
  readonly #gizmos = new Map<string, Gizmo<T>>();
  #items: readonly GizmoItem<T>[] = [];
  #hovered: { key: string; handle: Handle; group: number; position?: Vec3 } | undefined;
  #drag: Drag<T> | undefined;
  /** This controller holds the trackball's drags (`TrackballDrags`) */
  #holding = false;
  /** A release's commit in flight: no drag starts until it is done (it would start from what the commit replaces) */
  #committing = false;
  /** A cancel's revert waiting on a preview in flight (`#undo`): no drag starts until it is done (it would start from the
   * preview it reverts, and the revert land in it) */
  #undoing = false;
  #syncing = 0;
  #disposed = false;

  constructor(canvas: Canvas3D, binding: GizmoBinding<T>, options: GizmoControllerOptions = {}) {
    this.#canvas = canvas;
    this.#binding = binding;
    this.#options = {
      bodyColor: Color(0x8fb8ec),
      bodyAlpha: 0.14,
      glow: 0.35,
      lineWidth: 4,
      fadeSlides: true,
      pointSize: 12,
      ...options,
    };
    this.#set = new VisualSet(canvas);
    this.#subs.push(
      binding.items.subscribe((items) => {
        this.#items = items;
        void this.#sync();
      }),
      canvas.interaction.hover.subscribe(({ current, position }) => this.#hover(current, position)),
      canvas.input.drag.subscribe((e) => this.#onDrag(e.isStart, e.x, e.y, e.dx, e.dy)),
      canvas.input.interactionEnd.subscribe(() => void this.#end()),
      canvas.input.keyDown.subscribe((e) => {
        if (e.key === 'Escape' && this.#drag) this.cancel();
      }),
      canvas.didDraw.subscribe(() => this.#fadeOnCameraMove()),
    );
    if (binding.invalidated) {
      this.#subs.push(
        binding.invalidated.subscribe((key) => {
          if (this.#drag?.key === key) this.cancel();
        }),
      );
    }
  }

  /** Whether slides fade out toward the viewer (`GizmoControllerOptions.fadeSlides`), from now on */
  setFadeSlides(fadeSlides: boolean) {
    this.#options = { ...this.#options, fadeSlides };
    this.#fadeAll();
    this.#canvas.requestDraw();
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
    this.#unadorn(d);
    this.#fadeAll();
    this.#forgetGrab();
    this.#undo(d);
    if (!this.#hovered) this.#restoreBindings();
  }

  /** The binding puts a drag's item back: after a preview in flight (it would land after the revert otherwise) */
  #undo(d: Drag<T>) {
    const running = this.#running;
    const undo = () => {
      this.#binding.cancel(d.key, d.from);
      // the gizmos as the items are now (the drag's skipped its own, and they may have moved on under it)
      void this.#sync();
      this.#canvas.requestDraw();
    };
    if (running) {
      this.#undoing = true;
      void running.then(() => {
        this.#undoing = false;
        undo();
      });
    } else undo();
  }

  /** After a drag: the grab's point is the last hover's, stale until the pointer moves again (a new hover) */
  #forgetGrab() {
    if (this.#hovered?.position) {
      const { key, handle, group } = this.#hovered;
      this.#hovered = { key, handle, group };
    }
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
    // the item dragged went away (removed, another target): the drag with it
    const dragged = this.#drag?.key;
    if (dragged !== undefined && !items.some((i) => i.key === dragged)) this.cancel();
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
      const layout = this.#binding.layout(item);
      // not shown (now): the one drawn goes
      if (!layout) {
        if (g) {
          this.#remove(g);
          this.#gizmos.delete(item.key);
        }
        continue;
      }
      if (g && this.#binding.sameGeometry(g.item.value, item.value)) {
        this.#place(g, item, layout);
        continue;
      }
      const made = await this.#make(item, layout);
      if (version !== this.#syncing || this.#disposed) {
        this.#remove(made);
        return;
      }
      // a drag's own remake took its place meanwhile (`#show`): that one stays
      if (this.#gizmos.get(item.key) !== g || this.#drag?.key === item.key) {
        this.#remove(made);
        continue;
      }
      if (g) this.#remove(g);
      this.#gizmos.set(item.key, made);
    }
    this.#canvas.requestDraw();
  }

  async #make(item: GizmoItem<T>, layout: GizmoLayout): Promise<Gizmo<T>> {
    const { bodyColor } = this.#options;
    // the loci of every handle and hit part: its labels the tooltips', no geometry of its own
    const handleShape = Shape.create(
      item.label,
      {},
      Mesh.createEmpty(),
      (g) => layout.groups[g]?.color ?? bodyColor,
      () => 1,
      (g) => `${item.label}: ${layout.groups[g]?.label ?? ''}`,
      Math.max(layout.groups.length, 1),
    );
    const parts = layout.parts.map((part) => this.#drawPart(item, layout, part, handleShape));
    // above the scene, out of its bounds (`Representation.State.layer`)
    const helper = { markerActions: MarkerActions.Highlighting, layer: 'overlay' as const };
    const extras = (layout.extras ?? []).map((e) =>
      this.#set.create((graphics) => e.build(graphics, handleShape), e.props, e.data, helper),
    );
    const built = await Promise.all([...parts, ...extras]);
    const g: Gizmo<T> = {
      item,
      layout,
      all: built,
      parts: built.slice(0, parts.length),
      handles: built.filter((_, i) => i >= parts.length || layout.parts[i].role !== 'body'),
      shape: handleShape,
      faded: [],
      shown: layout.parts.map(() => false),
      front: '',
    };
    this.#setFrame(g);
    return g;
  }

  /**
   * A part as a shape representation in the overlay: a body its own loci, a handle or a hit area the handles'
   * (`handleShape`); a `ghost` (`#ghosts`) grey and not picked, its `state` given
   */
  #drawPart(
    item: GizmoItem<T>,
    layout: GizmoLayout,
    part: GizmoPart,
    handleShape: Shape,
    ghost?: Partial<Representation.State>,
  ): Promise<Representation.Any> {
    const { bodyColor, bodyAlpha, glow, lineWidth, pointSize } = this.#options;
    const occludedAlpha = { ...DefaultOccludedAlpha, ...layout.occludedAlpha };
    const groupCount = Math.max(layout.groups.length, 1);
    const body = part.role === 'body';
    const color = ghost
      ? () => GhostColor
      : body
        ? () => bodyColor
        : (g: number) => layout.groups[g]?.color ?? bodyColor;
    const state: Partial<Representation.State> = {
      // above the scene, out of its bounds (`Representation.State.layer`)
      markerActions: MarkerActions.Highlighting,
      layer: 'overlay',
      ...(body ? { pickable: false } : part.role === 'hit' ? { pickOnly: true } : {}),
      ...(part.when ? { visible: false } : {}),
      ...ghost,
    };
    const common = {
      alpha: body ? bodyAlpha : ghost ? GhostAlpha : (part.alpha ?? 1),
      occludedAlpha: body ? occludedAlpha.body : occludedAlpha.handles,
      ...(part.front && !ghost ? { clip: this.#frontClip(layout, false) } : {}),
    };
    const draw = <G extends Mesh | Lines | Points>(
      utils: GeometryUtils<G>,
      geometry: G,
      size: number,
      props: object,
    ) => {
      const shape = Shape.create(
        item.label,
        {},
        geometry,
        color,
        () => size,
        () => (body ? item.label : ''),
        body ? 1 : groupCount,
      );
      return this.#set.create(
        (graphics) =>
          ShapeRepresentation((_, s: Shape<G>) => s, utils, {
            graphics,
            ...(body ? {} : { lociShape: () => handleShape }),
          }),
        { ...common, ...props } as never,
        shape,
        state,
      );
    };
    const drawn = part.role === 'handle' && !ghost;
    const glowing = drawn ? { emissive: glow } : {};
    const screenSize = part.screenSize ?? layout.screenSize ?? 0;
    const g = part.geometry;
    switch (g.kind) {
      case 'mesh':
        return draw(Mesh.Utils, g, 1, { screenSize, ignoreLight: true, ...(drawn ? { emissive: glow } : {}) });
      // the thin parts, as wide on screen at any zoom
      case 'lines':
        return draw(Lines.Utils, g, lineWidth, { ...glowing, sizeFactor: 1, lineSizeAttenuation: false });
      // knobs, circles and discs facing the camera, as large on screen at any zoom
      case 'points': {
        const p = part.point;
        const perUnit = screenSize || 1;
        return draw(Points.Utils, g, p ? 2 * p.radius * perUnit : pointSize, {
          ...glowing,
          sizeFactor: 1,
          pointSizeAttenuation: false,
          pointStyle: p?.style ?? 'circle',
          pointRingWidth: p?.ringWidth ?? 1,
          pointDepthOffset: (p?.depth ?? 0) * perUnit,
          // a knob in front of the lines ending in it, drawn and picked (the default knobs'; a part's where it says so)
          pointFront: p ? !!p.front : true,
        });
      }
    }
  }

  /**
   * The front parts' clip: a plane through the gizmo's origin across the view, the half beyond it cut away (a plane
   * clips the side its normal, local +Y, points to); `whole`: the plane pushed far behind (`WholeOffset`)
   */
  #frontClip(layout: GizmoLayout, whole: boolean) {
    const s = this.#canvas.camera.state;
    const origin = Vec3.fromArray(Vec3(), layout.frame, 12);
    const away = Vec3.sub(Vec3(), s.mode === 'orthographic' ? s.target : origin, s.position);
    if (Vec3.magnitude(away) < 1e-9) Vec3.sub(away, s.target, s.position);
    Vec3.normalize(away, away);
    const position = whole ? Vec3.scaleAndAdd(Vec3(), origin, away, WholeOffset) : origin;
    const axis = Vec3.cross(Vec3(), Up, away);
    const cos = Math.max(-1, Math.min(1, Vec3.dot(Up, away)));
    const rotation =
      Vec3.magnitude(axis) < 1e-9
        ? { axis: Vec3.create(1, 0, 0), angle: cos > 0 ? 0 : 180 }
        : { axis: Vec3.normalize(axis, axis), angle: (Math.acos(cos) * 180) / Math.PI };
    const object = {
      type: 'plane' as const,
      invert: false,
      position,
      rotation,
      scale: Vec3.create(1, 1, 1),
      transform: Mat4.identity(),
    };
    return { variant: 'pixel' as const, objects: [object] };
  }

  /** The front parts cut as the camera and the gizmo are now (only where that changed) */
  #front(g: Gizmo<T>) {
    if (!g.layout.parts.some((p) => p.front)) return;
    const d = this.#drag;
    // a ring turned is drawn whole
    const clip = this.#frontClip(g.layout, d?.key === g.item.key && d.handle.kind === 'rotate');
    const { position: p, rotation: r } = clip.objects[0];
    const key = [p[0], p[1], p[2], r.axis[0], r.axis[1], r.axis[2], r.angle].map((x) => x.toPrecision(6)).join();
    if (key === g.front) return;
    g.front = key;
    g.layout.parts.forEach((part, i) => {
      if (part.front) void this.#set.update(g.parts[i], { clip } as never);
    });
  }

  /** The parts shown in a state (`GizmoPart.when`), as the hover and the drag are now */
  #reveal(g: Gizmo<T>) {
    const d = this.#drag;
    const h = this.#hovered;
    g.layout.parts.forEach((part, i) => {
      const when = part.when;
      if (!when) return;
      const shown = d
        ? d.key === g.item.key && !!when.drag?.includes(d.handle.kind)
        : h?.key === g.item.key && when.hover === h.group;
      if (g.shown[i] === shown) return;
      g.shown[i] = shown;
      this.#set.setState(g.parts[i], { visible: shown });
    });
  }

  /** A drag's ghosts (a move's) and guides (`Guided`), made once and removed when it ends (`#unadorn`) */
  async #adorn(d: Drag<T>, g: Gizmo<T>) {
    const kind = d.handle.kind;
    const made: Promise<Representation.Any>[] = Moving.has(kind) ? this.#ghosts(g, d) : [];
    const guides = Guided.has(kind) ? this.#guides(d, g) : undefined;
    const reprs = await Promise.all(guides ? [...made, guides] : made);
    if (d.adorned.done || this.#disposed) {
      for (const r of reprs) this.#set.remove(r);
      return;
    }
    d.adorned.reprs.push(...reprs);
    if (guides) {
      d.adorned.guides = reprs[reprs.length - 1];
      // moved meanwhile
      this.#guide(d);
    }
  }

  #unadorn(d: Drag<T>) {
    d.adorned.done = true;
    for (const r of d.adorned.reprs) this.#set.remove(r);
    d.adorned.reprs.length = 0;
    d.adorned.guides = undefined;
  }

  /** Grey copies of the dragged handle where the drag started: the parts drawn at rest, the other groups hidden */
  #ghosts(g: Gizmo<T>, d: Drag<T>): Promise<Representation.Any>[] {
    const { item, layout } = g;
    const shape = Shape.create(
      item.label,
      {},
      Mesh.createEmpty(),
      () => GhostColor,
      () => 1,
      () => '',
      Math.max(layout.groups.length, 1),
    );
    const others = layout.groups.map((_, i) => i).filter((i) => i !== d.group);
    const ids = OrderedSet.ofSortedArray(others);
    const transparency = others.length
      ? Transparency('group-loci', [{ loci: ShapeGroup.Loci(shape, [{ ids, instance: 0 }]), value: 1 }])
      : Transparency.Empty;
    const state = { pickable: false, transform: Mat4.clone(layout.frame), transparency };
    return layout.parts
      .filter((p) => p.role === 'handle' && !p.when)
      .map((p) => this.#drawPart(item, layout, p, shape, state));
  }

  /** The guides' representation: lines in the world, not picked, as wide on screen at any zoom */
  #guides(d: Drag<T>, g: Gizmo<T>) {
    const color = g.layout.groups[d.group]?.color ?? GuideColor;
    // each move's into the last one's geometry (`guideLines`' `into`)
    const getShape = (_: unknown, f: GuideFrame, __: unknown, last?: Shape<Lines>) =>
      Shape.create(
        'guides',
        {},
        guideLines(d.handle, f, last?.geometry),
        (group) => (group === 0 ? color : GuideColor),
        () => this.#options.lineWidth * 0.75,
        () => '',
        2,
      );
    return this.#set.addShape(
      getShape,
      Lines.Utils,
      this.#guideFrame(d),
      { sizeFactor: 1, lineSizeAttenuation: false, occludedAlpha: 1 },
      { layer: 'overlay', pickable: false },
    );
  }

  /** Where the guides go as the drag is now */
  #guideFrame(d: Drag<T>): GuideFrame {
    const origin = Vec3.fromArray(Vec3(), d.pressFrame, 12);
    const eye = this.#canvas.camera.state.position;
    const { right, up, perPixel, grabPoint } = d.frame;
    // a point on the screen at the origin's depth
    const atDepth = (x: number, y: number) => {
      const [cx, cy] = grabFrom(this.#canvas, origin, x, y);
      const p = Vec3.scaleAndAdd(Vec3(), origin, right, cx * perPixel);
      return Vec3.scaleAndAdd(p, p, up, -cy * perPixel);
    };
    const reach = Math.max(2 * Vec3.distance(eye, origin), 1);
    const depth = -Vec3.dot(Vec3.sub(Vec3(), origin, eye), Vec3.cross(Vec3(), right, up));
    const start = d.handle.kind === 'roll' ? atDepth(d.x0, d.y0) : grabPoint;
    // the start turned as the item is now: the frame now over the frame at the press
    const g = this.#gizmos.get(d.key);
    const now = g ? g.layout.frame : d.pressFrame;
    const turned = Mat4.mul(Mat4(), now, Mat4.invert(Mat4(), d.pressFrame));
    const current = start ? Vec3.transformMat4(Vec3(), start, turned) : undefined;
    return {
      frame: d.pressFrame,
      origin,
      reach,
      depth,
      cursor: atDepth(d.x, d.y),
      right,
      up,
      perPixel,
      ...(start ? { start } : {}),
      ...(current ? { current } : {}),
    };
  }

  #guide(d: Drag<T>) {
    const guides = d.adorned.guides;
    if (guides) void this.#set.update(guides, {}, this.#guideFrame(d));
  }

  #remove(g: Gizmo<T>) {
    for (const r of g.all) this.#set.remove(r);
  }

  #setFrame(g: Gizmo<T>) {
    for (const r of g.all) this.#set.setTransform(r, g.layout.frame);
    this.#fade(g);
  }

  /**
   * The handles' transparency per group: a slide or a plane faded as it turns edge-on to the view (not the one being
   * dragged: it stays under the cursor); during a drag of the gizmo the others hidden by kind (`Turning`, `Moving`) or
   * dimmed. Its front parts cut and its parts shown by state follow (`#front`, `#reveal`)
   */
  #fade(g: Gizmo<T>) {
    this.#front(g);
    this.#reveal(g);
    const d = this.#drag;
    const own = d?.key === g.item.key;
    const dragged = own ? d.group : -1;
    const kind = own ? d.handle.kind : undefined;
    const faded = g.layout.groups.map(({ handle }, i) => {
      if (i === dragged) return 0;
      // a turn hides every other handle, a move those that do not move (a turn's, a size's)
      if (kind && (Turning.has(kind) || (Moving.has(kind) && !Moving.has(handle.kind)))) return 1;
      let value = d ? Fade.dimmed : 0;
      if ((handle.kind === 'slide' && this.#options.fadeSlides) || handle.kind === 'plane') {
        const cos = this.#towardViewer(g, handle.kind === 'slide' ? handle.axis : handle.normal);
        const [a, b] = handle.kind === 'slide' ? Fade.slide : Fade.plane;
        const edgeOn = handle.kind === 'slide' ? smoothstep(a, b, cos) : 1 - smoothstep(a, b, cos);
        value = Math.max(value, edgeOn);
      }
      return Math.round(value / Fade.step) * Fade.step;
    });
    if (faded.length === g.faded.length && faded.every((v, i) => v === g.faded[i])) return;
    g.faded = faded;
    // one layer a value, over its groups
    const byValue = new Map<number, number[]>();
    faded.forEach((v, i) => {
      if (v > 0) byValue.set(v, [...(byValue.get(v) ?? []), i]);
    });
    const layers = [...byValue].map(([value, ids]) => ({
      loci: ShapeGroup.Loci(g.shape, [{ ids: OrderedSet.ofSortedArray(ids), instance: 0 }]),
      value,
    }));
    // on the handles' loci, which every part takes: a handle faded out is not picked by its hit area either
    const transparency = layers.length ? Transparency('group-loci', layers) : Transparency.Empty;
    for (const r of g.handles) this.#set.setState(r, { transparency });
  }

  /** |cos| of a local direction of `g` and the view */
  #towardViewer(g: Gizmo<T>, local: Vec3) {
    const s = this.#canvas.camera.state;
    const t = g.layout.frame;
    const view = Vec3();
    if (s.mode === 'orthographic') Vec3.sub(view, s.target, s.position);
    else Vec3.sub(view, Vec3.fromArray(view, t, 12), s.position);
    Vec3.normalize(view, view);
    return Math.abs(Vec3.dot(Vec3.normalize(tmpAxis, Vec3.transformDirection(tmpAxis, local, t)), view));
  }

  /** What a pick on `group` grabs: its handle, but for an unfaded slide seen end-on the move handle, which it covers */
  #grabbedGroup(g: Gizmo<T>, group: number) {
    const handle = g.layout.groups[group]?.handle;
    if (handle?.kind !== 'slide' || this.#options.fadeSlides) return group;
    if (this.#towardViewer(g, handle.axis) < Fade.slide[1]) return group;
    const move = g.layout.groups.findIndex((h) => h.handle.kind === 'move');
    return move >= 0 ? move : group;
  }

  #fadeAll() {
    for (const g of this.#gizmos.values()) this.#fade(g);
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

  /**
   * After a frame: where the camera moved, the handles faded as they now face it (one more frame, then still); the
   * gizmos facing it are turned by the vertex stage
   */
  readonly #lastEye = Vec3.create(Number.NaN, 0, 0);
  readonly #lastTarget = Vec3();
  #fadeOnCameraMove() {
    const { position: eye, target } = this.#canvas.camera.state;
    if (Vec3.exactEquals(eye, this.#lastEye) && Vec3.exactEquals(target, this.#lastTarget)) return;
    Vec3.copy(this.#lastEye, eye);
    Vec3.copy(this.#lastTarget, target);
    this.#fadeAll();
  }

  #hover(current: Representation.Loci, position?: Vec3) {
    if (this.#drag) return;
    // a pick on the handles' loci (an extra of other loci is drawn only)
    const loci = current.loci;
    const entry = ShapeGroup.isLoci(loci)
      ? [...this.#gizmos].find(([, g]) => loci.shape === g.shape && !!current.repr && g.handles.includes(current.repr))
      : undefined;
    if (!entry || !ShapeGroup.isLoci(loci) || !loci.groups.length) {
      if (this.#hovered) {
        this.#hovered = undefined;
        this.#restoreBindings();
        this.#revealAll();
      }
      return;
    }
    const [key, g] = entry;
    const group = this.#grabbedGroup(g, OrderedSet.start(loci.groups[0].ids));
    const handle = g.layout.groups[group]?.handle;
    if (!handle) return;
    this.#hovered = { key, handle, group, ...(position ? { position: Vec3.clone(position) } : {}) };
    this.#takeDrags();
    this.#revealAll();
  }

  #revealAll() {
    for (const g of this.#gizmos.values()) this.#reveal(g);
  }

  #onDrag(isStart: boolean, x: number, y: number, fx: number, fy: number) {
    if (isStart) {
      if (!this.#hovered || this.#drag || this.#committing || this.#undoing) return;
      const [x0, y0] = pressPoint(this.#canvas, x, y, fx, fy);
      this.#begin(this.#hovered, x0, y0);
    }
    const d = this.#drag;
    if (!d) return;
    // a turn, a roll and a ring by steps; their grab the last cursor from the origin
    const g = d.frame.grab;
    const stepFrame = g ? { ...d.frame, grab: [g[0] + d.x - d.x0, g[1] + d.y - d.y0] as const } : d.frame;
    const stepped = Stepped.has(d.handle.kind);
    // the camera moved (a zoom): the moves from here on in a frame made here (the steps go by angles, not by Å)
    if (!stepped && !isSameView(d.view, viewKey(this.#canvas, this.#viewNow))) this.#rebase(d);
    d.next = stepped
      ? this.#binding.apply(d.next, d.handle, stepFrame, x - d.x, y - d.y)
      : this.#binding.apply(d.base, d.handle, d.frame, x - d.bx, y - d.by);
    d.x = x;
    d.y = y;
    this.#guide(d);
    this.#live();
  }

  #begin(
    { key, handle, group, position: grabPoint }: { key: string; handle: Handle; group: number; position?: Vec3 },
    x0: number,
    y0: number,
  ) {
    const g = this.#gizmos.get(key);
    if (!g) return;
    const transform = g.layout.frame;
    const from = g.item.value;
    this.#binding.begin?.(key, from);
    const frame = { ...this.#dragFrameAt(transform, handle, x0, y0), ...(grabPoint ? { grabPoint } : {}) };
    const adorned = { reprs: [], done: false };
    const d: Drag<T> = {
      key,
      handle,
      group,
      from,
      frame,
      x0,
      y0,
      x: x0,
      y: y0,
      next: from,
      base: from,
      bx: x0,
      by: y0,
      view: viewKey(this.#canvas, []),
      pressFrame: transform,
      adorned,
    };
    this.#drag = d;
    this.dragging.next(key);
    this.#fadeAll();
    void this.#adorn(d, g);
  }

  /** A drag's frame for an item at `transform`, pressed at `x`, `y`: a handle with a pivot grabbed from it */
  #dragFrameAt(transform: Mat4, handle: Handle, x: number, y: number) {
    const at = Vec3.fromArray(Vec3(), transform, 12);
    // a handle with a pivot is grabbed from it (a rim's center), the rest from the item's origin
    const local = handle.kind === 'scale' || handle.kind === 'custom' ? handle.pivot : undefined;
    const pivot = local ? Vec3.transformMat4(Vec3(), local, transform) : at;
    return dragFrame(this.#canvas, at, pivot, x, y);
  }

  readonly #viewNow: number[] = [];

  /** The moves from the last one on in a frame made where the item is now, as the camera is now (`Drag.base`) */
  #rebase(d: Drag<T>) {
    const g = this.#gizmos.get(d.key);
    const transform = g ? g.layout.frame : d.pressFrame;
    const { grabPoint } = d.frame;
    d.frame = { ...this.#dragFrameAt(transform, d.handle, d.x, d.y), ...(grabPoint ? { grabPoint } : {}) };
    d.base = d.next;
    d.bx = d.x;
    d.by = d.y;
    viewKey(this.#canvas, d.view as number[]);
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
        if (g && this.#drag === d) {
          await this.#show(g, next);
          // the radius now where the item is now
          this.#guide(d);
        }
        this.#canvas.requestDraw();
      } while (this.#again);
    };
    // a preview that fails is the binding's to say; the drag goes on (its release commits or cancels as ever)
    this.#running = run()
      .catch((e) => console.warn('[gizmo] a preview failed:', e))
      .finally(() => {
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
    this.#unadorn(d);
    this.#fadeAll();
    this.#forgetGrab();
    if (!this.#hovered) this.#restoreBindings();
    // taken: the items follow (sync); not taken, or failed: the binding puts it back (this drag's item, not a later one's)
    this.#committing = true;
    let done = false;
    try {
      done = await this.#binding.commit(d.key, d.from, d.next);
    } catch (e) {
      console.warn('[gizmo] a commit failed:', e);
    } finally {
      this.#committing = false;
    }
    if (!done && !this.#disposed) this.#undo(d);
  }

  #takeDrags() {
    if (this.#holding) return;
    this.#holding = true;
    TrackballDrags.take(this.#canvas, this);
  }

  #restoreBindings() {
    if (!this.#holding || this.#drag) return;
    this.#holding = false;
    TrackballDrags.release(this.#canvas, this);
  }
}
