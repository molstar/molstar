/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type Placement, PlacePoint, type Hit as RenderHit } from '@molstar/seqx-core/render';
import { cancelFrameCallback, DomHost, requestFrameCallback } from '@molstar/ui/dom';
import type { Axis } from '@molstar/seqx-core/view';
import {
  type KeyChord,
  matchesChord,
  type ModifierState,
  modifiersOf,
  parseChord,
} from '@molstar/seqx-viewer/bindings';
import { changeKindOf } from './commands.js';
import {
  type HoverBehavior,
  HoverProvider,
  isDragSelect,
  type NavigateBehavior,
  NavigateProvider,
  type SelectBehavior,
  SelectProvider,
  type WindowBehavior,
  WindowProvider,
} from '@molstar/seqx-viewer/behaviors';
import {
  type Anchor,
  axisTotal,
  type Emphasized,
  type Hit,
  type HitTarget,
  hitTarget,
  isNothing,
  nearestOffsetOf,
  cellOffset,
  offsetOf,
} from '@molstar/seqx-viewer/models';
import type { SeqViewer } from './viewer.js';

// The viewer's input (_spec/README.md "Input"): pointer, wheel and key events in canvas coordinates (CSS px from the
// canvas's top left), DOM-free, so the mounted UI only translates events and a test drives the same calls. Every
// gesture goes through the binding table: what a drag, a wheel or a pinch does is the table's, the behaviors do it.
// Per pointer move nothing is allocated: one hit record and one modifiers record are reused.

/** The modifier keys of an event, and where it is. */
export interface KeyState {
  readonly shiftKey?: boolean;
  readonly metaKey?: boolean;
  readonly ctrlKey?: boolean;
  readonly altKey?: boolean;
}

/** A click's state: its count (`detail`: 2 for the second click of a double-click). */
export interface ClickState extends KeyState {
  readonly detail?: number;
}

/** How far a pointer moves before a press is a drag, not a click (CSS px). */
const DragSlop = 3;

/**
 * A drag-select held past a scrolling box's edge scrolls it toward the pointer at this many px a second for each px
 * past the edge, up to `EdgeScrollMax` (fixed: a reading speed near the edge, a page a second far past it; timed by the
 * clock, never by the frame count).
 */
const EdgeScrollRate = 8;
const EdgeScrollMax = 1200;

export interface WheelState extends KeyState {
  readonly deltaX: number;
  readonly deltaY: number;
  /** 0: pixels, 1: lines, 2: pages. */
  readonly deltaMode?: number;
}

const NoModifiers: ModifierState = {};

function platformIsMac(): boolean {
  return typeof navigator !== 'undefined' && /Mac|iPhone|iPad/.test(navigator.platform ?? '');
}

export class ViewerInput {
  /** Cmd is `mod` on macOS, Ctrl elsewhere. */
  mac = platformIsMac();
  private _space = false;
  private readonly _mods: ModifierState = {};
  private readonly _hit: { -readonly [K in keyof Hit]: Hit[K] } = {} as never;
  /** A point the placement found, reused. */
  private readonly _point = new PlacePoint();
  /** The pointer's last place while it rests over the canvas (EC-143's re-hit reads it). */
  private _x = Number.NaN;
  private _y = Number.NaN;
  private _drag:
    | {
        x: number;
        y: number;
        x0: number;
        y0: number;
        section: string;
        target: HitTarget;
        select: boolean;
        /** A window's drag (its part): the `window` behavior's, before the bindings. */
        window: boolean;
        moved: boolean;
      }
    | undefined;
  /** A drag-select's scroll past its box's edge: its frame callback and the time it last stepped. */
  private _edge: { frame: number; t: number } | undefined;
  /** A drag moved before the button went up: the click the platform sends after it is not one. */
  private _suppressClick = false;
  /** Space was held while a drag began: its release is no tap (no toggle). */
  private _spaceUsed = false;

  /** The configured keys of commands (`view.keys`), parsed once: a command listed there answers these alone. */
  private readonly _keys: ReadonlyMap<string, readonly KeyChord[]>;

  constructor(private readonly viewer: SeqViewer) {
    const keys = new Map<string, KeyChord[]>();
    for (const [id, chords] of Object.entries(viewer.config.view.keys))
      keys.set(
        id,
        chords.map(parseChord).filter((c): c is KeyChord => c !== undefined),
      );
    this._keys = keys;
  }

  private get hoverBehavior(): HoverBehavior | undefined {
    return this.viewer.behavior(HoverProvider);
  }

  private get navigateBehavior(): NavigateBehavior | undefined {
    return this.viewer.behavior(NavigateProvider);
  }

  private get selectBehavior(): SelectBehavior | undefined {
    return this.viewer.behavior(SelectProvider);
  }

  private get windowBehavior(): WindowBehavior | undefined {
    return this.viewer.behavior(WindowProvider);
  }

  /** Whether the pointer is over the canvas. */
  get resting(): boolean {
    return !Number.isNaN(this._x);
  }

  /** The section of the layout at a y (a pointer over its label, its tracks or the gaps within), by address. */
  sectionAt(y: number): string | undefined {
    const layout = this.viewer.frame?.layout;
    if (!layout) return undefined;
    for (const s of layout.sections) if (y >= s.y && y < s.y + s.height) return s.section.address;
    return undefined;
  }

  /**
   * The display unit under a page point on a section, from its placement in the last frame (_spec/README.md "Input":
   * the layout's inverse, never the input's arithmetic): the nearest when between lines or in a gutter.
   */
  offsetAt(section: string, x: number, y: number): number | undefined {
    const p = this.viewer.placementOf(section);
    if (!p || !p.strips.n) return undefined;
    const o = p.offsetAt(x, y);
    return o >= 0 ? o : undefined;
  }

  /**
   * Round a ring (a layout whose window wraps the origin): the display units the pointer swept from one point to the
   * next, read from the placement by angle, the shorter way round; undefined elsewhere.
   */
  private _arcUnits(section: string, x0: number, y0: number, x1: number, y1: number): number | undefined {
    if (!this.viewer.view.wraps(section)) return undefined;
    const a = this.offsetAt(section, x0, y0);
    const b = this.offsetAt(section, x1, y1);
    const v = this.viewer.view.get(section);
    const axis = this.viewer.view.resolved.value?.sections.find((s) => s.address === section)?.axis;
    const t = axis && axisTotal(axis);
    if (a === undefined || b === undefined || !v || !t) return undefined;
    const d = b - a;
    return d - Math.round(d / t) * t;
  }

  /** Whether a section scrolls inside its box and its content overflows it (a text taller than its box). */
  private _overflows(p: Placement | undefined): p is Placement {
    return !!p?.clip && p.contentHeight > p.box.height;
  }

  /**
   * The wheel over a section that scrolls inside its box (_spec/README.md "Input": a text): its content scrolls by the
   * wheel's vertical delta through its layout's arithmetic while it overflows the box; at its top or bottom, or with
   * nothing to scroll, the wheel is given back to the page (false: scroll chaining).
   */
  private _scrollBox(section: string, dy: number): boolean {
    const view = this.viewer.view;
    const before = view.get(section)?.start;
    if (dy === 0 || before === undefined || !view.panPx(section, 0, dy)) return false;
    if (view.get(section)?.start === before) return false;
    this.viewer.schedule();
    return true;
  }

  /**
   * A drag-select's pointer past its scrolling box's top or bottom (_spec/README.md "Input"): the box scrolls toward
   * it, faster the further past, timed by the clock; the selection extends to what comes under the held pointer.
   * Stops when the pointer comes back in, the content ends or the drag ends.
   */
  private _edgeScroll(): void {
    const d = this._drag;
    const p = d?.select ? this.viewer.placementOf(d.section) : undefined;
    const past = this._overflows(p) && d ? pastEdge(p, d.y) : 0;
    if (!past) return this._stopEdge();
    if (this._edge) return;
    const step = () => {
      const e = this._edge;
      const g = this._drag;
      if (!e || !g?.select || this.viewer.retired) return this._stopEdge();
      const pl = this.viewer.placementOf(g.section);
      const by = this._overflows(pl) ? pastEdge(pl, g.y) : 0;
      const now = DomHost.now();
      const dt = Math.max(0, now - e.t);
      e.t = now;
      if (!by) return this._stopEdge();
      const speed = Math.sign(by) * Math.min(EdgeScrollMax, EdgeScrollRate * Math.abs(by));
      const before = this.viewer.view.get(g.section)?.start;
      this.viewer.view.panPx(g.section, 0, (speed * dt) / 1000);
      // what the last frame placed under the held pointer joins the drag; the next step reads the frame this one asks
      this._dragSelect(g.section, g.x, g.y);
      if (this.viewer.view.get(g.section)?.start === before && dt > 0) return this._stopEdge();
      this.viewer.schedule();
      e.frame = requestFrameCallback(step);
    };
    this._edge = { frame: requestFrameCallback(step), t: DomHost.now() };
  }

  private _stopEdge(): void {
    if (this._edge) cancelFrameCallback(this._edge.frame);
    this._edge = undefined;
  }

  /** A drag-select to the pointer: the position under it in reading order (the placement's nearest unit). */
  private _dragSelect(section: string, x: number, y: number): void {
    this.selectBehavior?.dragTo(section, x - this.viewer.dataX, this._hitAt(x, y), this.offsetAt(section, x, y));
  }

  /** The hit under a point as the hover model takes it (the section is read from the layer), in the reused record. */
  private _hitAt(x: number, y: number): Hit | undefined {
    const h: RenderHit | undefined = this.viewer.hit(x, y);
    if (!h) return undefined;
    const out = this._hit as Record<string, unknown>;
    out.layer = h.layer;
    out.system = h.system;
    out.position = h.position;
    out.element = h.element;
    out.extent = h.extent;
    out.hull = h.hull;
    // a placement's own targets (_spec/README.md "Input"): a window with its part, a gutter's row
    out.target = h.target === 'window' || h.target === 'row' ? h.target : undefined;
    out.part = h.part;
    // a block's hit: its positions on the block's own axis (a section of blocks on several systems)
    out.block = h.block;
    return this._hit;
  }

  private _modifiers(e: KeyState): ModifierState {
    return modifiersOf(e, this.mac, this._space, this._mods);
  }

  /** The pointer moved to (x, y): a drag pans, and the hover follows through the bindings. */
  move(x: number, y: number, e: KeyState = {}): void {
    if (this.viewer.retired) return;
    this._x = x;
    this._y = y;
    const mods = this._modifiers(e);
    const d = this._drag;
    const hit = this._hitAt(x, y);
    if (d) {
      if (!d.moved && Math.max(Math.abs(x - d.x0), Math.abs(y - d.y0)) >= DragSlop) d.moved = true;
      if (d.window) {
        d.x = x;
        d.y = y;
        const at = this.offsetAt(d.section, x, y);
        if (d.moved && at !== undefined) this.windowBehavior?.dragTo(at);
      } else if (d.select) {
        d.x = x;
        d.y = y;
        // the range from where it began to the pointer, in reading order, read from the placement even past the data
        if (d.moved) {
          this.selectBehavior?.dragTo(d.section, x - this.viewer.dataX, hit, this.offsetAt(d.section, x, y));
          this._edgeScroll();
        }
      } else {
        const dx = x - d.x;
        const dy = y - d.y;
        const units = this._arcUnits(d.section, d.x, d.y, x, y);
        d.x = x;
        d.y = y;
        this.navigateBehavior?.drag(d.section, dx, mods, d.target, dy, units);
      }
    }
    this.hoverBehavior?.move(hit, mods);
  }

  /** The pointer left the canvas: the hover it owned goes, and no re-hit follows. */
  leave(): void {
    this._x = this._y = Number.NaN;
    this.hoverBehavior?.leave();
  }

  /**
   * A resting pointer hit again (EC-143): the view moved under it, so what it is over is not what it was. Once a frame
   * (the viewer's draw asks after a layout); true when the hover changed.
   */
  rehit(): boolean {
    if (this.viewer.retired || !this.resting || this._drag) return false;
    const before = this.viewer.hover.version.value;
    this.hoverBehavior?.move(this._hitAt(this._x, this._y), this._mods);
    return this.viewer.hover.version.value !== before;
  }

  /** The primary button went down at (x, y): a drag begins when the bindings give one here. */
  down(x: number, y: number, e: KeyState = {}): boolean {
    if (this.viewer.retired) return false;
    // a jump's motion stops where it is drawn: the input reads one picture
    this.viewer.stopMotion();
    // a new press: a click after it is its own
    this._suppressClick = false;
    const section = this.sectionAt(y);
    if (section === undefined) return false;
    const mods = this._modifiers(e);
    const hit = this._hitAt(x, y);
    const target = hitTarget(hit);
    // a window's part is dragged by the window behavior (_spec/README.md "Behaviors", Window), Space+drag still pans
    const win = this.windowBehavior;
    if (target === 'window' && hit?.part && win && !this._space) {
      const at = this.offsetAt(section, x, y);
      if (at === undefined || !win.begin(section, hit.part, at)) return false;
      this._drag = { x, y, x0: x, y0: y, section, target, select: false, window: true, moved: false };
      return true;
    }
    const action = this.viewer.bindings.resolve('drag', mods, target);
    if (action === undefined) return false;
    // a drag-select (of any mode) begins on what was hit; it selects once the pointer moves past the slop (else it is
    // a click)
    const selecting = isDragSelect(action);
    const select = selecting && !!this.selectBehavior?.beginDrag(section, hit, mods, target);
    if (selecting && !select) return false;
    if (this._space) this._spaceUsed = true;
    this._drag = { x, y, x0: x, y0: y, section, target, select, window: false, moved: false };
    return true;
  }

  /**
   * The button went up (or the gesture was cancelled). The platform calls it more than once for one press (`pointerup`,
   * then `lostpointercapture`): only the call that ends a drag touches what follows it.
   */
  up(): void {
    const d = this._drag;
    if (!d) return;
    this._stopEdge();
    if (d.window) this.windowBehavior?.end();
    if (d.select) this.selectBehavior?.endDrag();
    this._suppressClick = d.moved;
    this._drag = undefined;
  }

  /**
   * A click at (x, y) through the bindings (select, extend, toggle, clear): what was hit is the granularity. The click
   * that ends a drag, and the second click of a double-click, are not clicks. True when a binding took it.
   */
  click(x: number, y: number, e: ClickState = {}): boolean {
    if (this.viewer.retired) return false;
    // a jump's motion stops where it is drawn: the input reads one picture
    this.viewer.stopMotion();
    if (this._suppressClick) {
      this._suppressClick = false;
      return false;
    }
    if ((e.detail ?? 1) > 1) {
      this._restoreCleared();
      return false;
    }
    this._cleared = undefined;
    const select = this.selectBehavior;
    if (!select) return false;
    const section = this.sectionAt(y);
    // the track's "n filtered" chip: its click clears the layer's filter (one patch)
    const chip = this.viewer.hit(x, y);
    if (chip?.target === 'filter') {
      void this.viewer.execute('layer.filter', { address: chip.layer });
      return true;
    }
    const hit = this._hitAt(x, y);
    const sel = this.viewer.selection;
    const before = sel.value.value;
    const anchor = sel.anchor;
    // the behavior keeps nothing of the reused record but what it copies
    this._clicking = section;
    let taken: boolean;
    try {
      taken = select.click(
        section,
        hit ? { ...hit, ...(section !== undefined ? { section } : {}) } : undefined,
        this._modifiers(e),
      );
    } finally {
      this._clicking = undefined;
    }
    // a click that cleared the selection may be a double-click's first: kept until the next click says; so is a
    // Shift+click's extend (a Shift+double-click zooms out and selects nothing)
    const after = sel.value.value;
    if ((!isNothing(before) && isNothing(after)) || (e.shiftKey && after !== before))
      this._cleared = { value: before, anchor, after };
    return taken;
  }

  /**
   * What the first click of a double-click changed (clicking the only selected element cleared it; a Shift+click
   * extended it), until the next click: the selection before and after it.
   */
  private _cleared: { value: Emphasized; anchor: Anchor | undefined; after: Emphasized } | undefined;

  /** A double-click's second click (or its `dblclick`) gives back what its first click changed, while it stands. */
  private _restoreCleared(): void {
    const c = this._cleared;
    this._cleared = undefined;
    if (c && this.viewer.selection.value.value === c.after) this.viewer.selection.set(c.value, c.anchor);
  }

  get dragging(): boolean {
    return this._drag !== undefined;
  }

  /** The section a drag-select or a click is selecting in, while it does (a linked reveal leaves it be). */
  get selectingIn(): string | undefined {
    return this._drag?.select ? this._drag.section : this._clicking;
  }

  private _clicking: string | undefined;

  /**
   * A wheel at (x, y). A trackpad's pinch arrives as a wheel with Ctrl held and small pixel deltas: it is the `pinch`
   * gesture (`pointer.pinchFromCtrlWheel`), its zoom factor `exp(-deltaY · pinchSpeed)`; any other wheel is the `wheel`
   * gesture. True when a binding took it (the page should not scroll).
   */
  wheel(x: number, y: number, e: WheelState): boolean {
    if (this.viewer.retired) return false;
    // a jump's motion stops where it is drawn: the input reads one picture
    this.viewer.stopMotion();
    const section = this.sectionAt(y);
    const nav = this.navigateBehavior;
    if (section === undefined || !nav) return false;
    const { pointer } = this.viewer.config;
    // lines are `pointer.lineHeight` px; a page is the viewport's own extent on each axis
    const ux = e.deltaMode === 1 ? pointer.lineHeight : e.deltaMode === 2 ? this.viewer.width : 1;
    const uy = e.deltaMode === 1 ? pointer.lineHeight : e.deltaMode === 2 ? this.viewer.height : 1;
    const dx = e.deltaX * ux;
    const dy = e.deltaY * uy;
    const target = hitTarget(this._hitAt(x, y));
    // over the label column the zoom is about the section's first visible position, never before it
    const width = this.viewer.view.get(section)?.width ?? 0;
    const local = Math.min(Math.max(x - this.viewer.dataX, 0), width);
    // the unit under the pointer, from the placement (a text's in reading order): what a zoom keeps where it is
    const at = this.offsetAt(section, x, y);
    // on macOS Ctrl is no zoom modifier, so a Ctrl+wheel is a pinch whatever its delta; elsewhere a large delta is a mouse
    const pinch =
      e.ctrlKey && pointer.pinchFromCtrlWheel && (this.mac || (!e.deltaMode && Math.abs(dy) <= pointer.pinchMaxDelta));
    if (pinch) {
      // the pinch's own modifiers: none (the Ctrl is the gesture's, not a key held)
      const taken = this.viewer.bindings.resolve('pinch', NoModifiers, target) !== undefined;
      if (taken) nav.pinch(section, Math.exp(-dy * pointer.pinchSpeed), local, NoModifiers, target, at);
      this.viewer.schedule();
      return taken;
    }
    const mods = this._modifiers(e);
    // a section that scrolls inside its box (a text) takes the wheel only while it overflows, and gives it back at its
    // ends; Cmd/Ctrl+wheel zooms it as any section
    if (this.viewer.placementOf(section)?.clip && !e.ctrlKey && !e.metaKey) return this._scrollBox(section, dy);
    const taken = this.viewer.bindings.resolve('wheel', mods, target) !== undefined;
    if (taken) nav.wheel(section, dx, dy, local, mods, target, at);
    // the frame asked for after the behavior's own, so a wheel's deltas apply in the frame that draws them
    this.viewer.schedule();
    return taken;
  }

  /** A double-click: the bound command, with the element or position under the pointer as what to focus. */
  doubleClick(x: number, y: number, e: KeyState = {}): boolean {
    if (this.viewer.retired) return false;
    // a jump's motion stops where it is drawn: the input reads one picture
    this.viewer.stopMotion();
    this._restoreCleared();
    const section = this.sectionAt(y);
    if (section === undefined) return false;
    const mods = this._modifiers(e);
    const hit = this._hitAt(x, y);
    const action = this.viewer.bindings.resolve('double-click', mods, hitTarget(hit));
    if (action === undefined) return false;
    const focus = hit ? this._focusOf(section, hit) : undefined;
    // a step zoom is about the pointer (the unit the placement has under it), by the double-click's own factor
    const p = this.viewer.placementOf(section);
    const pt = this._point;
    const at = p?.at(x, y, pt) && (pt.kind === 'data' || pt.kind === 'gap') ? p.offsetAt(x, y) : -1;
    this.viewer.runCommand(action, {
      section,
      ...(focus ? { focus } : {}),
      ...(at >= 0 ? { at } : {}),
      step: this.viewer.config.navigate.doubleClickZoomStep,
    });
    return true;
  }

  /**
   * A mouse's side button (back or forward) at (x, y): the bound command (`view.back`, `view.forward` by default). True
   * when one ran, so the page does not go back in its own history.
   */
  sideButton(which: 'back' | 'forward', x: number, y: number, e: KeyState = {}): boolean {
    if (this.viewer.retired) return false;
    this._sideTaken = undefined;
    const section = this.sectionAt(y) ?? this.viewer.view.resolved.value?.sections[0]?.address;
    const action = this.viewer.bindings.resolve(`${which}-button`, this._modifiers(e), hitTarget(this._hitAt(x, y)));
    if (action === undefined) return false;
    // only a command that can act takes the button: with no step to go, the page keeps its own Back
    const ctx = { section };
    const can = this.viewer.commands.find((c) => c.id === action)?.when?.(ctx as never);
    if (can !== undefined && can !== true) return false;
    this._sideTaken = which;
    this.viewer.runCommand(action, ctx);
    return true;
  }

  /** Whether the side button's press was taken (its release is then not the page's), once. */
  sideButtonTaken(which: 'back' | 'forward'): boolean {
    const taken = this._sideTaken === which;
    this._sideTaken = undefined;
    return taken;
  }

  private _sideTaken: 'back' | 'forward' | undefined;

  /**
   * What to focus for a hit: the element's hull on the section (its extent carried along the layer's path: a domain on
   * the gene it is coded by spans its exons), else its extent where it is on the section's system, else the position's
   * cell, in the section's offsets (the cells' offsets, so a reversed axis keeps the whole hull).
   */
  private _focusOf(section: string, hit: Hit): { start: number; end: number } | undefined {
    const s = this.viewer.view.resolved.value?.sections.find((x) => x.address === section);
    if (!s) return undefined;
    // a section of blocks on several systems (its axis its own units): the block's axis, from the block's offset
    const block = hit.block !== undefined && s.axis.system !== hit.system ? s.blocks?.[hit.block] : undefined;
    const at = block ? (block.offset ?? 0) : 0;
    const r = this._focusOn(block ? block.axis : s.axis, hit, !block && this.viewer.view.wraps(section));
    return r && at ? { start: r.start + at, end: r.end + at } : r;
  }

  private _focusOn(axis: Axis, hit: Hit, arc: boolean): { start: number; end: number } | undefined {
    // an item's whole extent (a layer's hit is the piece hit: one row of an item over a cyclic origin)
    const extent = (hit.element && this.viewer.selection.extentOfHit(hit)) || hit.extent;
    const part = extent?.parts.find((p) => p.system === hit.system);
    // round a ring: the arc its pieces cover, the circle less the widest gap between them (over the origin)
    if (arc && part && part.ranges.length >= 4) {
      const r = coveringArc(axis, part.ranges);
      if (r) return r;
    }
    // the edges' offsets, in display order whatever the axis's direction
    if (hit.hull && hit.hull[1] > hit.hull[0]) {
      const a = nearestOffsetOf(axis, hit.hull[0]);
      const b = nearestOffsetOf(axis, hit.hull[1]);
      if (a !== undefined && b !== undefined && a !== b) return { start: Math.min(a, b), end: Math.max(a, b) };
    }
    if (part && part.ranges.length >= 2) {
      const a = offsetOf(axis, part.ranges[0]);
      const b = offsetOf(axis, part.ranges[part.ranges.length - 1]);
      return a !== undefined && b !== undefined ? { start: Math.min(a, b), end: Math.max(a, b) } : undefined;
    }
    if (hit.position === undefined) return undefined;
    const start = cellOffset(axis, Math.floor(hit.position));
    return start === undefined ? undefined : { start, end: start + 1 };
  }

  /** Whether a binding takes Space as a modifier: the page's own Space (a scroll) is then the viewer's while it has focus. */
  get spaceIsModifier(): boolean {
    return this.viewer.bindings.bindings.some((b) => Array.isArray(b.modifiers) && b.modifiers.includes('space'));
  }

  /** The Space key, held for a drag that pans on any target. */
  setSpace(held: boolean): void {
    if (held && !this._space) this._spaceUsed = false;
    this._space = held;
  }

  /**
   * Space released: a tap (no drag began while it was held) is the Space key's command (the cursor's toggle); a hold
   * that panned is not. True when a command ran.
   */
  spaceUp(e: KeyState = {}): boolean {
    const tap = this._space && !this._spaceUsed;
    this.setSpace(false);
    return tap
      ? this.key({ key: ' ', shiftKey: !!e.shiftKey, altKey: !!e.altKey, metaKey: !!e.metaKey, ctrlKey: !!e.ctrlKey })
      : false;
  }

  /**
   * A key: the first command whose shortcut it is and that can run now (Esc cancels a gesture's preview first, then
   * clears the selection), on the section under the pointer (else the first). True when one ran. Keys go to the
   * commands' shortcuts; the keymap of `@molstar/ui` replaces this when the commands join a manager.
   */
  key(e: KeyState & { key: string }): boolean {
    if (this.viewer.retired) return false;
    // Esc during a drag-select cancels it: the selection before it comes back, and the click after it is none
    // Esc during a window's drag puts the view back where the press found it
    if (e.key === 'Escape' && this._drag?.window) {
      this.windowBehavior?.cancel();
      this._drag = undefined;
      this._suppressClick = true;
      return true;
    }
    if (e.key === 'Escape' && this._drag?.select) {
      this._stopEdge();
      this.selectBehavior?.cancelDrag();
      this._drag = undefined;
      this._suppressClick = true;
      return true;
    }
    const section =
      (this.resting ? this.sectionAt(this._y) : undefined) ?? this.viewer.view.resolved.value?.sections[0]?.address;
    const ctx = { section };
    let matched = false;
    // a gesture in progress first (a preview's command), then the rest in the spec's order
    for (const previews of [true, false])
      for (const c of this.viewer.commands) {
        if ((changeKindOf(c.id) === 'preview') !== previews) continue;
        // a command's configured keys (`view.keys`), else its shortcut
        const keys = this._keys.get(c.id);
        if (
          keys ? !keys.some((k) => matchesChord(k, e, this.mac)) : !c.shortcut || !shortcutIs(c.shortcut, e, this.mac)
        )
          continue;
        const can = c.when?.(ctx as never);
        // the history's keys are the page's Back and Forward when the view has no step to go
        if (can !== undefined && can !== true) {
          if (!YieldWhenIdle.has(c.id)) matched = true;
          continue;
        }
        matched = true;
        this.viewer.runCommand(c.id, ctx);
        return true;
      }
    // a shortcut with modifiers is the viewer's even when nothing can run now, but the history's (Alt+Left is the
    // browser's Back when the view has no step back); a plain key nothing takes (Esc with nothing to cancel) stays the
    // host's
    return matched && !!(e.altKey || e.shiftKey || e.metaKey || e.ctrlKey);
  }
}

/** Commands whose keys and buttons are the page's when they cannot run: the view's history (the browser's Back). */
const YieldWhenIdle: ReadonlySet<string> = new Set(['view.back', 'view.forward']);

/** Whether a key event is a command's shortcut exactly (its key and every modifier). */
function shortcutIs(k: KeyChord, e: KeyState & { readonly key: string }, mac: boolean): boolean {
  return (
    k.key === e.key &&
    !!k.alt === !!e.altKey &&
    !!k.shift === !!e.shiftKey &&
    !!k.mod === (mac ? !!e.metaKey : !!e.ctrlKey)
  );
}

/** How far a point is past a box's top (negative) or bottom (positive), px; 0 within it. */
function pastEdge(p: Placement, y: number): number {
  const b = p.box;
  return y < b.y ? y - b.y : y > b.y + b.height ? y - (b.y + b.height) : 0;
}

/**
 * The arc of a cyclic axis that pieces cover (positions `[start, end, …]`): from the end of the widest gap between them
 * round to its start (past the total over the origin); undefined when a piece is not shown.
 */
function coveringArc(axis: Axis, ranges: ArrayLike<number>): { start: number; end: number } | undefined {
  const t = axisTotal(axis);
  if (!t) return undefined;
  const spans: [number, number][] = [];
  for (let i = 0; i + 1 < ranges.length; i += 2) {
    const a = offsetOf(axis, ranges[i]);
    const b = offsetOf(axis, ranges[i + 1]);
    if (a === undefined || b === undefined) return undefined;
    spans.push([Math.min(a, b), Math.max(a, b)]);
  }
  spans.sort((x, y) => x[0] - y[0]);
  let gap = -1;
  let at = 0;
  for (let k = 0; k < spans.length; k++) {
    const next = k + 1 < spans.length ? spans[k + 1][0] : spans[0][0] + t;
    const g = next - spans[k][1];
    if (g > gap) [gap, at] = [g, k];
  }
  const start = at + 1 < spans.length ? spans[at + 1][0] : spans[0][0];
  const end = spans[at][1];
  return { start, end: end > start ? end : end + t };
}
