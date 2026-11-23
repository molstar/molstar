/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node as SpecNode } from '@molstar/core/repr/tree';
import {
  BlockPlacement,
  floorExtent,
  type GridPlacement,
  type Placement,
  PlacePoint,
  type Rail,
  panOnly,
  RowsRail,
  SectionRail,
  railUnitAt,
  type Hit as RenderHit,
  type SectionLayout,
  type WindowPart,
} from '@molstar/seqx-core/render';
import type { HeaderActionEntry } from './header-entries.js';
import { type PartCursor, partCursor } from './part-cursor.js';
import { cancelFrameCallback, requestFrameCallback } from '@molstar/ui/dom/frame';
import { DomHost } from '@molstar/ui/dom/host';
import { type Location, rangesLocation } from '@molstar/seqx-core/query';
import { type Axis, blockSpan } from '@molstar/seqx-core/view';
import {
  CladeActions,
  type KeyChord,
  matchesChord,
  type ModifierState,
  modifiersOf,
  parseChord,
} from '@molstar/seqx-viewer/bindings';
import { changeKindOf } from './commands.js';
import type {
  HoverBehavior,
  HoverCardBehavior,
  KeysPanelBehavior,
  NavigateBehavior,
  NavigatorBehavior,
  RowsBehavior,
  RowsRailGesture,
  RowsRailWindow,
  SelectBehavior,
  WindowBehavior,
} from '@molstar/seqx-viewer/behaviors';
import { isDragSelect } from '@molstar/seqx-viewer/behaviors/select';
import { trackOfLayer } from '@molstar/seqx-viewer/behaviors/two-axes';
import { SectionRailGesture, SectionWheel, scrolledSection } from './section-scroll.js';
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
import { BridgeView } from '@molstar/seqx-viewer/models/bridge-hooks';

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
  /** The pointer's kind (a pointer event's): `touch` widens a window's hit. */
  readonly pointerType?: string;
}

/** A click's state: its count (`detail`: 2 for the second click of a double-click). */
export interface ClickState extends KeyState {
  readonly detail?: number;
}

/** How far a pointer moves before a press is a drag, not a click (CSS px). */
export const DragSlop = 3;

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
  /** A page point for a placement's `at`, which reads it from an array. */
  private readonly _at = new Float64Array(2);
  private _pageAt(x: number, y: number): Float64Array {
    this._at[0] = x;
    this._at[1] = y;
    return this._at;
  }
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
        /** A ruler rail's press (its part, outside it a click or a sweep): the `navigator` behavior's. */
        rail?: Rail;
        part?: WindowPart;
        /** A rows rail's press (upright: along the page's y): the rows rail's gesture. */
        rowsRail?: Rail;
      }
    | undefined;
  /** A rows rail's press that never moved: its click centres the rows there (outside the window). */
  private _rowsRailClick: { track: string; at: number; part: WindowPart; window: RowsRailWindow } | undefined;
  /** The rows rails' gestures (a rows band's scroll and zoom by its rail). */
  private _rowsRailGesture: RowsRailGesture | undefined;
  /** The rows rail's gesture: the `rows` behavior's (none without it: a rows rail then does not move). */
  private get rowsRailGesture(): RowsRailGesture | undefined {
    this._rowsRailGesture ??= this.viewer.behavior<RowsBehavior>('rows')?.railGesture();
    return this._rowsRailGesture;
  }
  /** A section rail's gesture: the viewer's own (a section past its cap scrolls in every composition). */
  private _sectionRailGesture: SectionRailGesture | undefined;
  /** The wheel over a section past its cap, coalesced to the frame. */
  private _sectionWheel: SectionWheel | undefined;
  // The two getters below make their objects in methods of their own: a closure in a getter makes V8 allocate the
  // getter's context (40 bytes) on every call, the fast path's too, once the slow path has run in optimized code.
  private get sectionWheel(): SectionWheel {
    return (this._sectionWheel ??= this._newSectionWheel());
  }
  private _newSectionWheel(): SectionWheel {
    return new SectionWheel(this.viewer.view, () => this.viewer.frame?.layout);
  }
  private get sectionRailGesture(): SectionRailGesture {
    return (this._sectionRailGesture ??= this._newSectionRailGesture());
  }
  private _newSectionRailGesture(): SectionRailGesture {
    return new SectionRailGesture(this.viewer.view, () => this.viewer.frame?.layout);
  }
  /** A rail's press that never moved: its click (the platform's, after the button went up) centres or not. */
  private _railClick: { section: string; at: number; part: WindowPart } | undefined;
  /** The pointer is a finger (the last pointer event's kind): a window's hit is wider. */
  private _touch = false;
  /** The window part under the pointer (a rail's, a ring's): it lights; copied, kept while it holds. */
  private _windowHit: RenderHit | undefined;
  /** The window being dragged: its rail is drawn grown wherever the pointer goes. */
  private _dragHit: RenderHit | undefined;
  /** What the last hit found, as the render gave it (a window's part lights from it). */
  private _renderHit: RenderHit | undefined;
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
    viewer.selection.value.subscribe(() => this._selectionChanges++);
  }

  /** How many times the selection changed: what a click cleared is given back only if nothing changed it since. */
  private _selectionChanges = 0;
  /** The pointer as the hover card is told it (reused). */
  private readonly _pointerAt = new Float64Array(2);

  private get hoverBehavior(): HoverBehavior | undefined {
    return this.viewer.behavior<HoverBehavior>('hover');
  }

  private get navigateBehavior(): NavigateBehavior | undefined {
    return this.viewer.behavior<NavigateBehavior>('navigate');
  }

  private get selectBehavior(): SelectBehavior | undefined {
    return this.viewer.behavior<SelectBehavior>('select');
  }

  private get windowBehavior(): WindowBehavior | undefined {
    return this.viewer.behavior<WindowBehavior>('window');
  }

  private get navigatorBehavior(): NavigatorBehavior | undefined {
    return this.viewer.behavior<NavigatorBehavior>('navigator');
  }

  /**
   * The hover card and the keys panel, by name and only when listed (no issue when not: they are no gesture's), so a
   * viewer without them bundles nothing of theirs.
   */
  private get cardBehavior(): HoverCardBehavior | undefined {
    return this.viewer.behaviors.get('hoverCard');
  }

  private get keysPanel(): KeysPanelBehavior | undefined {
    return this.viewer.behaviors.get('keysPanel');
  }

  /** The window part under the pointer, lit (a rail's, a ring's); undefined elsewhere. */
  get windowHit(): RenderHit | undefined {
    return this._windowHit;
  }

  /** The window being dragged (a rail's grows while its drag runs). */
  get dragHit(): RenderHit | undefined {
    return this._dragHit;
  }

  /**
   * The cursor over the canvas by the part (`partCursor`): a rail's grip resizes along the rail, its body grabs
   * (`grabbing` while dragged), outside it a click centres; a ring's window grabs.
   */
  get cursor(): PartCursor | undefined {
    const h = this._dragHit ?? this._windowHit;
    if (!h?.part) return undefined;
    if (h.navigator !== 'rail' && !upright(h.navigator)) return h === this._dragHit ? 'grabbing' : 'grab';
    return partCursor(h.part, h.navigator === 'rail' ? 'positions' : 'rows', h === this._dragHit);
  }

  /** The window under the pointer as lit: copied when it changes (its section, part, navigator or layer), a frame asked. */
  private _windowHover(hit: RenderHit | undefined): void {
    const h = hit?.target === 'window' && hit.part ? hit : undefined;
    const w = this._windowHit;
    if (h?.section === w?.section && h?.part === w?.part && h?.navigator === w?.navigator && h?.layer === w?.layer)
      return;
    this._windowHit = h && { ...h };
    this.viewer.schedule();
  }

  /** A section's rail in the last frame. */
  /** A rows rail of the last frame, by its track's address. */
  private _rowsRailOf(track: string): Rail | undefined {
    for (const s of this.viewer.frame?.layout.sections ?? []) {
      for (const r of s.rowsRails ?? []) if (r.key === track) return r;
      for (const r of s.lanesRails ?? []) if (r.key === track) return r;
      if (s.sectionRail?.key === track) return s.sectionRail;
    }
    return undefined;
  }

  private _railOf(section: string): Rail | undefined {
    return this._sectionOf(section)?.rail;
  }

  /** The type name of a section's layer by its key (a scoped binding's), in the last frame. */
  private _layerTypeOf(section: string, key: string): string | undefined {
    const s = this._sectionOf(section);
    if (!s) return undefined;
    for (const l of s.sectionLays) if (l.key === key) return l.type.name;
    for (const t of s.tracks) for (const l of t.layers) if (l.key === key) return l.type.name;
    for (const ls of s.lays) for (const l of ls) if (l.key === key) return l.type.name;
    return undefined;
  }

  /** A section's layout in the last frame, by address. */
  private _sectionOf(section: string): SectionLayout | undefined {
    const sections = this.viewer.frame?.layout.sections;
    if (!sections) return undefined;
    for (let i = 0; i < sections.length; i++) if (sections[i].section.address === section) return sections[i];
    return undefined;
  }

  /** Whether the pointer is over the canvas. */
  get resting(): boolean {
    return !Number.isNaN(this._x);
  }

  /** The section of the layout at a y (a pointer over its label, its tracks or the gaps within), by address. */
  sectionAt(y: number): string | undefined {
    const layout = this.viewer.frame?.layout;
    if (!layout) return undefined;
    // a pinned rail is the frame's foot's, past the sections (or over the foot of a section the frame clips)
    for (const s of layout.sections) {
      const rs = s.rails;
      for (let i = 0; rs && i < rs.length; i++)
        if (rs[i].pinned && y >= rs[i].top && y < rs[i].bottom) return s.section.address;
    }
    for (const s of layout.sections) if (y >= s.y && y < s.y + s.height) return s.section.address;
    return undefined;
  }

  /** The bridge whose band holds a y (the bridge chunk's): a drag or a wheel there moves both its sides. */
  private _band(y: number): string | undefined {
    const layout = this.viewer.frame?.layout;
    return layout && BridgeView.parts?.bandAt(layout, y);
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

  /** The display unit down a grid's y under a page point (its plot's second axis), undefined elsewhere. */
  offsetYAt(section: string, x: number, y: number): number | undefined {
    const p = this.viewer.placementOf(section);
    const pt = this._point;
    if (!p || !p.strips.n || !p.at(this._pageAt(x, y), pt) || pt.strip < 0) return undefined;
    const u = p.strips.ys[pt.strip]?.unitAt(pt.y);
    return u !== undefined && u >= 0 ? u : undefined;
  }

  /**
   * Round a ring (a layout whose window wraps the origin): the display units the pointer swept from one point to the
   * next, read from the placement by angle, the shorter way round; undefined elsewhere.
   */
  private _arcUnits(section: string, x0: number, y0: number, x1: number, y1: number): number | undefined {
    // the window follows the pointer round a ring (one drawing the whole circle); a linear view pans by px
    // a polar placement (`warpOf`): the rings layout read through the placement contract, never by its class, so a
    // composition that lists no rings bundles none
    if (!this.viewer.view.arcs(section) || !this.viewer.placementOf(section)?.warpOf) return undefined;
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
   * The wheel's vertical delta over a tracks section past its cap (_spec/README.md "Input", The one rail): its tracks
   * below the ruler scroll by it (the viewport's `top`, clamped to its content); at its top or bottom the wheel is the
   * page's (false: scroll chaining).
   */
  private _scrollSection(section: string, dy: number): boolean {
    if (!this.sectionWheel.wheel(section, dy)) return false;
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
    const h: RenderHit | undefined = this.viewer.hit(x, y, this._touch);
    this._renderHit = h;
    if (!h) return undefined;
    const out = this._hit as Record<string, unknown>;
    out.layer = h.layer;
    // a bridge's ribbon is in no section (its `section` the bridge's id): named, or its card finds no bridge
    out.section = h.bridge === undefined ? undefined : h.section;
    out.system = h.system;
    out.position = h.position;
    out.element = h.element;
    out.extent = h.extent;
    out.hull = h.hull;
    out.ySystem = h.ySystem;
    out.yPosition = h.yPosition;
    // the hit's own targets (_spec/README.md "Input"): a window with its part, a gutter's or a header's row, a cell of
    // two axes, a tree node, a summary row, a header's slot, state or line; the rest read from what was hit
    out.target = h.target === 'break' ? undefined : h.target;
    out.part = h.part;
    out.navigator = h.navigator;
    // a block's hit: its positions on the block's own axis (a section of blocks on several systems)
    out.block = h.block;
    // a grid over blocks' cell: its y block too (the rect drag reads its units on the hit strip's blocks)
    out.yBlock = h.yBlock;
    // what a layer names within the position, for the targets that read it (a rows column cell: its value, which the
    // tip reads whole; a title's field; the hidden-columns chip): another's (a logo's letter) would make each move within
    // its position a new hover (the card re-armed, the announcer speaking again)
    const t = h.target;
    out.state =
      t === 'row' ||
      t === 'column-title' ||
      t === 'columns-hidden' ||
      t === 'group' ||
      t === 'run' ||
      t === 'header-slot' ||
      t === 'header-state' ||
      (t === 'window' && h.count !== undefined)
        ? h.state
        : undefined;
    out.members = h.members;
    // a run's count (its blocks: "one of 2,841"; the ruler rail's place over a run too)
    out.count = t === 'run' || t === 'window' ? h.count : undefined;
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
    this._touch = e.pointerType === 'touch';
    const mods = this._modifiers(e);
    const d = this._drag;
    const hit = this._hitAt(x, y);
    this._windowHover(this._renderHit);
    if (d) {
      if (!d.moved && Math.max(Math.abs(x - d.x0), Math.abs(y - d.y0)) >= DragSlop) d.moved = true;
      if (d.rowsRail) {
        d.x = x;
        d.y = y;
        const at = railUnitAt(d.rowsRail, y);
        if (d.moved && d.rowsRail.hit.navigator === SectionRail) this.sectionRailGesture.dragTo(at);
        else if (d.moved) this.rowsRailGesture?.dragTo(at);
      } else if (d.rail) {
        d.x = x;
        d.y = y;
        if (d.moved) this.navigatorBehavior?.dragTo(railUnitAt(d.rail, x));
      } else if (d.window) {
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
    this.hoverBehavior?.move(this._hoverable(hit), mods);
    // the card is told where the pointer is (two numbers): it shows on rest, hides while a drag runs
    // in an array: two fractional numbers passed to a call that is not inlined are boxed (a move allocates nothing)
    const card = this.cardBehavior;
    if (card) {
      this._pointerAt[0] = x;
      this._pointerAt[1] = y;
      card.pointer(this._pointerAt, d !== undefined && d.moved);
    }
  }

  /** What the hover model (and so the card and the data hover) takes of a hit: a rail's is the rail's own, none. */
  private _hoverable(hit: Hit | undefined): Hit | undefined {
    if (hit?.target !== 'window') return hit;
    const n = hit.navigator;
    return n === 'rail' || upright(n) ? undefined : hit;
  }

  /** The pointer left the canvas: the hover it owned goes, and no re-hit follows. */
  leave(): void {
    this._x = this._y = Number.NaN;
    // a rail never stays lit (grown) without a pointer over it or its drag running
    this._windowHover(undefined);
    this.hoverBehavior?.leave();
    this.cardBehavior?.leave();
  }

  /**
   * A resting pointer hit again (EC-143): the view moved under it, so what it is over is not what it was. Once a frame
   * (the viewer's draw asks after a layout); true when the hover changed.
   */
  rehit(): boolean {
    if (this.viewer.retired || !this.resting || this._drag) return false;
    const before = this.viewer.hover.version.value;
    this.hoverBehavior?.move(this._hoverable(this._hitAt(this._x, this._y)), this._mods);
    return this.viewer.hover.version.value !== before;
  }

  /** The primary button went down at (x, y): a drag begins when the bindings give one here. */
  down(x: number, y: number, e: KeyState = {}): boolean {
    if (this.viewer.retired) return false;
    // a jump's motion stops where it is drawn: the input reads one picture
    this.viewer.stopMotion();
    // a press may begin a drag: the card goes until the hover rests again
    this.cardBehavior?.interrupt();
    // a new press: a click after it is its own
    this._suppressClick = false;
    this._railClick = undefined;
    this._rowsRailClick = undefined;
    this._touch = e.pointerType === 'touch';
    const inSection = this.sectionAt(y);
    const section = inSection ?? this._band(y);
    if (section === undefined) return false;
    const mods = this._modifiers(e);
    // a section's wheel still pending lands first: the press reads the model
    this._sectionWheel?.apply();
    const hit = this._hitAt(x, y);
    // over a bridge's band a drag moves both its sides and selects nothing (_spec/README.md "Bridges")
    const target = inSection === undefined ? 'empty' : hitTarget(hit);
    // a rows rail's part (_spec/README.md "Input", Navigators, The rows rail): its window scrolls, a grip zooms the rows
    const rowsRail =
      target === 'window' && upright(hit?.navigator) && hit?.layer && !this._space
        ? this._rowsRailOf(hit.layer)
        : undefined;
    if (rowsRail && hit?.part) {
      const at = railUnitAt(rowsRail, y);
      const pressed =
        rowsRail.hit.navigator === SectionRail
          ? this.sectionRailGesture.press(rowsRail.key, hit.part, at)
          : this.rowsRailGesture?.press(rowsRail.key, hit.part, at, rowsRail);
      if (!pressed) return false;
      this._drag = {
        x,
        y,
        x0: x,
        y0: y,
        section,
        target,
        select: false,
        window: true,
        moved: false,
        rowsRail,
        part: hit.part,
      };
      this._dragHit = this._renderHit && { ...this._renderHit };
      this.viewer.schedule();
      return true;
    }
    // a ruler rail's part (_spec/README.md "Input", Navigators): its window dragged, outside it a click or a sweep, by
    // the navigator behavior at the rail's scale
    const rail = target === 'window' && hit?.navigator === 'rail' && !this._space ? this._railOf(section) : undefined;
    if (rail && hit?.part) {
      const nav = this.navigatorBehavior;
      const at = railUnitAt(rail, x);
      if (!nav || !nav.press(section, hit.part, at, railDrawn(rail, this.viewer))) return false;
      this._drag = {
        x,
        y,
        x0: x,
        y0: y,
        section,
        target,
        select: false,
        window: true,
        moved: false,
        rail,
        part: hit.part,
      };
      this._dragHit = this._renderHit && { ...this._renderHit };
      this.viewer.schedule();
      return true;
    }
    // a window's part is dragged by the window behavior (_spec/README.md "Behaviors", Window), Space+drag still pans
    const win = this.windowBehavior;
    if (target === 'window' && hit?.part && win && !this._space) {
      const at = this.offsetAt(section, x, y);
      if (at === undefined || !win.begin(section, hit.part, at)) return false;
      this._drag = { x, y, x0: x, y0: y, section, target, select: false, window: true, moved: false };
      this._dragHit = this._renderHit && { ...this._renderHit };
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
    this.cardBehavior?.released();
    this._stopEdge();
    if (this._dragHit) {
      this._dragHit = undefined;
      this.viewer.schedule();
    }
    if (d.rowsRail) {
      const r = d.rowsRail;
      const g = r.hit.navigator === SectionRail ? this.sectionRailGesture : this.rowsRailGesture;
      // a press that never moved waits for its click: outside the window it centres the rows there
      if (!d.moved) {
        g?.cancel();
        this._rowsRailClick = {
          track: r.key,
          at: railUnitAt(r, d.y0),
          part: d.part ?? 'body',
          window: { start: r.start, end: r.end, total: r.total, width: r.width, hit: r.hit },
        };
      } else g?.release();
      this._suppressClick = true;
      this._drag = undefined;
      return;
    }
    if (d.rail) {
      const nav = this.navigatorBehavior;
      // a press that never moved waits for its click (the platform's, with its count): outside the window it centres
      if (!d.moved) {
        nav?.cancel();
        this._railClick = { section: d.section, at: railUnitAt(d.rail, d.x0), part: d.part ?? 'body' };
      } else nav?.release();
      // the click the platform sends after it is the rail's
      this._suppressClick = true;
      this._drag = undefined;
      return;
    }
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
    // what a click cleared is kept for that click's double-click only: any click but a double-click's second drops it
    if ((e.detail ?? 1) <= 1) this._cleared = undefined;
    // a rail's click: outside its window it centres there (a jump, its own step); its second a double-click's
    const rr = this._rowsRailClick;
    if (rr) {
      this._rowsRailClick = undefined;
      this._suppressClick = false;
      if (rr.window.hit?.navigator === SectionRail)
        this.sectionRailGesture.click(rr.track, rr.part, rr.at, e.detail ?? 1);
      else this.rowsRailGesture?.click(rr.track, rr.part, rr.at, rr.window, e.detail ?? 1);
      return true;
    }
    const rc = this._railClick;
    if (rc) {
      this._railClick = undefined;
      this._suppressClick = false;
      this.navigatorBehavior?.click(rc.section, rc.part === 'outside' ? rc.at : undefined, e.detail ?? 1);
      return true;
    }
    // a rows layer's column title: the rows ordered by its field (one patch), again the other direction; before the
    // double-click's drop, so a fast second click flips too (the order read from the spec as it is now)
    const title = this.viewer.hit(x, y);
    if (title?.target === 'column-title' && title.state !== undefined) {
      let n: SpecNode | undefined = this.viewer.spec.root.value;
      for (const id of title.layer.split('/').filter(Boolean)) n = n?.children?.find((c) => c.id === id);
      const now = n?.params?.order as { field?: unknown; direction?: unknown } | undefined;
      const here = typeof now === 'object' && now?.field === title.state;
      // ascending, descending, then none: the third click clears this column's sort
      if (here && now!.direction === 'descending')
        void this.viewer.execute('rows.unsort', { track: trackOfLayer(title.layer) });
      else
        void this.viewer.execute('rows.order', {
          track: trackOfLayer(title.layer),
          order: { field: title.state, direction: here ? 'descending' : 'ascending' },
        });
      this._suppressClick = false;
      return true;
    }
    if (this._suppressClick) {
      this._suppressClick = false;
      return false;
    }
    if ((e.detail ?? 1) > 1) {
      this._restoreCleared();
      return false;
    }
    const chip = this.viewer.hit(x, y);
    // the hidden-columns chip: its hover names them, a click does nothing
    if (chip?.target === 'columns-hidden') return true;
    // a narrow frame's tree toggle: the tree opened (session state, `rows.tree`)
    if (chip?.target === 'tree-toggle') {
      void this.viewer.execute('rows.tree', { track: trackOfLayer(chip.layer) });
      return true;
    }
    // a header's parts are never a selection's (_spec/README.md "Input"): the `headers` behavior's (by name: a viewer
    // without it bundles none of it); without it a state listed as an entry with the layers it is about runs that
    // entry's command over them (the "n filtered" chip clears its layers' filters), the other states are facts
    const t = chip?.target;
    if (chip && (t === 'header-state' || t === 'header-slot' || t === 'header-label' || t === 'section-header')) {
      const headers = this.viewer.behaviors.get('headers');
      if (headers) headers.click(x, y, this._touch);
      else if (t === 'header-state') {
        const command = this.viewer.registry.find(
          (r): r is HeaderActionEntry =>
            r.kind === 'header-action' && (r as HeaderActionEntry).entry.name === chip.state,
        )?.entry.command;
        const st = this.viewer.frame?.layout.sections
          .flatMap((s) => s.tracks)
          .find((x) => x.track.address === chip.layer)
          ?.header?.states.find((s) => s.name === chip.state);
        if (command && st?.addresses) void this.viewer.execute(command, { addresses: st.addresses });
      }
      return true;
    }
    // the rest is a selection's (none without the select behavior)
    const select = this.selectBehavior;
    if (!select) return false;
    const section = this.sectionAt(y);
    // a tree node's click: its bound command (toggles the clade's fold by default), unless it is a selection's
    if (chip?.target === 'tree-node' && chip.element) {
      const action = this.viewer.bindings.resolve('click', this._modifiers(e), 'tree-node');
      if (action !== undefined && !action.startsWith('selection.')) {
        void this.viewer.execute(action, { track: trackOfLayer(chip.layer), clade: chip.element.key });
        return true;
      }
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
      this._cleared = { value: before, anchor, changes: this._selectionChanges };
    return taken;
  }

  /**
   * What the first click of a double-click changed (clicking the only selected element cleared it; a Shift+click
   * extended it), until the next click: the selection before and after it.
   */
  private _cleared: { value: Emphasized; anchor: Anchor | undefined; changes: number } | undefined;

  /** A double-click's second click (or its `dblclick`) gives back what its first click changed, while it stands. */
  private _restoreCleared(): void {
    const c = this._cleared;
    this._cleared = undefined;
    if (c && this._selectionChanges === c.changes) this.viewer.selection.set(c.value, c.anchor);
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
   * gesture. True when a binding took it, and always over a plot's box (the wheel over a plot is the plot's: at an
   * axis's end it stops there and the page does not scroll); elsewhere the page should not scroll when true.
   */
  wheel(x: number, y: number, e: WheelState): boolean {
    if (this.viewer.retired) return false;
    // a jump's motion stops where it is drawn: the input reads one picture
    this.viewer.stopMotion();
    this.cardBehavior?.interrupt();
    const section = this.sectionAt(y) ?? this._band(y);
    const nav = this.navigateBehavior;
    if (section === undefined || !nav) return false;
    const { pointer } = this.viewer.config;
    // lines are `pointer.lineHeight` px; a page is the viewport's own extent on each axis
    const ux = e.deltaMode === 1 ? pointer.lineHeight : e.deltaMode === 2 ? this.viewer.width : 1;
    const uy = e.deltaMode === 1 ? pointer.lineHeight : e.deltaMode === 2 ? this.viewer.height : 1;
    const dx = e.deltaX * ux;
    const dy = e.deltaY * uy;
    // a plain vertical wheel over a section past its cap and over none of its scrolling bands: the section's at once
    // (the chain's second step, no hit read: a wheel step allocates nothing)
    const plain = !e.ctrlKey && !e.metaKey && !e.shiftKey && !e.altKey && Math.abs(dy) > Math.abs(dx);
    const capped = plain ? scrolledSection(this.viewer.frame?.layout, section) : undefined;
    if (capped && capped.layout.name !== 'text' && !overBand(capped, x, y)) {
      if (this._scrollSection(section, dy)) return true;
      return this._overPlot(section, x, y);
    }
    const target = hitTarget(this._hitAt(x, y));
    // over the label column the zoom is about the section's first visible position, never before it
    const width = this.viewer.view.get(section)?.width ?? 0;
    const local = Math.min(Math.max(x - this.viewer.dataX, 0), width);
    // the unit under the pointer, from the placement (a text's in reading order): what a zoom keeps where it is
    const at = this.offsetAt(section, x, y);
    const atY = this.offsetYAt(section, x, y);
    // on macOS Ctrl is no zoom modifier, so a Ctrl+wheel is a pinch whatever its delta; elsewhere a large delta is a mouse
    const pinch =
      e.ctrlKey && pointer.pinchFromCtrlWheel && (this.mac || (!e.deltaMode && Math.abs(dy) <= pointer.pinchMaxDelta));
    const plot = this._overPlot(section, x, y);
    if (pinch) {
      // the pinch's own modifiers: none (the Ctrl is the gesture's, not a key held)
      const taken = this.viewer.bindings.resolve('pinch', NoModifiers, target) !== undefined;
      if (taken) nav.pinch(section, Math.exp(-dy * pointer.pinchSpeed), local, NoModifiers, target, at, atY);
      this.viewer.schedule();
      return taken || plot;
    }
    const mods = this._modifiers(e);
    // the wheel's chain (_spec/README.md "Input", The one rail): its vertical delta scrolls a rows band or a lanes band
    // under the pointer first, then the section past its cap (a text's lines, a tracks section's tracks), then the page
    const overflowing = scrolledSection(this.viewer.frame?.layout, section) !== undefined;
    const textBox = this.viewer.placementOf(section)?.clip === true && this.viewer.view.layoutOf(section) === 'text';
    // a text takes no other wheel (Cmd/Ctrl+wheel zooms it as any section)
    const vertical = !mods.mod && !mods.shift && (textBox ? !e.ctrlKey && !e.metaKey : Math.abs(dy) > Math.abs(dx));
    // over a rows band (its cells, its header) the wheel's vertical delta scrolls its rows, the rest pans the columns
    const rowsHit =
      this._hit.ySystem !== undefined || target === 'row' || target === 'summary-row' || target === 'tree-node';
    const rows = rowsHit && !mods.mod && !mods.shift ? this.viewer.behavior<RowsBehavior>('rows') : undefined;
    // only over a rows layer's own band: a `row` hit of another layer (a gutter's) pans as any
    const over = rows && this._hit.layer !== undefined ? rows.lay(trackOfLayer(this._hit.layer)) : undefined;
    if (rows && over?.key === this._hit.layer && this._hit.layer !== undefined && Math.abs(dy) >= Math.abs(dx)) {
      const scrolled = rows.wheel(trackOfLayer(this._hit.layer), dy);
      this.viewer.schedule();
      if (scrolled || !overflowing) return scrolled;
      return this._scrollSection(section, dy);
    }
    // over a ranges layer whose lanes overflow its band the wheel's vertical delta scrolls them (the rows behavior's); at
    // their end it goes on to the section, then the page
    const lanesLayer = this._hitAt(x, y)?.layer;
    if (lanesLayer !== undefined && vertical && !textBox) {
      if (this.viewer.behavior<RowsBehavior>('rows')?.wheelLanes(lanesLayer, dy)) {
        this.viewer.schedule();
        return true;
      }
    }
    if (textBox && !e.ctrlKey && !e.metaKey) return this._scrollBox(section, dy);
    if (overflowing && vertical) return this._scrollSection(section, dy);
    const taken = this.viewer.bindings.resolve('wheel', mods, target) !== undefined;
    if (taken) nav.wheel(section, dx, dy, local, mods, target, at, atY);
    // the frame asked for after the behavior's own, so a wheel's deltas apply in the frame that draws them
    this.viewer.schedule();
    return taken || plot;
  }

  /**
   * Whether a point is over a grid section's plot: its data box (x's extent, the blocks' gaps included) and the well's
   * height round its cells; never the label column or the y blocks' labels in it.
   */
  private _overPlot(section: string, x: number, y: number): boolean {
    // the grid layout by its registered name, never its class (a composition that lists no grid bundles none)
    const s = this._sectionOf(section);
    if (s?.layout.name !== 'grid' || !s.placement.strips.n) return false;
    const p = s.placement as GridPlacement;
    const r = p.xRow;
    const a = p.strips.axes[0];
    const x0 = r ? r.x0 : a.x0;
    const x1 = r ? r.x0 + r.width : a.x0 + a.width;
    const b = p.box;
    return x >= x0 && x < x1 && y >= b.y && y < b.y + b.height;
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
    // on a rows rail: its window resets the rows' zoom; outside it nothing (its first click centred, or not)
    if (hit?.target === 'window' && upright(hit.navigator)) {
      if (hit.navigator === RowsRail && hit.part !== 'outside' && hit.layer !== undefined)
        this.rowsRailGesture?.resetZoom(hit.layer);
      return true;
    }
    // on a ruler rail: the navigator's double-click (fit, or fitted a window `navigator.zoomPx` px wide there)
    const rail = hit?.target === 'window' && hit.navigator === 'rail' ? this._railOf(section) : undefined;
    const nav = this.navigatorBehavior;
    if (rail && nav)
      return nav.doubleClick(
        section,
        railUnitAt(rail, x),
        (this.viewer.config.navigator.zoomPx / Math.max(1e-9, rail.width)) * rail.total,
      );
    // over a plot every point is a cell's, a gap between blocks too: it zooms as a cell does, never fits
    const p = this.viewer.placementOf(section);
    const inPlot = this._overPlot(section, x, y);
    const target = inPlot ? 'cell' : hitTarget(hit);
    // a binding scoped to the hit's layer type (an extension's) wins there: run with the hit
    const type = hit?.layer === undefined ? undefined : this._layerTypeOf(section, hit.layer);
    const own = type && this.viewer.bindings.resolveIn(type, 'double-click', mods, target);
    if (own !== undefined && hit) {
      this.viewer.runCommand(own, { section, hit: { ...hit } });
      return true;
    }
    const action = this.viewer.bindings.resolve('double-click', mods, target);
    if (action === undefined) return false;
    // a clade's command on a summary row (its element the clade): the clade it names
    if (hit?.target === 'summary-row' && hit.element && hit.layer !== undefined && CladeActions.includes(action)) {
      void this.viewer.execute(action, { track: trackOfLayer(hit.layer), clade: hit.element.key });
      return true;
    }
    // a step zoom is about the pointer (the unit the placement has under it), by the double-click's own factor
    const pt = this._point;
    const at = p?.at(this._pageAt(x, y), pt) && (pt.kind === 'data' || pt.kind === 'gap') ? p.offsetAt(x, y) : -1;
    // a grid's cell (or a gap in its plot): its focus zooms a step about the pointer on both axes (`cell`), never a
    // jump to its two residues; on a grid whose y follows x the point keeps its place on x, y following
    const cell = inPlot && at >= 0;
    const atY = cell ? this.offsetYAt(section, x, y) : undefined;
    // on a block with nothing under the point (its label, its band's ground) or only a position (a ruler, letters, a
    // value: no item to focus): fit it, and from that fit the whole (`navigator.doubleClick: toggle-fit`, one level
    // deeper; never a zoom to one base of a chromosome)
    const positional = !!hit && !hit.element && !hit.extent && !hit.hull;
    const hull = hit ? this._focusOf(section, hit) : undefined;
    const block =
      !cell && ((!hull && action === 'view.fit') || positional) ? this._blockFit(section, p, x, y) : undefined;
    if (block === 'whole' && positional) {
      this.viewer.runCommand('view.fit', { section });
      return true;
    }
    if (block && block !== 'whole') {
      this.viewer.runCommand('view.focus', { section, focus: block });
      return true;
    }
    // a position alone (a ruler's) where a unit is under a pixel: one zoom step about the pointer, never the closest
    // zoom to a base (chr17 at the fit to 20 bases); a unit that is drawn (a letter) still focuses
    const shown = this.viewer.frame?.layout.sections.find((x) => x.section.address === section);
    const scale = shown?.layout.name === 'tracks' ? Math.abs(shown.axis?.scale ?? 1) : 1;
    if (positional && !cell && !block && at >= 0 && scale < 1 && action === 'view.focus') {
      this.viewer.runCommand('view.zoom-in', {
        section,
        at,
        step: this.viewer.config.navigate.doubleClickZoomStep,
      });
      return true;
    }
    const focus = hull ?? (cell ? { start: Math.floor(at), end: Math.floor(at) + 1 } : undefined);
    const location = hull && hit ? this._locationOf(hit) : undefined;
    this.viewer.runCommand(action, {
      section,
      ...(focus ? { focus } : {}),
      ...(location ? { location } : {}),
      ...(at >= 0 ? { at } : {}),
      ...(cell ? { cell: true } : {}),
      ...(atY !== undefined ? { atY } : {}),
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
    const r = this._focusOn(block ? block.axis : s.axis, hit, !block && this.viewer.view.arcs(section));
    return r && at ? { start: r.start + at, end: r.end + at } : r;
  }

  /**
   * The units of the block under a point of a row of blocks (its strip, or a run's block), `whole` when the window
   * shows it already (then a fit of the whole); undefined elsewhere.
   */
  private _blockFit(
    section: string,
    p: Placement | undefined,
    x: number,
    y: number,
  ): { start: number; end: number } | 'whole' | undefined {
    const s = this.viewer.view.resolved.value?.sections.find((v) => v.address === section);
    if (!s?.blocks || !(p instanceof BlockPlacement) || p.stacked) return undefined;
    const pt = this._point;
    if (!p.at(this._pageAt(x, y), pt)) return undefined;
    const b = pt.strip >= 0 ? p.strips.block[pt.strip] : p.row.blockAt(x);
    const block = b >= 0 ? s.blocks[b] : undefined;
    if (!block || block.offset === undefined) return undefined;
    const start = block.offset;
    const end = start + blockSpan(block);
    const r = p.row;
    // already this block (within a unit either side): the whole, as a fit
    if (Math.abs(r.u0 - start) <= 1 && Math.abs(r.u1 - end) <= 1) return 'whole';
    return { start, end };
  }

  /** What a focused hit names: the element's whole extent, else its position (the view's last jump says it). */
  private _locationOf(hit: Hit): Location | undefined {
    const extent = (hit.element && this.viewer.selection.extentOfHit(hit)) || hit.extent;
    if (extent) return extent;
    const p = hit.position;
    return p === undefined ? undefined : rangesLocation(hit.system, [Math.floor(p), Math.floor(p) + 1]);
  }

  private _focusOn(axis: Axis, hit: Hit, arc: boolean): { start: number; end: number } | undefined {
    // an item's whole extent (a layer's hit is the piece hit: one row of an item over a cyclic origin)
    const extent = (hit.element && this.viewer.selection.extentOfHit(hit)) || hit.extent;
    const part = extent?.parts.find((p) => p.system === hit.system);
    // round a ring: the arc its pieces cover, the circle less the widest gap between them (over the origin)
    if ((arc || axis.cyclic) && part && part.ranges.length >= 4) {
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
    // Esc hides the card (whatever else it does)
    if (e.key === 'Escape') this.cardBehavior?.interrupt();
    // Esc during a drag-select cancels it: the selection before it comes back, and the click after it is none
    // Esc during a window's drag puts the view back where the press found it
    if (e.key === 'Escape' && this._drag?.window) {
      if (this._drag.rowsRail?.hit.navigator === SectionRail) this.sectionRailGesture.cancel();
      else if (this._drag.rowsRail) this.rowsRailGesture?.cancel();
      else if (this._drag.rail) this.navigatorBehavior?.cancel();
      else this.windowBehavior?.cancel();
      this._dragHit = undefined;
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
    // a running drag's Esc is its cancel (above); then an open keys panel's: it closes, and does nothing else
    if (e.key === 'Escape' && this.keysPanel?.close()) return true;
    const section =
      (this.resting ? this.sectionAt(this._y) : undefined) ?? this.viewer.view.resolved.value?.sections[0]?.address;
    // the layer under the pointer: the lanes band a page key scrolls
    const ctx = { section, layer: this.resting ? this._hitAt(this._x, this._y)?.layer : undefined };
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
        // a key that moved the keyboard cursor: the card shows at its cell, after the frame that revealed it
        if (c.id.startsWith('cursor.')) this.cardBehavior?.cursorMoved();
        return true;
      }
    // a shortcut with modifiers is the viewer's even when nothing can run now, but the history's (Alt+Left is the
    // browser's Back when the view has no step back); a plain key nothing takes (Esc with nothing to cancel) stays the
    // host's
    return matched && !!(e.altKey || e.shiftKey || e.metaKey || e.ctrlKey);
  }
}

/**
 * Commands whose keys and buttons are the page's when they cannot run: the view's history (the browser's Back), the
 * cursor's first and last row with no cursor on rows (Cmd+↑ ↓ scroll the page to its ends).
 */
const YieldWhenIdle: ReadonlySet<string> = new Set([
  'view.back',
  'view.forward',
  'cursor.first-row',
  'cursor.last-row',
]);

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

/** Whether a point is over a band of a section that scrolls on its own (a lanes or a rows band with its rail). */
function overBand(s: SectionLayout, x: number, y: number): boolean {
  return overRails(s.lanesRails, x, y) || overRails(s.rowsRails, x, y);
}

/** Whether a point is beside one of upright rails along it, left of its right side (its band's rows). */
function overRails(list: readonly Rail[] | undefined, x: number, y: number): boolean {
  if (list)
    for (let i = 0; i < list.length; i++) {
      const r = list[i];
      if (y >= r.x - r.inset && y < r.x + r.width + r.inset && x < r.bottom) return true;
    }
  return false;
}

/** An upright rail's navigator: a rows band's, or a pan-only one (a lanes band's, a section's): the one gesture. */
function upright(n: string | undefined): boolean {
  return n === RowsRail || panOnly(n);
}

/** A rail's window as drawn, in display units, when the floor drew it wider than the truth (else undefined). */
function railDrawn(rail: Rail, viewer: SeqViewer): readonly [number, number] | undefined {
  const k = rail.width / Math.max(1e-9, rail.total);
  const a = rail.x + rail.start * k;
  const b = rail.x + Math.min(rail.end, rail.total) * k;
  const out = RailEnds;
  out[0] = a;
  out[1] = b;
  out[2] = rail.x;
  out[3] = rail.x + rail.width;
  floorExtent(out, 0, viewer.context.renderConfig.navigator.window);
  if (out[1] - out[0] <= b - a + 1e-9) return undefined;
  return [railUnitAt(rail, out[0]), railUnitAt(rail, out[1])];
}

const RailEnds = new Float64Array(4);
