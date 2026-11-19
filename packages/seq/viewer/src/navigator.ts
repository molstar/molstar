/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ReactiveModel, type ReadonlyState, State } from '@molstar/core/reactive';
import type { Node } from '@molstar/core/repr/tree';
import type { Want } from '@molstar/seqx-core/data';
import {
  type NavigatorFrame,
  type NavigatorHit,
  NavigatorRenderer,
  NavWindowFlags,
  navigatorBox,
  type WindowPart,
  windowText,
} from '@molstar/seqx-core/render';
import { axisTotal, type ResolvedNavigator, type ResolvedView, resolveNavigatorNode } from '@molstar/seqx-core/view';
import {
  NavigateProvider,
  type NavigatorBehavior,
  NavigatorProvider,
  windowIdle,
} from '@molstar/seqx-viewer/behaviors';
import type { ModifierState } from '@molstar/seqx-viewer/bindings';
import type { ViewSnapshot } from '@molstar/seqx-viewer/models';
import { DragSlop } from './input.js';
import type { SeqViewer } from './viewer.js';

// The navigator's model (_spec/README.md "Models", NavigatorModel): DOM-free, one a navigator a viewer shows (a
// declared node of the view by its id, or a host's definition, `viewer.navigator(def)`). It holds the navigator's
// picture (its own renderer over the viewer's store, theme, image cache and held domains), the bound windows' frames as
// the viewer's last frame placed them, the lit part, idle per window and axis, the gesture in progress, the focused
// slider and the readout; the canvas view and an export read it. Its pointer and keys are in its own px, carried to the
// bound section's display units through the picture, and run through the `navigator` behavior. A hover stays here (a
// readout, never the viewer's hover).

/** A host's navigator: the params of a `navigator` node with its id, not saved in the spec. */
export interface NavigatorDef {
  readonly id: string;
  readonly of: string | readonly string[];
  readonly axes?: 'positions' | 'rows' | 'both';
  readonly content?: 'axis' | 'tracks' | 'bands' | 'rows';
  readonly tracks?: readonly string[];
  readonly bands?: string;
  readonly system?: string;
  readonly place?: 'top' | 'bottom' | 'left' | 'right';
  readonly size?: number | readonly [number, number];
}

/** A window as its slider says it (its `aria-valuetext`): set when the view settles, never a frame. */
export interface NavigatorValue {
  readonly section: string;
  readonly text: string;
  readonly disabled: boolean;
}

/** The pointer's state for a gesture: touch widens the hits. */
export interface NavigatorPointer {
  readonly touch?: boolean;
}

/** A press in progress: where, what it hit, whether it moved past the slop; a rows press's band and its step. */
interface Press {
  readonly x0: number;
  readonly y0: number;
  readonly hit: NavigatorHit;
  moved: boolean;
  /** The press drags the positions (through the `navigator` behavior). */
  readonly positions: boolean;
  /** The press drags the rows: their band's track, its scroll at the press, the viewports then (a far drag's step). */
  readonly rows?: { readonly track: string; readonly top: number; readonly before: ViewSnapshot };
}

export class NavigatorModel extends ReactiveModel {
  /** Bumped when what the navigator draws changes (hover, a gesture, its size): its view redraws. */
  private readonly _version = new State(0);
  readonly version: ReadonlyState<number> = this._version.asReadonly;
  /** The readout under the pointer (a position, a row, a band's name), or none. */
  private readonly _readout = new State<string | undefined>(undefined);
  readonly readout: ReadonlyState<string | undefined> = this._readout.asReadonly;
  /** Each window's slider value, set when the view settles. */
  private readonly _values = new State<readonly NavigatorValue[]>([]);
  readonly values: ReadonlyState<readonly NavigatorValue[]> = this._values.asReadonly;
  /** The slider that has the focus (roving), by its window's index; a press on a window focuses its slider. */
  private readonly _focused = new State(0);
  readonly focused: ReadonlyState<number> = this._focused.asReadonly;
  /** Its size, CSS px: the host's (`resize`), each the form's box for the view where the host gives none. */
  width = 0;
  height = 0;
  private hostWidth = 0;
  private hostHeight = 0;
  private renderer: NavigatorRenderer | undefined;
  private rendererFor: unknown[] = [];
  private readonly wanter;
  private lastFrame: NavigatorFrame | undefined;
  private lit: { section: string; part: WindowPart } | undefined;
  private pointerX = Number.NaN;
  private pointerY = Number.NaN;
  /** The section the pointer last acted on: its window on top. */
  private top: string | undefined;
  private readonly moves: Record<string, number> = {};
  private press: Press | undefined;
  /** The part has the keyboard focus (the host says: `setFocus`): a move made elsewhere then announces nothing. */
  private hasFocus = false;
  /** The view rested since the last draw; this part's own gesture or key moved it. */
  private restDue = false;
  private ownMove = false;

  constructor(
    readonly viewer: SeqViewer,
    readonly id: string,
    private readonly def?: NavigatorDef,
  ) {
    super();
    this.wanter = viewer.store.store.wanter();
    this.registerPermanentDispose(() => this.wanter.dispose());
    this.listenPermanently(viewer.drawn, () => this.changed());
    this.listenPermanently(viewer.rested, () => {
      this.restDue = true;
      this.changed();
    });
    this.listenPermanently(viewer.tokensVersion, () => this.renderer?.setTokens(viewer.tokens));
  }

  /** The navigator as resolved over the viewer's view; undefined when the view has none of this id. */
  get navigator(): ResolvedNavigator | undefined {
    const view = this.viewer.view.resolved.value;
    if (!view) return undefined;
    if (!this.def) return view.navigators.find((n) => n.id === this.id || n.address === this.id);
    return this.hostNavigator(view);
  }

  private hosted: [ResolvedView, ResolvedNavigator] | undefined;

  private hostNavigator(view: ResolvedView): ResolvedNavigator {
    if (this.hosted?.[0] === view) return this.hosted[1];
    const { id, ...params } = this.def!;
    const node: Node = { kind: 'navigator', id, params: params as never };
    const n = resolveNavigatorNode(
      this.viewer.spec.model.value,
      view.sections,
      node,
      `${view.address}/navigators/${id}`,
    );
    this.hosted = [view, n];
    return n;
  }

  /** Whether this is a host's navigator (`viewer.navigator(def)`, never saved). */
  get isHosted(): boolean {
    return this.def !== undefined;
  }

  /** Sets its size (a host's element): a side of 0 is the form's box's for the viewer's width and height. */
  resize(width: number, height = 0): void {
    if (width === this.hostWidth && height === this.hostHeight) return;
    this.hostWidth = Math.max(0, width);
    this.hostHeight = Math.max(0, height);
    this.place();
    this.changed();
  }

  /** Its size from the host's and the form's box. */
  private place(): void {
    const nav = this.navigator;
    if (!nav) return;
    const b = navigatorBox(nav, this.hostWidth || this.viewer.width, this.viewer.height);
    this.width = this.hostWidth || b.width;
    this.height = this.hostHeight || b.height;
  }

  private changed(): void {
    this._version.set(this._version.value + 1);
  }

  private rendererOf(nav: ResolvedNavigator): NavigatorRenderer {
    const model = this.viewer.spec.model.value;
    const view = this.viewer.view.resolved.value!;
    const viewRenderer = this.viewer.renderer;
    const f = this.rendererFor;
    if (!this.renderer || f[0] !== model || f[1] !== view || f[2] !== nav || f[3] !== viewRenderer) {
      this.renderer = new NavigatorRenderer({
        model,
        view,
        navigator: nav,
        data: this.viewer.store.store,
        tokens: this.viewer.tokens,
        text: this.viewer.context.text,
        config: this.viewer.context.renderConfig,
        direction: this.viewer.config.navigator.direction,
        storeSettings: this.viewer.store.store.settings,
        // the view's image budget and held domains: a value has one color in both
        ...(viewRenderer ? { images: viewRenderer.images, domains: viewRenderer.layoutState } : {}),
      });
      this.rendererFor = [model, view, nav, viewRenderer];
    }
    return this.renderer;
  }

  /** What a gesture can move: each bound window's axes (`NavWindowFlags`), 0 idle. */
  private updateMoves(nav: ResolvedNavigator): void {
    const behavior = this.viewer.behaviors.has(NavigatorProvider.name);
    const view = this.viewer.view;
    for (const a of nav.of) {
      let m = 0;
      if (behavior && view.canMove(a)) {
        if (nav.axes !== 'rows' && !windowIdle(view, a, this.viewer.config)) m |= NavWindowFlags.x;
        if (nav.axes !== 'positions') {
          const w = this.lastFrame?.windows.find((x) => x.section === a);
          if (w && w.rows > w.r1 - w.r0 + 1e-9) m |= NavWindowFlags.y;
        }
      }
      this.moves[a] = m;
    }
  }

  /**
   * Draws a frame of the navigator over the viewer's last frame (its picture kept while its key holds): undefined
   * when the navigator is hidden, unknown, or the viewer has drawn nothing.
   */
  draw(): NavigatorFrame | undefined {
    const nav = this.navigator;
    const layout = this.viewer.frame?.layout;
    if (!nav || nav.hidden || !layout) return undefined;
    this.place();
    if (!(this.width > 0) || !(this.height > 0)) return undefined;
    const r = this.rendererOf(nav);
    this.updateMoves(nav);
    const sweep = this.behavior?.sweep.value;
    const f = r.frame({
      size: { width: this.width, height: this.height },
      dpr: this.viewer.dpr,
      settled: this.viewer.store.store.settled.value,
      overlay: this.viewer.store.overlay(
        this.viewer.spec.model.value,
        this.viewer.view.resolved.value!,
        this.viewer.context,
      ),
      view: layout,
      moves: this.moves,
      ...(this.lit ? { lit: this.lit } : {}),
      ...(this.top ? { top: this.top } : {}),
      pointerX: this.pointerX,
      pointerY: this.pointerY,
      emphasis: this.viewer.emphasis,
      cvdSafe: this.viewer.legends.cvdSafe.value,
      ...(sweep && this.press ? { sweep: this.sweepPx(sweep) } : {}),
    });
    if (f.rebuilt) this.wanter.want(f.wants as Want[]);
    this.lastFrame = f;
    // the sliders' values: first drawn, and when the view rested (a move made elsewhere while the part has the focus
    // announces nothing: its value waits for the focus to leave)
    if (
      !this.press &&
      (this._values.value.length !== f.windows.length || (this.restDue && (this.ownMove || !this.hasFocus)))
    ) {
      this.settle(nav, f);
      this.ownMove = false;
    }
    this.restDue = false;
    return f;
  }

  /** The navigator's last frame. */
  get frame(): NavigatorFrame | undefined {
    return this.lastFrame;
  }

  /** The issues the navigator met (its bindings, what it left out, its picture's). */
  get issues() {
    return this.renderer?.allIssues ?? [];
  }

  /** The part took or lost the keyboard focus (its host's element): losing it, its sliders say where the view is. */
  setFocus(focused: boolean): void {
    if (this.hasFocus === focused) return;
    this.hasFocus = focused;
    const nav = this.navigator;
    if (!focused && nav && this.lastFrame) this.settle(nav, this.lastFrame);
  }

  /** The sliders' values (never a frame): compared, set only when they change. */
  private settle(nav: ResolvedNavigator, f: NavigatorFrame): void {
    const model = this.viewer.spec.model.value;
    const prev = this._values.value;
    let same = prev.length === f.windows.length;
    const next: NavigatorValue[] = [];
    for (let i = 0; i < f.windows.length; i++) {
      const w = f.windows[i];
      const v: NavigatorValue = {
        section: w.section,
        text: w.shown ? windowText(model, nav, w) : '',
        disabled: !w.shown || !this.moves[w.section],
      };
      if (same && (prev[i].text !== v.text || prev[i].disabled !== v.disabled || prev[i].section !== v.section))
        same = false;
      next.push(v);
    }
    if (!same) this._values.set(next);
  }

  private get behavior(): NavigatorBehavior | undefined {
    return this.viewer.behavior(NavigatorProvider);
  }

  /** A display unit range of the bound section as x on the navigator (a sweep's). */
  private sweepPx(s: readonly [number, number]): { x0: number; x1: number } {
    const r = this.renderer!;
    const section = this.press!.hit.section;
    return { x0: r.xOfSectionUnit(section, s[0]), x1: r.xOfSectionUnit(section, s[1]) };
  }

  /** What is under a point, idle windows' parts read as `outside` (an idle window is no target). */
  hit(x: number, y: number, touch = false): NavigatorHit | undefined {
    const h = this.renderer?.hit(x, y, touch, this.top);
    if (!h) return undefined;
    const m = this.moves[h.section] ?? 0;
    const along = h.axis === 'rows' ? NavWindowFlags.y : NavWindowFlags.x;
    if (!m) return { ...h, part: 'outside' };
    if (h.part !== 'outside' && h.part !== 'body' && !(m & along)) return { ...h, part: 'body' };
    return h;
  }

  /** The cursor over the navigator: `grab` over a window that can move, `grabbing` during a drag. */
  get cursor(): 'grab' | 'grabbing' | 'default' | 'pointer' {
    if (this.press) return this.press.hit.part === 'outside' ? 'default' : 'grabbing';
    const l = this.lit;
    if (!l) return 'default';
    return l.part === 'outside' ? 'pointer' : 'grab';
  }

  /** The pointer moved over the navigator: the lit part and the readout follow; a drag moves its window. */
  move(x: number, y: number, p: NavigatorPointer = {}): void {
    this.pointerX = x;
    this.pointerY = y;
    const d = this.press;
    if (d) {
      if (!d.moved && Math.max(Math.abs(x - d.x0), Math.abs(y - d.y0)) >= DragSlop) d.moved = true;
      if (d.moved) this.dragTo(d, x, y);
      this.changed();
      return;
    }
    const h = this.hit(x, y, p.touch);
    this.lit = h && (this.moves[h.section] ?? 0) ? { section: h.section, part: h.part } : undefined;
    this._readout.set(h ? this.readoutOf(h) : undefined);
    this.changed();
  }

  private readoutOf(h: NavigatorHit): string {
    const nav = this.navigator;
    const name = (nav?.system && this.viewer.spec.model.value.systems.get(nav.system)?.name) ?? '';
    const parts: string[] = [];
    if (h.position !== undefined) parts.push(`${name}:${(h.position + 1).toLocaleString('en-US')}`);
    if (h.row !== undefined) parts.push(`row ${(Math.floor(h.row) + 1).toLocaleString('en-US')}`);
    if (h.band !== undefined) parts.push(h.band);
    return parts.join(' · ');
  }

  /** The pointer left: no lit part, no readout, no hairline. */
  leave(): void {
    this.pointerX = this.pointerY = Number.NaN;
    this.lit = undefined;
    this._readout.set(undefined);
    this.changed();
  }

  /** The bound rows track a rows window moves (its section's first rows track). */
  private rowsTrack(section: string): string | undefined {
    const s = this.viewer.view.resolved.value?.sections.find((x) => x.address === section);
    return s?.tracks.find((t) => t.layers.some((l) => l.type === 'rows'))?.address;
  }

  /** The window of a bound section in the last frame. */
  private windowOf(section: string) {
    return this.lastFrame?.windows.find((v) => v.section === section);
  }

  /** The rows band's scroll clamped to what it can show (its rows less those in view). */
  private clampTop(section: string, top: number): number {
    const w = this.windowOf(section);
    const most = w ? Math.max(0, w.rows - (w.r1 - w.r0)) : Number.POSITIVE_INFINITY;
    return Math.max(0, Math.min(most, top));
  }

  /**
   * The primary button went down: true when a gesture began (the host captures the pointer). A press on a window
   * focuses its slider. A miniature's body drags both axes; elsewhere the axis of the part hit.
   */
  down(x: number, y: number, p: NavigatorPointer = {}): boolean {
    const b = this.behavior;
    const nav = this.navigator;
    const h = this.hit(x, y, p.touch);
    if (!b || !nav || !h || !(this.moves[h.section] ?? 0)) return false;
    this.top = h.section;
    const i = nav.of.indexOf(h.section);
    if (i >= 0 && i !== this._focused.value) this._focused.set(i);
    const m = this.moves[h.section] ?? 0;
    const inside = h.part !== 'outside';
    // what the press drags: a rows part its rows, a positions part its positions, a miniature's body both; outside the
    // window a positions press (a click or a sweep), and on rows a press for its click
    const both = nav.axes === 'both' && h.part === 'body';
    let rows: Press['rows'];
    if (inside && (h.axis === 'rows' || both) && m & NavWindowFlags.y) {
      const track = this.rowsTrack(h.section);
      if (track) rows = { track, top: this.viewer.view.band(track)?.bandTop ?? 0, before: this.viewer.view.snapshot() };
    }
    const dragPositions =
      nav.axes !== 'rows' && !!(m & NavWindowFlags.x) && (!inside || h.axis === 'positions' || both);
    const positions = dragPositions && this.pressPositions(b, h, x);
    if (!positions && !rows && (inside || nav.axes === 'positions')) return false;
    this.press = { x0: x, y0: y, hit: h, moved: false, positions, ...(rows ? { rows } : {}) };
    this.changed();
    return true;
  }

  private pressPositions(b: NavigatorBehavior, h: NavigatorHit, x: number): boolean {
    const r = this.renderer!;
    const at = r.sectionUnitAt(h.section, x);
    if (at === undefined) return false;
    const w = this.windowOf(h.section);
    const drawn =
      w && w.dx1 - w.dx0 > w.x1 - w.x0 + 1e-9
        ? ([r.sectionUnitAt(h.section, w.dx0) ?? at, r.sectionUnitAt(h.section, w.dx1) ?? at] as const)
        : undefined;
    return b.press(h.section, h.part, at, drawn && drawn[0] <= drawn[1] ? drawn : undefined);
  }

  private dragTo(d: Press, x: number, y: number): void {
    const r = this.renderer!;
    if (d.rows && d.hit.part !== 'outside') {
      const row = r.rowAt(y);
      const row0 = r.rowAt(d.y0);
      if (row !== undefined && row0 !== undefined)
        this.viewer.view.setBand(d.rows.track, { bandTop: this.clampTop(d.hit.section, d.rows.top + row - row0) });
    }
    if (!d.positions) return;
    const at = r.sectionUnitAt(d.hit.section, x);
    if (at !== undefined) this.behavior?.dragTo(at);
  }

  /**
   * The button went up: a drag ends (a sweep zooms, a far drag pushes its step: a rows drag farther than the rows in
   * view pushes the viewports at its press); a press that never moved waits for its click (`click`).
   */
  up(): void {
    const d = this.press;
    if (!d) return;
    this.press = undefined;
    const b = this.behavior;
    if (!d.moved) {
      if (d.positions) b?.cancel();
    } else {
      if (d.positions) b?.release();
      if (d.rows) {
        const w = this.windowOf(d.hit.section);
        const now = this.viewer.view.band(d.rows.track)?.bandTop ?? d.rows.top;
        const step = this.viewer.config.navigator.dragStep;
        const far = w !== undefined && Math.abs(now - d.rows.top) > w.r1 - w.r0;
        if (step === 'always' || (step === 'far' && far)) this.viewer.view.record(d.rows.before);
      }
      this.ownMove = true;
    }
    this.lastPress = d;
    this.changed();
  }

  /** The last press that ended: its click (the platform's, after the button went up) reads it. */
  private lastPress: Press | undefined;

  /**
   * A click (the platform's, after a press that never moved; `detail` its count): outside the window it centres
   * there, both axes on a miniature (a jump, its own step); on the window nothing. A second click is a
   * double-click's: nothing (`doubleClick` follows).
   */
  click(x: number, y: number, detail = 1): boolean {
    const d = this.lastPress;
    this.lastPress = undefined;
    const b = this.behavior;
    const nav = this.navigator;
    if (!d || d.moved || !b || !nav || detail > 1) return false;
    if (d.hit.part !== 'outside') {
      b.click(d.hit.section, undefined);
      return false;
    }
    this.ownMove = true;
    let moved = false;
    if (nav.axes !== 'positions') {
      const track = this.rowsTrack(d.hit.section);
      const w = this.windowOf(d.hit.section);
      const row = this.renderer?.rowAt(y);
      if (track && w && row !== undefined)
        moved =
          this.viewer.view.jump(() => {
            const before = this.viewer.view.band(track)?.bandTop;
            this.viewer.view.setBand(track, { bandTop: this.clampTop(d.hit.section, row - (w.r1 - w.r0) / 2) });
            return this.viewer.view.band(track)?.bandTop !== before;
          }) === 'moved';
      if (nav.axes === 'rows') return moved;
    }
    const at = this.renderer?.sectionUnitAt(d.hit.section, x);
    return (at !== undefined && b.click(d.hit.section, at, detail)) || moved;
  }

  /** A double-click at a point: fit, or fitted, a window `navigator.zoomPx` navigator px wide there. */
  doubleClick(x: number, y: number): boolean {
    const b = this.behavior;
    const nav = this.navigator;
    const h = this.hit(x, y);
    if (!b || !nav || !h || nav.axes === 'rows' || !(this.moves[h.section] ?? 0)) return false;
    const r = this.renderer!;
    const at = r.sectionUnitAt(h.section, x);
    const half = this.viewer.config.navigator.zoomPx / 2;
    const a = r.sectionUnitAt(h.section, x - half);
    const c = r.sectionUnitAt(h.section, x + half);
    if (at === undefined || a === undefined || c === undefined) return false;
    this.ownMove = true;
    return b.doubleClick(h.section, at, Math.abs(c - a));
  }

  /**
   * The wheel over the navigator: as over the bound section (the same deltas in its px), a zoom about the window's
   * centre, never the pointer (a position here may be far outside the view). True when taken.
   */
  wheel(x: number, y: number, dx: number, dy: number, modifiers: ModifierState): boolean {
    const nav = this.navigator;
    const h = this.hit(x, y);
    const section = h?.section ?? nav?.of[0];
    const v = section && this.viewer.view.get(section);
    const navigate = this.viewer.behavior(NavigateProvider);
    if (!section || !v || !navigate) return false;
    navigate.wheel(section, dx, dy, v.width / 2, modifiers, undefined as never, (v.start + v.end) / 2);
    this.ownMove = true;
    return true;
  }

  /** Esc during a drag or a sweep: cancelled, the view back where the press found it. */
  cancel(): boolean {
    const d = this.press;
    if (!d) return false;
    this.press = undefined;
    if (d.rows) this.viewer.view.setBand(d.rows.track, { bandTop: d.rows.top });
    if (d.positions) this.behavior?.cancel();
    this.changed();
    return true;
  }

  /**
   * A key while the part has the focus (its window sliders, roving: PageDown and PageUp move the focus to the next
   * and the previous window): arrows pan (Up and Down on rows), Shift+arrows resize about the centre, Home and End go
   * to the ends (jumps), Enter centres on the readout's position, Esc cancels; the view's keys act on the focused
   * window's section (Alt+Left and Alt+Right whether it can move or not: they walk the history). True when taken.
   */
  key(key: string, modifiers: { readonly shiftKey?: boolean; readonly altKey?: boolean } = {}): boolean {
    const nav = this.navigator;
    const section = nav?.of[this._focused.value] ?? nav?.of[0];
    if (!nav || !section) return false;
    if (key === 'Escape') return this.cancel();
    if (key === 'PageDown' || key === 'PageUp') return this.focusNext(key === 'PageDown' ? 1 : -1);
    const view = this.viewer.view;
    if (modifiers.altKey && (key === 'ArrowLeft' || key === 'ArrowRight'))
      return key === 'ArrowLeft' ? view.back() : view.forward();
    if (!(this.moves[section] ?? 0)) return false;
    const b = this.behavior;
    const run = (id: string) => this.viewer.runCommand(id, { section }) !== '';
    const v = view.get(section);
    if (!v || !b) return false;
    this.ownMove = true;
    if (nav.axes !== 'positions' && (key === 'ArrowUp' || key === 'ArrowDown')) {
      const track = this.rowsTrack(section);
      const w = this.windowOf(section);
      if (!track || !w) return false;
      const step = Math.max(1, (w.r1 - w.r0) * this.viewer.config.navigate.keyPanFraction);
      const top = (view.band(track)?.bandTop ?? 0) + (key === 'ArrowUp' ? -step : step);
      view.setBand(track, { bandTop: this.clampTop(section, top) });
      return true;
    }
    switch (key) {
      case 'ArrowLeft':
      case 'ArrowRight': {
        const back = key === 'ArrowLeft';
        if (modifiers.shiftKey) return run(back ? 'window.shrink' : 'window.grow');
        return run(back ? 'window.pan-back' : 'window.pan-on');
      }
      case 'Home':
      case 'End': {
        const w = v.end - v.start;
        const axis = view.resolved.value?.sections.find((s) => s.address === section)?.axis;
        const total = axis && axisTotal(axis);
        if (total === undefined) return false;
        return b.zoomTo(section, key === 'Home' ? 0 : total - w, key === 'Home' ? w : total);
      }
      case 'Enter': {
        const at = Number.isFinite(this.pointerX) ? this.renderer?.sectionUnitAt(section, this.pointerX) : undefined;
        return at !== undefined && b.center(section, at);
      }
      case '=':
      case '+':
        return view.jump(() => view.zoom(section, this.viewer.config.navigate.doubleClickZoomStep)) === 'moved';
      case '-':
        return view.jump(() => view.zoom(section, 1 / this.viewer.config.navigate.doubleClickZoomStep)) === 'moved';
      case '0':
        return view.jump(() => view.fit(section)) === 'moved';
      default:
        return false;
    }
  }

  /** Moves the slider focus to the next or previous window (roving: one tab stop for the part); true when it moved. */
  focusNext(step = 1): boolean {
    const n = this.navigator?.of.length ?? 0;
    if (n < 2) return false;
    this._focused.set((((this._focused.value + step) % n) + n) % n);
    // the newly focused slider says its value
    const nav = this.navigator;
    if (nav && this.lastFrame) this.settle(nav, this.lastFrame);
    return true;
  }
}
