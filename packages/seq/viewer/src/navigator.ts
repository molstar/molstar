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
import { type ResolvedNavigator, type ResolvedView, resolveNavigatorNode } from '@molstar/seqx-core/view';
import {
  NavigateProvider,
  type NavigatorBehavior,
  NavigatorProvider,
  windowIdle,
} from '@molstar/seqx-viewer/behaviors';
import type { ModifierState } from '@molstar/seqx-viewer/bindings';
import type { SeqViewer } from './viewer.js';

// The navigator's model (_spec/README.md "Models", NavigatorModel): DOM-free, one a navigator a viewer shows (a
// declared node of the view by its id, or a host's definition, `viewer.navigator(def)`). It holds the navigator's
// picture (its own renderer over the viewer's store and theme), the bound windows' frames as the viewer's last frame
// placed them, the lit part, idle per window and axis, the gesture in progress and the readout; the canvas view and an
// export read it. Its pointer and keys are in its own px, carried to the bound section's display units through the
// picture, and run through the `navigator` behavior. A hover stays here (a readout, never the viewer's hover).

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

/** A window as its slider says it (its `aria-valuetext`), set when the viewport settles, never a frame. */
export interface NavigatorValue {
  readonly section: string;
  readonly text: string;
  readonly disabled: boolean;
}

/** The pointer's state for a gesture: touch widens the hits. */
export interface NavigatorPointer {
  readonly touch?: boolean;
}

/** The drag slop, CSS px: a press that moves less is a click. */
const Slop = 3;

export class NavigatorModel extends ReactiveModel {
  /** Bumped when what the navigator draws changes (hover, a gesture, its size): its view redraws. */
  private readonly _version = new State(0);
  readonly version: ReadonlyState<number> = this._version.asReadonly;
  /** The readout under the pointer (a position, a row, a band's name), or none. */
  private readonly _readout = new State<string | undefined>(undefined);
  readonly readout: ReadonlyState<string | undefined> = this._readout.asReadonly;
  /** Each window's slider value, set on settle. */
  private readonly _values = new State<readonly NavigatorValue[]>([]);
  readonly values: ReadonlyState<readonly NavigatorValue[]> = this._values.asReadonly;
  /** The slider that has the focus (roving), by its window's index. */
  focused = 0;
  /** Its size, CSS px: the host's (`resize`), each the form's box for the view where the host gives none. */
  width = 0;
  height = 0;
  private hostWidth = 0;
  private hostHeight = 0;
  private renderer: NavigatorRenderer | undefined;
  private rendererFor: [unknown, unknown, unknown] = [undefined, undefined, undefined];
  private readonly wanter;
  private lastFrame: NavigatorFrame | undefined;
  private lit: { section: string; part: WindowPart } | undefined;
  private pointerX = Number.NaN;
  private pointerY = Number.NaN;
  /** The section the pointer last acted on: its window on top. */
  private top: string | undefined;
  private readonly moves: Record<string, number> = {};
  private drag:
    | { x0: number; y0: number; hit: NavigatorHit; moved: boolean; rows?: { track: string; top: number } }
    | undefined;

  constructor(
    readonly viewer: SeqViewer,
    readonly id: string,
    private readonly def?: NavigatorDef,
  ) {
    super();
    this.wanter = viewer.store.store.wanter();
    this.registerPermanentDispose(() => this.wanter.dispose());
    this.listenPermanently(viewer.drawn, () => this.changed());
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
    const f = this.rendererFor;
    if (!this.renderer || f[0] !== model || f[1] !== view || f[2] !== nav) {
      this.renderer = new NavigatorRenderer({
        model,
        view,
        navigator: nav,
        data: this.viewer.store.store,
        tokens: this.viewer.tokens,
        text: this.viewer.context.text,
        config: this.viewer.context.renderConfig,
        direction: this.viewer.config.navigator.direction,
      });
      this.rendererFor = [model, view, nav];
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
      ...(sweep && this.drag ? { sweep: this.sweepPx(sweep) } : {}),
    });
    if (f.rebuilt) this.wanter.want(f.wants as Want[]);
    this.lastFrame = f;
    if (!this.drag) this.settle(nav, f);
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

  /** The sliders' values, said when the windows settle (never a frame): compared, set only when they change. */
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
    return { x0: r.xOfSectionUnit(this.drag!.hit.section, s[0]), x1: r.xOfSectionUnit(this.drag!.hit.section, s[1]) };
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
    if (this.drag) return this.drag.hit.part === 'outside' ? 'default' : 'grabbing';
    const l = this.lit;
    if (!l) return 'default';
    return l.part === 'outside' ? 'pointer' : 'grab';
  }

  /** The pointer moved over the navigator: the lit part and the readout follow; a drag moves its window. */
  move(x: number, y: number, p: NavigatorPointer = {}): void {
    this.pointerX = x;
    this.pointerY = y;
    const d = this.drag;
    if (d) {
      if (!d.moved && Math.max(Math.abs(x - d.x0), Math.abs(y - d.y0)) >= Slop) d.moved = true;
      if (d.moved) this.dragTo(x, y);
      this.changed();
      return;
    }
    const h = this.hit(x, y, p.touch);
    const lit = h && (this.moves[h.section] ?? 0) ? { section: h.section, part: h.part } : undefined;
    if (lit?.section !== this.lit?.section || lit?.part !== this.lit?.part) this.lit = lit;
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

  /** The primary button went down: true when a gesture began (the host captures the pointer). */
  down(x: number, y: number, p: NavigatorPointer = {}): boolean {
    const b = this.behavior;
    const h = this.hit(x, y, p.touch);
    if (!b || !h || !(this.moves[h.section] ?? 0)) return false;
    this.top = h.section;
    const r = this.renderer!;
    if (h.axis === 'rows' && h.part !== 'outside') {
      const track = this.rowsTrack(h.section);
      if (!track) return false;
      this.drag = {
        x0: x,
        y0: y,
        hit: h,
        moved: false,
        rows: { track, top: this.viewer.view.band(track)?.bandTop ?? 0 },
      };
      // a miniature's body pans both axes: the positions through the window behavior too
      if (this.navigator?.axes === 'both' && h.part === 'body') this.pressPositions(b, h, x);
      this.changed();
      return true;
    }
    if (this.navigator?.axes === 'rows') {
      this.drag = { x0: x, y0: y, hit: h, moved: false };
      return true;
    }
    if (!this.pressPositions(b, h, x)) return false;
    this.drag = { x0: x, y0: y, hit: h, moved: false };
    void r;
    this.changed();
    return true;
  }

  private pressPositions(b: NavigatorBehavior, h: NavigatorHit, x: number): boolean {
    const r = this.renderer!;
    const at = r.sectionUnitAt(h.section, x);
    if (at === undefined) return false;
    const w = this.lastFrame?.windows.find((v) => v.section === h.section);
    const drawn =
      w && w.dx1 - w.dx0 > w.x1 - w.x0 + 1e-9
        ? ([r.sectionUnitAt(h.section, w.dx0) ?? at, r.sectionUnitAt(h.section, w.dx1) ?? at] as const)
        : undefined;
    return b.press(h.section, h.part, at, drawn && drawn[0] <= drawn[1] ? drawn : undefined);
  }

  private dragTo(x: number, y: number): void {
    const d = this.drag!;
    const r = this.renderer!;
    if (d.rows) {
      const row = r.rowAt(y);
      const row0 = r.rowAt(d.y0);
      if (row !== undefined && row0 !== undefined)
        this.viewer.view.setBand(d.rows.track, { bandTop: Math.max(0, d.rows.top + row - row0) });
    }
    if (this.navigator?.axes === 'rows') return;
    if (d.rows && !(this.navigator?.axes === 'both' && d.hit.part === 'body')) return;
    const at = r.sectionUnitAt(d.hit.section, x);
    if (at !== undefined) this.behavior?.dragTo(at);
  }

  /** The button went up: a drag ends (a sweep zooms, a far drag pushes its step); a press that never moved is a click. */
  up(): void {
    const d = this.drag;
    if (!d) return;
    this.drag = undefined;
    const b = this.behavior;
    if (!d.moved) {
      // a click outside the window centres there; on the window it is a press for a drag: nothing
      if (this.navigator?.axes !== 'rows') b?.cancel();
      if (d.hit.part === 'outside') this.clickAt(d.hit, d.x0, d.y0);
    } else if (!d.rows || (this.navigator?.axes === 'both' && d.hit.part === 'body')) b?.release();
    this.changed();
  }

  private clickAt(h: NavigatorHit, x: number, y: number): void {
    const b = this.behavior;
    if (!b) return;
    if (h.axis === 'rows' || this.navigator?.axes === 'rows') {
      const track = this.rowsTrack(h.section);
      const w = this.lastFrame?.windows.find((v) => v.section === h.section);
      const row = this.renderer?.rowAt(y);
      if (track && w && row !== undefined)
        this.viewer.view.jump(() => {
          const before = this.viewer.view.band(track)?.bandTop;
          this.viewer.view.setBand(track, { bandTop: Math.max(0, row - (w.r1 - w.r0) / 2) });
          return this.viewer.view.band(track)?.bandTop !== before;
        });
      if (this.navigator?.axes === 'rows') return;
    }
    const at = this.renderer?.sectionUnitAt(h.section, x);
    if (at !== undefined) b.click(h.section, at);
  }

  /** A double-click at a point: fit, or fitted, a window `navigator.zoomPx` navigator px wide there. */
  doubleClick(x: number, _y: number): boolean {
    const b = this.behavior;
    const nav = this.navigator;
    const h = this.hit(x, _y);
    if (!b || !nav || !h || nav.axes === 'rows' || !(this.moves[h.section] ?? 0)) return false;
    const r = this.renderer!;
    const at = r.sectionUnitAt(h.section, x);
    const half = this.viewer.config.navigator.zoomPx / 2;
    const a = r.sectionUnitAt(h.section, x - half);
    const c = r.sectionUnitAt(h.section, x + half);
    if (at === undefined || a === undefined || c === undefined) return false;
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
    return true;
  }

  /** Esc during a drag or a sweep: cancelled, the view back where the press found it. */
  cancel(): boolean {
    const d = this.drag;
    if (!d) return false;
    this.drag = undefined;
    if (d.rows) this.viewer.view.setBand(d.rows.track, { bandTop: d.rows.top });
    this.behavior?.cancel();
    this.changed();
    return true;
  }

  /**
   * A key while the part has the focus (its window slider, roving): arrows pan (Up and Down on rows), Shift+arrows
   * resize about the centre, Home and End go to the ends (jumps), Enter centres on the readout's position, Esc cancels;
   * the view's keys act on the focused window's section. True when taken.
   */
  key(key: string, modifiers: { readonly shiftKey?: boolean; readonly altKey?: boolean } = {}): boolean {
    const nav = this.navigator;
    const section = nav?.of[this.focused] ?? nav?.of[0];
    if (!nav || !section) return false;
    if (key === 'Escape') return this.cancel();
    if (!(this.moves[section] ?? 0)) return false;
    const view = this.viewer.view;
    const b = this.behavior;
    const run = (id: string) => this.viewer.runCommand(id, { section }) !== '';
    const v = view.get(section);
    if (!v || !b) return false;
    if (nav.axes === 'rows' && (key === 'ArrowUp' || key === 'ArrowDown')) {
      const track = this.rowsTrack(section);
      const w = this.lastFrame?.windows[this.focused];
      if (!track || !w) return false;
      const step = Math.max(1, (w.r1 - w.r0) * this.viewer.config.navigate.keyPanFraction);
      view.setBand(track, {
        bandTop: Math.max(0, (view.band(track)?.bandTop ?? 0) + (key === 'ArrowUp' ? -step : step)),
      });
      return true;
    }
    switch (key) {
      case 'ArrowLeft':
      case 'ArrowRight': {
        const back = key === 'ArrowLeft';
        if (modifiers.altKey) return back ? view.back() : view.forward();
        if (modifiers.shiftKey) return run(back ? 'window.shrink' : 'window.grow');
        return run(back ? 'window.pan-back' : 'window.pan-on');
      }
      case 'Home':
      case 'End': {
        const w = v.end - v.start;
        const total = this.total(section);
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
      case 'Tab':
        return false;
      default:
        return false;
    }
  }

  private total(section: string): number | undefined {
    const r = this.viewer.view.resolved.value?.sections.find((s) => s.address === section)?.axis.ranges;
    if (!r) return undefined;
    let t = 0;
    for (let k = 0; k + 1 < r.length; k += 2) t += r[k + 1] - r[k];
    return t;
  }

  /** Moves the slider focus to the next or previous window (roving: one tab stop for the part). */
  focusNext(step = 1): void {
    const n = this.navigator?.of.length ?? 0;
    if (n) this.focused = (((this.focused + step) % n) + n) % n;
  }
}
