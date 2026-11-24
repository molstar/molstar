/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ReactiveModel, type ReadonlyState, State } from '@molstar/core/reactive';
import { type Issue, issueOf, type Node, type Value } from '@molstar/core/repr/tree';
import type { Want } from '@molstar/seqx-core/data';
import { NavWindowFlags, overviewWidth, type WindowPart } from '@molstar/seqx-core/render';
import {
  type NavigatorFrame,
  type NavigatorFrameInput,
  type NavigatorHit,
  NavigatorRenderer,
  navigatorBox,
  type NavigatorWindow,
  windowText,
} from '@molstar/seqx-core/render/navigator';
import type { SquareNavigator } from '@molstar/seqx-core/render/navigator-square';
import { isSquareNavigator } from '@molstar/seqx-core/render/navigator-square-kind';
import {
  axisTotal,
  type ResolvedNavigator,
  type ResolvedView,
  resolveNavigatorNode,
  ViewIssues,
} from '@molstar/seqx-core/view';
import type {
  NavigateBehavior,
  NavigatorBehavior,
  RowsBehavior,
  RowsRailGesture,
} from '@molstar/seqx-viewer/behaviors';
// the module, not the barrel: the behaviors' providers make this model (an import cycle through the barrel)
import { windowIdle } from '@molstar/seqx-viewer/behaviors/window';
import type { ModifierState } from '@molstar/seqx-viewer/bindings';
import type { ViewSnapshot } from '@molstar/seqx-viewer/models';
import type { StaticModel } from '@molstar/seqx-core/resolve';
import { blockOfUnit, blockSwitchRange } from './block-switch.js';
import { DragSlop } from './input.js';
import { setNarrowed } from './narrowed.js';
import { type PartCursor, partCursor } from './part-cursor.js';
import type { SquareSetGesture } from './square-set.js';
import type { SeqViewer } from './viewer.js';

// The navigator's model (_spec/README.md "Models", NavigatorModel): DOM-free, one a navigator a viewer shows (a
// declared node of the view by its id, or a host's definition, `viewer.navigator(def)`). It holds the navigator's
// picture (its own renderer over the viewer's store, theme, image cache and held domains), the bound windows' frames as
// the viewer's last frame placed them, the lit part, idle per window and axis, the gesture in progress, the focused
// slider and the readout; the canvas view and an export read it. Its pointer and keys are in its own px, carried to the
// bound section's display units through the picture, and run through the `navigator` behavior. A hover stays here (a
// readout, never the viewer's hover).

/**
 * The readout's line, CSS px (the `msui-text-xs` leading): the first line of every navigator part, above its picture,
 * right-aligned, in the same place whatever the window or the pointer does (_spec/README.md "Composition", Navigators,
 * The readout). A side part's picture is the form's box less it, so the part is the view's height.
 */
export const ReadoutLinePx = 16;

/** The square's chunk (`square-chunk.ts`: core's square, the set's gestures). */
type SquareModule = typeof import('./square-chunk.js');

/**
 * The square's door: how its chunk loads, set by the navigators part (`navigators.ts`), never imported here (a dynamic
 * import is a chunk of every bundle that reaches its module, the minimal viewer's too). Why lazy: a viewer whose views
 * have no square pays nothing for one (_spec/README.md "Composition", Navigators, The square).
 */
export const SquareDoor: { load?: () => Promise<SquareModule> } = {};

/** The square's chunk once landed. */
let Square: SquareModule | undefined;
let squareLoading: Promise<void> | undefined;

/**
 * Loads the square's chunk (once; a failed load is tried again by the next ask): a view showing a square waits for it
 * (resolve's lazy kinds, the navigators part's `loaded`), a host's square draws once it landed.
 */
export function loadSquare(): Promise<void> {
  const load = SquareDoor.load;
  if (!load) return Promise.reject(new Error('No square here: the navigators part is not composed'));
  squareLoading ??= load().then(
    (m) => {
      Square = m;
    },
    (e) => {
      squareLoading = undefined;
      throw e;
    },
  );
  return squareLoading;
}

/** Whether the square's chunk landed. */
export const squareLoaded = (): boolean => Square !== undefined;

/** A square's renderer (a type guard: the class is the chunk's, never imported here). */
const isSquare = (r: NavigatorRenderer | SquareNavigator | undefined): r is SquareNavigator =>
  (r as { square?: boolean } | undefined)?.square === true;

/** A host's navigator: the params of a `navigator` node with its id, not saved in the spec. */
export interface NavigatorDef {
  readonly id: string;
  readonly of: string | readonly string[];
  readonly axes?: 'positions' | 'rows' | 'both';
  /** Its picture: a built-in's name or a registered `<id>:<name>` (the axis when the composition lacks it). */
  readonly content?: string;
  readonly tracks?: readonly string[];
  readonly bands?: string;
  readonly system?: string;
  readonly place?: 'top' | 'bottom' | 'left' | 'right';
  readonly size?: number | readonly [number, number];
  /** Its picture's own params, plain (a host's or an extension's picture reads them as written in a spec). */
  readonly params?: Readonly<Record<string, Value>>;
}

/** A readout in parts: the name (a system's, absent when the picture says it), the value, the rest (` · row 47`). */
export interface ReadoutParts {
  readonly name: string;
  readonly value: string;
  readonly rest: string;
}

/** A window as its slider says it (its `aria-valuetext`): set when the view settles, never a frame. */
export interface NavigatorValue {
  readonly section: string;
  readonly text: string;
  readonly disabled: boolean;
}

/** The pointer's state for a gesture: touch widens the hits. */
/** A click's modifiers: Cmd (Mac) or Ctrl adds or removes on a square over a set. */
export interface ClickModifiers {
  readonly metaKey?: boolean;
  readonly ctrlKey?: boolean;
}

/** A key's modifiers. */
export interface KeyModifiers extends ClickModifiers {
  readonly shiftKey?: boolean;
  readonly altKey?: boolean;
}

/** No pointer state (a mouse): `move`'s default, one object (a move allocates nothing). */
const NoPointer: NavigatorPointer = Object.freeze({});

/** What a drag over a set that grips no box says: only a run's box drags (a scattered set changes by its cells). */
const ScatteredHint: ReadoutParts = {
  name: '',
  value: 'Only a run’s box drags',
  rest: ' · click cells, Cmd or Ctrl to add, Home for all',
};

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
  /** The press is on a rows grip: the rows rail's gesture zooms them about the other end (`RowsRailGesture`). */
  readonly grip?: RowsRailGesture;
  /** The press is on a square over a set: a drag of its run's box (`SquareSetGesture`), else a click's. */
  readonly set?: true;
}

export class NavigatorModel extends ReactiveModel {
  /** Bumped when what the navigator draws changes (hover, a gesture, its size): its view redraws. */
  private readonly _version = new State(0);
  readonly version: ReadonlyState<number> = this._version.asReadonly;
  /** The readout under the pointer (a position, a row, a band's name), or none. */
  private readonly _readout = new State<string | undefined>(undefined);
  readonly readout: ReadonlyState<string | undefined> = this._readout.asReadonly;
  private readonly _readoutParts = new State<ReadoutParts | undefined>(undefined);
  /** The readout in parts, for a host that has little room: the name gives way, the value never. */
  readonly readoutParts: ReadonlyState<ReadoutParts | undefined> = this._readoutParts.asReadonly;
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
  /** The width it shares with the canvas beside it (`setRoom`), 0 unknown. */
  private room = 0;
  /** Its renderer: a picture's, or the square's (two axes over blocks × blocks). */
  private renderer: NavigatorRenderer | SquareNavigator | undefined;
  private rendererFor: unknown[] = [];
  private readonly wanter;
  private lastFrame: NavigatorFrame | undefined;
  private lit: { section: string; part: WindowPart; axis: 'positions' | 'rows' } | undefined;
  private pointerX = Number.NaN;
  private pointerY = Number.NaN;
  /** The section the pointer last acted on: its window on top. */
  private top: string | undefined;
  private readonly moves: Record<string, number> = {};
  private press: Press | undefined;
  /** The square over a set: its box's drag, its clicks and keys (each change one `section.include`), the hover's cursor. */
  private _setGesture: SquareSetGesture | undefined;
  /** Made on the first square over a set (the chunk landed: a set's square is the chunk's). */
  private get setGesture(): SquareSetGesture {
    this._setGesture ??= new Square!.SquareSetGesture((include) => this.writeSet(include));
    return this._setGesture;
  }
  private setHover: PartCursor = 'default';
  /** A key drove the square over a set since the last press or focus: its keyboard's cell is drawn. */
  private keyed = false;
  /** The press over a set that grips no box said so (`ScatteredHint`). */
  private hinted = false;
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
  /** A host's picture the composition lacks, said as a declared navigator's is (the view's code). */
  private pictureIssue: Issue | undefined;

  private hostNavigator(view: ResolvedView): ResolvedNavigator {
    if (this.hosted?.[0] === view) return this.hosted[1];
    const { id, params: own, ...params } = this.def!;
    const node: Node = { kind: 'navigator', id, params: { ...own, ...params } as never };
    const n = resolveNavigatorNode(
      this.viewer.spec.model.value,
      view.sections,
      node,
      `${view.address}/navigators/${id}`,
      this.viewer.context.kinds,
    );
    this.hosted = [view, n];
    const content = this.def!.content;
    this.pictureIssue =
      content !== undefined && content !== n.content
        ? issueOf(ViewIssues, 'seq.view.navigator.picture', `No picture \`${content}\` here: drawn as the axis`, {
            address: n.address,
            hint: ViewIssues.codes['seq.view.navigator.picture'].hint,
          })
        : undefined;
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

  /**
   * Sets the width it shares with the canvas beside it (the default `ui.navigators`' row, less the other side parts):
   * a vertical overview of the form's size narrows and hides in a narrow frame by the plot it leaves (core's "Render",
   * A narrow frame: `overviewWidth`).
   */
  setRoom(width: number): void {
    if (width === this.room) return;
    this.room = Math.max(0, width);
    this.changed();
  }

  /** The box the picture is drawn in: the host's and the form's (a zoom strip draws below it). */
  private boxHeight = 0;

  /**
   * Its size from the host's and the form's box, and its zoom strip's height (one navigator: the pair's): the form's
   * box and the strip below it; a height the host fixes holds both (the strip out of it, the context the rest).
   */
  private place(): void {
    const nav = this.navigator;
    if (!nav) return;
    const b = navigatorBox(nav, this.hostWidth || this.viewer.width, this.viewer.height, this.viewer.context.pictures);
    const zoom = this.renderer?.zoomHeight ?? 0;
    this.width = this.hostWidth || b.width;
    // a vertical overview at its form's size beside the canvas gives way to the plot (A narrow frame)
    const narrow = this.viewer.frame?.layout.narrow;
    if (!this.hostWidth && this.room > 0 && narrow && nav.axes === 'rows' && nav.size === undefined)
      if (nav.place === 'left' || nav.place === 'right')
        this.width = overviewWidth(this.room - narrow.columnsPx, b.width, this.viewer.context.renderConfig);
    // hidden for want of width: as a hidden navigator, its rows rail back under `unlessOverview` (the next frame's)
    if (setNarrowed(this.viewer, this.id, b.width > 0 && this.width === 0)) this.viewer.schedule();
    // a side part's picture gives the readout's line its room: the part is the form's box, the view's height
    const side = nav.place === 'left' || nav.place === 'right' ? ReadoutLinePx : 0;
    this.boxHeight = this.hostHeight ? Math.max(1, this.hostHeight - zoom) : Math.max(1, b.height - side);
    this.height = this.hostHeight || this.boxHeight + zoom;
  }

  private changed(): void {
    this._version.set(this._version.value + 1);
  }

  private rendererOf(nav: ResolvedNavigator): NavigatorRenderer | SquareNavigator | undefined {
    const model = this.viewer.spec.model.value;
    const view = this.viewer.view.resolved.value!;
    const viewRenderer = this.viewer.renderer;
    const f = this.rendererFor;
    if (!this.renderer || f[0] !== model || f[1] !== view || f[2] !== nav || f[3] !== viewRenderer) {
      const square = isSquareNavigator(nav, view);
      // a square before its chunk landed: nothing drawn, drawn when it lands
      if (square && !Square) {
        loadSquare().then(
          () => {
            this.changed();
            this.viewer.schedule();
          },
          () => undefined,
        );
        return undefined;
      }
      this.renderer = square
        ? new Square!.SquareNavigator({
            model,
            view,
            navigator: nav,
            tokens: this.viewer.tokens,
            config: this.viewer.context.renderConfig,
            // a picture's cells over the lattice: read from the viewer's store with the composition's types
            data: this.viewer.store.store,
            text: this.viewer.context.text,
            types: this.viewer.context.types,
            layouts: this.viewer.context.layouts,
            pictures: this.viewer.context.pictures,
            kinds: this.viewer.context.kinds,
            // the square it replaces: its cells kept when made of the same (an edit elsewhere, the ring's `include`)
            ...(isSquare(this.renderer) ? { previous: this.renderer } : {}),
          })
        : this.pictureRenderer(nav, model, view);
      this.rendererFor = [model, view, nav, viewRenderer];
    }
    return this.renderer;
  }

  /** A picture's renderer over the viewer's store, theme, image cache and held domains. */
  private pictureRenderer(nav: ResolvedNavigator, model: StaticModel, view: ResolvedView): NavigatorRenderer {
    const viewRenderer = this.viewer.renderer;
    return new NavigatorRenderer({
      model,
      view,
      navigator: nav,
      data: this.viewer.store.store,
      tokens: this.viewer.tokens,
      text: this.viewer.context.text,
      config: this.viewer.context.renderConfig,
      types: this.viewer.context.types,
      layouts: this.viewer.context.layouts,
      // the composition's pictures and layer kinds: a picture's track of an extension's type drawn as that type
      pictures: this.viewer.context.pictures,
      kinds: this.viewer.context.kinds,
      sizings: this.viewer.context.sizings,
      ...(this.viewer.context.trackSizing ? { trackSizing: this.viewer.context.trackSizing } : {}),
      direction: this.viewer.config.navigator.direction,
      zoomStrip: this.viewer.config.navigator.zoomStrip,
      storeSettings: this.viewer.store.store.settings,
      // the view's image budget and held domains: a value has one color in both
      ...(viewRenderer ? { images: viewRenderer.images, domains: viewRenderer.layoutState } : {}),
    });
  }

  /** What a gesture can move: each bound window's axes (`NavWindowFlags`), 0 idle. */
  private updateMoves(nav: ResolvedNavigator): void {
    const behavior = this.viewer.behaviors.has('navigator');
    const view = this.viewer.view;
    // (loops, not `find`: a frame allocates no closure)
    for (let k = 0; k < nav.of.length; k++) {
      const a = nav.of[k];
      let m = 0;
      if (behavior && view.canMove(a)) {
        if (nav.axes !== 'rows' && !windowIdle(view, a, this.viewer.config)) m |= NavWindowFlags.x;
        const ws = nav.axes !== 'positions' ? this.lastFrame?.windows : undefined;
        for (let i = 0; ws && i < ws.length; i++) {
          const w = ws[i];
          if (w.section !== a) continue;
          if (w.rows > w.r1 - w.r0 + 1e-9) m |= NavWindowFlags.y;
          break;
        }
      }
      if (this.moves[a] !== m) this.moves[a] = m;
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
    if (!r) return undefined;
    this.updateMoves(nav);
    let f = this.drawFrame(r);
    let rebuilt = f.rebuilt;
    // a zoom strip came or went in this frame: the part's size with it (a height the host fixed drawn again at once:
    // the context gives the strip its room)
    this.place();
    if (f.lists && f.lists.height !== this.height) {
      f = this.drawFrame(r);
      rebuilt ||= f.rebuilt;
      this.place();
    }
    if (rebuilt) this.wanter.want(f.wants as Want[]);
    this.lastFrame = f;
    this.publishReadout();
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

  /** A frame of the picture at the box's size, over the viewer's last frame. */
  private drawFrame(r: NavigatorRenderer | SquareNavigator): NavigatorFrame {
    const layout = this.viewer.frame!.layout;
    const sweep = this.behavior?.sweep.value;
    // over a set the keyboard's cell is drawn while the part has the focus and the keys drive it (a pointer's press
    // hides it: a cursor no key moved is no part of the window)
    if (isSquare(r) && r.setMode) {
      const g = this.setGesture;
      const keyed = this.hasFocus && this.keyed;
      if (keyed) g.ensureKeyCell(r);
      r.setKeyCell(keyed ? g.keyX : -1, keyed ? g.keyY : -1);
    }
    // one input, its fields set a frame (a drag's frame allocates nothing; the renderers read its values, never keep it)
    const input = this.frameInput;
    input.size.width = this.width;
    input.size.height = this.boxHeight;
    input.dpr = this.viewer.dpr;
    input.settled = this.viewer.store.store.settled.value;
    input.overlay = this.viewer.store.overlay(
      this.viewer.spec.model.value,
      this.viewer.view.resolved.value!,
      this.viewer.context,
    );
    input.view = layout;
    input.moves = this.moves;
    input.lit = this.lit;
    input.top = this.top;
    input.pointerX = this.pointerX;
    input.pointerY = this.pointerY;
    input.emphasis = this.viewer.emphasis;
    input.cvdSafe = this.viewer.legends.cvdSafe.value;
    input.sweep = sweep && this.press ? this.sweepPx(sweep) : undefined;
    return r.frame(input as NavigatorFrameInput);
  }

  /** The frame's input, reused (`drawFrame`): an absent part undefined. */
  private readonly frameInput: { -readonly [K in keyof NavigatorFrameInput]: NavigatorFrameInput[K] | undefined } & {
    size: { width: number; height: number };
  } = {
    size: { width: 0, height: 0 },
    dpr: 1,
    settled: false,
    overlay: undefined,
    view: undefined,
    moves: undefined,
    lit: undefined,
    top: undefined,
    pointerX: Number.NaN,
    pointerY: Number.NaN,
    emphasis: undefined,
    cvdSafe: false,
    sweep: undefined,
  };

  /** The navigator's last frame. */
  get frame(): NavigatorFrame | undefined {
    return this.lastFrame;
  }

  /** The issues the navigator met (its bindings, what it left out, its picture's). */
  get issues(): readonly Issue[] {
    const own = this.renderer?.allIssues ?? [];
    return this.pictureIssue ? [this.pictureIssue, ...own] : own;
  }

  /**
   * x on the navigator of positions `[start, end)` of its system, `[x0, x1]` in its display order (a reversed axis's,
   * collapsed introns' nearest), after a `draw`; undefined before, or where none is shown. What a host drawing its own
   * navigator places its marks with.
   */
  xRange(start: number, end: number): readonly [number, number] | undefined {
    return this.renderer?.xRangeOf(start, end);
  }

  /** The position of the navigator's system at x (its own px), clamped to its picture; undefined before a `draw`. */
  positionAt(x: number): number | undefined {
    return this.renderer?.positionAtX(x);
  }

  /** The part took or lost the keyboard focus (its host's element): losing it, its sliders say where the view is. */
  setFocus(focused: boolean): void {
    if (this.hasFocus === focused) return;
    this.hasFocus = focused;
    this.keyed = false;
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
    return this.viewer.behavior<NavigatorBehavior>('navigator');
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
    // over a set the parts are the run's box's, whatever the axis's window does
    if (!h || h.context || this.setSquare) return h;
    const m = this.moves[h.section] ?? 0;
    const along = h.axis === 'rows' ? NavWindowFlags.y : NavWindowFlags.x;
    if (!m) return { ...h, part: 'outside' };
    if (h.part !== 'outside' && h.part !== 'body' && !(m & along)) return { ...h, part: 'body' };
    return h;
  }

  /** The cursor over the navigator, by the part and the gesture (`partCursor`), held through a drag. */
  get cursor(): PartCursor {
    // a press over a set that grips no box (a scattered set's cell, a label) drags nothing: said, and the cursor says so
    if (this.press?.set)
      return this.setGesture.dragging ? this.setGesture.cursor : this.press.moved ? 'not-allowed' : 'pointer';
    if (this.press) return partCursor(this.press.hit.part, this.press.hit.axis, true);
    if (this.setSquare) return this.setHover;
    const l = this.lit;
    return l ? partCursor(l.part, l.axis, false) : 'default';
  }

  /** The square over a set (one section on both axes with a switch), else undefined. */
  private get setSquare(): SquareNavigator | undefined {
    const r = this.renderer;
    return isSquare(r) && r.setMode ? r : undefined;
  }

  /** The pointer moved over the navigator: the lit part and the readout follow; a drag moves its window. */
  move(x: number, y: number, p: NavigatorPointer = NoPointer): void {
    this.pointerX = x;
    this.pointerY = y;
    const d = this.press;
    if (d) {
      if (!d.moved && Math.max(Math.abs(x - d.x0), Math.abs(y - d.y0)) >= DragSlop) d.moved = true;
      const set = d.set ? this.setSquare : undefined;
      if (d.moved && set) {
        set.pointer[0] = x;
        set.pointer[1] = y;
        this.setGesture.drag(set);
        if (!this.setGesture.dragging && !this.hinted) {
          // once a press: a scattered set (or a label) has no box to drag
          this.hinted = true;
          this.say(ScatteredHint);
        }
      } else if (d.moved) this.dragTo(d, x, y);
      this.changed();
      return;
    }
    const set = this.setSquare;
    if (set) {
      const g = set.gripAt(x, y, p.touch);
      this.setHover = g
        ? Square!.setCursor(g.grip, g.axis)
        : set.cellAt(x, y) || set.labelAt(x, y) >= 0
          ? 'pointer'
          : 'default';
      this.say(Square!.setReadout(set, x, y));
      this.changed();
      return;
    }
    const h = this.hit(x, y, p.touch);
    // the context of a zoom strip is no target: nothing lit, the default cursor
    this.lit =
      h && !h.context && (this.moves[h.section] ?? 0) ? { section: h.section, part: h.part, axis: h.axis } : undefined;
    this.say(h ? this.readoutOf(h) : undefined);
    this.changed();
  }

  /**
   * Over a window, what it shows (its positions, its rows); over a zoom strip's mark, the region it marks; elsewhere
   * the position, row and band under the pointer.
   */
  private readoutOf(h: NavigatorHit): ReadoutParts {
    if (h.context === 'mark') return { name: '', value: this.renderer?.markText() ?? '', rest: '' };
    const nav = this.navigator;
    const w = h.part !== 'outside' ? this.windowOf(h.section) : undefined;
    if (nav && w?.shown) return { name: '', value: windowText(this.viewer.spec.model.value, nav, w), rest: '' };
    const systems = this.viewer.spec.model.value.systems;
    const name = (nav?.system && systems.get(nav.system)?.name) ?? '';
    const rest: string[] = [];
    let value = '';
    let named = '';
    if (h.position !== undefined && h.system !== undefined) {
      // over blocks: the block's system and its position there (`chr17: 7,668,402`)
      value = (Math.floor(h.position) + 1).toLocaleString('en-US');
      named = systems.get(h.system)?.name ?? '';
    } else if (h.position !== undefined) {
      const own = this.renderer?.readout(h.position);
      value = own ?? (h.position + 1).toLocaleString('en-US');
      if (own === undefined) named = name;
    }
    if (h.row !== undefined) rest.push(`row ${(Math.floor(h.row) + 1).toLocaleString('en-US')}`);
    if (h.band !== undefined) rest.push(h.band);
    // no position: the first of the rest is the value
    if (h.position === undefined && rest.length) value = rest.shift() as string;
    return { name: named, value, rest: rest.map((r) => ` \u00b7 ${r}`).join('') };
  }

  /** What the pointer is over (`say`), and whether the readout follows the window (it moved since the pointer did). */
  private hover: ReadoutParts | undefined;
  private following = false;
  /** The window the readout last read (its ends, as numbers: compared each frame, its text made when they move). */
  private readParts: ReadoutParts | undefined;
  private readP0 = Number.NaN;
  private readP1 = Number.NaN;
  private readR0 = Number.NaN;
  private readR1 = Number.NaN;
  private readS0 = '';
  private readS1 = '';

  /** The pointer's readout (a hover's, a press on the context's): shown unless the window moved since. */
  private say(parts: ReadoutParts | undefined): void {
    this.hover = parts;
    this.following = false;
    this.publishReadout();
  }

  /**
   * The readout is one line in one place (_spec/README.md "Composition", Navigators, The readout) and says what is
   * drawn: while a window moves (a drag, a wheel, a key, a jump, a move made elsewhere) the window's region as this
   * frame draws it, until the pointer moves again; else what is under the pointer, else, at rest, the window's region
   * (the one last acted on). Called with each frame and with each hover.
   */
  private publishReadout(): void {
    const nav = this.navigator;
    const f = this.lastFrame;
    // the window last acted on, else the first shown (an indexed loop: no closure a frame)
    let w: NavigatorWindow | undefined;
    if (nav && f) {
      const ws = f.windows;
      for (let i = 0; i < ws.length && !w; i++) if (ws[i].section === this.top) w = ws[i];
      for (let i = 0; i < ws.length && !w; i++) if (ws[i].shown) w = ws[i];
    }
    let own: ReadoutParts | undefined;
    if (nav && w?.shown) {
      if (
        w.p0 !== this.readP0 ||
        w.p1 !== this.readP1 ||
        w.r0 !== this.readR0 ||
        w.r1 !== this.readR1 ||
        w.s0 !== this.readS0 ||
        w.s1 !== this.readS1 ||
        !this.readParts
      ) {
        // a window that moved by itself (not just drawn again) takes the readout from the pointer
        if (this.readParts) this.following = true;
        this.readP0 = w.p0;
        this.readP1 = w.p1;
        this.readR0 = w.r0;
        this.readR1 = w.r1;
        this.readS0 = w.s0;
        this.readS1 = w.s1;
        this.readParts = this.windowParts(nav, w);
      }
      own = this.readParts;
    } else this.readParts = undefined;
    if (this.press?.moved) this.following = true;
    // a press that grips no box says its hint while it runs (no window moves)
    const parts = this.press && this.hinted ? this.hover : this.following ? (own ?? this.hover) : (this.hover ?? own);
    const prev = this._readoutParts.value;
    if (parts === prev) return;
    if (prev && parts && prev.name === parts.name && prev.value === parts.value && prev.rest === parts.rest) return;
    this._readoutParts.set(parts);
    this._readout.set(parts && `${parts.name ? `${parts.name}: ` : ''}${parts.value}${parts.rest}`);
  }

  /**
   * A window as the readout says it: its region, short where the part is narrow. A rows-only or vertical overview: the
   * rows alone (`401–440 / 10,000`); the miniature (two axes): the positions then the rows (`1–300 · 1–12 / 10,000`);
   * the square over blocks (two axes, a block's range): the block and its range without the assembly (`chr3: 1–181,000,000`).
   */
  private windowParts(nav: ResolvedNavigator, w: NavigatorWindow): ReadoutParts {
    const n = (x: number) => Math.round(x).toLocaleString('en-US');
    const rows = `${n(w.r0 + 1)}–${n(w.r1)} / ${n(w.rows)}`;
    const text = windowText(this.viewer.spec.model.value, nav, w);
    const side = nav.place === 'left' || nav.place === 'right';
    let value = text;
    if (nav.axes === 'rows' || (nav.axes === 'both' && side && !nav.blocks)) value = rows;
    else if (nav.axes === 'both' && nav.blocks) value = text.replace(/ of [^,]*$/, '');
    else if (nav.axes === 'both') value = `${text.replace(/^[^:]*:\s?/, '').replace(/ of [^,]*,.*$/, '')} · ${rows}`;
    return { name: '', value, rest: '' };
  }

  /** The pointer left: no lit part, no readout, no hairline. */
  leave(): void {
    this.pointerX = this.pointerY = Number.NaN;
    this.lit = undefined;
    this.say(undefined);
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
    // over a set: a press on the run's box drags it, on a cell or a label it is a click's
    const set = this.setSquare;
    if (set && b && nav) {
      const label = set.labelAt(x, y);
      if (!h && label < 0) return false;
      const hit: NavigatorHit = h ?? { section: nav.of[0], part: 'outside', axis: 'positions' };
      this.setGesture.press(set, x, y, p.touch);
      this.keyed = false;
      this.hinted = false;
      this.press = { x0: x, y0: y, hit, moved: false, positions: false, set: true };
      this.changed();
      return true;
    }
    // over blocks the context is the genome: a press there is a click's (a jump to the block, a switch)
    if (h?.context && nav?.blocks && b) {
      this.press = { x0: x, y0: y, hit: h, moved: false, positions: false };
      return true;
    }
    // a press on a zoom strip's context is no gesture: it names what is under it, as a hover (a tap's readout)
    if (h?.context) {
      this.say(this.readoutOf(h));
      this.changed();
      return false;
    }
    // a switch takes a click on any block outside the window, idle or not (a ring's window is idle unless followed)
    if (b && nav && h && nav.switch === 'include' && nav.blocks && (h.part === 'outside' || !this.moves[h.section])) {
      this.press = { x0: x, y0: y, hit: { ...h, part: 'outside' }, moved: false, positions: false };
      return true;
    }
    if (!b || !nav || !h || !(this.moves[h.section] ?? 0)) return false;
    this.top = h.section;
    const i = nav.of.indexOf(h.section);
    if (i >= 0 && i !== this._focused.value) this._focused.set(i);
    const m = this.moves[h.section] ?? 0;
    const inside = h.part !== 'outside';
    // what the press drags: a rows part its rows, a positions part its positions, a miniature's body both; outside the
    // window a positions press (a click or a sweep), and on rows a press for its click
    const both = nav.axes === 'both' && h.part === 'body';
    // a rows grip (the window's top or bottom bar) zooms the rows about the other end, as the rows rail's grip does
    const grip = h.axis === 'rows' && (h.part === 'start' || h.part === 'end') && m & NavWindowFlags.y;
    const g = grip ? this.pressGrip(h, y) : undefined;
    if (g) {
      this.press = { x0: x, y0: y, hit: h, moved: false, positions: false, grip: g };
      this.changed();
      return true;
    }
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

  /** A press on a rows grip at `y`: the rows behavior's rail gesture over the bound band (none without it). */
  private pressGrip(h: NavigatorHit, y: number): RowsRailGesture | undefined {
    const rowsBehavior = this.viewer.behavior<RowsBehavior>('rows');
    const track = this.rowsTrack(h.section);
    const w = this.windowOf(h.section);
    const l = track && rowsBehavior?.lay(track);
    const at = this.renderer?.rowAt(y);
    if (!rowsBehavior || !track || !w || !l || at === undefined) return undefined;
    const g = rowsBehavior.railGesture();
    const window = { start: w.r0, end: w.r1, total: w.rows, width: Math.max(0, l.box.height - l.pinnedHeight) };
    return g.press(track, h.part, at, window) ? g : undefined;
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
    if (d.grip) {
      const row = r.rowAt(y);
      if (row !== undefined) d.grip.dragTo(row);
      return;
    }
    if (d.rows && d.hit.part !== 'outside') {
      const row = r.rowAt(y);
      const row0 = r.rowAt(d.y0);
      if (row !== undefined && row0 !== undefined)
        this.viewer.view.setBand(d.rows.track, { bandTop: this.clampTop(d.hit.section, d.rows.top + row - row0) });
    }
    if (!d.positions) return;
    // the square's sweep (from a cell off the window): whole cells, the drawn rectangle and the release's zoom alike
    if (isSquare(r) && d.hit.part === 'outside') {
      const snapped = r.snapSweep(d.hit.section, d.x0, x);
      if (snapped) {
        // the end the pointer is at is the one away from the press's cell
        const forward = x >= d.x0;
        this.behavior?.dragTo(forward ? snapped[1] : snapped[0], forward ? snapped[0] : snapped[1]);
        return;
      }
    }
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
    const set = d.set && this.setSquare;
    if (set) {
      // a drag of the run's box writes its run once (one patch, one undo step); a press that never moved is a click's
      if (d.moved) this.setGesture.release(set);
      else this.setGesture.cancel(set);
    } else if (d.grip) {
      if (d.moved) d.grip.release();
      else d.grip.cancel();
      if (d.moved) this.ownMove = true;
    } else if (!d.moved) {
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

  /** A square's click chain (its double-click's clicks): the view before the first, its history's count, the last click's time. */
  private chain: { before: ViewSnapshot; back0: number; t: number } | undefined;

  /** The last press that ended: its click (the platform's, after the button went up) reads it. */
  private lastPress: Press | undefined;

  /**
   * A click (the platform's, after a press that never moved; `detail` its count): outside the window it centres
   * there, both axes on a miniature (a jump, its own step); on the window nothing. A second click is a
   * double-click's: nothing (`doubleClick` follows).
   */
  click(x: number, y: number, detail = 1, modifiers: ClickModifiers = {}): boolean {
    const d = this.lastPress;
    this.lastPress = undefined;
    const b = this.behavior;
    const nav = this.navigator;
    const square = isSquare(this.renderer) ? this.renderer : undefined;
    // a square with a switch: the clicks of one double-click are one gesture with its `dblclick` (it replaces their steps)
    if (square && nav?.switch === 'include') {
      const view = this.viewer.view;
      if (detail === 1 || !this.chain)
        this.chain = { before: view.snapshot({ includes: true }), back0: view.history.value.back, t: 0 };
      this.chain.t = performance.now();
    }
    // the square with a switch: every click on a cell shows its pair at once, a quick second one too (a platform's
    // double-click count never turns a click into a select first: _spec/README.md, The square over a set)
    if (!d || d.moved || !b || !nav || (detail > 1 && !(square && nav.switch === 'include'))) return false;
    // over a set: a cell its pair, a label its block; Cmd or Ctrl adds or removes them
    if (square?.setMode) {
      const toggle = !!(modifiers.metaKey || modifiers.ctrlKey);
      const label = square.labelAt(x, y);
      if (label >= 0) return this.setGesture.label(square, label, toggle);
      const cell = square.cellAt(x, y);
      return !!cell && this.setGesture.cell(square, cell.x, cell.y, toggle);
    }
    // `navigator.switch: include`: a click on a block shows it alone on the bound section (on what it already shows
    // alone: no change); the window, where it has one, keeps its gestures
    const ctx = !!d.hit.context;
    // the square: a click on a cell shows its two blocks with a switch, else centres the section on the x block
    if (square && nav.switch === 'include') return this.switchToPair(square, d.hit.section, x, y);
    // the square: a click on a cell outside the window centres both axes there (one step; one section on both axes:
    // its window on the x block's place); on the window nothing
    if (square) return d.hit.part === 'outside' && this.centreSquare(square, b, x, y);
    // a click that never moved switches, on the window too (a drag of it keeps its gestures)
    if (nav.switch === 'include' && nav.blocks) return this.switchTo(d.hit.section, x, ctx);
    // over blocks, a click on the context (the genome above a zoom strip) centres there: a jump to that block
    if (ctx) {
      const at = this.renderer?.sectionUnitAt(d.hit.section, x, true);
      this.ownMove = true;
      if (at === undefined) return false;
      return this.switchBlock(b, d.hit.section, at) ?? b.click(d.hit.section, at, detail);
    }
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
          this.viewer.view.jump(
            () => {
              const before = this.viewer.view.band(track)?.bandTop;
              this.viewer.view.setBand(track, { bandTop: this.clampTop(d.hit.section, row - (w.r1 - w.r0) / 2) });
              return this.viewer.view.band(track)?.bandTop !== before;
            },
            true,
            'move',
          ) === 'moved';
      if (nav.axes === 'rows') return moved;
    }
    const at = this.renderer?.sectionUnitAt(d.hit.section, x);
    if (at === undefined) return moved;
    // over blocks a click on another block switches to it by the rule (whole from whole, else the zoom kept)
    const switched = nav.blocks ? this.switchBlock(b, d.hit.section, at) : undefined;
    return (switched ?? b.click(d.hit.section, at, detail)) || moved;
  }

  /**
   * A click at unit `at` over blocks on a block other than the view's (the one holding its centre): a jump by the
   * switching rule (`blockSwitchRange`); undefined on the view's own block (a click centres there).
   */
  private switchBlock(b: NavigatorBehavior, section: string, at: number): boolean | undefined {
    const s = this.viewer.view.resolved.value?.sections.find((x) => x.address === section);
    const blocks = s?.blocks;
    const v = this.viewer.view.get(section);
    if (!blocks || blocks.length < 2 || !v || s.params.arrange === 'stack') return undefined;
    const to = blockOfUnit(blocks, at);
    if (to === blockOfUnit(blocks, (v.start + v.end) / 2)) return undefined;
    const r = blockSwitchRange(blocks, v, to, v.width > 0 ? (v.end - v.start) / v.width : 0);
    this.ownMove = true;
    return b.zoomTo(section, r[0], r[1]);
  }

  /** The square's click: the x section centred on the cell's x place, a second section on its y place (one step). */
  private centreSquare(square: SquareNavigator, b: NavigatorBehavior, x: number, y: number): boolean {
    const [sx, sy] = this.navigator!.of;
    const ux = square.sectionUnitAt(sx, x);
    const uy = sy !== undefined && sy !== sx ? square.sectionUnitAtY(sy, y) : undefined;
    this.ownMove = true;
    // each axis by the switching rule on another block (whole from whole, else the zoom kept), else centred there
    const go = (s: string, u: number, record: boolean) => {
      const v = this.viewer.view.get(s);
      const blocks = this.viewer.view.resolved.value?.sections.find((x) => x.address === s)?.blocks;
      if (v && blocks && blocks.length > 1) {
        const to = blockOfUnit(blocks, u);
        if (to !== blockOfUnit(blocks, (v.start + v.end) / 2)) {
          const r = blockSwitchRange(blocks, v, to, v.width > 0 ? (v.end - v.start) / v.width : 0);
          return b.zoomTo(s, r[0], r[1], record);
        }
      }
      return b.center(s, u, record);
    };
    const movedX = ux !== undefined && go(sx, ux, true);
    const movedY = uy !== undefined && go(sy!, uy, !movedX);
    return movedX || movedY;
  }

  /**
   * A double-click at a point: on a window over zoomed rows, the rows' zoom reset (one jump, the rows rail's
   * `resetZoom`); else fit, or fitted, a window `navigator.zoomPx` navigator px wide there.
   */
  doubleClick(x: number, y: number): boolean {
    const b = this.behavior;
    const nav = this.navigator;
    // a square with a switch: anywhere on it the whole genome (the clicks before it showed pairs; _spec/README.md,
    // The square over a set), one step
    if (nav?.switch === 'include' && isSquare(this.renderer)) {
      // the steps its own clicks pushed (they showed the cell's pair) are replaced by its one: Back returns the window
      // before the gesture
      const c = this.chain;
      this.chain = undefined;
      const fresh = c !== undefined && performance.now() - c.t < 600;
      if (fresh) this.viewer.view.dropSteps(this.viewer.view.history.value.back - c.back0);
      return this.showAllOnSquare(true, fresh ? c.before : undefined);
    }
    const h = this.hit(x, y);
    // over blocks a double-click opens a chromosome from the fit too (its window idle there), on the context too
    const blocks = !!nav?.blocks && !!h && this.viewer.view.canMove(h.section);
    if (!nav || !h || (!blocks && (h.context || !(this.moves[h.section] ?? 0)))) return false;
    const track = nav.axes !== 'positions' && h.part !== 'outside' ? this.rowsTrack(h.section) : undefined;
    if (track && this.viewer.config.rows.zoom.enabled && this.viewer.view.band(track)?.rowScale !== undefined) {
      const reset = !!this.viewer.behavior<RowsBehavior>('rows')?.railGesture().resetZoom(track);
      if (reset) this.ownMove = true;
      return reset;
    }
    if (!b || nav.axes === 'rows') return false;
    const r = this.renderer!;
    // over blocks (`navigator.doubleClick: toggle-fit`, one level deeper): the block under the pointer fitted, and
    // from that fit the whole
    // the square: the cell's two blocks fitted (one section on both axes: their hull), again the whole
    if (isSquare(r)) {
      const pair = r.pairAt(h.section, x, y);
      const v = this.viewer.view.get(h.section);
      if (!pair || !v) return false;
      this.ownMove = true;
      const view = this.viewer.view;
      const shown = Math.abs(v.start - pair.start) <= 1 && Math.abs(v.end - pair.end) <= 1;
      return shown
        ? view.jump(() => view.fit(h.section), true, 'fit') === 'moved'
        : b.zoomTo(h.section, pair.start, pair.end);
    }
    if (nav.blocks) {
      const block = r.blockAt(h.section, x, !!h.context);
      const v = this.viewer.view.get(h.section);
      if (!block || !v || !Number.isFinite(block.start)) return false;
      this.ownMove = true;
      const shown = Math.abs(v.start - block.start) <= 1 && Math.abs(v.end - block.end) <= 1;
      const view = this.viewer.view;
      return shown
        ? view.jump(() => view.fit(h.section), true, 'fit') === 'moved'
        : b.zoomTo(h.section, block.start, block.end);
    }
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
    const navigate = this.viewer.behavior<NavigateBehavior>('navigate');
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
    d.grip?.cancel();
    const set = d.set && this.setSquare;
    if (set) this.setGesture.cancel(set);
    if (d.rows) this.viewer.view.setBand(d.rows.track, { bandTop: d.rows.top });
    if (d.positions) this.behavior?.cancel();
    this.changed();
    return true;
  }

  /**
   * A key while the part has the focus (its window sliders, roving: PageDown and PageUp move the focus to the next
   * and the previous window): arrows pan (Up and Down on rows), Shift+arrows resize about the centre, Home and End go
   * to the ends (jumps), Enter centres on the readout's position, Esc cancels; the view's keys act on the focused
   * window's section (Alt+Left and Alt+Right whether it can move or not: they walk the history). Over a set the keys
   * are the set's (`SquareSetGesture.key`), Home showing all. True when taken.
   */
  key(key: string, modifiers: KeyModifiers = {}): boolean {
    const nav = this.navigator;
    const section = nav?.of[this._focused.value] ?? nav?.of[0];
    if (!nav || !section) return false;
    if (key === 'Escape') return this.cancel();
    const set = this.setSquare;
    // the set takes the arrows and Enter; PageDown, PageUp, End and the rest are the slider's (a block's step: `block.next`)
    if (set && !modifiers.altKey && (key === 'Enter' || key.startsWith('Arrow'))) {
      const toggle = !!(modifiers.metaKey || modifiers.ctrlKey);
      const taken = this.setGesture.key(set, key, !!modifiers.shiftKey, toggle);
      if (taken && !this.keyed) {
        this.keyed = true;
        this.changed();
      }
      // the keyboard's cell said (a move of the run says its run through the slider)
      if (taken) this.say(Square!.setCellReadout(set, this.setGesture.keyX, this.setGesture.keyY));
      return taken;
    }
    // over blocks the slider's large step is a block (`block.next`, `block.previous`); Home with a switch shows all
    if (nav.blocks && (key === 'PageDown' || key === 'PageUp'))
      return this.runIf(key === 'PageDown' ? 'block.next' : 'block.previous', { section });
    if (nav.switch === 'include' && (key === 'Home' || (key === '0' && isSquare(this.renderer))))
      return this.showAllOnSquare(key === '0');
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
        return (
          view.jump(() => view.zoom(section, this.viewer.config.navigate.doubleClickZoomStep), true, 'zoom') === 'moved'
        );
      case '-':
        return (
          view.jump(() => view.zoom(section, 1 / this.viewer.config.navigate.doubleClickZoomStep), true, 'zoom') ===
          'moved'
        );
      case '0':
        return view.jump(() => view.fit(section), true, 'fit') === 'moved';
      default:
        return false;
    }
  }

  /**
   * `navigator.switch: include`: the bound section shows the block under x alone (`section.include` with its name, one
   * undo step); on what it already shows alone, nothing. True when the command ran.
   */
  private switchTo(section: string, x: number, context = false): boolean {
    const block = this.renderer?.blockAt(section, x, context);
    if (!block) return false;
    const name = this.viewer.spec.model.value.systems.get(block.system)?.name;
    if (!name) return false;
    const node = this.viewer.view.resolved.value?.sections.find((s) => s.address === section);
    const shown = node?.blocks;
    if (shown?.length === 1 && shown[0].axis.system === block.system) return false;
    return this.runIf('section.include', { address: section, include: [name] });
  }

  /** The square with a switch: the bound section shows the cell's two blocks (one when on the diagonal). */
  private switchToPair(square: SquareNavigator, section: string, x: number, y: number): boolean {
    const cell = square.cellAt(x, y);
    if (!cell) return false;
    // in the lattice's order (the ring's), whichever axis the click came by
    const names = [...new Set([Math.min(cell.x, cell.y), Math.max(cell.x, cell.y)])]
      .map((k) => square.blockName(k))
      .filter(Boolean);
    return names.length > 0 && this.runIf('section.include', { address: section, include: names });
  }

  /**
   * `navigator.switch: include`: every block shown again on the bound section (its `include` removed; a click on the
   * navigator's title, Home on its focus). True when the command ran.
   */
  /** A square over a set writes its bound section's `include` (undefined: every block), one patch. */
  private writeSet(include: readonly string[] | undefined): boolean {
    const section = this.navigator?.of[0];
    if (!section) return false;
    return this.runIf('section.include', include ? { address: section, include } : { address: section });
  }

  /** The square's way home: every block shown, nothing when they all are (no step); one step otherwise. */
  private showAllOnSquare(guard = true, before?: ViewSnapshot): boolean {
    const set = this.setSquare;
    if (guard && set) {
      let all = true;
      for (let k = 0; k < set.blockCount && all; k++) all = set.shows(k);
      if (all) return false;
    }
    return this.showAll(before);
  }

  showAll(before?: ViewSnapshot): boolean {
    const nav = this.navigator;
    const section = nav?.of[0];
    if (!nav || nav.switch !== 'include' || !section) return false;
    return this.runIf('section.include', { address: section }, before);
  }

  /** Runs a command when it can run now; true when it ran. */
  private runIf(id: string, context: unknown, taken?: ViewSnapshot): boolean {
    const c = this.viewer.commands.find((x) => x.id === id);
    const can = c?.when?.(context as never);
    if (!c || (can !== undefined && can !== true)) return false;
    // a change of the blocks shown is one step of the view's history, pushed once it landed (a drag's on release:
    // never a frame's): Back shows the blocks it left
    const view = this.viewer.view;
    const before = id === 'section.include' ? (taken ?? view.snapshot({ includes: true })) : undefined;
    const done = this.viewer.execute(id, context);
    if (before) {
      if (typeof done === 'string') view.record(before);
      else void done.then(() => view.record(before));
    }
    return true;
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
