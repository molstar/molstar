/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type Issue, issueOf } from '@molstar/core/repr/tree';
import type { StaticModel } from '@molstar/seqx-core/resolve';
import {
  blockSpan,
  placeOfUnit,
  type ResolvedNavigator,
  type ResolvedSection,
  type ResolvedView,
  sameBlocks,
  unitOfPlace,
} from '@molstar/seqx-core/view';
import { DefaultRenderConfig, type RenderConfig } from './config.js';
import { Anchor, DrawLists, EmphasisPasses, NoColor, type PassName, type PassWriter, Round } from './draw-lists.js';
import { RenderIssues } from './issues.js';
import { bundledWidth } from './text.js';
import type { SeqTokens } from './tokens.js';
import type { ViewLayout } from './tracks.js';
import type { NavigatorFrame, NavigatorFrameInput, NavigatorHit } from './navigator.js';
import { NavigatorWindow } from './navigator.js';
import { SquarePictureCells, type SquarePictureOptions } from './navigator-square-cells.js';
import { latticeAssembly, latticeOf, sectionTotal, sharedPrefix, shows, washRuns } from './navigator-square-lattice.js';
import { drawShownSet, gripOf, type SquareGripHit, type SquareRun, shownRun } from './navigator-square-set.js';
import type { SquareLattice } from './layer-type.js';

export { washRuns } from './navigator-square-lattice.js';
export type { SquareGrip, SquareGripAxis, SquareGripHit, SquareRun } from './navigator-square-set.js';

// The square (_spec/README.md "Render", Navigators, The square): a navigator on two axes over blocks × blocks, its
// lattice the assembly's chromosome pairs (after the assembly's `include`, never a section's), the cells the bound
// sections show plain and the rest washed, the window the x and y windows as one rectangle. What the cells hold (pairs
// by block, a grid's tiles by block pair) is region reads': the picture here is the lattice and what static data gives.

export { isSquareNavigator } from './navigator-square-kind.js';

/** A cell of the square: its x and y blocks (indices in the lattice), their systems. */
export interface SquareCell {
  readonly x: number;
  readonly y: number;
  readonly xSystem: string;
  readonly ySystem: string;
}

/** The square's options; a picture's cells read the store and the composition's types (`SquarePictureOptions`). */
export interface SquareNavigatorOptions extends SquarePictureOptions {
  readonly model: StaticModel;
  readonly view: ResolvedView;
  readonly navigator: ResolvedNavigator;
  readonly tokens: SeqTokens;
  readonly config?: RenderConfig;
}

/** No wants (a square without a picture): one array, never one a frame. */
const NoWants: readonly never[] = [];
/** The keyboard's cell's inset in its cell and its corners' radius (CSS px). */
const KeyInset = 2;
/** The lattice's inset from the box (CSS px): its outer lines are seen. */
const Inset = 1;
/** The x labels' band above the lattice (CSS px); the y labels' column is as wide as the longest (at most a fifth). */
const LabelBand = 12;
/** The labels' font: a block's name where its cell has room for it, else none. */
const LabelFont = '10px sans';
const LabelPad = 2;
/** The alpha (of 255) of a sweep's tint: the window's own. */
const SweepTint = 0x2e;
const AnchorMiddle: number = Anchor.middle;
const AnchorEnd: number = Anchor.end;
/** The set's paint (`drawShownSet`): its ink, stroke width and grips' size. */
const SetStyle = new Float64Array(3);
/** What a window's move redraws (the base kept). */
const WindowPasses: readonly PassName[] = EmphasisPasses;

/**
 * Draws a square: its lattice over the blocks of the assembly on both axes (its lines past `navigator.square.maxCells`
 * pairs not drawn, said once), the cells the bound sections do not show washed, the windows as one rectangle; hits a
 * cell. Its draw lists are its own, made again on a size, a token or a window change.
 */
export class SquareNavigator {
  /** A square (a caller's type guard: the class is in a lazy chunk, so `instanceof` would load it). */
  readonly square = true;
  readonly navigator: ResolvedNavigator;
  readonly config: RenderConfig;
  private tokens: SeqTokens;
  /** The lattice: the assembly's blocks as one section (its offsets the units of both axes). */
  private readonly lattice: ResolvedSection | undefined;
  private readonly total: number;
  private readonly lists = new DrawLists();
  private readonly issues: Issue[] = [];
  private readonly windows: NavigatorWindow[];
  private box = { x: 0, y: 0, width: 0, height: 0 };
  /** The lattice's box in the navigator's (below the x labels' band, right of the y labels' column). */
  private grid = { x: 0, y: 0, width: 0, height: 0 };
  /** Which lattice blocks each axis's bound section shows (by index). */
  private readonly shownX: boolean[] = [];
  private readonly shownY: boolean[] = [];
  /** Each lattice block's name (its system's), and the widest's width in the labels' font. */
  private readonly names: string[];
  private labelWidth: number;
  /** The washed cells as rectangles of lattice block indices (i0, i1, j0, j1 each), made once a layout. */
  private readonly washRuns: Int32Array;
  /** What the base passes (the wash, the lattice, the labels) were drawn for: drawn again on a change only. */
  private baseW = Number.NaN;
  private baseH = Number.NaN;
  private basePx = Number.NaN;
  private baseTokens: SeqTokens | undefined;
  /** The picture's cells (a registered picture's layer type with `square`), else none: the lattice alone. */
  private readonly cells: SquarePictureCells | undefined;
  /** The lattice as the cells see it. */
  private readonly latticeView: SquareLattice;
  /**
   * One section on both axes with a switch (`navigator.switch: include`): the window is the blocks it shows (a run on
   * the diagonal, or scattered cells), its gestures the viewer's (_spec/README.md "Render", Navigators, The square).
   */
  readonly setMode: boolean;
  /** The shown run `[i0, i1)` of lattice blocks (`i0` -1: scattered), and a gesture's preview of another (-1: none). */
  readonly run: SquareRun = { i0: -1, i1: -1 };
  private readonly preview: SquareRun = { i0: -1, i1: -1 };
  /** The bound section's view shows its whole axis (else its window is drawn too: a followed ring's, EC-110). */
  private wholeView = true;
  private readonly gripHit: SquareGripHit = { grip: 'body', axis: 'both' };
  /** What each window was last placed from (its section, its units, the lattice's box) and whether it showed. */
  private readonly placed = [0, 1].map(() => ({
    section: undefined as unknown,
    u0: Number.NaN,
    u1: Number.NaN,
    box0: Number.NaN,
    box1: Number.NaN,
    box2: Number.NaN,
    box3: Number.NaN,
    shown: false,
    whole: true,
  }));
  private readonly edgeFn = (k: number, horizontal: boolean) => this.edgePx(k, horizontal);
  private readonly frameRecord: {
    lists: DrawLists;
    rebuilt: boolean;
    base: boolean;
    wants: readonly never[];
    windows: readonly NavigatorWindow[];
    issues: readonly Issue[];
  };

  constructor(readonly options: SquareNavigatorOptions) {
    this.navigator = options.navigator;
    this.config = options.config ?? DefaultRenderConfig;
    this.tokens = options.tokens;
    this.lattice = latticeOf(options.model, options.view, options.navigator);
    let total = 0;
    for (const b of this.lattice?.blocks ?? []) total = Math.max(total, (b.offset ?? 0) + blockSpan(b));
    this.total = total;
    this.windows = this.navigator.of.map(() => new NavigatorWindow());
    const blocks = this.lattice?.blocks ?? [];
    const sections = this.navigator.of.map((a) => options.view.sections.find((s) => s.address === a));
    for (let k = 0; k < blocks.length; k++) {
      this.shownX.push(shows(sections[0], blocks[k].axis.system));
      this.shownY.push(shows(sections[1], blocks[k].axis.system));
    }
    this.names = blocks.map((b) => options.model.systems.get(b.axis.system)?.name ?? b.label);
    this.labelWidth = 0;
    for (const n of this.names) this.labelWidth = Math.max(this.labelWidth, bundledWidth(LabelFont, n));
    // the wash and the lattice within `maxCells` pairs (at most as many rectangles), else neither: never n² a frame
    const max = this.config.navigator.square.maxCells;
    const within = blocks.length * blocks.length <= max;
    this.washRuns = within ? washRuns(this.shownX, this.shownY) : new Int32Array(0);
    this.setMode =
      this.navigator.switch === 'include' && this.navigator.of[0] === this.navigator.of[1] && blocks.length > 0;
    if (this.setMode) shownRun(this.shownX, this.run);
    const of = this.lattice ? latticeAssembly(options.model, options.view, options.navigator) : undefined;
    this.cells =
      of !== undefined
        ? SquarePictureCells.of(options.model, options.view, this.navigator, of, this.tokens, this.config, options)
        : undefined;
    const self = this;
    this.latticeView = {
      systems: blocks.map((b) => b.axis.system),
      names: this.names,
      edge: this.edgeFn,
      get tokens() {
        return self.tokens;
      },
    };
    if (!within)
      this.issues.push(
        issueOf(
          RenderIssues,
          'seq.render.navigator.square',
          `The square has ${blocks.length.toLocaleString('en-US')} blocks a side, past ${max.toLocaleString('en-US')} cells: its lattice and wash are not drawn`,
          { address: this.navigator.address, hint: RenderIssues.codes['seq.render.navigator.square'].hint },
        ),
      );
    this.frameRecord = {
      lists: this.lists,
      rebuilt: true,
      base: true,
      wants: [],
      windows: this.windows,
      issues: this.issues,
    };
  }

  setTokens(tokens: SeqTokens): void {
    this.tokens = tokens;
    this.cells?.setTokens(tokens);
  }

  /** The square adds no zoom strip. */
  get zoomHeight(): number {
    return 0;
  }

  get allIssues(): readonly Issue[] {
    return this.issues;
  }

  /** x (or y) of a lattice unit in the box: each block's cell as long as the block (sized by length). */
  private at(u: number, horizontal: boolean): number {
    const b = this.grid;
    const len = horizontal ? b.width : b.height;
    return (horizontal ? b.x : b.y) + (this.total > 0 ? (u / this.total) * len : 0);
  }

  /** The lattice unit at x (or y), clamped. */
  private unitAt(p: number, horizontal: boolean): number {
    const b = this.grid;
    const len = Math.max(1e-9, horizontal ? b.width : b.height);
    const f = (p - (horizontal ? b.x : b.y)) / len;
    return Math.max(0, Math.min(this.total, f * this.total));
  }

  /** The box and the lattice's box in it: the x labels' band above, the y labels' column left (none when tiny). */
  private place(width: number, height: number): void {
    this.box = { x: 0, y: 0, width, height };
    const labels = width >= 60 && height >= 60;
    const left = labels ? Math.min(width / 5, Math.ceil(this.labelWidth) + 2 * LabelPad) : 0;
    const top = labels ? LabelBand : 0;
    this.grid = {
      x: left + Inset,
      y: top + Inset,
      width: Math.max(0, width - left - 2 * Inset),
      height: Math.max(0, height - top - 2 * Inset),
    };
  }

  /** The lattice block at a unit (the last starting at or before it). */
  private blockOf(u: number): number {
    const blocks = this.lattice?.blocks ?? [];
    let k = 0;
    while (k + 1 < blocks.length && (blocks[k + 1].offset ?? 0) <= u) k++;
    return k;
  }

  /** The cell under a point; undefined off the lattice. */
  cellAt(x: number, y: number): SquareCell | undefined {
    const blocks = this.lattice?.blocks;
    const b = this.grid;
    if (!blocks?.length || x < b.x || y < b.y || x >= b.x + b.width || y >= b.y + b.height) return undefined;
    const i = this.blockOf(this.unitAt(x, true));
    const j = this.blockOf(this.unitAt(y, false));
    return { x: i, y: j, xSystem: blocks[i].axis.system, ySystem: blocks[j].axis.system };
  }

  /** A bound section's window on the lattice's units (identity where the block lists are equal, else by places). */
  private windowOn(s: ResolvedSection, start: number, end: number): readonly [number, number] | undefined {
    const lattice = this.lattice;
    if (!lattice) return undefined;
    if (sameBlocks(s, lattice)) return [start, end];
    const a = placeOfUnit(s, start, 'start');
    const b = placeOfUnit(s, end, 'end');
    const ua = a && unitOfPlace(lattice, a);
    const ub = b && unitOfPlace(lattice, b);
    return ua !== undefined && ub !== undefined ? [Math.min(ua, ub), Math.max(ua, ub)] : undefined;
  }

  frame(input: NavigatorFrameInput): NavigatorFrame {
    const { width, height } = input.size;
    const lists = this.lists;
    const c = this.tokens.colors;
    const px = 1 / (input.dpr ?? 1);
    // the base (the cells, the wash, the lattice, the labels) drawn once a size and tokens, and when the cells change; a
    // window move redraws the window alone (nothing made: the key is numbers)
    const sized =
      width !== this.baseW || height !== this.baseH || px !== this.basePx || this.tokens !== this.baseTokens;
    if (sized) {
      this.baseW = width;
      this.baseH = height;
      this.basePx = px;
      this.baseTokens = this.tokens;
      this.place(width, height);
    }
    const cells = this.cells;
    // asked first with no number passed (a call that is not inlined boxes one): a window's move lays nothing out
    const settled = input.settled ?? false;
    const changed =
      !!cells &&
      cells.due(settled, sized) &&
      cells.layout(this.grid.width, input.dpr ?? 1, settled, input.overlay, sized);
    const rebuilt = sized || changed;
    if (rebuilt) {
      lists.reset();
      lists.width = width;
      lists.height = height;
      this.drawBase(px);
    } else lists.reset(WindowPasses);
    this.frameRecord.rebuilt = rebuilt;
    this.frameRecord.base = rebuilt;
    this.frameRecord.wants = (cells?.wants ?? NoWants) as never[];
    this.placeWindows(input.view);
    const out = lists.writer('window');
    const w = this.windows[0];
    const wy = this.windows[1];
    if (this.setMode) {
      const p = this.preview;
      const shown = p.i0 >= 0 ? p : this.run;
      // its paint through an array (a number passed to a call that is not inlined is boxed)
      const style = SetStyle;
      style[0] = c['seq-window'];
      style[1] = Math.max(px, 1.5);
      style[2] = shown.i0 >= 0 ? this.config.navigator.window.handlePx : 0;
      drawShownSet(out, this.shownX, shown, this.edgeFn, style, this.maxCells);
      if (this.keyI >= 0 && this.keyJ >= 0 && this.keyI < this.blockCount && this.keyJ < this.blockCount) {
        // inset and rounded in the accent: a cursor, never read as the window's outline
        const g = out.args;
        g[0] = this.edgePx(this.keyI, true) + KeyInset;
        g[1] = this.edgePx(this.keyJ, false) + KeyInset;
        g[2] = Math.max(px, this.edgePx(this.keyI + 1, true) - g[0] - KeyInset);
        g[3] = Math.max(px, this.edgePx(this.keyJ + 1, false) - g[1] - KeyInset);
        g[4] = NoColor;
        g[5] = c['msui-accent'];
        g[6] = Math.max(px, 1.5);
        g[7] = KeyInset;
        g[8] = Round.both;
        out.rectArgs();
      }
    }
    // over a set the window is drawn too where the view does not show the whole section (a followed ring's, EC-110)
    if (w.shown && (!this.setMode || !this.wholeView)) {
      const y0 = wy.shown ? wy.y0 : w.y0;
      const y1 = wy.shown ? wy.y1 : w.y1;
      const g = out.args;
      g[0] = w.x0;
      g[1] = y0;
      g[2] = Math.max(px, w.x1 - w.x0);
      g[3] = Math.max(px, y1 - y0);
      g[4] = NoColor;
      g[5] = c['seq-window'];
      g[6] = Math.max(px, 1.5);
      g[7] = 0;
      g[8] = 0;
      out.rectArgs();
    }
    const sw = input.sweep;
    if (sw && w.shown) this.drawSweep(out, sw, wy.shown ? wy : w, px, c['seq-window']);
    return this.frameRecord;
  }

  /**
   * A sweep in progress (_spec/README.md "Render", Navigators, The square): the rectangle its release zooms to (whole cells, `snapSweep`), tinted and
   * outlined as the window is. Along x the swept range; along y the same range where one section is on both axes (its
   * window the diagonal), else the y section's window as it is.
   */
  private drawSweep(out: PassWriter, sw: { x0: number; x1: number }, wy: NavigatorWindow, px: number, ink: number) {
    const x0 = Math.min(sw.x0, sw.x1);
    const x1 = Math.max(sw.x0, sw.x1);
    if (!(x1 - x0 > 0)) return;
    const one = this.navigator.of[1] === this.navigator.of[0];
    const y0 = one ? this.at(this.unitAt(x0, true), false) : wy.y0;
    const y1 = one ? this.at(this.unitAt(x1, true), false) : wy.y1;
    const g = out.args;
    g[0] = x0;
    g[1] = y0;
    g[2] = Math.max(px, x1 - x0);
    g[3] = Math.max(px, y1 - y0);
    g[4] = ((ink & 0xffffff00) | SweepTint) >>> 0;
    g[5] = ink;
    g[6] = Math.max(px, 1.5);
    g[7] = 0;
    g[8] = 0;
    out.rectArgs();
  }

  /**
   * The base passes: the background, the washed cells (by runs of blocks: a column run of blocks the x section does not
   * show, a row run the y section does not show over the columns it does; at most `maxCells` rectangles, else none),
   * the lattice's lines at the blocks' edges (not past `maxCells` pairs) and each block's name where its cell has room.
   */
  private drawBase(px: number): void {
    const lists = this.lists;
    const c = this.tokens.colors;
    const blocks = this.lattice?.blocks ?? [];
    const g = this.grid;
    lists.writer('background').fill(0, 0, this.box.width, this.box.height, c['seq-track']);
    const data = lists.writer('data');
    // the picture's cells under the wash and the lines (always drawn: never bounded by `maxCells`)
    this.cells?.draw(data, this.latticeView);
    // nothing to show yet: the pending hatch over the lattice's box, the cells' words in its middle (`Loading…`, a failure)
    const note = this.cells?.note();
    if (note) {
      const pending = c['seq-pending'];
      for (let x = g.x; x < g.x + g.width; x += 8) data.fill(x, g.y, Math.min(4, g.x + g.width - x), g.height, pending);
      lists
        .writer('labels')
        .text(
          note,
          g.x + g.width / 2,
          g.y + g.height / 2 + 4,
          AnchorMiddle,
          LabelFont,
          c['msui-text-muted'],
          bundledWidth(LabelFont, note),
        );
    }
    if (blocks.length * blocks.length > this.config.navigator.square.maxCells) return;
    const edge = (k: number) => (k < blocks.length ? (blocks[k].offset ?? 0) : this.total);
    const runs = this.washRuns;
    for (let r = 0; r + 3 < runs.length; r += 4) {
      const x0 = this.at(edge(runs[r]), true);
      const x1 = this.at(edge(runs[r + 1]), true);
      const y0 = this.at(edge(runs[r + 2]), false);
      const y1 = this.at(edge(runs[r + 3]), false);
      data.fill(x0, y0, x1 - x0, y1 - y0, c['seq-wash']);
    }
    for (let k = 0; k <= blocks.length; k++) {
      const x = this.at(edge(k), true);
      const y = this.at(edge(k), false);
      data.fill(x - px / 2, g.y, px, g.height, c['msui-text-muted']);
      data.fill(g.x, y - px / 2, g.width, px, c['msui-text-muted']);
    }
    // the names: along x in the band above, along y in the column left (right-aligned), where the cell has room
    if (g.y <= Inset && g.x <= Inset) return;
    const text = lists.writer('labels');
    const ink = c['msui-text-muted'];
    // along x a name too long for its column drops the prefix every name shares (`chr7` → `7`)
    const prefix = sharedPrefix(this.names);
    for (let k = 0; k < blocks.length; k++) {
      const name = this.names[k];
      const w = bundledWidth(LabelFont, name);
      const x0 = this.at(edge(k), true);
      const x1 = this.at(edge(k + 1), true);
      const short = name.slice(prefix.length);
      const ws = bundledWidth(LabelFont, short);
      const top = x1 - x0 >= w + LabelPad ? name : short && x1 - x0 >= ws + LabelPad ? short : '';
      if (g.y > Inset && top)
        text.text(top, (x0 + x1) / 2, g.y - 3, AnchorMiddle, LabelFont, ink, top === name ? w : ws);
      const y0 = this.at(edge(k), false);
      const y1 = this.at(edge(k + 1), false);
      if (g.x > Inset && y1 - y0 >= 10 && w <= g.x - Inset - LabelPad)
        text.text(name, g.x - LabelPad, (y0 + y1) / 2 + 3.5, AnchorEnd, LabelFont, ink, w);
    }
  }

  /** The bound sections' windows from the view's last layout: the first along x, the second along y. */
  private placeWindows(view: ViewLayout | undefined): void {
    const nav = this.navigator;
    for (let k = 0; k < this.windows.length; k++) {
      const w = this.windows[k];
      w.section = nav.of[k];
      w.shown = false;
      w.axes = 0;
      // one section on both axes: its one window, drawn on the diagonal
      if (k === 1 && nav.of[1] === nav.of[0]) continue;
      // (a loop, not `find`: a frame allocates no closure)
      let s: ViewLayout['sections'][number] | undefined;
      const sections = view?.sections;
      for (let i = 0; sections && i < sections.length && !s; i++)
        if (sections[i].section.address === nav.of[k]) s = sections[i];
      if (!s?.axis) continue;
      // the same window over the same lattice as last placed: its numbers and places kept (a gesture's preview, a
      // move of another window: nothing made)
      const last = this.placed[k];
      const box = this.grid;
      if (
        last.section === s.section &&
        last.u0 === s.axis.u0 &&
        last.u1 === s.axis.u1 &&
        last.box0 === box.x &&
        last.box1 === box.y &&
        last.box2 === box.width &&
        last.box3 === box.height
      ) {
        if (k === 0) this.wholeView = last.whole;
        w.shown = last.shown;
        w.axes = last.shown ? 1 : 0;
        continue;
      }
      last.section = s.section;
      last.u0 = s.axis.u0;
      last.u1 = s.axis.u1;
      last.box0 = box.x;
      last.box1 = box.y;
      last.box2 = box.width;
      last.box3 = box.height;
      last.shown = false;
      last.whole = s.axis.u0 <= 1e-6 && s.axis.u1 >= sectionTotal(s.section) - 1e-6;
      if (k === 0) this.wholeView = last.whole;
      const r = this.windowOn(s.section, s.axis.u0, s.axis.u1);
      if (!r) continue;
      const horizontal = k === 0;
      const a = this.at(r[0], horizontal);
      const b = this.at(r[1], horizontal);
      if (horizontal) {
        w.x0 = w.dx0 = a;
        w.x1 = w.dx1 = b;
        w.lo = this.grid.x;
        w.hi = this.grid.x + this.grid.width;
        // on the diagonal until a y window says otherwise
        w.y0 = w.dy0 = this.at(r[0], false);
        w.y1 = w.dy1 = this.at(r[1], false);
      } else {
        w.y0 = w.dy0 = a;
        w.y1 = w.dy1 = b;
      }
      // its ends as places (what its title and slider say: `chr3: 1–198,295,559 of Toy38`)
      const lattice = this.lattice!;
      const pa = placeOfUnit(lattice, r[0], 'start');
      const pb = placeOfUnit(lattice, r[1], 'end');
      w.s0 = pa?.system ?? '';
      w.s1 = pb?.system ?? '';
      w.p0 = pa?.position ?? 0;
      w.p1 = pb?.position ?? 0;
      w.axes = 1;
      w.shown = true;
      last.shown = true;
    }
  }

  /** The lattice's blocks a side. */
  get blockCount(): number {
    return this.lattice?.blocks?.length ?? 0;
  }

  /** A lattice block's name (its system's). */
  blockName(k: number): string {
    return this.names[k] ?? '';
  }

  /** A lattice block's system. */
  blockSystem(k: number): string {
    return this.lattice?.blocks?.[k]?.axis.system ?? '';
  }

  /** Whether the bound section shows lattice block `k` (along x: with one section on both axes, the set's). */
  shows(k: number): boolean {
    return this.shownX[k] === true;
  }

  /** The px of lattice block `k`'s start along x (or y); `k` the count: the lattice's end. */
  edgePx(k: number, horizontal: boolean): number {
    const blocks = this.lattice?.blocks;
    const u = !blocks || k >= blocks.length ? this.total : (blocks[Math.max(0, k)].offset ?? 0);
    return this.at(u, horizontal);
  }

  private get maxCells(): number {
    return this.config.navigator.square.maxCells;
  }

  /** What a cell's picture says of it (`37 links`), undefined without one. */
  cellText(i: number, j: number): string | undefined {
    return this.cells?.text(i, j, this.latticeView);
  }

  /** The keyboard's cell (lattice indices; -1 none), outlined in the accent while the square has the focus. */
  setKeyCell(i: number, j: number): void {
    this.keyI = i;
    this.keyJ = j;
  }

  private keyI = -1;
  private keyJ = -1;

  /** A gesture's run drawn instead of the shown one (`[i0, i1)` of lattice blocks); `i0` -1 clears it. */
  setPreview(i0: number, i1: number): void {
    this.preview.i0 = i0;
    this.preview.i1 = i0 < 0 ? -1 : i1;
  }

  /**
   * The grip of the shown run's box under a point (its body, its start or end edge or corner), within
   * `navigator.window.hitPx` / 2 of an edge; undefined off the box, over a scattered set, or not a set. The object
   * returned is reused.
   */
  gripAt(x: number, y: number, touch = false): SquareGripHit | undefined {
    const r = this.run;
    if (!this.setMode || r.i0 < 0) return undefined;
    const cfg = this.config.navigator.window;
    const slop = (touch ? cfg.touchHitPx : cfg.hitPx) / 2;
    const x0 = this.edgePx(r.i0, true);
    const x1 = this.edgePx(r.i1, true);
    const y0 = this.edgePx(r.i0, false);
    const y1 = this.edgePx(r.i1, false);
    return gripOf(x, y, x0, x1, y0, y1, slop, this.gripHit) ? this.gripHit : undefined;
  }

  /**
   * The lattice edge nearest a point (0 to the count) along x, y, or the diagonal (`both`: the mean of its two units):
   * where a resized run's edge snaps (whole blocks).
   */
  edgeAt(x: number, y: number, axis: 'x' | 'y' | 'both'): number {
    this.pointer[0] = x;
    this.pointer[1] = y;
    return this.edgeAtPointer(axis);
  }

  /** The lattice block under a point along x, y or the diagonal (a run's move: the block the pointer is over). */
  blockAtPoint(x: number, y: number, axis: 'x' | 'y' | 'both'): number {
    this.pointer[0] = x;
    this.pointer[1] = y;
    return this.blockAtPointer(axis);
  }

  /**
   * A gesture's point, set by the caller before `blockAtPointer` and `edgeAtPointer` (a number passed to a call that is
   * not inlined is boxed: a drag's moves read it here, nothing allocated).
   */
  readonly pointer = new Float64Array(2);

  /** `blockAtPoint` of `pointer`. */
  blockAtPointer(axis: 'x' | 'y' | 'both'): number {
    return this.blockOf(this.pointerUnit(axis));
  }

  /** `edgeAt` of `pointer`. */
  edgeAtPointer(axis: 'x' | 'y' | 'both'): number {
    const blocks = this.lattice?.blocks;
    if (!blocks?.length) return 0;
    const u = this.pointerUnit(axis);
    let best = 0;
    let d = Number.POSITIVE_INFINITY;
    for (let k = 0; k <= blocks.length; k++) {
      const dk = Math.abs((k < blocks.length ? (blocks[k].offset ?? 0) : this.total) - u);
      if (dk < d) {
        d = dk;
        best = k;
      }
    }
    return best;
  }

  /** `pointer`'s lattice unit along x, y or the diagonal (read from the array: nothing boxed). */
  private pointerUnit(axis: 'x' | 'y' | 'both'): number {
    const g = this.grid;
    const p = this.pointer;
    const ux = Math.max(0, Math.min(this.total, ((p[0] - g.x) / Math.max(1e-9, g.width)) * this.total));
    const uy = Math.max(0, Math.min(this.total, ((p[1] - g.y) / Math.max(1e-9, g.height)) * this.total));
    return axis === 'x' ? ux : axis === 'y' ? uy : (ux + uy) / 2;
  }

  /** The block whose name is under a point (the x labels' band above the lattice, the y labels' column left); -1 else. */
  labelAt(x: number, y: number): number {
    const blocks = this.lattice?.blocks;
    const g = this.grid;
    if (!blocks?.length) return -1;
    if (y < g.y && y >= 0 && x >= g.x && x < g.x + g.width) return this.blockOf(this.unitAt(x, true));
    if (x < g.x && x >= 0 && y >= g.y && y < g.y + g.height) return this.blockOf(this.unitAt(y, false));
    return -1;
  }

  /** Whether a point is inside the window's rectangle (both axes). */
  private inWindow(x: number, y: number): boolean {
    const w = this.windows[0];
    const wy = this.windows[1];
    const y0 = wy?.shown ? wy.y0 : w.y0;
    const y1 = wy?.shown ? wy.y1 : w.y1;
    return w.shown && x >= w.x0 && x <= w.x1 && y >= y0 && y <= y1;
  }

  /** A bound section's unit at y (the second; one section on both axes: the first), through the lattice's places. */
  sectionUnitAtY(section: string, y: number): number | undefined {
    const s = this.options.view.sections.find((v) => v.address === section);
    const lattice = this.lattice;
    if (!s || !lattice) return undefined;
    const u = this.unitAt(y, false);
    if (sameBlocks(s, lattice)) return u;
    const at = placeOfUnit(lattice, u, 'start');
    return at ? unitOfPlace(s, at) : undefined;
  }

  /** A point: the cell's x block and position on it (the first bound section's axis), never a part but the window's body. */
  hit(x: number, y: number): NavigatorHit | undefined {
    const cell = this.cellAt(x, y);
    const section = this.navigator.of[0];
    if (!cell || section === undefined) return undefined;
    const blocks = this.lattice!.blocks!;
    const u = this.unitAt(x, true);
    const block = blocks[cell.x];
    const position = Math.max(0, (u - (block.offset ?? 0)) / (block.scale ?? 1)) + (block.axis.ranges?.[0] ?? 0);
    // over a set: the run's grip under the point (its body, an edge or corner), never the scattered cells
    const grip = this.setMode ? this.gripAt(x, y) : undefined;
    const inside = this.setMode ? grip !== undefined : this.inWindow(x, y);
    return {
      section,
      part: !inside ? 'outside' : (grip?.grip ?? 'body'),
      axis: 'positions',
      position,
      system: cell.xSystem,
      block: cell.x,
      ySystem: cell.ySystem,
      yBlock: cell.y,
    };
  }

  /**
   * A sweep from x `xa` to x `xb` snapped to whole cells (the lattice blocks under its ends, both included), as the
   * bound section's units `[start, end]`: what the sweep draws and its release fits (_spec/README.md "Render",
   * Navigators, The square). Undefined when the section shows none of it.
   */
  snapSweep(section: string, xa: number, xb: number): readonly [number, number] | undefined {
    const s = this.options.view.sections.find((v) => v.address === section);
    const lattice = this.lattice;
    const blocks = lattice?.blocks;
    if (!s || !lattice || !blocks?.length) return undefined;
    const lo = this.blockOf(this.unitAt(Math.min(xa, xb), true));
    const hi = this.blockOf(this.unitAt(Math.max(xa, xb), true));
    const u0 = blocks[lo].offset ?? 0;
    const u1 = hi + 1 < blocks.length ? (blocks[hi + 1].offset ?? 0) : this.total;
    if (sameBlocks(s, lattice)) return [u0, u1];
    const a = placeOfUnit(lattice, u0, 'start');
    const b = placeOfUnit(lattice, u1, 'end');
    const ua = a && unitOfPlace(s, a);
    const ub = b && unitOfPlace(s, b);
    return ua !== undefined && ub !== undefined ? [Math.min(ua, ub), Math.max(ua, ub)] : undefined;
  }

  /** The first bound section's unit at x (a click centres it there), through the places of the lattice. */
  sectionUnitAt(section: string, x: number, _context = false): number | undefined {
    const s = this.options.view.sections.find((v) => v.address === section);
    const lattice = this.lattice;
    if (!s || !lattice) return undefined;
    const u = this.unitAt(x, true);
    if (sameBlocks(s, lattice)) return u;
    const at = placeOfUnit(lattice, u, 'start');
    return at ? unitOfPlace(s, at) : undefined;
  }

  xOfSectionUnit(section: string, u: number): number {
    const s = this.options.view.sections.find((v) => v.address === section);
    const r = s && this.windowOn(s, u, u);
    return r ? this.at(r[0], true) : Number.NaN;
  }

  /** The first bound section's units of the x block under x: what a double-click fits. */
  blockAt(
    section: string,
    x: number,
    _context = false,
  ): { start: number; end: number; system: string; block: number } | undefined {
    const cell = this.cellAt(x, this.grid.y);
    const s = this.options.view.sections.find((v) => v.address === section);
    const lattice = this.lattice;
    if (!cell || !s || !lattice) return undefined;
    const b = lattice.blocks![cell.x];
    const o = b.offset ?? 0;
    const r = sameBlocks(s, lattice) ? ([o, o + blockSpan(b)] as const) : this.sectionRange(s, b.axis.system);
    return {
      start: r ? r[0] : Number.NaN,
      end: r ? r[1] : Number.NaN,
      system: b.axis.system,
      block: cell.x,
    };
  }

  /**
   * The cell under a point as a bound section's units (one section on both axes: the hull of its two blocks, what a
   * double-click fits) and its two systems; undefined off the lattice or where the section shows neither.
   */
  pairAt(section: string, x: number, y: number): { start: number; end: number; cell: SquareCell } | undefined {
    const cell = this.cellAt(x, y);
    const s = this.options.view.sections.find((v) => v.address === section);
    if (!cell || !s) return undefined;
    const a = this.sectionRange(s, cell.xSystem);
    const b = this.sectionRange(s, cell.ySystem);
    if (!a && !b) return undefined;
    const lo = Math.min(a?.[0] ?? Infinity, b?.[0] ?? Infinity);
    const hi = Math.max(a?.[1] ?? -Infinity, b?.[1] ?? -Infinity);
    return { start: lo, end: hi, cell };
  }

  /** A bound section's units of the block showing a system (its first), undefined when it shows none. */
  private sectionRange(s: ResolvedSection, system: string): readonly [number, number] | undefined {
    const b = s.blocks?.find((x) => x.axis.system === system);
    return b ? [b.offset ?? 0, (b.offset ?? 0) + blockSpan(b)] : undefined;
  }

  markText(): string | undefined {
    return undefined;
  }

  xRangeOf(_start: number, _end: number): readonly [number, number] | undefined {
    return undefined;
  }

  positionAtX(_x: number): number | undefined {
    return undefined;
  }

  readout(_position: number): string | undefined {
    return undefined;
  }

  rowAt(_y: number): number | undefined {
    return undefined;
  }
}
