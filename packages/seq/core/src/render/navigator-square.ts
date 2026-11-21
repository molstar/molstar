/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type Issue, issueOf, type Node, type Value } from '@molstar/core/repr/tree';
import type { StaticModel } from '@molstar/seqx-core/resolve';
import {
  blockSpan,
  placeOfUnit,
  type ResolvedNavigator,
  type ResolvedSection,
  type ResolvedView,
  resolveViewNode,
  sameBlocks,
  unitOfPlace,
} from '@molstar/seqx-core/view';
import { DefaultRenderConfig, type RenderConfig } from './config.js';
import { DrawLists, NoColor } from './draw-lists.js';
import { RenderIssues } from './issues.js';
import type { SeqTokens } from './tokens.js';
import type { ViewLayout } from './tracks.js';
import type { NavigatorFrame, NavigatorFrameInput, NavigatorHit } from './navigator.js';
import { NavigatorWindow } from './navigator.js';

// The square (_spec/README.md "Render", Navigators, The square): a navigator on two axes over blocks × blocks, its
// lattice the assembly's chromosome pairs (after the assembly's `include`, never a section's), the cells the bound
// sections show plain and the rest washed, the window the x and y windows as one rectangle. What the cells hold (pairs
// by block, a grid's tiles by block pair) is region reads': the picture here is the lattice and what static data gives.

/** Whether a navigator is a square: two axes over blocks, bound to two sections of blocks (or one section twice). */
export function isSquareNavigator(nav: ResolvedNavigator, view: ResolvedView): boolean {
  if (!nav.blocks || nav.axes !== 'both' || nav.of.length !== 2) return false;
  return nav.of.every((a) => view.sections.find((s) => s.address === a)?.blocks?.length);
}

/** A cell of the square: its x and y blocks (indices in the lattice), their systems. */
export interface SquareCell {
  readonly x: number;
  readonly y: number;
  readonly xSystem: string;
  readonly ySystem: string;
}

export interface SquareNavigatorOptions {
  readonly model: StaticModel;
  readonly view: ResolvedView;
  readonly navigator: ResolvedNavigator;
  readonly tokens: SeqTokens;
  readonly config?: RenderConfig;
}

/** The lattice's inset from the box (CSS px): its outer lines are seen. */
const Inset = 1;

/**
 * Draws a square: its lattice over the blocks of the assembly on both axes (its lines past `navigator.square.maxCells`
 * pairs not drawn, said once), the cells the bound sections do not show washed, the windows as one rectangle; hits a
 * cell. Its draw lists are its own, made again on a size, a token or a window change.
 */
export class SquareNavigator {
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
  /** Which lattice blocks each axis's bound section shows (by index). */
  private readonly shownX: boolean[] = [];
  private readonly shownY: boolean[] = [];
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
    if (blocks.length * blocks.length > this.config.navigator.square.maxCells)
      this.issues.push(
        issueOf(RenderIssues, 'seq.render.navigator.layer', 'The square has too many cells: its lattice is not drawn', {
          address: this.navigator.address,
          hint: RenderIssues.codes['seq.render.navigator.layer'].hint,
        }),
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
  }

  /** The square adds no zoom strip. */
  get zoomHeight(): number {
    return 0;
  }

  get allIssues(): readonly Issue[] {
    return this.issues;
  }

  /** x (or y) of a lattice unit in the box. */
  private at(u: number, horizontal: boolean): number {
    const b = this.box;
    const len = (horizontal ? b.width : b.height) - 2 * Inset;
    return (horizontal ? b.x : b.y) + Inset + (this.total > 0 ? (u / this.total) * len : 0);
  }

  /** The lattice unit at x (or y), clamped. */
  private unitAt(p: number, horizontal: boolean): number {
    const b = this.box;
    const len = Math.max(1e-9, (horizontal ? b.width : b.height) - 2 * Inset);
    const f = (p - (horizontal ? b.x : b.y) - Inset) / len;
    return Math.max(0, Math.min(this.total, f * this.total));
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
    const b = this.box;
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
    this.box = { x: 0, y: 0, width, height };
    const lists = this.lists;
    lists.reset();
    lists.width = width;
    lists.height = height;
    const c = this.tokens.colors;
    const blocks = this.lattice?.blocks ?? [];
    const px = 1 / (input.dpr ?? 1);
    const bg = lists.writer('background');
    bg.fill(0, 0, width, height, c['seq-track']);
    const data = lists.writer('data');
    // the cells a bound section does not show: washed
    for (let i = 0; i < blocks.length; i++)
      for (let j = 0; j < blocks.length; j++) {
        if (this.shownX[i] && this.shownY[j]) continue;
        const o = (k: number) => blocks[k].offset ?? 0;
        const x0 = this.at(o(i), true);
        const x1 = this.at(o(i) + blockSpan(blocks[i]), true);
        const y0 = this.at(o(j), false);
        const y1 = this.at(o(j) + blockSpan(blocks[j]), false);
        data.fill(x0, y0, x1 - x0, y1 - y0, c['seq-wash']);
      }
    // the lattice's lines at each block's edges (not past `maxCells`)
    if (blocks.length * blocks.length <= this.config.navigator.square.maxCells)
      for (let k = 0; k <= blocks.length; k++) {
        const u = k < blocks.length ? (blocks[k].offset ?? 0) : this.total;
        const x = this.at(u, true);
        const y = this.at(u, false);
        data.fill(x - px / 2, this.box.y + Inset, px, height - 2 * Inset, c['msui-text-muted']);
        data.fill(this.box.x + Inset, y - px / 2, width - 2 * Inset, px, c['msui-text-muted']);
      }
    this.placeWindows(input.view);
    const out = lists.writer('window');
    const w = this.windows[0];
    const wy = this.windows[1];
    if (w.shown) {
      const y0 = wy.shown ? wy.y0 : w.x0 - this.box.x + this.box.y;
      const y1 = wy.shown ? wy.y1 : w.x1 - this.box.x + this.box.y;
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
    return this.frameRecord;
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
      const s = view?.sections.find((x) => x.section.address === nav.of[k]);
      if (!s?.axis) continue;
      const r = this.windowOn(s.section, s.axis.u0, s.axis.u1);
      if (!r) continue;
      const horizontal = k === 0;
      const a = this.at(r[0], horizontal);
      const b = this.at(r[1], horizontal);
      if (horizontal) {
        w.x0 = w.dx0 = a;
        w.x1 = w.dx1 = b;
        w.lo = this.box.x;
        w.hi = this.box.x + this.box.width;
      } else {
        w.y0 = w.dy0 = a;
        w.y1 = w.dy1 = b;
      }
      w.axes = 1;
      w.shown = true;
    }
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
    const w = this.windows[0];
    const inside = w.shown && x >= w.x0 && x <= w.x1;
    return {
      section,
      part: inside ? 'body' : 'outside',
      axis: 'positions',
      position,
      system: cell.xSystem,
      block: cell.x,
      ySystem: cell.ySystem,
      yBlock: cell.y,
    };
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
    const cell = this.cellAt(x, this.box.y + Inset);
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

/** Whether a bound section shows a system (one of its blocks is on it). */
function shows(s: ResolvedSection | undefined, system: string): boolean {
  return !!s?.blocks?.some((b) => b.axis.system === system);
}

/** The lattice: a section of the navigator's assembly (its `include` alone), else of the first bound section's `of`. */
function latticeOf(model: StaticModel, view: ResolvedView, nav: ResolvedNavigator): ResolvedSection | undefined {
  const first = view.sections.find((s) => s.address === nav.of[0]);
  let of: Value | undefined = nav.assemblies
    ? nav.assemblies.length === 1
      ? nav.assemblies[0]
      : (nav.assemblies as Value)
    : undefined;
  if (of === undefined && first) {
    let n: Node | undefined = model.spec;
    for (const id of first.address.split('/').filter(Boolean)) n = n?.children?.find((c) => c.id === id);
    of = n?.params?.of as Value | undefined;
  }
  if (of === undefined) return undefined;
  const node: Node = {
    kind: 'view',
    id: nav.id,
    params: { layout: 'tracks' },
    children: [{ kind: 'section', id: 'square', params: { of } }],
  };
  return resolveViewNode(model, node, nav.address).sections[0];
}
