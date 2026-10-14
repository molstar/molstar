/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// What the channels light (seqspec:interaction#channels, the marks and the channel table; CP2): a `SectionMarks` per
// section, each channel's `MarkSet` rebuilt by the projector only when that channel or the view changes, never on a
// layout change or another channel's change. The hover lights its element closed by its closure (pointer and keyboard
// `containing`, external `element` by default), an item its rows; the selection what was stored (closed at click time);
// the highlight channels their union per style, each closed by its own closure. In a section of an alignment's rows,
// the selection's and highlights' row residues and cells (a gap's too) are also rectangles. One per surface (CS1);
// `sync` before reading. The hover's rebuild allocates nothing once buffers have room.

import type { ResolvedView } from "@molstar/seq-core/format";
import {
  type AlignedLoci,
  Extents,
  CoordinateSystemKind as K,
  type Loci,
  OrderedSet,
  type SeqSet,
} from "@molstar/seq-core/model";
import type { Assembled, ViewModel } from "./graph.js";
import type { HighlightManager, HighlightStyle } from "./highlight.js";
import { HoverKind } from "./hit.js";

/** Read once: a member read in a hover-rate path allocates under the test runner's module wrapper. */
const HEADER = HoverKind.header;

import type { HoverManager, HoverState } from "./hover.js";
import type { SplitFeatures } from "./layout.js";
import { MarkBit, type MarkSet, type RowMarks, SectionMarks, type ShownRuns } from "./marks.js";
import { type Coordinates, isCoordinates, Projector } from "./projection.js";
import type { Selection } from "./selection.js";
import type { SeqState, Unsubscribable } from "./state.js";
import type { SystemRegistry } from "./systems.js";

/** The managers a surface's channels read (a context is one). */
export interface ChannelSources {
  readonly graph: { readonly assembled: SeqState<Assembled | null> };
  readonly hover: HoverManager;
  readonly selection: { readonly current: SeqState<Selection> };
  readonly highlight: HighlightManager;
  readonly systems: SystemRegistry;
  /** A section's axis, for the hover's primary; null: not laid out (the primary is then none). */
  readonly layout?: { axis(section: string): ShownRuns | null };
}

export const Channel = { hover: 1, selection: 2, highlight: 4 } as const;
const ALL = 7;
const UNIT = K.unit;
type Target = Loci | Coordinates;
const PARTS = ["seq", "boundaries", "rows", "columns", "items", "aligned"] as const;
const STYLES: readonly [HighlightStyle, number][] = [
  ["highlight", MarkBit.highlight],
  ["accent", MarkBit.accent],
  ["band", MarkBit.band],
];

export class ChannelMarks {
  /** This surface's projector over the shown view and set. */
  readonly projector = new Projector();
  /** Rebuilds per channel, each over every section (counters for tests). */
  readonly rebuilds = { hover: 0, selection: 0, highlight: 0 };
  private readonly src: ChannelSources;
  private readonly sections: SectionMarks[] = [];
  private view: ResolvedView | null = null;
  private model: ViewModel | null = null;
  /** The model whose split rows' units the sections hold (`ownUnits`). */
  private ownFor: ViewModel | null = null;
  private set: SeqSet | null = null;
  private dirty = ALL;
  private readonly ext = Extents.create(16);
  /** `ext` less a rows section's own rows and cells (its rectangles, not columns); `sys` copies an entry's system. */
  private readonly pos = Extents.create(16);
  /** An extent less the runs on its section's other own units (`only`). */
  private readonly mine = Extents.create(16);
  private readonly sys = { kind: 0, id: -1, row: -1 };
  private readonly subs: Unsubscribable[];

  constructor(src: ChannelSources) {
    this.src = src;
    this.projector.systems = src.systems;
    const on = (bits: number) => () => {
      this.dirty |= bits;
    };
    this.subs = [
      src.graph.assembled.subscribe((a) => {
        this.view = a?.view?.resolved ?? null;
        this.model = a?.view ?? null;
        this.set = a?.set ?? null;
        this.projector.update(this.view, this.set);
        this.dirty = ALL;
      }),
      src.hover.live.version.subscribe(on(Channel.hover)),
      src.selection.current.subscribe(on(Channel.selection)),
      src.highlight.version.subscribe(on(Channel.highlight)),
      src.systems.version.subscribe(on(ALL)),
    ];
  }

  /** The hover these marks light (an export's own unless it draws the hover): what a hover's drawing reads. */
  get hovered(): HoverState {
    return this.src.hover.live.current;
  }

  /** Section `s`'s marks (the view's section index), its changed channels rebuilt first; null: no such section. */
  marks(s: number): SectionMarks | null {
    this.sync();
    return s >= 0 && s < this.count ? this.sections[s] : null;
  }

  /** Rebuilds the channels that changed since the last call. */
  sync(): void {
    const d = this.dirty;
    if (d === 0) return;
    this.dirty = 0;
    while (this.sections.length < this.count) this.sections.push(new SectionMarks());
    if (this.ownFor !== this.model) this.ownUnits();
    if (d & Channel.hover) this.hover();
    if (d & Channel.selection) this.selection();
    if (d & Channel.highlight) this.highlights();
  }

  private get count(): number {
    return this.view ? this.view.sections.count : 0;
  }

  private hover() {
    const h = this.src.hover.live.current;
    const p = this.projector;
    const ext = this.ext;
    ext.count = 0;
    if (h.external) p.fromLoci(h.external, ext);
    else p.fromHover(h, ext);
    // An item's rows and a header's set close nothing (seqspec:interaction#channels).
    if (h.closure === "containing" && h.kind !== HoverKind.item && h.kind !== HEADER) p.close(ext);
    for (let s = 0; s < this.count; s++) {
      const m = this.sections[s];
      p.project(ext, s, m.hover);
      for (let k = 0; k < m.ownUnits.length; k++) p.project(this.only(ext, m.ownUnits, k), s, m.own[k].hover);
      m.hoverLayer = !h.external && h.kind === HoverKind.item ? h.layer : -1;
      m.hoverAnnotation = h.item.annotation;
      m.hoverItem = h.item.item;
      const v = this.view as ResolvedView;
      // The hover's own hit (seqspec:interaction#hover): in its own section only, never an external hover.
      const own = !h.external && h.kind !== HoverKind.none && h.section === s;
      m.hoverTrack = own ? h.track : -1;
      m.hoverBlock = own ? h.block : -1;
      m.hoverIndex = h.index;
      m.hoverState = h.state;
      m.hoverPart = h.part;
      m.hoverPairUnit = own ? h.pairUnit : -1;
      m.hoverPairIndex = h.pairIndex;
      m.hoverPairBlock = own ? h.pairBlock : -1;
      m.layerTrack = v.layers.track;
      const axis = this.src.layout?.axis((this.set as SeqSet).strings.get(v.sections.ref[s])) ?? null;
      if (axis) m.updatePrimary(axis, v.sections.blockCount[s]);
    }
    this.rebuilds.hover++;
  }

  private selection() {
    const sel = this.src.selection.current.value;
    const targets: Loci[] = [];
    for (const f of PARTS) if (sel[f]) targets.push(sel[f] as Loci);
    this.light(targets, [], MarkBit.selection, (m) => m.selection);
    this.rebuilds.selection++;
  }

  private highlights() {
    const all = [...this.src.highlight.channels.value.values()];
    // A channel of one section (a section's `highlights`) lights there only: each section then from its own list.
    const scoped = all.some((e) => e.section !== null);
    const strings = (this.set as SeqSet | null)?.strings;
    for (const [style, bit] of STYLES) {
      const of = (m: SectionMarks) =>
        bit === MarkBit.highlight ? m.highlight : bit === MarkBit.accent ? m.accent : m.band;
      for (let s = scoped ? 0 : -1; s < (scoped ? this.count : 0); s++) {
        const ref = s >= 0 && strings ? strings.get((this.view as ResolvedView).sections.ref[s]) : null;
        // Channels closed by containment first, then closed once together; the others after, as given.
        const mine = all.filter((e) => e.style === style && (e.section === null || e.section === ref));
        const closed = mine.filter((e) => e.closure === "containing").map((e) => e.target);
        const given = mine.filter((e) => e.closure !== "containing").map((e) => e.target);
        this.light(given, closed, bit, of, s);
      }
    }
    this.rebuilds.highlight++;
  }

  /**
   * One channel's targets into every section's `MarkSet` and rectangles under `bit` (`only` ≥ 0: that section's alone):
   * `closed` closed by containment, then `given` as they are; aligned loci last, so the rows' rectangles take their
   * cells as given (a gap's too).
   */
  private light(
    given: readonly Target[],
    closed: readonly Target[],
    bit: number,
    of: (m: SectionMarks) => MarkSet,
    only = -1,
  ) {
    const p = this.projector;
    const ext = this.ext;
    ext.count = 0;
    for (const t of closed) p.fromLoci(t, ext);
    if (closed.length) p.close(ext);
    const aligned = given.filter((t): t is AlignedLoci => !isCoordinates(t) && t.kind === "aligned-loci");
    for (const t of given) if (isCoordinates(t) || t.kind !== "aligned-loci") p.fromLoci(t, ext);
    const rows = ext.count;
    for (const t of aligned) p.fromLoci(t, ext);
    for (let s = Math.max(only, 0); s < (only >= 0 ? only + 1 : this.count); s++) {
      const m = this.sections[s];
      const pe = this.positions(ext, rows, this.rowsOf(s));
      p.project(pe, s, of(m));
      for (let k = 0; k < m.ownUnits.length; k++) p.project(this.only(pe, m.ownUnits, k), s, of(m.own[k]));
      m.rows.clear(bit);
      const n = ext.count;
      ext.count = rows; // the row residues before the aligned loci's
      p.rows(ext, s, bit, m.rows);
      ext.count = n;
      for (const l of aligned) this.cells(l, s, bit, m.rows);
    }
  }

  /**
   * Each section's split rows on units of their own (`SectionMarks.ownUnits`): a split row's one unit
   * (`SplitRows.unit`) that no block of the section shows. Per view.
   */
  private ownUnits() {
    const v = this.view;
    this.ownFor = this.model;
    for (let s = 0; s < this.count; s++) {
      const units: number[] = [];
      for (let t = 0; v && t < v.tracks.count; t++) {
        const k = v.tracks.section[t] === s ? v.tracks.split[t] : -1;
        const split = k >= 0 ? (this.model?.prepared[k]?.state as SplitFeatures | undefined)?.split : null;
        for (const u of split?.unit ?? []) if (u >= 0 && !units.includes(u) && !this.shows(s, u)) units.push(u);
      }
      const m = this.sections[s];
      m.ownUnits = Int32Array.from(units);
      while (m.own.length < units.length) m.own.push(new SectionMarks());
    }
  }
  private shows(s: number, u: number): boolean {
    const v = this.view as ResolvedView;
    const b0 = v.sections.blockStart[s];
    for (let b = b0; b < b0 + v.sections.blockCount[s]; b++) if (v.blocks.unit[b] === u) return true;
    return false;
  }
  /** `ext` less its runs on `units` other than `units[k]`: what a split row on that unit lights. Hover rate. */
  private only(ext: Extents, units: Int32Array, k: number): Extents {
    const out = this.mine;
    const sys = this.sys;
    out.count = 0;
    for (let i = 0; i < ext.count; i++) {
      if (ext.kind[i] === UNIT && ext.id[i] !== units[k] && units.indexOf(ext.id[i]) >= 0) continue;
      sys.kind = ext.kind[i];
      sys.id = ext.id[i];
      sys.row = ext.row[i];
      Extents.add(out, sys, ext.start[i], ext.stop[i]);
    }
    return out;
  }

  /** The alignment whose rows section `s` shows (a track with lanes; its first block's alignment); -1: none. */
  private rowsOf(s: number): number {
    const v = this.view as ResolvedView;
    for (let t = 0; t < v.tracks.count; t++)
      if (v.tracks.section[t] === s && v.tracks.laneCount[t] > 0) return v.blocks.alignment[v.sections.blockStart[s]];
    return -1;
  }
  /**
   * What lights positions in a section of alignment a's rows: `ext` less a's row residues and cells (entries from
   * `rows` on are the aligned loci's), which are its rectangles (never whole columns). Edge rate.
   */
  private positions(ext: Extents, rows: number, a: number): Extents {
    if (a < 0) return ext;
    const out = this.pos;
    const sys = this.sys;
    out.count = 0;
    for (let i = 0; i < ext.count; i++) {
      if (ext.id[i] === a && (ext.kind[i] === K.row || (i >= rows && ext.kind[i] === K.columns))) continue;
      sys.kind = ext.kind[i];
      sys.id = ext.id[i];
      sys.row = ext.row[i];
      Extents.add(out, sys, ext.start[i], ext.stop[i]);
    }
    return out;
  }

  /** An aligned loci's cells, gaps included, as rectangles on the blocks of section `s` that show its alignment. */
  private cells(l: AlignedLoci, s: number, bit: number, out: RowMarks) {
    const v = this.view as ResolvedView;
    const b0 = v.sections.blockStart[s];
    for (let b = 0; l.set === this.set && b < v.sections.blockCount[s]; b++) {
      if (v.blocks.alignment[b0 + b] !== l.alignment) continue;
      const [lo, hi] = [v.blocks.start[b0 + b], v.blocks.stop[b0 + b]];
      const runs = (cols: OrderedSet, row: (c0: number, c1: number) => void) => {
        for (let i = 0, n = OrderedSet.size(cols); i < n; ) {
          const c0 = OrderedSet.getAt(cols, i);
          let c1 = c0 + 1;
          while (++i < n && OrderedSet.getAt(cols, i) === c1) c1++;
          if (Math.max(c0, lo) < Math.min(c1, hi)) row(Math.max(c0, lo), Math.min(c1, hi));
        }
      };
      const el = l.elements;
      // A band: by column run, then row, so consecutive rows join one rectangle.
      if ("rows" in el)
        runs(el.columns, (c0, c1) => OrderedSet.forEach(el.rows, (r) => out.add(bit, b, r, c0, c1), null));
      else for (const e of el) runs(e.columns, (c0, c1) => out.add(bit, b, e.row, c0, c1));
    }
  }

  dispose(): void {
    for (const s of this.subs) s.unsubscribe();
    this.projector.update(null, null);
    this.view = this.set = null;
  }
}
