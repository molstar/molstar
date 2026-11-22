/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { issueOf } from '@molstar/core/repr/tree';
import { formatGenomic } from '@molstar/core/util/number';
import {
  findNumbering,
  type Numbering,
  type NumberingOf,
  parseLabel,
  type QuerySystem,
} from '@molstar/seqx-core/query';
import type { NumberingEntry, StaticModel } from '@molstar/seqx-core/resolve';
import type { ResolvedLayer } from '@molstar/seqx-core/view';
import type { AxisTransform } from './axis.js';
import { fontSize } from './config.js';
import { Anchor, drawPass, grown, NoColor, type PassName, type PassWriter, passDraws } from './draw-lists.js';
import { RenderIssues } from './issues.js';
import { type Box, type Emphasis, type Hit, type LayerType, type LayoutContext, sectionId } from './layer-type.js';
import { readSourceNumbering } from './numberings.js';
import { mutedInk, type SeqTokens, withAlpha } from './tokens.js';

// the anchor a number read once (see ranges-draw.ts)
const AnchorMiddle: number = Anchor.middle;
const AnchorStart: number = Anchor.start;

// The `ruler` layer type (_spec/README.md "Render", Layer types): the section's axis labelled by a numbering. Ticks at
// the labels that are multiples of a step of 1, 2 or 5 × 10ⁿ (the least whose px apart reach the label gap); insertion
// labels (100A, 111.1) are minor ticks, labeled where they clear what is drawn; a label that does not fit is skipped,
// never shrunk.
// It labels the positions where they are, so a reversed axis counts down and a collapsed one labels its pieces. Where a
// numbering has no label (past its coverage: a light chain's constant domain under Kabat) `unnumbered` says what is
// there (axis.md, "Where a numbering does not reach"): `mark` (default) a muted dashed line over the run, one muted note
// at its start and no ticks; `fallback` the run labelled by `fallbackNumbering` (ordinal when absent) in the muted
// style. A run is only the head and the tail past the numbering's labelled extent (its first and last labelled
// position, found once a numbering): a position without a label between them stays blank, as it was.

/** The ruler's own position loop covers this many positions a piece at most; past it labels are looked up by number. */
const LoopLimit = 4096;
/** Insertion labels are minor ticks from this many px a position (a tick at each is a comb below it). */
const MinorScale = 4;
/** The labels' font: tabular figures (`mono`). */
const RulerFont = '10px mono';
const MajorHeight = 6;
const MinorHeight = 3;
/** The ruler's text ink per token table: computed once a table, never per layout. */
const RulerInks = new WeakMap<SeqTokens, number>();
function rulerInk(tokens: SeqTokens): number {
  let ink = RulerInks.get(tokens);
  if (ink === undefined) {
    ink = mutedInk(tokens, tokens.colors['seq-ruler'], [tokens.colors['seq-track'], tokens.colors['msui-bg']]);
    RulerInks.set(tokens, ink);
  }
  return ink;
}
/** The muted style of an unnumbered run: its line and fallback ticks, its note and fallback labels (alpha of the ruler's). */
const MutedLine = 0.55;
const Dash = 3;
/** The note over an unnumbered run, its inset from the run's start and the least gap to a label. */
const Note = 'not numbered';
const NoteInset = 4;
const NoteGap = 4;
const LineArgs = new Float64Array(4);

export interface RulerPrep {
  /** The numbering id asked for (`ordinal` when absent). */
  readonly id: string;
  /**
   * The system's numbering of that id; undefined: a source's (`source`, read at layout) or unknown (reported, and
   * `ordinal` is used).
   */
  readonly numbering?: NumberingOf;
  /** A numbering the system declares whose labels are a source's rows: read from the store at layout. */
  readonly source?: { readonly model: StaticModel; readonly entry: NumberingEntry };
  readonly system?: QuerySystem;
  /** `units: bases`: an `ordinal` ruler labels in bp, kb or Mb by its step (_spec/README.md "Render", Layer types). */
  readonly bases?: boolean;
  /** `unnumbered`: `mark` (the default, and for anything else) or `fallback`. */
  readonly unnumbered?: 'mark' | 'fallback';
  /** `fallback`: the numbering that labels an unnumbered run (the system's `ordinal` when `fallbackNumbering` is absent). */
  readonly fallback?: NumberingOf;
  /** A `fallbackNumbering` the system does not have: reported, and `mark` is used. */
  readonly unknownFallback?: string;
  /** The note's text with the track label it was made for (made once a label, at layout: a patch of the label remakes it). */
  noteLabel?: string;
  note?: string;
}

export interface RulerLay {
  key: string;
  system: string;
  section: string;
  box: Box;
  axis: AxisTransform;
  tokens: SeqTokens;
  font: string;
  fontSize: number;
  /** The ticks in view: their center x, whether major, and their label (empty: none drawn) with its measured width. */
  n: number;
  x: Float32Array;
  major: Uint8Array;
  labels: string[];
  widths: Float32Array;
  /** Whether the section prints a label at each group's end (`group_labels`): the ruler leaves those positions to it. */
  grouped: '' | 'above' | 'inline';
  /** Labels drawn: indices into the ticks, with the room taken. */
  shown: Uint8Array;
  /** Whether a tick (and its label) stands in an unnumbered run, in the fallback's muted style. */
  muted: Uint8Array;
  /** The unnumbered runs in view under `mark`: their left and right x (clipped to the box), pairs. */
  runs: Float32Array;
  nruns: number;
  /** The note over each run: its x (left edge, in `noteX`, NaN when it does not fit), measured width. */
  noteX: Float32Array;
  noteWidth: number;
  /** The ruler's color and its muted one, a tick's and a label's (by `muted`), made at layout: a draw allocates nothing. */
  tickColors: Float64Array;
  inkColors: Float64Array;
  /** The note's text: `not numbered`, with ` by <numbering>` where the track's label does not say it. */
  note: string;
}

export const RulerType: LayerType<RulerPrep, RulerLay> = {
  name: 'ruler',

  prepare(layer: ResolvedLayer, model: StaticModel): RulerPrep {
    const p = layer.params.numbering;
    const id = typeof p === 'string' && p ? p : 'ordinal';
    const system = layer.data ? model.system(layer.data.system) : undefined;
    const numbering = system ? findNumbering(system, id) : undefined;
    const entry = !numbering && system ? model.numberingOf(system.address, id) : undefined;
    const fallbackId = layer.params.fallbackNumbering;
    const named = typeof fallbackId === 'string' && fallbackId ? fallbackId : undefined;
    const fallback = system ? findNumbering(system, named ?? 'ordinal') : undefined;
    return {
      id,
      ...(numbering ? { numbering } : {}),
      ...(entry?.source ? { source: { model, entry } } : {}),
      ...(system ? { system } : {}),
      bases: layer.params.units === 'bases',
      unnumbered: layer.params.unnumbered === 'fallback' ? 'fallback' : 'mark',
      ...(fallback ? { fallback } : {}),
      ...(named && !fallback ? { unknownFallback: named } : {}),
    };
  },

  layout(prep, ctx: LayoutContext, data, _overlay, prev?: RulerLay): RulerLay {
    const { axis, box, config, tokens } = ctx;
    const lay: RulerLay = prev ?? {
      key: ctx.key,
      system: axis.system,
      section: sectionId(ctx.section),
      box,
      axis,
      tokens,
      font: RulerFont,
      fontSize: fontSize(RulerFont),
      n: 0,
      x: new Float32Array(64),
      major: new Uint8Array(64),
      labels: [],
      widths: new Float32Array(64),
      shown: new Uint8Array(64),
      muted: new Uint8Array(64),
      runs: new Float32Array(8),
      nruns: 0,
      noteX: new Float32Array(4),
      noteWidth: 0,
      tickColors: new Float64Array(2),
      inkColors: new Float64Array(2),
      note: Note,
      grouped: '',
    };
    lay.key = ctx.key;
    lay.system = axis.system;
    lay.section = sectionId(ctx.section);
    lay.box = box;
    lay.axis = axis;
    lay.tokens = tokens;
    lay.font = RulerFont;
    lay.fontSize = fontSize(RulerFont);
    lay.n = 0;
    lay.nruns = 0;
    lay.tickColors[0] = tokens.colors['seq-ruler'];
    lay.tickColors[1] = withAlpha(tokens.colors['seq-ruler'], MutedLine);
    // tick labels and the note hold 4.5:1 on the ruler's well and the page (`mutedInk`); a quieter label is the same ink
    const ink = rulerInk(tokens);
    lay.inkColors[0] = ink;
    lay.inkColors[1] = ink;
    lay.note = Note;
    const gl = axis.gsize > 0 ? ctx.section.params.group_labels : undefined;
    lay.grouped = axis.gsize > 0 && gl !== 'none' ? (gl === 'inline' ? 'inline' : 'above') : '';
    // a source's labels, once every tile of its system is in: no ticks while they come (pending, never ordinals); a
    // blank one is reported and labels by ordinals, as an unknown one does
    const read = prep.source
      ? readSourceNumbering(prep.source.model, prep.source.entry, data, {
          maxPositions: config.maxNumberingPositions,
          want: (w) => ctx.want(w),
        })
      : undefined;
    const fromSource = read?.status === 'held' ? read.numbering : undefined;
    // only where the group labels print this ruler's numbering (the text's: its first layer naming one, else ordinals)
    if (
      lay.grouped &&
      (prep.numbering || fromSource ? prep.id : 'ordinal') !== textNumberingId(ctx.section, prep.system)
    )
      lay.grouped = '';
    if (read?.status === 'pending') {
      fit(lay, ctx);
      return lay;
    }
    if (read?.status === 'held') for (const i of read.issues) ctx.report(i);
    else if (read?.status === 'blank') ctx.report(read.issue);
    else if (!prep.numbering && prep.id !== 'ordinal' && prep.id !== 'index')
      ctx.report(
        issueOf(RenderIssues, 'seq.render.numbering', `The system has no numbering ${prep.id}`, {
          address: ctx.key,
          hint: RenderIssues.codes['seq.render.numbering'].hint,
        }),
      );
    const numbering =
      prep.numbering ??
      (fromSource ? { numbering: fromSource, shift: 0 } : undefined) ??
      (prep.system ? findNumbering(prep.system, 'ordinal') : undefined);
    if (prep.unknownFallback !== undefined && prep.unnumbered === 'fallback')
      ctx.report(
        issueOf(
          RenderIssues,
          'seq.render.numbering.fallback',
          `The system has no numbering ${prep.unknownFallback} for the unnumbered positions: they are marked`,
          { address: ctx.key, hint: RenderIssues.codes['seq.render.numbering.fallback'].hint },
        ),
      );
    // a block's ruler in a row of many (a genome's chromosomes, _spec/README.md "Render", Blocks): ticks as close as half
    // its width, so a block a few labels wide still counts its own positions (its labels still clear each other)
    const gap = ctx.section.blocks
      ? Math.min(config.ruler.labelGap, Math.max(axis.width / 2, 24))
      : config.ruler.labelGap;
    const step = labelStep(axis.scale, gap);
    const ranges = axis.visible([]);
    const bases = prep.bases === true;
    // a numbering of its own (not a built-in, which labels every position) may not reach the whole system: the head
    // and the tail past its labelled extent are unnumbered
    const own =
      numbering && numbering.numbering.id !== 'ordinal' && numbering.numbering.id !== 'index' ? numbering : undefined;
    const fallback = prep.unnumbered === 'fallback' ? prep.fallback : undefined;
    const extent = own ? labelledExtent(own.numbering) : undefined;
    // positions of the system: `[lo, hi)` is labelled; none labelled: the whole system is unnumbered
    const lo = own && extent ? extent[0] + own.shift : Number.POSITIVE_INFINITY;
    const hi = own && extent ? extent[1] + 1 + own.shift : Number.POSITIVE_INFINITY;
    const unnumbered = (from: number, to: number) => {
      if (to <= from) return;
      if (fallback) tickRange(lay, fallback, from, to, step, false, true);
      else {
        lay.note = noteFor(prep, ctx);
        unnumberedRun(lay, from, to);
      }
    };
    for (let r = 0; r + 1 < ranges.length; r += 2) {
      const a = ranges[r];
      const b = ranges[r + 1];
      if (!own) {
        tickRange(lay, numbering, a, b, step, bases, false);
        continue;
      }
      unnumbered(a, Math.min(b, lo));
      tickRange(lay, numbering, Math.max(a, lo), Math.min(b, hi), step, bases, false);
      unnumbered(Math.max(a, hi), b);
    }
    fit(lay, ctx);
    return lay;
  },

  draw(lay, _pass: PassName, out: PassWriter, emphasis?: Emphasis) {
    drawPass(RulerDraws, lay, out, emphasis);
  },

  hit(lay, at): Hit | undefined {
    const x = at[0];
    const y = at[1];
    const { box } = lay;
    if (y < box.y || y >= box.y + box.height) return undefined;
    const position = lay.axis.positionAt(x);
    if (position < 0) return undefined;
    return { layer: lay.key, section: lay.section, position, system: lay.system };
  },
};

/** The numbering id a text section's numbers print: its first layer naming one the system has, else `ordinal`. */
function textNumberingId(section: LayoutContext['section'], system: QuerySystem | undefined): string {
  for (const t of section.tracks)
    for (const l of t.layers) {
      const id = l.params.numbering;
      if (typeof id !== 'string' || !id || id === 'ordinal') continue;
      return system && findNumbering(system, id) ? id : 'ordinal';
    }
  return 'ordinal';
}

/** The baseline of a ruler's labels in a band at `y` of `height` (what a text's line numbers share). */
export function rulerBaseline(y: number, height: number): number {
  return y + height - MajorHeight - 3;
}

/** The least step of 1, 2, 5 × 10ⁿ whose ticks are at least `gap` px apart at `scale` px a position. */
export function labelStep(scale: number, gap: number): number {
  const want = gap / Math.max(scale, 1e-9);
  const e = 10 ** Math.floor(Math.log10(want));
  for (const m of [1, 2, 5, 10]) if (m * e >= want) return m * e;
  return 10 * e;
}

/** The next tick's index, the arrays grown to hold it. */
function nextTick(lay: RulerLay): number {
  const j = lay.n++;
  if (j >= lay.x.length) {
    lay.x = grown(lay.x, j + 1);
    lay.major = grown(lay.major, j + 1);
    lay.widths = grown(lay.widths, j + 1);
    lay.shown = grown(lay.shown, j + 1);
    lay.muted = grown(lay.muted, j + 1);
  }
  return j;
}

/** The first and last labelled index of a numbering, once (a held numbering is one object); undefined: none labelled. */
const Extents = new WeakMap<Numbering, readonly [number, number] | null>();

function labelledExtent(n: Numbering): readonly [number, number] | undefined {
  let e = Extents.get(n);
  if (e === undefined) {
    let first = 0;
    while (first < n.length && n.labelAt(first) === undefined) first++;
    let last = n.length - 1;
    while (last > first && n.labelAt(last) === undefined) last--;
    e = first < n.length ? [first, last] : null;
    Extents.set(n, e);
  }
  return e ?? undefined;
}

/** The note over an unnumbered run: names the numbering where the ruler's track label does not (made once a label). */
function noteFor(prep: RulerPrep, ctx: LayoutContext): string {
  let label = '';
  for (const t of ctx.section.tracks) if (t.layers.some((l) => l.address === ctx.key)) label = t.label;
  if (prep.note === undefined || prep.noteLabel !== label) {
    prep.noteLabel = label;
    prep.note = label.toLowerCase().includes(prep.id.toLowerCase()) ? Note : `${Note} by ${prep.id}`;
  }
  return prep.note;
}

/** An unnumbered run of the positions `[a, b)` in view: its x extent clipped to the box (a pair in `runs`). */
function unnumberedRun(lay: RulerLay, a: number, b: number) {
  const axis = lay.axis;
  if (!axis.cell(a)) return;
  let x0 = Math.min(axis.sx0, axis.sx1);
  let x1 = Math.max(axis.sx0, axis.sx1);
  if (!axis.cell(b - 1)) return;
  x0 = Math.min(x0, axis.sx0, axis.sx1);
  x1 = Math.max(x1, axis.sx0, axis.sx1);
  x0 = Math.max(x0, lay.box.x);
  x1 = Math.min(x1, lay.box.x + lay.box.width);
  if (x1 <= x0) return;
  const j = 2 * lay.nruns++;
  if (j + 2 > lay.runs.length) lay.runs = grown(lay.runs, j + 2);
  lay.runs[j] = x0;
  lay.runs[j + 1] = x1;
}

/** A tick at position `p` with a label, when its cell is in view. */
function tickAt(lay: RulerLay, p: number, label: string, major: boolean, muted: boolean) {
  const axis = lay.axis;
  if (!axis.cell(p)) return;
  const x = (axis.sx0 + axis.sx1) / 2;
  if (x < lay.box.x || x > lay.box.x + lay.box.width) return;
  const j = nextTick(lay);
  lay.x[j] = x;
  lay.major[j] = major ? 1 : 0;
  lay.muted[j] = muted ? 1 : 0;
  // a group's end number stands over its last cell: the position carries that label, the tick none (text layout, Numbers)
  lay.labels[j] = lay.grouped && axis.endsGroup(axis.unitAt(x), lay.grouped === 'inline') ? '' : label;
}

/**
 * A label in bases (`formatGenomic`) by its step, then its number, kept between layouts: a pan's frame lays out the
 * same ticks again and formats none of them (formatting one parses a format and builds its string: 16 KB a frame of
 * p53's chromosome). Dropped whole past `GenomicLabelsCap` labels.
 */
const GenomicLabels = new Map<number, Map<number, string>>();
const GenomicLabelsCap = 4096;
let genomicLabelCount = 0;

function genomicLabel(n: number, step: number): string {
  let byN = GenomicLabels.get(step);
  if (!byN) {
    byN = new Map();
    GenomicLabels.set(step, byN);
  }
  let label = byN.get(n);
  if (label === undefined) {
    if (genomicLabelCount >= GenomicLabelsCap) {
      for (const m of GenomicLabels.values()) m.clear();
      genomicLabelCount = 0;
    }
    label = formatGenomic(n, step);
    byN.set(n, label);
    genomicLabelCount++;
  }
  return label;
}

/**
 * The ticks of the positions `[a, b)` in view (one piece): the labels that are multiples of the step, the insertions;
 * with `bases`, an `ordinal` ruler's labels in bp, kb or Mb by the step (another numbering keeps its own labels).
 */
function tickRange(
  lay: RulerLay,
  of: NumberingOf | undefined,
  a: number,
  b: number,
  step: number,
  bases: boolean,
  muted: boolean,
) {
  if (!of || b <= a) return;
  const { numbering, shift } = of;
  const builtIn = numbering.id === 'ordinal' || numbering.id === 'index';
  if (b - a > LoopLimit || builtIn) {
    // labels by number: every multiple of the step between the piece's first and last label
    const l0 = parseLabel(numbering.labelAt(a - shift) ?? '').int;
    const l1 = parseLabel(numbering.labelAt(b - 1 - shift) ?? '').int;
    if (l0 === undefined || l1 === undefined) return;
    const lo = Math.min(l0, l1);
    const hi = Math.max(l0, l1);
    if ((hi - lo) / step > 4 * LoopLimit) return;
    const genomic = bases && numbering.id === 'ordinal';
    for (let n = Math.ceil(lo / step) * step; n <= hi; n += step) {
      const label = genomic ? genomicLabel(n, step) : String(n);
      const i = builtIn ? n - (numbering.id === 'index' ? 0 : 1) : numbering.indexOf(label);
      const p = i + shift;
      if (i >= 0 && p >= a && p < b) tickAt(lay, p, label, true, muted);
    }
    return;
  }
  // an insertion code (100A, 111.1) a minor tick from MinorScale, never in the step's count: the step holds across
  // insertions, and an insertion never stands in for a round value the numbering lacks (an IMGT gap)
  const minors = lay.axis.scale >= MinorScale;
  for (let p = a; p < b; p++) {
    const label = numbering.labelAt(p - shift);
    if (label === undefined) continue;
    const parsed = parseLabel(label);
    if (parsed.int === undefined) tickAt(lay, p, label, true, muted);
    else if (parsed.kind === 0) {
      if (parsed.int % step === 0) tickAt(lay, p, label, true, muted);
    } else if (minors) tickAt(lay, p, label, false, muted);
  }
}

/**
 * Which labels are drawn: majors first, then minors, each only where it clears every one drawn and the data box. The
 * minors in priority order (one rule, the order decides who wins): the nearest a drawn major first, then the shorter
 * label, then by position, so the first insertion after a round value is labelled before the third, and a number
 * stands between every pair of majors where one fits. Decided at layout, never on a hover: no jitter.
 * @internal
 */
export function fit(lay: RulerLay, ctx: Pick<LayoutContext, 'config' | 'text'>) {
  // the gap a field, never an argument: a fractional number passed to a call is boxed
  FitGap[0] = Math.max(4, ctx.config.ruler.labelGap / 8);
  TakenN = 0;
  if (Taken.length < 2 * lay.n) Taken = new Float64Array(Math.max(2 * lay.n, 2 * Taken.length));
  // majors first, in tick order; the minors gathered (scratch kept between layouts: a layout allocates none of it)
  if (Minors.length < lay.n) Minors = new Int32Array(Math.max(lay.n, 2 * Minors.length));
  if (Near.length < lay.n) Near = new Float64Array(Math.max(lay.n, 2 * Near.length));
  const minors = Minors;
  const near = Near;
  let m = 0;
  for (let j = 0; j < lay.n; j++) {
    lay.shown[j] = 0;
    if (!lay.labels[j]) continue;
    lay.widths[j] = ctx.text.measure(lay.font, lay.labels[j]);
    if (lay.major[j]) placeLabel(lay, j);
    else minors[m++] = j;
  }
  if (m) {
    // each minor's distance to the nearest drawn major: one sweep each way over the ticks in order
    nearestMajors(lay, near);
    // by priority (nearest a major, then shorter, then by position): an insertion sort, the minors in view are few
    for (let a = 1; a < m; a++) {
      const j = minors[a];
      let b = a - 1;
      while (b >= 0 && before(lay, near, j, minors[b])) {
        minors[b + 1] = minors[b];
        b--;
      }
      minors[b + 1] = j;
    }
    for (let k = 0; k < m; k++) placeLabel(lay, minors[k]);
  }
  // one note at the start of each unnumbered run, past any label reaching into it, only where the run holds it
  if (lay.nruns) {
    lay.noteWidth = ctx.text.measure(lay.font, lay.note);
    if (lay.noteX.length < lay.nruns) lay.noteX = new Float32Array(lay.nruns);
    for (let r = 0; r < lay.nruns; r++)
      lay.noteX[r] = noteStart(Taken, lay.runs[2 * r], lay.runs[2 * r + 1], lay.noteWidth, TakenN);
  }
}

/** The labels drawn by a fit (start, end pairs), the minors in priority order, their distance to a major: scratch. */
let Taken = new Float64Array(128);
let TakenN = 0;
let Minors = new Int32Array(64);
let Near = new Float64Array(64);

/**
 * Each minor tick's distance to the nearest drawn major (`Infinity` with none), one sweep each way over the ticks in
 * order (the last major's index kept, never its x).
 */
function nearestMajors(lay: RulerLay, near: Float64Array): void {
  let last = -1;
  for (let j = 0; j < lay.n; j++) {
    if (lay.major[j] && lay.shown[j]) last = j;
    else near[j] = last >= 0 ? Math.abs(lay.x[j] - lay.x[last]) : Number.POSITIVE_INFINITY;
  }
  last = -1;
  for (let j = lay.n - 1; j >= 0; j--) {
    if (lay.major[j] && lay.shown[j]) last = j;
    else if (last >= 0) near[j] = Math.min(near[j], Math.abs(lay.x[j] - lay.x[last]));
  }
}

/** Whether minor `a` goes before minor `b`: nearer a drawn major, then shorter, then first. */
function before(lay: RulerLay, near: Float64Array, a: number, b: number): boolean {
  if (near[a] !== near[b]) return near[a] < near[b];
  if (lay.widths[a] !== lay.widths[b]) return lay.widths[a] < lay.widths[b];
  return a < b;
}

/** The gap a label keeps from those drawn (`fit`'s). */
const FitGap = new Float64Array(1);

/** Draws tick `j`'s label where it clears every label drawn (`Taken`) by the gap and stays within the data box. */
function placeLabel(lay: RulerLay, j: number): void {
  const gap = FitGap[0];
  const w = lay.widths[j];
  const a = lay.x[j] - w / 2;
  const b = lay.x[j] + w / 2;
  if (a < lay.box.x || b > lay.box.x + lay.box.width) return;
  for (let k = 0; k + 1 < TakenN; k += 2) if (!(b + gap <= Taken[k] || a >= Taken[k + 1] + gap)) return;
  Taken[TakenN++] = a;
  Taken[TakenN++] = b;
  lay.shown[j] = 1;
}

/**
 * Where a note of `width` starts in the run `[x0, x1]`: its inset from the run's start, moved past every label
 * (`taken`: start, end pairs in any order, a reversed axis's descending) that reaches it, NaN where it then does not fit.
 */
export function noteStart(taken: ArrayLike<number>, x0: number, x1: number, width: number, n = taken.length): number {
  let x = x0 + NoteInset;
  // a move past one label may land on another: until none reaches it (each move passes a label for good)
  for (let moved = true; moved; ) {
    moved = false;
    for (let k = 0; k + 1 < n; k += 2)
      if (taken[k + 1] + NoteGap > x && taken[k] < x + width + NoteGap && taken[k] < x1) {
        x = taken[k + 1] + NoteGap;
        moved = true;
      }
  }
  return x + width <= x1 - NoteInset ? x : Number.NaN;
}

/** Ticks on the box's bottom edge (toward the data) and labels above them; the dashed line over an unnumbered run. */
function drawTicks(lay: RulerLay, out: PassWriter) {
  const { box } = lay;
  const px = 1 / lay.axis.dpr;
  const g = out.args;
  for (let r = 0; r < lay.nruns; r++) {
    LineArgs[0] = lay.runs[2 * r];
    LineArgs[1] = box.y + box.height - px / 2;
    LineArgs[2] = lay.runs[2 * r + 1];
    LineArgs[3] = LineArgs[1];
    g[4] = NoColor;
    g[5] = lay.tickColors[1];
    g[6] = px;
    g[7] = Dash;
    out.pathArgs(LineArgs, 2, false);
  }
  for (let j = 0; j < lay.n; j++) {
    const h = lay.major[j] ? MajorHeight : MinorHeight;
    g[0] = Math.round(lay.x[j] * lay.axis.dpr) / lay.axis.dpr - px / 2;
    g[1] = box.y + box.height - h;
    g[2] = px;
    g[3] = h;
    g[4] = lay.tickColors[lay.muted[j]];
    g[5] = NoColor;
    g[6] = 0;
    g[7] = 0;
    g[8] = 0;
    out.rectArgs();
  }
}

function drawLabels(lay: RulerLay, out: PassWriter, _emphasis: Emphasis | undefined) {
  const { box } = lay;
  const g = out.args;
  for (let j = 0; j < lay.n; j++) {
    if (!lay.shown[j]) continue;
    g[0] = lay.x[j];
    g[1] = rulerBaseline(box.y, box.height);
    g[2] = lay.widths[j];
    g[3] = Number.POSITIVE_INFINITY;
    g[4] = lay.inkColors[lay.muted[j]];
    out.textArgs(lay.labels[j], AnchorMiddle, lay.font);
  }
  for (let r = 0; r < lay.nruns; r++) {
    const x = lay.noteX[r];
    if (Number.isNaN(x)) continue;
    g[0] = x;
    g[1] = rulerBaseline(box.y, box.height);
    g[2] = lay.noteWidth;
    g[3] = Number.POSITIVE_INFINITY;
    g[4] = lay.inkColors[1];
    out.textArgs(lay.note, AnchorStart, lay.font);
  }
}

const RulerDraws = passDraws<RulerLay, Emphasis | undefined>({ data: drawTicks, labels: drawLabels });
