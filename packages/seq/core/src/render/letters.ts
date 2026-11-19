/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { issueOf } from '@molstar/core/repr/tree';
import type { LettersTable, Want } from '@molstar/seqx-core/data';
import { rangesLocation } from '@molstar/seqx-core/query';
import {
  codonPairs,
  complementCode,
  type GeneticCode,
  projectBoundary,
  projectRanges,
  type StaticModel,
} from '@molstar/seqx-core/resolve';
import type { ResolvedLayer } from '@molstar/seqx-core/view';
import type { AxisTransform } from './axis.js';
import { type DensityStep, densityStep, fontSize, type RenderConfig } from './config.js';
import { Anchor, drawPass, grown, NoColor, type PassName, type PassWriter, passDraws } from './draw-lists.js';
import { RenderIssues } from './issues.js';
import { dataRegions, type LayerSource, layerSource, mergeNear, placeholderSpans, readTiles } from './layer-data.js';
import {
  type Box,
  type Emphasis,
  type Hit,
  type LayerType,
  type LayoutContext,
  type Overlay,
  sectionId,
  type StoreReader,
} from './layer-type.js';
import {
  Black,
  type Channel,
  contrastRgba,
  inkAgainst,
  lettersDefault,
  type ThemeBody,
  White,
} from '@molstar/seqx-core/theme/index';
import { channelForms, channelLegends, type LayerChannels, PlainDatum } from './channels.js';
import { drawFailed, drawSpans, hoverCell } from './marks.js';
import { inkOn, type SeqTokens } from './tokens.js';

// the anchor a number read once: a draw reads no property of an imported object (that site's feedback, cleared, made
// a warm draw deoptimize and allocate without end once the 150,000-piece scenario had run first)
const AnchorMiddle: number = Anchor.middle;

// The `letters` layer type (_spec/README.md "Render", Layer types): a system's letters, a cell a position. By density
// step: letters on their colored cells, cells only, a tick a position, then one bar a piece in view (presence needs no
// data). Letters are centered in their cells, one size a row, ink by the fill. Another system's letters are drawn
// through a path that pairs positions one to one (a letter its image's cell; a position with no letter a gap mark, a
// letter with no position a caret) or a residue to a codon (its fine places as parts, the letter in each part with
// room, three-letter codes where the narrowest part in view has room, `code`); any other path is refused. A nucleotide
// read against its strand on screen shows its complement (`complement: auto`). A cell's fill is its `color` channel
// (the alphabet's token scheme by default: `classes`, `nucleotide-classes`), its letter's ink the `ink` channel (`auto`:
// follows its fill at 4.5:1; a themed ink held at 3:1 against its fill).

/** A letters layer's alphabet (the system's); the rows and the logo read the same. */
export type Alphabet = 'protein' | 'dna' | 'rna';

/** The color channel's default body and environment by alphabet, made once (a layout makes none to read a kept channel). */
export const LettersColorFallback: Record<Alphabet, ThemeBody> = {
  protein: { kind: 'scheme', params: { name: lettersDefault('protein') } },
  dna: { kind: 'scheme', params: { name: lettersDefault('dna') } },
  rna: { kind: 'scheme', params: { name: lettersDefault('rna') } },
};
const LettersEnv: Record<Alphabet, { alphabet: Alphabet }> = {
  protein: { alphabet: 'protein' },
  dna: { alphabet: 'dna' },
  rna: { alphabet: 'rna' },
};

/** How a letters layer reaches its section: its own letters, a one-to-one path, a codon step, or none drawn. */
type Through = 'own' | 'one' | 'codon' | 'refused';

export interface LettersPrep {
  readonly src?: LayerSource;
  readonly alphabet: 'protein' | 'dna' | 'rna';
  /** Its channels: `color` (cells) and `ink` (letters). */
  readonly channels: LayerChannels;
  readonly through: Through;
  readonly complement: 'auto' | 'always' | 'never';
  readonly code: 'auto' | 'one' | 'three';
  /** Across a codon step: the section's own letters (each codon's status) and the map's genetic code. */
  readonly sectionSrc?: LayerSource;
  readonly genetic?: GeneticCode;
}

export interface LettersLay {
  key: string;
  system: string;
  /** The section's id (hits name it). */
  section: string;
  box: Box;
  axis: AxisTransform;
  step: DensityStep;
  /** How it reaches its section, and the data's system (a residue's extent is on it). */
  through: Through;
  dataSystem: string;
  src?: LayerSource;
  /**
   * Cells in view: position (the data's: a residue a codon wide is a cell a part), letter code (0: not loaded), fill,
   * ink, snapped x extents, whether the letter is drawn (its part has room).
   */
  n: number;
  pos: Int32Array;
  code: Uint16Array;
  /** The letter codes laid out (ASCII), 1 where one was: what its legend lists. */
  seen: Uint8Array;
  fill: Uint32Array;
  ink: Uint32Array;
  x0: Float32Array;
  x1: Float32Array;
  show: Uint8Array;
  font: string;
  fontSize: number;
  /** Three-letter codes (a residue a codon wide whose narrowest part in view has room). */
  three: boolean;
  /** The letters' baseline, made at layout (the draw passes it through the writer's `args`). */
  glyphY: number;
  /** `x0, x1, cut` triples: the bars of the coarse steps, the pending and failed placeholders, the gap marks. */
  bars: number[];
  pending: number[];
  failed: number[];
  gaps: number[];
  /** A caret's x a letter with no position (an insertion, one to one). */
  carets: number[];
  /**
   * `x0, x1` pairs: the parts of residues their codon does not translate to (EC-16), marked under the cell in the
   * conflict ink; the cell keeps its residue's fill, so a lone conflict reads as a mark, never a bar.
   */
  conflicts: number[];
  tokens: SeqTokens;
  hoverWidth: number;
  hoverOffset: number;
  /** What the hover's squircle reads (its rule: `squircleRatio`, `squircleMax`, `squircleMin`). */
  config: RenderConfig;
  /** The color channel this layout read (its legend), the datum it read a cell through, a themed ink by fill. */
  color?: Channel;
  datum: PlainDatum;
  inks: Map<number, number>;
  /** A chip's label and its font (a one-position block): drawn as a pill over its one cell; '' on any other strip. */
  chip: string;
  chipFont: string;
  /** The chip label's baseline, its measured width and its ink, made at layout. */
  chipY: number;
  chipWidth: number;
  chipInk: number;
}

/** The letter code at a position of a table, or 0. */
function codeAt(t: LettersTable, p: number): number {
  const i = p - t.start;
  if (i < 0 || i >= t.letters.length) return 0;
  if (typeof t.letters === 'string') return t.letters.charCodeAt(i);
  return t.alphabet ? t.alphabet.charCodeAt(t.letters[i]) || 0 : 0;
}

/** The letter at a position of sorted tables, or 0. */
function letterIn(tables: readonly LettersTable[], p: number): number {
  for (let k = 0; k < tables.length; k++) {
    const c = codeAt(tables[k], p);
    if (c) return c;
  }
  return 0;
}

/** Three-letter codes by one-letter code (IUPAC), the stop as `Ter`. */
const ThreeLetter: Readonly<Record<string, string>> = {
  A: 'Ala',
  R: 'Arg',
  N: 'Asn',
  D: 'Asp',
  C: 'Cys',
  Q: 'Gln',
  E: 'Glu',
  G: 'Gly',
  H: 'His',
  I: 'Ile',
  L: 'Leu',
  K: 'Lys',
  M: 'Met',
  F: 'Phe',
  P: 'Pro',
  S: 'Ser',
  T: 'Thr',
  W: 'Trp',
  Y: 'Tyr',
  V: 'Val',
  U: 'Sec',
  O: 'Pyl',
  B: 'Asx',
  Z: 'Glx',
  J: 'Xle',
  X: 'Xaa',
  '*': 'Ter',
};
// dense (no holes: a read of one allocates nothing), `Xaa` for a code with none
const ThreeByCode: string[] = Array.from({ length: 128 }, () => 'Xaa');
for (const [k, v] of Object.entries(ThreeLetter)) {
  ThreeByCode[k.charCodeAt(0)] = v;
  ThreeByCode[k.toLowerCase().charCodeAt(0)] = v;
}

/**
 * How a path places letters: one to one when no section position takes two data positions (identity, offset,
 * segments, pairwise, an alignment row), a residue to its codon across a codon step from residues; else refused.
 */
function throughOf(src: LayerSource | undefined): Through {
  const p = src?.path;
  if (!p) return 'own';
  if (p.step === 1) return 'codon';
  if (p.step === -1) return 'refused';
  const runs = p.linear.runs();
  if (!runs) return 'one';
  const targets: number[] = [];
  for (const r of runs) targets.push(r.dir === 1 ? r.to : r.to - r.length + 1, r.length);
  const order = targets.map((_, i) => i).filter((i) => i % 2 === 0);
  order.sort((a, b) => targets[a] - targets[b]);
  for (let i = 1; i < order.length; i++)
    if (targets[order[i]] < targets[order[i - 1]] + targets[order[i - 1] + 1]) return 'refused';
  return 'one';
}

const Refused =
  'Letters through this path are not drawn: a nucleotide on a protein’s axis is a third of a residue and an intron or UTR nucleotide has no place there; a many-to-one path gives a cell several letters';

export const LettersType: LayerType<LettersPrep, LettersLay> = {
  name: 'letters',

  prepare(layer: ResolvedLayer, model: StaticModel): LettersPrep {
    const src = layerSource(model, layer);
    const system = src ? model.systems.get(src.dataSystem) : undefined;
    const through = throughOf(src);
    const complement = layer.params.complement;
    const code = layer.params.code;
    let sectionSrc: LayerSource | undefined;
    let genetic: GeneticCode | undefined;
    if (through === 'codon' && src?.path) {
      const to = src.path.to;
      sectionSrc = layerSource(model, {
        ...layer,
        data: { kind: 'address', address: to, system: to, shape: 'letters' },
        placed: true,
        path: undefined,
      } as ResolvedLayer);
      for (const e of src.path.edges) if (e.step !== 0) genetic = model.maps.get(e.map)?.code;
    }
    return {
      ...(src ? { src } : {}),
      alphabet: system?.alphabet ?? 'protein',
      channels: channelForms(layer.theme),
      through,
      complement: complement === 'always' || complement === 'never' ? complement : 'auto',
      code: code === 'one' || code === 'three' ? code : 'auto',
      ...(sectionSrc ? { sectionSrc } : {}),
      ...(genetic ? { genetic } : {}),
    };
  },

  layout(prep, ctx: LayoutContext, data: StoreReader, _overlay: Overlay, prev?: LettersLay): LettersLay {
    const { axis, box, config, tokens } = ctx;
    const lay: LettersLay = prev ?? {
      key: ctx.key,
      system: axis.system,
      section: sectionId(ctx.section),
      box,
      axis,
      step: 'cells',
      through: prep.through,
      dataSystem: prep.src?.dataSystem ?? axis.system,
      n: 0,
      pos: new Int32Array(64),
      code: new Uint16Array(64),
      seen: new Uint8Array(128),
      fill: new Uint32Array(64),
      ink: new Uint32Array(64),
      x0: new Float32Array(64),
      x1: new Float32Array(64),
      show: new Uint8Array(64),
      font: '',
      fontSize: 0,
      three: false,
      glyphY: 0.5,
      bars: [],
      pending: [],
      failed: [],
      gaps: [],
      carets: [],
      conflicts: [],
      tokens,
      hoverWidth: config.outline.hover,
      hoverOffset: config.emphasis.hoverOffset,
      config,
      datum: new PlainDatum(),
      inks: new Map(),
      chip: '',
      chipFont: '',
      chipY: 0,
      chipWidth: 0,
      chipInk: 0,
    };
    lay.key = ctx.key;
    lay.system = axis.system;
    lay.section = sectionId(ctx.section);
    lay.through = prep.through;
    lay.dataSystem = prep.src?.dataSystem ?? axis.system;
    if (prep.src) lay.src = prep.src;
    lay.hoverWidth = config.outline.hover;
    lay.hoverOffset = config.emphasis.hoverOffset;
    lay.config = config;
    lay.box = box;
    lay.axis = axis;
    lay.tokens = tokens;
    lay.chip = ctx.chip ?? '';
    lay.chipFont = config.labelFont;
    lay.chipY = lay.chip ? box.y + box.height / 2 + fontSize(config.labelFont) * 0.35 : 0;
    lay.chipWidth = lay.chip ? ctx.text.measure(config.labelFont, lay.chip) : 0;
    lay.n = 0;
    lay.seen.fill(0);
    lay.three = false;
    lay.bars.length = 0;
    lay.pending.length = 0;
    lay.failed.length = 0;
    lay.gaps.length = 0;
    lay.carets.length = 0;
    lay.conflicts.length = 0;
    const src = prep.src;
    if (!src) return lay;
    if (prep.through === 'refused') {
      ctx.report(
        issueOf(RenderIssues, 'seq.render.layer.letters', Refused, {
          address: ctx.key,
          hint: RenderIssues.codes['seq.render.layer.letters'].hint,
        }),
      );
      return lay;
    }
    // the layer's density follows its own cell: a residue a codon wide is three positions
    const width = prep.through === 'codon' ? 3 : 1;
    const step = densityStep(axis.scale * width, axis.dpr, config);
    lay.step = step;
    const vis = axis.visible([]);
    // the positions in view exactly (walked one by one), read merged across near gaps (one read over short pieces)
    const regions = prep.through === 'own' ? vis : dataRegions(src, axis, false, false);
    if (step === 'bars' || step === 'density') {
      if (prep.through === 'own') for (let i = 0; i + 1 < vis.length; i += 2) axis.spans(vis[i], vis[i + 1], lay.bars);
      else placeholderSpans(src, axis, regions, lay.bars);
      return lay;
    }

    const color = prep.channels.get('color', LettersColorFallback[prep.alphabet], ctx, LettersEnv[prep.alphabet]);
    lay.color = color;
    const inkForm = prep.channels.forms.ink;
    const inkChannel = inkForm === undefined || inkForm === 'auto' ? undefined : prep.channels.get('ink', 'auto', ctx);
    if (lay.inks.size > 4096) lay.inks.clear();
    const datum = lay.datum;
    const missing = tokens.colors['seq-missing'];
    const read = readTiles(ctx, data, src, mergeNear(regions), 1);
    const tables = read.entries.map((e) => e.tables.letters).filter((t): t is LettersTable => !!t);
    tables.sort((a, b) => a.start - b.start);
    // the section's own letters, for each codon's status: nothing said before they are held
    const own =
      prep.sectionSrc && step !== 'dots'
        ? readTiles(ctx, data, prep.sectionSrc, mergeNear(vis), 1)
            .entries.map((e) => e.tables.letters)
            .filter((t): t is LettersTable => !!t)
        : undefined;
    let count = 0;
    for (let i = 0; i + 1 < regions.length; i += 2) count += (regions[i + 1] - regions[i]) * (width === 3 ? 3 : 1);
    lay.pos = grown(lay.pos, count);
    lay.code = grown(lay.code, count);
    lay.fill = grown(lay.fill, count);
    lay.ink = grown(lay.ink, count);
    lay.x0 = grown(lay.x0, count);
    lay.x1 = grown(lay.x1, count);
    lay.show = grown(lay.show, count);
    const nucleotides = prep.alphabet !== 'protein';
    const rna = prep.alphabet === 'rna';
    const onOwn = prep.through === 'own';
    const parts = onOwn ? OwnCell : Parts;
    const dpr = axis.dpr;
    for (let i = 0; i + 1 < regions.length; i += 2) {
      for (let p = regions[i]; p < regions[i + 1]; p++) {
        let code = letterIn(tables, p);
        if (!code) continue;
        // where the letter lands: its own cell (set in place: a length set to 0 lets an array's storage go, made
        // again a letter), or its image's parts through the path
        if (onOwn) {
          parts[0] = p;
          parts[1] = p + 1;
        } else projectRanges(src.path!, OnePosition(p), parts);
        if (!parts.length) continue;
        // a nucleotide read against its strand on screen shows its complement: the axis's direction times its run's
        if (nucleotides && prep.complement !== 'never') {
          const run = prep.through === 'one' ? runDir(src, p) : 1;
          if (prep.complement === 'always' || axis.dir * run === -1) code = complementCode(code, rna);
        }
        if (code < 128) lay.seen[code] = 1;
        const ch = String.fromCharCode(code);
        datum.letter = ch;
        datum.start = p;
        datum.end = p + 1;
        const c = color.constant ?? color.at(datum);
        const fill = typeof c === 'number' ? c : missing;
        // a residue its codon does not translate to, once the section's letters are held (EC-16): marked
        const bad = !!own?.length && prep.through === 'codon' && !codonOk(src, own, p, ch, prep.genetic);
        const ink = inkOf(inkChannel, datum, fill, tokens, lay.inks);
        for (let k = 0; k + 1 < parts.length; k += 2) {
          axis.piecesIn(parts[k], parts[k + 1]);
          const ka = axis.ka;
          const kb = axis.kb;
          for (let a = ka; a <= kb; a++) {
            if (!axis.spanIn(a, parts[k], parts[k + 1])) continue;
            const j = lay.n++;
            if (j >= lay.pos.length) {
              lay.pos = grown(lay.pos, j + 1);
              lay.code = grown(lay.code, j + 1);
              lay.fill = grown(lay.fill, j + 1);
              lay.ink = grown(lay.ink, j + 1);
              lay.x0 = grown(lay.x0, j + 1);
              lay.x1 = grown(lay.x1, j + 1);
              lay.show = grown(lay.show, j + 1);
            }
            lay.pos[j] = p;
            lay.code[j] = code;
            lay.fill[j] = fill;
            lay.ink[j] = ink;
            // snapped here, not by a call: a number passed to a call that is not inlined is boxed
            lay.x0[j] = Math.round(axis.sx0 * dpr) / dpr;
            lay.x1[j] = Math.round(axis.sx1 * dpr) / dpr;
            if (bad) lay.conflicts.push(lay.x0[j], lay.x1[j]);
          }
        }
      }
    }
    if (prep.through === 'one') caretsAt(src, axis, regions, box, lay.carets);
    placeholderSpans(src, axis, read.pending, lay.pending);
    placeholderSpans(src, axis, read.failed, lay.failed);
    // one to one: a section position in view that no letter lands on is an empty cell (the gap mark)
    if (prep.through === 'one') gapMarks(src, axis, vis, lay.gaps);

    // a chip's label ink: the table's, else black or white, at 4.5:1 against its pill (a label is text, not a letter)
    if (lay.chip && lay.n) lay.chipInk = chipInk(lay.fill[0], tokens);
    // one size a row: the largest of the scale whose advance fits the cell and whose height fits the row
    lay.fontSize = 0;
    if (step === 'letters') {
      const cell = axis.scale * width;
      for (const s of config.letterSizes)
        if (s * config.monoAdvance <= cell - 2 && s <= box.height * 0.75) lay.fontSize = s;
      // the step starts where the smallest size fits the cell: a short row still shows its letters
      if (!lay.fontSize) lay.fontSize = config.letterSizes[0] ?? 9;
      lay.font = `${lay.fontSize}px mono`;
      lay.glyphY = box.y + box.height / 2 + lay.fontSize * 0.35;
      // each part names its residue where it has room (EC-117); three-letter codes where the narrowest part in view
      // has room for them (`code: auto`), or as the writer says
      const advance = lay.fontSize * config.monoAdvance;
      let narrowest = Infinity;
      for (let j = 0; j < lay.n; j++) {
        const w = lay.x1[j] - lay.x0[j];
        lay.show[j] = w >= advance + 1 ? 1 : 0;
        if (w < narrowest) narrowest = w;
      }
      if (prep.through === 'codon' && prep.alphabet === 'protein')
        lay.three = prep.code === 'three' || (prep.code === 'auto' && narrowest >= 3 * advance + 2);
      if (lay.three) for (let j = 0; j < lay.n; j++) lay.show[j] = lay.x1[j] - lay.x0[j] >= 3 * advance + 1 ? 1 : 0;
    }
    return lay;
  },

  draw(lay, pass: PassName, out: PassWriter, emphasis?: Emphasis) {
    drawPass(LettersDraws, lay, out, emphasis);
  },

  legend(lay, options) {
    // the letters drawn in view: the key lists only what the data holds (all of its letters before any is drawn)
    let painted = '';
    for (let c = 0; c < lay.seen.length; c++) if (lay.seen[c]) painted += String.fromCharCode(c);
    return channelLegends([lay.color], {
      ...options,
      layer: lay.key,
      ...(painted ? { painted: painted.toUpperCase() } : {}),
    });
  },

  hit(lay, at): Hit | undefined {
    const x = at[0];
    const y = at[1];
    const { box } = lay;
    if (y < box.y || y >= box.y + box.height) return undefined;
    const position = lay.axis.positionAt(x);
    if (position < 0) return undefined;
    if (lay.through !== 'codon' && lay.through !== 'one')
      return { layer: lay.key, section: lay.section, position, system: lay.system };
    // a letter through a path: its data element on `extent`, its parts on the section's system as `hull`
    for (let j = 0; j < lay.n; j++) {
      if (x < lay.x0[j] || x >= lay.x1[j]) continue;
      const p = lay.pos[j];
      const parts = projectRanges(lay.src!.path!, [p, p + 1]);
      let lo = Infinity;
      let hi = -Infinity;
      for (let k = 0; k + 1 < parts.length; k += 2) {
        lo = Math.min(lo, parts[k]);
        hi = Math.max(hi, parts[k + 1]);
      }
      return {
        layer: lay.key,
        section: lay.section,
        position,
        system: lay.system,
        extent: rangesLocation(lay.dataSystem, [p, p + 1]),
        ...(lo < hi ? { hull: [lo, hi] as const } : {}),
      };
    }
    return { layer: lay.key, section: lay.section, position, system: lay.system };
  },
};

/** The parts of one letter's image (kept: a layout runs to its end). */
const Parts: number[] = [];
/** A letter's own cell, `[p, p + 1]`. */
const OwnCell = [0, 1];
const One = [0, 1];
function OnePosition(p: number): number[] {
  One[0] = p;
  One[1] = p + 1;
  return One;
}

/**
 * Letters with no position beside the letters in view (EC-29): where a region in view ends (or starts) at data
 * positions the path places nowhere, a caret at their boundary, when it is in view; one a boundary. Read from the
 * regions' edges, so how the reads are merged never changes them.
 */
function caretsAt(src: LayerSource, axis: AxisTransform, regions: readonly number[], box: Box, out: number[]) {
  const path = src.path!;
  for (let i = 0; i < regions.length; i++) {
    // the position after a region's end (odd i) or before its start (even i)
    const p = i % 2 ? regions[i] : regions[i] - 1;
    if (p < 0 || (src.length !== undefined && p >= src.length)) continue;
    if (projectRanges(path, OnePosition(p), Parts).length) continue;
    const b = projectBoundary(path, i % 2 ? p : p + 1);
    if (b === undefined) continue;
    const x = axis.x(b);
    if (x >= box.x && x <= box.x + box.width && !out.includes(x)) out.push(x);
  }
}

/**
 * The letters a letters layer shows at section positions, in the order given (_spec/README.md "Render", Layer types:
 * what the viewer's `selection.copyLetters` copies): complemented as drawn on an axis of direction `dir`, or with
 * `stored` the letters as stored. A position no letter lands on (a gap through a path) gives none; letters across a
 * codon step are not read (none). What is not held yet is in `wants`: the text is then empty (load it and ask again).
 */
export function shownLetters(
  prep: LettersPrep,
  dir: 1 | -1,
  data: StoreReader,
  positions: readonly number[],
  stored = false,
): { readonly text: string; readonly wants: readonly Want[] } {
  const src = prep.src;
  if (!src || (prep.through !== 'own' && prep.through !== 'one') || !positions.length) return { text: '', wants: [] };
  // each section position's data position (one to one through the path read back), −1 for none
  const at = new Float64Array(positions.length);
  let lo = Infinity;
  let hi = -Infinity;
  for (let i = 0; i < positions.length; i++) {
    let d = positions[i];
    if (prep.through === 'one') {
      const r = src.back ? projectRanges(src.back, OnePosition(d), Parts) : [];
      d = r.length ? r[0] : -1;
    }
    at[i] = d;
    if (d < 0) continue;
    if (d < lo) lo = d;
    if (d > hi) hi = d;
  }
  if (lo > hi) return { text: '', wants: [] };
  const want: Want = { source: src.source, system: src.system, start: lo, end: hi + 1, level: 1, priority: 'view' };
  if (data.coverage(src.source) === undefined) return { text: '', wants: [want] };
  const keys = data.keysOf(want);
  // none listed past the store's bound on a want's tiles: not held, as a tile not loaded
  if (!keys.length && data.tooMany?.(want)) return { text: '', wants: [want] };
  const tables: LettersTable[] = [];
  for (const key of keys) {
    const e = data.get(key);
    if (e?.status === 'loaded') {
      if (e.tables.letters) tables.push(e.tables.letters);
    } else if (e?.status !== 'empty' && e?.status !== 'failed') return { text: '', wants: [want] };
  }
  tables.sort((a, b) => a.start - b.start);
  const nucleotides = prep.alphabet !== 'protein';
  const rna = prep.alphabet === 'rna';
  let text = '';
  for (let i = 0; i < at.length; i++) {
    const d = at[i];
    let code = d < 0 ? 0 : letterIn(tables, d);
    if (!code) continue;
    // a nucleotide read against its strand on screen shows its complement, as the layout draws it
    if (!stored && nucleotides && prep.complement !== 'never') {
      const run = prep.through === 'one' ? runDir(src, d) : 1;
      if (prep.complement === 'always' || dir * run === -1) code = complementCode(code, rna);
    }
    text += String.fromCharCode(code);
  }
  return { text, wants: [] };
}

/** The direction of the run a one-to-one path carries a position along. */
function runDir(src: LayerSource, p: number): 1 | -1 {
  const r = src.path!.linear.forward(p, p + 1)[0];
  return r ? r.dir : 1;
}

/** Whether residue `k` pairs with its codon in the section's letters (the codon read on its strand). */
function codonOk(
  src: LayerSource,
  own: readonly LettersTable[],
  k: number,
  residue: string,
  genetic: GeneticCode | undefined,
): boolean {
  let codon = '';
  for (const r of src.path!.linear.forward(3 * k, 3 * k + 3))
    for (let i = 0; i < r.length; i++) {
      const c = letterIn(own, r.to + r.dir * i);
      if (!c) return true;
      codon += String.fromCharCode(r.dir === 1 ? c : complementCode(c));
    }
  if (codon.length !== 3) return true;
  return genetic ? codonPairs(residue, codon, k === 0, genetic) : codonPairs(residue, codon, k === 0);
}

/** The section positions in view no letter lands on, through a one-to-one path: `x0, x1, cut` triples into `out`. */
function gapMarks(src: LayerSource, axis: AxisTransform, vis: readonly number[], out: number[]) {
  const back = src.back!;
  for (let i = 0; i + 1 < vis.length; i += 2) {
    let s = vis[i];
    for (let p = vis[i]; p <= vis[i + 1]; p++) {
      const mapped = p < vis[i + 1] && back.linear.forward(p, p + 1).length > 0;
      if (mapped || p === vis[i + 1]) {
        if (p > s) axis.spans(s, p, out);
        s = p + 1;
      }
    }
  }
}

// Each pass draws in a function of its own, called through a table (`drawPass`: nothing allocated once warm).

/** The data pass: bars, placeholders, gap marks, cells or ticks, the letters on them, carets. */
function drawCells(lay: LettersLay, out: PassWriter, _emphasis: Emphasis | undefined) {
  const { box, tokens } = lay;
  drawSpans(out, lay.bars, box.y, box.height, tokens.colors['seq-sequence']);
  drawSpans(out, lay.pending, box.y, box.height, tokens.colors['seq-pending']);
  drawFailed(out, lay.failed, box.y, box.height, tokens);
  // a path's gap marks and carets, and three-letter codes, each in a function of its own: this one stays small enough to
  // have its calls inlined (a call that is not boxes the numbers it is passed); every number through `args`
  if (lay.gaps.length || lay.carets.length) drawPathMarks(lay, out);
  if (lay.chip && lay.n) {
    drawChip(lay, out);
    return;
  }
  // a tick a position at `dots`, else the cell
  const tick = Math.max(2, Math.round(box.height * 0.3));
  const dots = lay.step === 'dots';
  const y = dots ? box.y + (box.height - tick) / 2 : box.y;
  const h = dots ? tick : box.height;
  const px = 1 / lay.axis.dpr;
  // a letter cell's hairline (`seq-cell-edge`: dark only, none in light), so dim dark cells read apart from the page
  const edge = lay.step === 'letters' ? tokens.colors['seq-cell-edge'] : NoColor;
  const g = out.args;
  for (let j = 0; j < lay.n; j++) {
    g[0] = lay.x0[j];
    g[1] = y;
    g[2] = Math.max(lay.x1[j] - lay.x0[j], px);
    g[3] = h;
    g[4] = lay.fill[j];
    g[5] = NoColor;
    g[6] = 0;
    g[7] = 0;
    g[8] = 0;
    out.rectArgs();
    if (edge !== NoColor) {
      g[5] = edge;
      g[6] = px;
      cellEdge(out);
    }
  }
  if (lay.conflicts.length) drawConflicts(lay, out);
  if (lay.step !== 'letters' || !(lay.fontSize > 0) || !lay.n) return;
  if (lay.three) {
    drawThree(lay, out);
    return;
  }
  g[1] = lay.glyphY;
  out.beginGlyphsArgs(lay.font, NoColor);
  for (let j = 0; j < lay.n; j++) {
    if (!lay.show[j]) continue;
    g[0] = (lay.x0[j] + lay.x1[j]) / 2;
    g[4] = lay.ink[j];
    out.letterArgs(lay.code[j]);
  }
  out.endGlyphs();
}

/**
 * A letter cell's hairline (`seq-cell-edge`) over the cell `out.args` holds (x, y, w, h; its color in `[5]`, a device
 * pixel in `[6]`: read from the array, a number passed to a call that is not inlined is boxed): drawn inside the cell,
 * a stroke a device pixel wide on the rect inset by half of one, so it covers whole device pixels of its own cell and a
 * neighbour's fill never paints over it (a stroke on the edge itself straddles two pixels). Shared with the rows layer.
 */
export function cellEdge(out: PassWriter) {
  const g = out.args;
  const px = g[6];
  if (g[2] <= 2 * px || g[3] <= 2 * px) return;
  g[0] += px / 2;
  g[1] += px / 2;
  g[2] -= px;
  g[3] -= px;
  g[4] = NoColor;
  out.rectArgs();
}

/** A one-to-one path's gap marks (a thin mark across the middle of an empty cell) and carets (a letter with no position). */
function drawPathMarks(lay: LettersLay, out: PassWriter) {
  const { box, tokens } = lay;
  const px = 1 / lay.axis.dpr;
  const g = out.args;
  for (let k = 0; k < 2; k++) {
    const gaps = k === 0;
    const list = gaps ? lay.gaps : lay.carets;
    const step = gaps ? 3 : 1;
    for (let i = 0; i + step - 1 < list.length; i += step) {
      g[0] = gaps ? list[i] : list[i] - px;
      g[1] = gaps ? box.y + box.height / 2 - 0.5 : box.y;
      g[2] = gaps ? Math.max(list[i + 1] - list[i], 1) : 2 * px;
      g[3] = gaps ? 1 : box.height;
      g[4] = tokens.colors[gaps ? 'seq-missing' : 'msui-text'];
      g[5] = NoColor;
      g[6] = 0;
      g[7] = 0;
      g[8] = 0;
      out.rectArgs();
    }
  }
}

/** A conflict's mark: a bar along the bottom of the cell, a fifth of its height (at least 2 px), in the conflict ink. */
function drawConflicts(lay: LettersLay, out: PassWriter) {
  const { box, tokens } = lay;
  const h = Math.max(2, Math.round(box.height * 0.2));
  const g = out.args;
  const c = lay.conflicts;
  for (let i = 0; i + 1 < c.length; i += 2) {
    g[0] = c[i];
    g[1] = box.y + box.height - h;
    g[2] = Math.max(c[i + 1] - c[i], 1);
    g[3] = h;
    g[4] = tokens.colors['seq-conflict'];
    g[5] = NoColor;
    g[6] = 0;
    g[7] = 0;
    g[8] = 0;
    out.rectArgs();
  }
}

/**
 * A chip's pill (_spec/README.md "Render", Blocks): its one cell's fill as a pill a pixel inside the cell, its label
 * centred on it in the ink that holds 4.5:1 against the fill, faded at the pill's width (the label is elided to its
 * slot: it fits); the letter cells' hairline (`seq-cell-edge`, dark only) round it, as round a cell.
 */
function drawChip(lay: LettersLay, out: PassWriter) {
  const { box } = lay;
  const g = out.args;
  const h = Math.max(box.height - 2, 1);
  const edge = lay.tokens.colors['seq-cell-edge'];
  g[0] = lay.x0[0] + 1;
  g[1] = box.y + 1;
  g[2] = lay.x1[0] - lay.x0[0] - 2;
  g[3] = h;
  g[4] = lay.fill[0];
  g[5] = edge;
  g[6] = edge === NoColor ? 0 : 1 / lay.axis.dpr;
  g[7] = h / 2;
  g[8] = 3;
  out.rectArgs();
  g[0] = (lay.x0[0] + lay.x1[0]) / 2;
  g[1] = lay.chipY;
  g[2] = lay.chipWidth;
  g[3] = lay.x1[0] - lay.x0[0];
  g[4] = lay.chipInk;
  out.textArgs(lay.chip, AnchorMiddle, lay.chipFont);
}

/** A chip's label ink against its fill: the tokens' ink when it holds 4.5:1, else black or white, the higher. */
function chipInk(fill: number, tokens: SeqTokens): number {
  const ink = inkOn(fill, tokens);
  if (contrastRgba(ink, fill) >= 4.5) return ink;
  return contrastRgba(Black, fill) >= contrastRgba(White, fill) ? Black : White;
}

/** Three-letter codes, centred in their parts. */
function drawThree(lay: LettersLay, out: PassWriter) {
  const g = out.args;
  for (let j = 0; j < lay.n; j++) {
    if (!lay.show[j]) continue;
    g[0] = (lay.x0[j] + lay.x1[j]) / 2;
    g[1] = lay.glyphY;
    g[2] = 0;
    g[3] = lay.x1[j] - lay.x0[j];
    g[4] = lay.ink[j];
    out.textArgs(ThreeByCode[lay.code[j] & 127], AnchorMiddle, lay.font);
  }
}

/**
 * The hover pass: a letter through a path outlines each part of its element (a residue's codon, split by an intron
 * too); a chip, its pill; else the hovered position's cell. Each a squircle (_spec/README.md "Render", Shapes: every
 * emphasis a squircle), concentric: a part's ends round where the element really ends (its outermost parts' outer
 * sides), square where an intron cuts it; square below `squircleMin`; a chip's round as its pill.
 */
function drawHover(lay: LettersLay, out: PassWriter, emphasis: Emphasis | undefined) {
  const h = emphasis?.hover;
  if (!h || h.layer !== lay.key) return;
  if (lay.chip !== '' && lay.n > 0) {
    drawChipHover(lay, out, emphasis);
    return;
  }
  const part = lay.through === 'own' ? undefined : h.extent?.parts[0];
  if (!part || part.system !== lay.dataSystem || !part.ranges.length) {
    hoverCell(lay, out, emphasis);
    return;
  }
  const p = part.ranges[0];
  const d = lay.hoverWidth / 2 + lay.hoverOffset;
  const g = out.args;
  const { box } = lay;
  // the element's outer sides: its parts' leftmost and rightmost x (the sides between parts are cut)
  let lo = Infinity;
  let hi = -Infinity;
  for (let j = 0; j < lay.n; j++) {
    if (lay.pos[j] !== p) continue;
    if (lay.x0[j] < lo) lo = lay.x0[j];
    if (lay.x1[j] > hi) hi = lay.x1[j];
  }
  const c = lay.config;
  const bh = box.height;
  for (let j = 0; j < lay.n; j++) {
    if (lay.pos[j] !== p) continue;
    const w = lay.x1[j] - lay.x0[j];
    const r =
      w < c.squircleMin || bh < c.squircleMin ? 0 : Math.min(bh * c.squircleRatio, c.squircleMax, w / 2, bh / 2);
    const round = (lay.x0[j] <= lo ? 1 : 0) | (lay.x1[j] >= hi ? 2 : 0);
    g[0] = lay.x0[j] - d;
    g[1] = box.y - d;
    g[2] = w + 2 * d;
    g[3] = bh + 2 * d;
    g[4] = NoColor;
    g[5] = lay.tokens.colors['seq-hover'];
    g[6] = lay.hoverWidth;
    g[7] = r > 0 && round !== 0 ? r + d : 0;
    g[8] = round;
    out.rectArgs();
  }
}

/** A chip's hover: outside its pill by the offset, concentric with it (the pill's radius plus the gap to it). */
function drawChipHover(lay: LettersLay, out: PassWriter, emphasis: Emphasis | undefined) {
  const h = emphasis!.hover!;
  if (h.system !== lay.axis.system || h.position === undefined || !lay.axis.cell(h.position)) return;
  const { box } = lay;
  // the pill is a pixel inside the cell (`drawChip`); the outline `d` outside the cell
  const d = lay.hoverWidth / 2 + lay.hoverOffset;
  const ph = Math.max(box.height - 2, 1);
  const g = out.args;
  g[0] = lay.x0[0] - d;
  g[1] = box.y - d;
  g[2] = lay.x1[0] - lay.x0[0] + 2 * d;
  g[3] = box.height + 2 * d;
  g[4] = NoColor;
  g[5] = lay.tokens.colors['seq-hover'];
  g[6] = lay.hoverWidth;
  g[7] = ph / 2 + 1 + d;
  g[8] = 3;
  out.rectArgs();
}

const LettersDraws = passDraws<LettersLay, Emphasis | undefined>({ data: drawCells, hover: drawHover });

/**
 * A letter's ink: `auto` follows its fill (the table's inks at 4.5:1); a themed ink is held at 3:1 against its fill,
 * memoized by the pair (a fill's ink is the same every cell).
 */
function inkOf(
  channel: Channel | undefined,
  d: PlainDatum,
  fill: number,
  tokens: SeqTokens,
  memo: Map<number, number>,
): number {
  if (!channel) return inkOn(fill, tokens);
  const c = channel.constant ?? channel.at(d);
  if (typeof c !== 'number') return inkOn(fill, tokens);
  // the pair as one key: a fill's RGB and an ink's RGB fit 48 bits
  const key = (fill >>> 8) * 0x1000000 + (c >>> 8);
  let ink = memo.get(key);
  if (ink === undefined) memo.set(key, (ink = inkAgainst(c, fill)));
  return ink;
}
