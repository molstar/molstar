/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The letters layer on what is not letters, and on letters whose component is not their standard one (_spec/README.md
// "Render", Layer types: letters, Molecules): a molecule's residues labeled by their component codes and its atoms by
// their element symbols, a cell a position; a modified position's mark (`modified`: `mark`, a dot at the cell's top
// right in its letter's ink, a bar over a component of several positions, a tick at the top of a strip below the
// letters; `code`, the component's code in place of the letter where three letters fit, the mark elsewhere; `tint`,
// the cell outlined in `seq-modified`; `none`). No value changes a cell's width (compact by default).

import { projectRanges, type ModifiedIndex, type StaticModel } from '@molstar/seqx-core/resolve';
import type { ResolvedLayer } from '@molstar/seqx-core/view';
import type { AxisTransform } from './axis.js';
import type { RenderConfig } from './config.js';
import { Anchor, grown, NoColor, type PassName, type PassWriter } from './draw-lists.js';
import type { LayerSource } from './layer-data.js';
import type { Emphasis, LayerType } from './layer-type.js';
import { type LettersLay, type LettersPrep, LettersType } from './letters.js';
import { chipInk } from './letters-draw.js';
import { inkOn, type SeqTokens } from './tokens.js';

const AnchorMiddle: number = Anchor.middle;

/** How a letters layer draws a position whose component is not standard. */
export type ModifiedStyle = 'mark' | 'code' | 'tint' | 'none';

/** The `modified` param read: `mark` when absent or unknown. */
export function modifiedStyle(v: unknown): ModifiedStyle {
  return v === 'code' || v === 'tint' || v === 'none' ? v : 'mark';
}

/**
 * A molecule's labels a position (its residues' component codes, or its atoms' element symbols), on the system a
 * letters layer reads: a molecule's residues, a unit of a molecule, an atoms system. Undefined for a sequence.
 */
export interface MoleculeLabels {
  readonly kind: 'residues' | 'atoms';
  readonly labels: readonly string[];
  /** Its molecule's class, name and first component: the fields a chip's color reads (`class`, `name`, `component`). */
  readonly class: string;
  readonly name: string;
  readonly component: string;
}

const Labels = new WeakMap<StaticModel, Map<string, MoleculeLabels | null>>();

/** The datum a molecule's color is read at (one, reused: a layout makes none): its fields, no letter or value. */
const Datum = {
  of: undefined as MoleculeLabels | undefined,
  field(name: string): string | undefined {
    const m = this.of;
    return !m
      ? undefined
      : name === 'class'
        ? m.class
        : name === 'name'
          ? m.name
          : name === 'component'
            ? m.component
            : undefined;
  },
  value: Number.NaN,
  letter: '',
  start: 0,
  end: 1,
};
const moleculeDatum = (m: MoleculeLabels) => {
  Datum.of = m;
  return Datum;
};

/** A system's molecule labels (made once a model and system), or undefined when it is no molecule's. */
export function moleculeLabels(model: StaticModel, system: string): MoleculeLabels | undefined {
  let byModel = Labels.get(model);
  if (!byModel) Labels.set(model, (byModel = new Map()));
  let out = byModel.get(system);
  if (out === undefined) byModel.set(system, (out = makeLabels(model, system) ?? null));
  return out ?? undefined;
}

function makeLabels(model: StaticModel, system: string): MoleculeLabels | undefined {
  const s = model.systems.get(system);
  const m = model.molecular;
  const fields = (address: string) => {
    const unit = m.units.get(address);
    const molecule = m.molecules.get(unit ? unit.of : address);
    return { molecule, unit };
  };
  if (s?.kind === 'atoms') {
    const a = m.atoms.get(system);
    if (!a) return undefined;
    const n = a.offsets[a.offsets.length - 1];
    const { molecule } = fields(a.of);
    return {
      kind: 'atoms',
      labels: Array.from({ length: n }, (_, i) => a.elements.get(i) ?? ''),
      class: molecule?.class ?? 'ligand',
      name: molecule?.name ?? '',
      component: molecule?.residues?.[0] ?? '',
    };
  }
  const { molecule, unit } = fields(system);
  const residues = molecule?.residues;
  if (!molecule || !residues) return undefined;
  const start = unit?.start ?? 0;
  const end = unit?.end ?? residues.length;
  return {
    kind: 'residues',
    labels: residues.slice(start, end),
    class: molecule.class,
    name: molecule.name ?? '',
    component: residues[start] ?? '',
  };
}

/** What a letters layout keeps for molecules and modified positions. */
export interface MoleculeLay {
  /** A molecule's labels by position (undefined on letters). */
  labels: readonly string[] | undefined;
  /** Modified positions in view: `x0, x1, row` triples (the row into `modComponents`). */
  mods: number[];
  modStyle: ModifiedStyle;
  modComponents: readonly string[];
  /** Whether a row covers several positions (a bar), by row. */
  modLong: Uint8Array | undefined;
}

/** The parts of one row's image (kept: a layout runs to its end). */
const Parts: number[] = [];
const Range = [0, 1];
const Spans: number[] = [];

/**
 * The modified positions in view as `x0, x1, row` triples into `out`: each row's positions on the data's system, on
 * the section's own or through a one-to-one path; a row over several positions one span a piece.
 */
export function modifiedSpans(
  index: ModifiedIndex,
  src: LayerSource | undefined,
  axis: AxisTransform,
  out: number[],
): void {
  out.length = 0;
  const path = src?.path;
  for (let i = 0; i < index.position.length; i++) {
    const s = index.position[i];
    const e = s + index.length[i];
    Spans.length = 0;
    if (!path) axis.spans(s, e, Spans);
    else {
      Range[0] = s;
      Range[1] = e;
      projectRanges(path, Range, Parts);
      for (let k = 0; k + 1 < Parts.length; k += 2) axis.spans(Parts[k], Parts[k + 1], Spans);
    }
    // alternatives at one position: one mark (the first row's)
    if (i > 0 && index.position[i - 1] === s) continue;
    for (let k = 0; k + 2 < Spans.length; k += 3) out.push(Spans[k], Spans[k + 1], i);
  }
}

/** Whether each row covers several positions (a bar), or undefined when none does. */
export function longRows(index: ModifiedIndex): Uint8Array | undefined {
  let any = false;
  const out = new Uint8Array(index.length.length);
  for (let i = 0; i < out.length; i++)
    if (index.length[i] > 1) {
      out[i] = 1;
      any = true;
    }
  return any ? out : undefined;
}

/** The cells a molecule's letters layout fills, the arrays it writes and how its labels fit. */
interface Cells {
  n: number;
  pos: Int32Array;
  code: Uint16Array;
  fill: Uint32Array;
  ink: Uint32Array;
  x0: Float32Array;
  x1: Float32Array;
  show: Uint8Array;
  font: string;
  fontSize: number;
  glyphY: number;
}

/**
 * A molecule's cells in view (a cell a position, its fill `fill`), each labeled where its label fits at the largest
 * letter size the cell and the row allow: residue codes, element symbols.
 */
export function layMoleculeCells(
  lay: Cells,
  labels: readonly string[],
  axis: AxisTransform,
  vis: readonly number[],
  box: { readonly y: number; readonly height: number },
  config: RenderConfig,
  tokens: SeqTokens,
  fill: number,
  letters: boolean,
): void {
  let count = 0;
  for (let i = 0; i + 1 < vis.length; i += 2) count += vis[i + 1] - vis[i];
  lay.pos = grown(lay.pos, count);
  lay.code = grown(lay.code, count);
  lay.fill = grown(lay.fill, count);
  lay.ink = grown(lay.ink, count);
  lay.x0 = grown(lay.x0, count);
  lay.x1 = grown(lay.x1, count);
  lay.show = grown(lay.show, count);
  const ink = inkOn(fill, tokens);
  const dpr = axis.dpr;
  let narrowest = Infinity;
  for (let i = 0; i + 1 < vis.length; i += 2)
    for (let p = vis[i]; p < vis[i + 1]; p++) {
      // its cell in each piece in view, as the letters' own loop reads it (`piecesIn`, `spanIn`)
      axis.piecesIn(p, p + 1);
      const ka = axis.ka;
      const kb = axis.kb;
      let k = ka;
      while (k <= kb && !axis.spanIn(k, p, p + 1)) k++;
      if (k > kb) continue;
      const j = lay.n++;
      if (j >= lay.pos.length) break;
      lay.pos[j] = p;
      lay.code[j] = 0;
      lay.fill[j] = fill;
      lay.ink[j] = ink;
      lay.x0[j] = Math.round(axis.sx0 * dpr) / dpr;
      lay.x1[j] = Math.round(axis.sx1 * dpr) / dpr;
      lay.show[j] = 0;
      narrowest = Math.min(narrowest, lay.x1[j] - lay.x0[j]);
    }
  lay.fontSize = 0;
  if (!letters || !lay.n) return;
  for (const s of config.letterSizes)
    if (s * config.monoAdvance <= narrowest - 2 && s <= box.height * 0.75) lay.fontSize = s;
  if (!lay.fontSize) lay.fontSize = config.letterSizes[0] ?? 9;
  lay.font = `${lay.fontSize}px mono`;
  lay.glyphY = box.y + box.height / 2 + lay.fontSize * 0.35;
  const advance = lay.fontSize * config.monoAdvance;
  for (let j = 0; j < lay.n; j++) {
    const label = labels[lay.pos[j]] ?? '';
    lay.show[j] = label && label.length * advance + 1 <= lay.x1[j] - lay.x0[j] ? 1 : 0;
  }
}

/** A molecule's labels drawn centred in their cells, where they fit (`show`). */
export function drawMoleculeLabels(lay: Cells & MoleculeLay, out: PassWriter): void {
  const labels = lay.labels!;
  const g = out.args;
  for (let j = 0; j < lay.n; j++) {
    if (!lay.show[j]) continue;
    g[0] = (lay.x0[j] + lay.x1[j]) / 2;
    g[1] = lay.glyphY;
    g[2] = 0;
    g[3] = lay.x1[j] - lay.x0[j];
    g[4] = lay.ink[j];
    out.textArgs(labels[lay.pos[j]] ?? '', AnchorMiddle, lay.font);
  }
}

/**
 * Before the letters are drawn: `code` hides the letter of a cell that shows its component's code instead (the cell
 * has room for three letters), so the two are never drawn over each other.
 */
export function hideCoded(lay: Cells & MoleculeLay & { step: string }, config: RenderConfig): void {
  if (lay.modStyle !== 'code' || lay.step !== 'letters' || !(lay.fontSize > 0)) return;
  const room = 3 * lay.fontSize * config.monoAdvance + 1;
  const m = lay.mods;
  for (let i = 0; i + 2 < m.length; i += 3) {
    if (m[i + 1] - m[i] < room || lay.modLong?.[m[i + 2]]) continue;
    for (let j = 0; j < lay.n; j++) if (lay.x0[j] >= m[i] - 0.5 && lay.x1[j] <= m[i + 1] + 0.5) lay.show[j] = 0;
  }
}

/**
 * A modified cell's ink: its letter's, from the cell under mark `i`'s middle (the page's text where none is laid out).
 * The mark by its index: a number passed to a call that is not inlined is boxed (a warm draw allocates nothing).
 */
function inkUnder(lay: Cells & MoleculeLay, i: number, tokens: SeqTokens): number {
  const x = (lay.mods[i] + lay.mods[i + 1]) / 2;
  for (let j = 0; j < lay.n; j++) if (x >= lay.x0[j] && x < lay.x1[j]) return lay.ink[j];
  return tokens.colors['msui-text'];
}

/**
 * The modified positions' marks over the cells (after the letters): a dot at a cell's top right, a bar over a long
 * component, a tick at the strip's top below the letters' sizes; `code` the component's code where it has room;
 * `tint` the cell outlined in `seq-modified`. Every position in view gets one (EC-117's codes where they fit).
 */
export function drawModified(
  lay: Cells & MoleculeLay & { step: string; box: { y: number; height: number }; axis: AxisTransform },
  out: PassWriter,
  tokens: SeqTokens,
  config: RenderConfig,
): void {
  const m = lay.mods;
  const g = out.args;
  const { box } = lay;
  const px = 1 / lay.axis.dpr;
  const size = config.letters.modifiedMark.size;
  const underline = config.letters.modifiedMark.shape === 'underline';
  const fine = lay.step === 'letters' || lay.step === 'cells';
  const room = 3 * lay.fontSize * config.monoAdvance + 1;
  for (let i = 0; i + 2 < m.length; i += 3) {
    const x0 = m[i];
    const x1 = m[i + 1];
    const row = m[i + 2];
    const w = x1 - x0;
    if (lay.modStyle === 'tint' && fine) {
      g[0] = x0 + px;
      g[1] = box.y + px;
      g[2] = Math.max(w - 2 * px, px);
      g[3] = Math.max(box.height - 2 * px, px);
      g[4] = NoColor;
      g[5] = tokens.colors['seq-modified'];
      g[6] = Math.max(1.5, 2 * px);
      g[7] = 0;
      g[8] = 0;
      out.rectArgs();
      continue;
    }
    const long = !!lay.modLong?.[row];
    if (lay.modStyle === 'code' && lay.step === 'letters' && lay.fontSize > 0 && !long && w >= room) {
      g[0] = (x0 + x1) / 2;
      g[1] = lay.glyphY;
      g[2] = 0;
      g[3] = w;
      g[4] = inkUnder(lay, i, tokens);
      out.textArgs(lay.modComponents[row] ?? '', AnchorMiddle, lay.font);
      continue;
    }
    const s = Math.min(size, Math.max(1, box.height / 3));
    g[5] = NoColor;
    g[6] = 0;
    if (!fine) {
      // below the cells: a tick at the strip's top, the mark surviving the zoom out
      g[0] = x0;
      g[1] = box.y;
      g[2] = Math.max(w, px);
      g[3] = s;
      g[4] = tokens.colors['msui-text'];
      g[7] = 0;
      g[8] = 0;
    } else if (long || underline) {
      // a bar along the cells' top (a component over several positions), or the configured underline
      g[0] = x0 + 1;
      g[1] = underline ? box.y + box.height - s - 1 : box.y + 1;
      g[2] = Math.max(w - 2, px);
      g[3] = Math.max(1, s * 0.67);
      g[4] = inkUnder(lay, i, tokens);
      g[7] = 0;
      g[8] = 0;
    } else {
      // a dot at the cell's top right, inside it, never over the letter's middle
      const d = Math.min(s, Math.max(1, w / 3));
      g[0] = x1 - d - Math.min(1, w / 6);
      g[1] = box.y + 1;
      g[2] = d;
      g[3] = d;
      g[4] = inkUnder(lay, i, tokens);
      g[7] = d / 2;
      g[8] = 3;
    }
    out.rectArgs();
  }
}

/**
 * The letters layer with molecules (_spec/README.md "Render", Layer types: letters, Molecules): the core letters, and
 * on a molecule's residues or atoms its labels, on a modified position its mark. A composition registers it in the
 * letters' place (`letters` by name): the default viewer does, an export loads it when its model holds molecules or
 * modified residues (`renderSvg`); a minimal composition need not, and a view of plain letters never reads it.
 */
export const MoleculeLettersType: LayerType<LettersPrep, LettersLay> = {
  ...LettersType,

  prepare(layer: ResolvedLayer, model: StaticModel): LettersPrep {
    const prep = LettersType.prepare(layer, model);
    const dataSystem = prep.src?.dataSystem ?? layer.data?.address;
    const molecule = dataSystem === undefined ? undefined : moleculeLabels(model, dataSystem);
    const modified = modifiedStyle(layer.params.modified);
    // a modified position's mark: on the section's own letters or through a one-to-one path (never a codon's parts)
    const marks =
      modified !== 'none' && dataSystem !== undefined && (prep.through === 'own' || prep.through === 'one')
        ? model.modified(dataSystem)
        : undefined;
    const longMarks = marks && longRows(marks);
    return {
      ...prep,
      ...(molecule ? { molecule } : {}),
      modified,
      ...(marks?.position.length ? { marks } : {}),
      ...(longMarks ? { longMarks } : {}),
    };
  },

  layout(prep, ctx, data, overlay, prev) {
    const lay = LettersType.layout(prep, ctx, data, overlay, prev);
    lay.modStyle = prep.modified ?? 'mark';
    lay.modComponents = prep.marks?.component ?? [];
    lay.modLong = prep.longMarks;
    const { axis, box, config, tokens } = ctx;
    // every modified position in view has its mark, whatever the step (a tick below the cells)
    if (prep.marks && prep.src) modifiedSpans(prep.marks, prep.src, axis, lay.mods);
    if (prep.molecule) {
      // a molecule's residues or atoms: a cell a position, labeled by its code or element where it fits
      lay.labels = prep.molecule.labels;
      if (lay.step === 'bars' || lay.step === 'density') return lay;
      const c = prep.channels.forms.color === undefined ? undefined : prep.channels.get('color', 'token:sequence', ctx);
      // a chip colored on its own: a theme over the molecule's fields (its class, name, component), once a block
      const own = c && c.constant === undefined ? c.at(moleculeDatum(prep.molecule)) : undefined;
      const fill =
        typeof c?.constant === 'number' ? c.constant : typeof own === 'number' ? own : tokens.colors['seq-sequence'];
      lay.n = 0;
      lay.pending.length = 0;
      const vis = axis.visible([]);
      layMoleculeCells(lay, prep.molecule.labels, axis, vis, box, config, tokens, fill, lay.step === 'letters');
      if (lay.chip && lay.n) lay.chipInk = chipInk(lay.fill[0], tokens);
      return lay;
    }
    if (lay.mods.length) hideCoded(lay, config);
    return lay;
  },

  draw(lay, pass: PassName, out: PassWriter, emphasis?: Emphasis) {
    LettersType.draw(lay, pass, out, emphasis);
    if (pass !== 'data') return;
    if (lay.labels && lay.step === 'letters' && lay.fontSize > 0 && lay.n && !lay.chip) drawMoleculeLabels(lay, out);
    if (lay.mods.length && !(lay.chip && lay.n)) drawModified(lay, out, lay.tokens, lay.config);
  },
};
