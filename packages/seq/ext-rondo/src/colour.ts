/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type LayoutContext, type PairsRead, type PairsSource, theme } from '@molstar/seqx-core/extension';
import { grow } from './geometry.js';
import type { RondoLay, RondoPrep } from './lay.js';

// What a ribbons layer reads besides its pairs (_spec/README.md, Colour, Kinds of links): its colour (the blocks' hues,
// a fixed colour, a field's a link), its rows' weights (`width`: the count, or a field summed) and the other ribbons
// layers of its section (a bin's card counts each). Kept buffers: a layout makes them once a size.

/** How the layer is coloured (_spec/README.md, Colour): its blocks' hues, one colour, or a field's colour a link. */
export const ColorMode = { system: 0, fixed: 1, field: 2 } as const;

/** This layout's rows' colours (a field colour's) and weights, as `colorsOf` and `weightsOf` left them. */
export const rowColors = (): Uint32Array => RowColor;
export const rowWeights = (): Float64Array => RowW;
/** A layout's colour classes start again (a ribbon's key reads its colour's class). */
export const clearClasses = () => ColorClass.clear();

/** The section's lays by its resolved section (a new resolve, a new list): each layout adds its own once. */
const Peers = new WeakMap<object, RondoLay[]>();

export function peersOf(lay: RondoLay, ctx: LayoutContext) {
  let list = Peers.get(ctx.section);
  if (!list) Peers.set(ctx.section, (list = []));
  // one lay a layer: a lay laid out again for the same key replaces the one before it
  const at = list.findIndex((l) => l === lay || l.key === lay.key);
  if (at < 0) list.push(lay);
  else list[at] = lay;
  lay.peers = list;
}

/** Each row's colour (a field colour), its ribbon's colour class by colour, made a layout; kept buffers. */
let RowColor = new Uint32Array(256);
const ColorClass = new Map<number, number>();
export const classOf = (c: number) => {
  let k = ColorClass.get(c);
  if (k === undefined) ColorClass.set(c, (k = Math.min(63, ColorClass.size)));
  return k;
};

/** The colour channel's default where the layer writes one it cannot read (a scheme that does not fit): the hues. */
const SystemColor = theme.field('system', { palette: theme.categorical({ colors: 'categories', repeat: true }) });
const PairEnds = new Set(['a_start', 'a_end', 'b_start', 'b_end']);

/**
 * The layer's colour (_spec/README.md, Colour): absent or by `system`, the blocks' hues; a fixed colour, every polygon's;
 * another field (`type`) read link by link: at the item level each ribbon its links' colour (a stem pair split by
 * colour), at the bins level the layer's one colour, its value when every link read shares it, else the ranges' ink and
 * `rondo.ribbons.color` said (a bin sums links of several values). The `opacity` channel's uniform value scales the
 * fill.
 */
export function colorsOf(lay: RondoLay, prep: RondoPrep, ctx: LayoutContext, read: PairsRead) {
  lay.mode = ColorMode.system;
  lay.channel = undefined;
  lay.one = -1;
  lay.byRow = false;
  const channels = prep.source.channels;
  const o = channels.get('opacity', 1, ctx).constant;
  lay.opacity = typeof o === 'number' && o >= 0 ? Math.min(1, o) : 1;
  if (!prep.colored) return;
  const d = read.data!;
  const domain = ctx.section.domains?.system;
  const c = channels.get('color', SystemColor, ctx, {
    numeric: (n: string) => PairEnds.has(n) || !!d.extra?.[n]?.some((v) => typeof v === 'number'),
    categories: (f: string) => (f === 'system' ? domain : undefined),
  });
  lay.channel = c;
  if (c.field === 'system') return;
  if (typeof c.constant === 'number') {
    lay.mode = ColorMode.fixed;
    lay.one = c.constant;
    return;
  }
  lay.mode = ColorMode.field;
  const datum = (read as PairsSource).datum;
  const ink = ctx.tokens.colors['seq-range'];
  RowColor = grow(RowColor, d.n);
  let first = -1;
  let mixed = false;
  for (let k = 0; k < d.n; k++) {
    if (!read.shown[k]) continue;
    const v = datum ? c.at(datum.set(k)) : undefined;
    const col = typeof v === 'number' ? v >>> 0 : ink;
    RowColor[k] = col;
    if (first < 0) first = col;
    else if (col !== first) mixed = true;
  }
  lay.one = first < 0 || mixed ? ink : first;
  lay.byRow = lay.items && mixed;
  if (mixed && !lay.items)
    ctx.report({
      code: 'rondo.ribbons.color',
      severity: 'info',
      message: `A colour by ${c.field ?? 'a field'} is read link by link: the bins sum links of several values, drawn in one colour`,
      address: ctx.key,
    } as never);
}

/** Each row's weight (`width`): its count, or the field summed (a held bin's sum, an item's value), kept. */
let RowW = new Float64Array(256);

export function weightsOf(lay: RondoLay, prep: RondoPrep, ctx: LayoutContext, read: PairsRead) {
  const d = read.data!;
  const w = read.weight;
  const col = prep.width === 'count' ? undefined : d.extra?.[prep.width];
  const numeric = !!col?.some((v) => typeof v === 'number');
  if (prep.width !== 'count' && !numeric)
    ctx.report({
      code: 'rondo.ribbons.width',
      severity: 'info',
      message: `Ribbons wide by count: the links read carry no numeric ${prep.width}`,
      address: ctx.key,
    } as never);
  RowW = grow(RowW, d.n);
  for (let k = 0; k < d.n; k++) RowW[k] = numeric ? Number(col![k]) || 0 : w ? w[k] : 1;
}
