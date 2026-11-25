/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import {
  type LayerType,
  type NavigatorPicture,
  type NavigatorPictureContext,
  PairsReader,
  type PairsSource,
  type SquareLattice,
} from '@molstar/seqx-core/extension';
import { HoleIssue } from './hole.js';

// The square's picture (_spec/README.md, The square navigator): `rondo:square` names one `rondo:cells` track over the
// bound ring's ribbons' data; `rondo:cells` reads that data's coarsest bins over the assembly (its section a hidden
// one the square lays out) and sums them by block pair: a cell a block pair coloured by its count on a square-root
// ramp (`seq-ramp-*`), the diagonal (links within a chromosome) included, mirrored for unsided data (one assembly),
// sided data (two) each pair in its own cell. It draws nothing on a strip: its cells are the square's.

/** The read's bins: as coarse as the data has (one want a block), never items. */
const CellRead = { gate: 0, binPx: Number.MAX_SAFE_INTEGER };

/** What the cells read, made once a resolved layer. */
export interface CellsPrep {
  readonly source: PairsSource;
  /** The pairs' two sides on different assemblies (each pair in its own cell, not mirrored). */
  readonly sided: boolean;
}

/** The cells laid out: the read's systems, its counts by system pair (`n × n`, row `a`), the largest. */
export class CellsLay {
  readonly reader = new PairsReader();
  systems: readonly string[] = [];
  /** A system's index among the read's. */
  readonly index = new Map<string, number>();
  n = 0;
  counts = new Float64Array(0);
  max = 0;
  /** No bins read yet, or a tile still awaited. */
  pending = true;
  /** The source's address (its failure's issue is found by it), and the store read through (the viewer's it is). */
  source: string | undefined;
  store: object | undefined;
}

function countCells(lay: CellsLay, prep: CellsPrep, settled: boolean): void {
  const read = lay.reader.last;
  const d = read?.data;
  lay.pending = !d ? !settled || lay.reader.waiting : lay.reader.waiting;
  if (!read || !d) return;
  const systems = read.endUnits ?? (read.aSystem === read.bSystem ? [read.aSystem] : [read.aSystem, read.bSystem]);
  if (systems !== lay.systems) {
    lay.systems = systems;
    lay.index.clear();
    for (let k = 0; k < systems.length; k++) lay.index.set(systems[k], k);
  }
  const n = systems.length;
  if (lay.counts.length < n * n) lay.counts = new Float64Array(n * n);
  lay.n = n;
  lay.counts.fill(0, 0, n * n);
  const w = read.weight;
  const aSys = 0;
  const bSys = n > 1 && !d.bUnit ? 1 : 0;
  let max = 0;
  for (let k = 0; k < d.n; k++) {
    if (!read.shown[k]) continue;
    const a = d.aUnit ? d.aUnit[k] : aSys;
    const b = d.bUnit ? d.bUnit[k] : bSys;
    if (a < 0 || b < 0) continue;
    const c = w ? w[k] : 1;
    const ab = (lay.counts[a * n + b] += c);
    if (ab > max) max = ab;
    if (!prep.sided && a !== b) lay.counts[b * n + a] += c;
  }
  lay.max = max;
}

/** A count's colour on the square-root ramp (five stops of `seq-ramp-*`). */
function rampOf(lattice: SquareLattice, c: number, max: number): number {
  const v = Math.sqrt(c / max);
  const k = Math.min(4, Math.floor(v * 5));
  const t = lattice.tokens.colors;
  return k === 0
    ? t['seq-ramp-0']
    : k === 1
      ? t['seq-ramp-1']
      : k === 2
        ? t['seq-ramp-2']
        : k === 3
          ? t['seq-ramp-3']
          : t['seq-ramp-4'];
}

/** The count of the lattice's cell `(i, j)`, 0 when either block is not read. */
function countAt(lay: CellsLay, lattice: SquareLattice, i: number, j: number): number {
  const a = lay.index.get(lattice.systems[i]);
  const b = lay.index.get(lattice.systems[j]);
  return a === undefined || b === undefined ? 0 : lay.counts[a * lay.n + b];
}

/** `rondo:cells`: the square's cells over pairs, by block pair. */
export const CellsType: LayerType<CellsPrep, CellsLay> = {
  name: 'rondo:cells',
  shape: 'pairs',
  scope: 'section',
  prepare(layer, model) {
    const source = PairsReader.source(layer, model);
    return { source, sided: source.aSystem !== source.bSystem && !source.endUnits };
  },
  layout(prep, ctx, data, _overlay, prev) {
    const lay = prev ?? new CellsLay();
    lay.reader.read(prep.source, ctx, data, CellRead);
    countCells(lay, prep, ctx.settled);
    lay.source = prep.source.region?.src.source;
    lay.store = data.origin ?? data;
    return lay;
  },
  // nothing on a strip: its cells are the square's
  draw() {},
  hit: () => undefined,
  square: {
    draw(lay, out, lattice) {
      if (!(lay.max > 0)) return;
      const n = lattice.systems.length;
      for (let i = 0; i < n; i++) {
        const x0 = lattice.edge(i, true);
        const x1 = lattice.edge(i + 1, true);
        for (let j = 0; j < n; j++) {
          const c = countAt(lay, lattice, i, j);
          if (!(c > 0)) continue;
          const y0 = lattice.edge(j, false);
          out.fill(x0, y0, x1 - x0, lattice.edge(j + 1, false) - y0, rampOf(lattice, c, lay.max));
        }
      }
    },
    text(lay, i, j, lattice) {
      if (lay.pending && !lay.max) return 'links loading';
      const c = Math.round(countAt(lay, lattice, i, j));
      return `${c.toLocaleString('en-US')} link${c === 1 ? '' : 's'}${i === j ? ' within' : ''}`;
    },
    pending: (lay) => lay.pending,
    // links pending: the hatch and `Loading…`; failed: the issue the source raised (as the ring's hole says it)
    note: (lay) =>
      lay.reader.failed
        ? (HoleIssue.of?.(undefined, lay.source ?? '', lay.store) ?? 'Links failed')
        : lay.pending && !lay.max
          ? 'Loading…'
          : undefined,
  },
};

/** A spec node (the picture's context's). */
type Node = NonNullable<NavigatorPictureContext['boundNode']>;

/** The first layer of a section node's tracks that `pick` takes. */
function layerIn(section: Node | undefined, pick: (l: Node) => boolean): Node | undefined {
  for (const t of section?.children ?? [])
    for (const l of t.children ?? []) if (l.kind === 'layer' && pick(l)) return l;
  return undefined;
}

/**
 * `rondo:square`: the square's picture, one `rondo:cells` track over the bound section's ribbons' data (its
 * `of_layer` param names the layer by id; else the first `rondo:ribbons` layer, else the first `pairs` one).
 */
export const RondoSquare: NavigatorPicture = {
  name: 'rondo:square',
  axes: ['both'],
  size: [220, 220],
  picture(ctx) {
    const named = ctx.navigator.params?.of_layer;
    const layer =
      (typeof named === 'string' ? layerIn(ctx.boundNode, (l) => l.id === named) : undefined) ??
      layerIn(ctx.boundNode, (l) => l.params?.type === 'rondo:ribbons') ??
      layerIn(ctx.boundNode, (l) => l.params?.type === 'pairs');
    const data = layer?.params?.data;
    if (data === undefined) return { tracks: [] };
    return {
      tracks: [
        {
          kind: 'track',
          id: 'cells',
          params: { label: '', height: 8 },
          children: [{ kind: 'layer', id: 'cells', params: { type: 'rondo:cells', data } }],
        },
      ],
    };
  },
};
