/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { RangesLay } from '@molstar/seqx-core/render';
import { P53Views } from '@molstar/seqx-examples';
import { nodeFixtureAssets } from '@molstar/seqx-examples/node-assets';
import { Domains, GeneExons, GeneFile } from '@molstar/seqx-examples/p53/data';
import { GeneOffset, Transcripts } from '@molstar/seqx-examples/p53/region';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { afterEach, beforeAll, describe, expect, it, vi } from 'vitest';
import { checkMarks, openExample } from './draw-checks.js';

// Drawing checked against the data (the p53 example views, at the fit view and across a tile edge): the oracle is the
// example's own tables (`Domains`, `GeneExons`, `Transcripts`), never the layout: every row of them in the window has
// its marks (a transcript's every exon, a domain's every piece cut into exons), and every mark drawn is hit. Canvas,
// SVG and hits can agree and all be wrong; the tables cannot.

const DelayMs = 300;
const ViewId = '/views/main';

let viewer: SeqViewer | undefined;
let gene = '';

beforeAll(async () => {
  gene = await nodeFixtureAssets().reader.readText(GeneFile, new AbortController().signal);
});

afterEach(() => {
  viewer?.retire();
  viewer = undefined;
  vi.useRealTimers();
});

const view = (id: string) => P53Views.find((v) => v.id === id)!;

async function open(id: string): Promise<SeqViewer> {
  viewer = await openExample(view(id), { gene, delayMs: DelayMs });
  return viewer;
}

/** A row of the data: an extent on the data's system, and the pieces it is drawn in on the section's. */
interface Row {
  readonly id: string;
  readonly start: number;
  readonly end: number;
  readonly pieces: readonly (readonly [number, number])[];
}

/** A domain's residues `[start, end)` cut into the exons they fall in (coding order, minus strand: genomic reversed). */
function cut(
  exons: readonly (readonly [number, number])[],
  strand: '+' | '-',
  start: number,
  end: number,
): [number, number][] {
  // each coding nucleotide's genomic position, in coding order
  const coding: number[] = [];
  const ordered = strand === '+' ? exons : [...exons].reverse();
  for (const [a, b] of ordered) {
    if (strand === '+') for (let p = a; p < b; p++) coding.push(p);
    else for (let p = b - 1; p >= a; p--) coding.push(p);
  }
  const hit = coding.slice(3 * start, 3 * end).sort((x, y) => x - y);
  const pieces: [number, number][] = [];
  for (const p of hit) {
    const last = pieces[pieces.length - 1];
    // a run continues within one exon only: positions an intron apart are two pieces
    if (last && last[1] === p) last[1] = p + 1;
    else pieces.push([p, p + 1]);
  }
  return pieces;
}

const domainRows = (pieces: (start: number, end: number) => [number, number][], limit = Infinity): Row[] =>
  Domains.id.flatMap((id, i) =>
    Domains.start[i] < limit
      ? [
          {
            id,
            start: Domains.start[i],
            end: Math.min(Domains.end[i], limit),
            pieces: pieces(Domains.start[i], Math.min(Domains.end[i], limit)),
          },
        ]
      : [],
  );

/**
 * The marks the layer draws for a row: the hit marks of the item that has the row's extent on the data's system (an
 * extent is one row's in these tables). The layer holds every item loaded; its hit marks are those in the window.
 */
function marksOf(lay: RangesLay, row: Row): [number, number][] | undefined {
  for (let i = 0; i < lay.items; i++) {
    if (lay.starts[i] !== row.start || lay.ends[i] !== row.end) continue;
    const out: [number, number][] = [];
    for (let k = 0; k < lay.hn; k++) if (lay.hitem[k] === i) out.push([lay.hx0[k], lay.hx1[k]]);
    return out;
  }
  return undefined;
}

/**
 * Every row with a piece in `[from, to)` has an item in the layer, drawn in as many marks as it has pieces in the
 * window, each inside the layer's box; no mark is drawn that no row accounts for. Returns what is wrong.
 */
function checkRows(lay: RangesLay, rows: readonly Row[], from: number, to: number): string[] {
  const wrong: string[] = [];
  const inWindow = (r: Row) => r.pieces.filter(([a, b]) => b > from && a < to);
  const expected = rows.filter((r) => inWindow(r).length);
  for (const r of expected) {
    const marks = marksOf(lay, r);
    if (!marks) {
      wrong.push(`${lay.key}: ${r.id} [${r.start}, ${r.end}) has no item in view`);
      continue;
    }
    if (marks.length !== inWindow(r).length)
      wrong.push(`${lay.key}: ${r.id} has ${marks.length} marks, its pieces in view are ${inWindow(r).length}`);
    for (const [x0, x1] of marks)
      if (!(x1 > x0) || x1 < lay.box.x || x0 > lay.box.x + lay.box.width)
        wrong.push(`${lay.key}: ${r.id} mark ${x0}–${x1} is outside its box`);
  }
  const total = expected.reduce((n, r) => n + inWindow(r).length, 0);
  if (lay.hn !== total) wrong.push(`${lay.key}: ${lay.hn} marks drawn, the data's pieces in view are ${total}`);
  return wrong;
}

const layers = (v: SeqViewer) => v.frame!.layout.sections.flatMap((s) => s.tracks.flatMap((t) => t.layers));
const layerIn = (v: SeqViewer, section: string, track: string) => {
  const t = v.frame!.layout.sections
    .find((s) => s.section.address === `${ViewId}/${section}`)
    ?.tracks.find((x) => x.track.address.endsWith(`/${track}`));
  return t?.layers.find((l) => l.type.name === 'ranges')?.lay as RangesLay | undefined;
};

describe('Tracks: drawing against the data', () => {
  const rows = domainRows((a, b) => [[a, b]]);

  it('at the fit view, every domain has its mark and every mark is hit', async () => {
    const v = await open('tracks');
    const lay = layerIn(v, 'protein', 'domains')!;
    expect(checkRows(lay, rows, 0, 393)).toEqual([]);
    expect(checkMarks(v).wrong).toEqual([]);
  });

  it('zoomed across the edge between two domains (the proline-rich and DNA-binding ones overlap its window)', async () => {
    const v = await open('tracks');
    v.view.setRange(`${ViewId}/protein`, 55, 100);
    v.draw();
    const lay = layerIn(v, 'protein', 'domains')!;
    expect(checkRows(lay, rows, 55, 100)).toEqual([]);
    expect(lay.hn).toBe(3);
    expect(checkMarks(v).wrong).toEqual([]);
  });
});

describe('Protein on its gene: drawing against the data', () => {
  const exons = GeneExons.start.map((s, i) => [s, GeneExons.end[i]] as const);
  const protein = domainRows((a, b) => [[a, b]]);
  const onGene = domainRows((a, b) => cut(exons, '-', a, b));

  it('at the fit view, every domain has its mark on the protein and each of its exon pieces on the gene', async () => {
    const v = await open('protein-on-gene');
    expect(checkRows(layerIn(v, 'protein', 'domains')!, protein, 0, 393)).toEqual([]);
    const gene = layerIn(v, 'gene', 'domains')!;
    const w = v.view.get(`${ViewId}/gene`)!;
    expect(checkRows(gene, onGene, Math.floor(w.start), Math.ceil(w.end))).toEqual([]);
    expect(checkMarks(v).wrong).toEqual([]);
  });

  it('zoomed across an intron edge, the domain cut by it has a piece each side and each is hit', async () => {
    const v = await open('protein-on-gene');
    // minus strand: the coding starts in the last exon (residues 0 to 154), the next holds 154 to 294. The gene's
    // section follows the protein's: a window over residue 154 is over the intron between them
    v.view.setRange(`${ViewId}/protein`, 144, 164);
    v.draw();
    const w = v.view.get(`${ViewId}/gene`)!;
    const lay = layerIn(v, 'gene', 'domains')!;
    expect(checkRows(lay, onGene, Math.floor(w.start), Math.ceil(w.end))).toEqual([]);
    // the DNA-binding domain (93 to 292) is in two pieces here, one each side of the intron: the check is not vacuous
    const dbd = onGene.find((r) => r.id === 'DBD')!;
    expect(dbd.pieces.filter(([a, b]) => b > w.start && a < w.end).length).toBe(2);
    expect(marksOf(lay, dbd)?.length).toBe(2);
    expect(checkMarks(v).wrong).toEqual([]);
  });
});

describe('Gene by region: drawing against the data', () => {
  const exonsOf = (id: string) => Transcripts[id].exons;
  const transcriptRows = (id: string): Row[] =>
    exonsOf(id).map(([a, b]) => ({ id, start: a, end: b, pieces: [[a, b]] }));

  /** A transcript's domains on the region: the protein's, cut to its length, each in its exon pieces. */
  const domainsOn = (id: string): Row[] => {
    const t = Transcripts[id];
    const length = Math.floor(exonsOf(id).reduce((n, [a, b]) => n + b - a, 0) / 3);
    return domainRows((a, b) => cut(exonsOf(id), t.strand, a, b), length);
  };

  const transcriptLayer = (v: SeqViewer, id: string) => {
    const t = v.frame!.layout.sections[0].tracks.find((x) => x.track.label === `Transcript ${id}`);
    return t?.layers.find((l) => l.type.name === 'ranges')?.lay as RangesLay | undefined;
  };
  const domainsLayer = (v: SeqViewer, id: string) => {
    const t = v.frame!.layout.sections[0].tracks.find((x) => x.track.label === `Domains ${id}-domains`);
    return t?.layers.find((l) => l.type.name === 'ranges')?.lay as RangesLay | undefined;
  };

  async function check(from: number, to: number, ids: readonly string[]) {
    const v = await open('gene-by-region');
    v.view.setRange(`${ViewId}/region`, from, to);
    await vi.advanceTimersByTimeAsync(20 * DelayMs);
    v.draw();
    const wrong: string[] = [];
    for (const id of ids) {
      const t = transcriptLayer(v, id);
      if (!t) wrong.push(`no track for ${id}`);
      else wrong.push(...checkRows(t, transcriptRows(id), from, to));
      if (Transcripts[id].domains) {
        const d = domainsLayer(v, id);
        if (!d) wrong.push(`no domains track for ${id}`);
        else wrong.push(...checkRows(d, domainsOn(id), from, to));
      }
    }
    expect(wrong).toEqual([]);
    expect(checkMarks(v).wrong).toEqual([]);
    return v;
  }

  it('at the fit view (the whole region), every transcript has its exons and its domains their pieces', async () => {
    const v = await check(0, 6000, ['T1', 'T2', 'T3', 'W1']);
    expect(layers(v).filter((l) => l.type.name === 'ranges').length).toBeGreaterThan(4);
  });

  it('zoomed across TP53’s tile edge (2500), the end of the first exon and the start of the second are drawn', async () => {
    const v = await check(GeneOffset + 350, GeneOffset + 550, ['T1', 'T2', 'T3']);
    // the first exon's end and the second's start, in T1 and T2 (T3 has only the last exon: none here)
    expect(transcriptLayer(v, 'T1')!.hn).toBe(2);
    expect(transcriptLayer(v, 'T2')!.hn).toBe(2);
    expect(transcriptLayer(v, 'T3')!.hn).toBe(0);
  });

  it('zoomed across WRAP53’s tile edge (4000), its two exons are drawn each side of it', async () => {
    const v = await check(3850, 4350, ['W1']);
    expect(transcriptLayer(v, 'W1')!.hn).toBe(2);
  });
});
