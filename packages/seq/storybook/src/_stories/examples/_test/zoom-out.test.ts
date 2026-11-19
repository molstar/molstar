/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { RangesLay } from '@molstar/seqx-core/render';
import { P53Views } from '@molstar/seqx-examples';
import { nodeFixtureAssets } from '@molstar/seqx-examples/node-assets';
import { GeneFile } from '@molstar/seqx-examples/p53/data';
import { regionSeqText } from '@molstar/seqx-examples/p53/region';
import { DefaultSeqViewerSpec } from '@molstar/seqx-viewer/default-spec';
import { createSeqViewer, type SeqViewer } from '@molstar/seqx-viewer/viewer';
import { bestOf, itWithin } from '@molstar/test-harness/timing';
import { afterEach, beforeAll, describe, expect, it, vi } from 'vitest';

// Zooming out on Protein on its gene, headless on a fake clock: a double-click on a domain on the following gene
// focuses it (its leader moves), a second goes back, both sections together; Shift+double-click and the keys zoom out;
// the gene keeps following. The performance representative: a back step (the command, then the frame) under 4 ms.

const Protein = '/views/main/protein';
const Gene = '/views/main/gene';

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

/** `clock`: the bound of the performance case reads the real `performance`, its timers faked alone. */
async function make(clock = true, id = 'protein-on-gene'): Promise<SeqViewer> {
  vi.useFakeTimers({
    toFake: clock ? ['setTimeout', 'clearTimeout', 'Date', 'performance'] : ['setTimeout', 'clearTimeout'],
  });
  const view = P53Views.find((v) => v.id === id)!;
  const text = id === 'gene-by-region' ? regionSeqText() : gene;
  const base = DefaultSeqViewerSpec();
  viewer = createSeqViewer({
    ...base,
    ui: undefined,
    registry: [...(base.registry ?? []), ...(view.registry?.({ delayMs: 300 }) ?? [])],
    ports: {
      assets: () => ({ reader: { readText: async () => text, readBytes: async () => new TextEncoder().encode(text) } }),
    },
  });
  viewer.resize(720);
  void viewer.load(view.build());
  await vi.advanceTimersByTimeAsync(6000);
  viewer.draw();
  return viewer;
}

/** A domain's marks on a section's domains track (the first such track, or the one whose address has `track`). */
function domainMarks(v: SeqViewer, section: string, key: string, track = '/domains') {
  v.draw();
  const s = v.frame!.layout.sections.find((x) => x.section.address === section)!;
  const t = s.tracks.find((x) => x.track.address.includes(track))!;
  const lay = t.layers[0].lay as RangesLay;
  const marks: [number, number][] = [];
  for (let j = 0; j < lay.hn; j++) if (lay.keys[lay.hitem[j]] === key) marks.push([lay.hx0[j], lay.hx1[j]]);
  return { marks, box: lay.box };
}

/** The centre of a domain's first mark on a section's domains track, in canvas px. */
function domainAt(v: SeqViewer, section: string, key: string) {
  v.draw();
  const s = v.frame!.layout.sections.find((x) => x.section.address === section)!;
  const lay = s.tracks.find((t) => t.track.address.endsWith('/domains'))!.layers[0].lay as RangesLay;
  for (let j = 0; j < lay.hn; j++)
    if (lay.keys[lay.hitem[j]] === key)
      return [(lay.hx0[j] + lay.hx1[j]) / 2, lay.box.y + (lay.hlane[j] + 0.5) * lay.laneHeight] as const;
  throw new Error(`no ${key} on ${section}`);
}

const rangeOf = (v: SeqViewer, s: string) => {
  const r = v.view.get(s)!;
  return [r.start, r.end];
};

describe('a double-click on an item focuses its whole extent on the section clicked', () => {
  it('Gene by region: T1’s DNA-binding domain, carried onto the region through its codon map from a nested source', async () => {
    const v = await make(true, 'gene-by-region');
    const region = v.view.resolved.value!.sections[0].address;
    const before = domainMarks(v, region, 'DBD', 'T1-domains');
    expect(before.marks.length).toBeGreaterThan(0);
    const [x0, x1] = before.marks[0];
    const y = v.frame!.layout.sections[0].tracks.find((t) => t.track.address.includes('T1-domains'))!;
    const lay = y.layers[0].lay as RangesLay;
    const j = Array.from({ length: lay.hn }, (_, k) => k).find((k) => lay.keys[lay.hitem[k]] === 'DBD')!;
    v.input.doubleClick((x0 + x1) / 2, lay.box.y + (lay.hlane[j] + 0.5) * lay.laneHeight);
    await vi.advanceTimersByTimeAsync(1000);
    // the domain's pieces (its codons over three exons): all in view, and they fill it from edge to edge
    const after = domainMarks(v, region, 'DBD', 'T1-domains');
    const lo = Math.min(...after.marks.map((m) => m[0]));
    const hi = Math.max(...after.marks.map((m) => m[1]));
    expect(lo).toBeGreaterThanOrEqual(after.box.x - 1);
    expect(hi).toBeLessThanOrEqual(after.box.x + after.box.width + 1);
    expect(hi - lo).toBeGreaterThan(0.95 * after.box.width);
    const r = rangeOf(v, region);
    // 199 residues are 597 nucleotides, and the introns between them
    expect(r[1] - r[0]).toBeGreaterThan(597);
  });

  it('Protein on its gene: the domain on the following gene focuses its hull there, the protein the domain', async () => {
    const v = await make();
    v.input.doubleClick(...domainAt(v, Gene, 'DBD'));
    await vi.advanceTimersByTimeAsync(1000);
    const after = domainMarks(v, Gene, 'DBD');
    const lo = Math.min(...after.marks.map((m) => m[0]));
    const hi = Math.max(...after.marks.map((m) => m[1]));
    expect(hi - lo).toBeGreaterThan(0.95 * after.box.width);
    // the protein leads: it shows the domain, residues 93 to 292
    const p = rangeOf(v, Protein);
    expect(p[0]).toBeCloseTo(93, 0);
    expect(p[1]).toBeCloseTo(292, 0);
  });
});

describe('zooming out on Protein on its gene', () => {
  it('a double-click on a domain on the gene focuses it through the leader; on it again, back: both sections', async () => {
    const v = await make();
    const protein = rangeOf(v, Protein);
    const geneRange = rangeOf(v, Gene);
    v.input.doubleClick(...domainAt(v, Gene, 'TET'));
    expect(rangeOf(v, Protein)).not.toEqual(protein);
    expect(v.view.isFollowing(Gene)).toBe(true);
    await vi.advanceTimersByTimeAsync(500);
    v.input.doubleClick(...domainAt(v, Gene, 'TET'));
    expect(rangeOf(v, Protein)).toEqual(protein);
    expect(rangeOf(v, Gene)).toEqual(geneRange);
    // forward, then Shift+double-click zooms out about the pointer; - more; 0 fits; Alt+Left back
    expect(v.input.key({ key: 'ArrowRight', altKey: true })).toBe(true);
    const focused = rangeOf(v, Protein);
    await vi.advanceTimersByTimeAsync(500);
    v.input.doubleClick(...domainAt(v, Protein, 'TET'), { shiftKey: true });
    const out = rangeOf(v, Protein);
    expect(out[1] - out[0]).toBeCloseTo(2 * (focused[1] - focused[0]));
    v.input.key({ key: '-' });
    expect(rangeOf(v, Protein)[1] - rangeOf(v, Protein)[0]).toBeGreaterThan(out[1] - out[0]);
    v.input.key({ key: '0' });
    expect(rangeOf(v, Protein)).toEqual(protein);
    v.input.key({ key: 'ArrowLeft', altKey: true });
    expect(rangeOf(v, Protein)[1] - rangeOf(v, Protein)[0]).toBeGreaterThan(out[1] - out[0]);
  });

  itWithin('a back step (the command, then the frame) under 4 ms', async (within) => {
    const v = await make(false);
    // a history to walk: focuses, back and forward in turn
    for (let i = 0; i < 36; i++)
      v.runCommand('view.focus', { section: Protein, focus: { start: 10 * i, end: 10 * i + 40 } });
    let forward = false;
    const step = () => {
      v.runCommand(forward ? 'view.forward' : 'view.back', { section: Protein });
      if (!(forward ? v.view.canForward : v.view.canBack)) forward = !forward;
      v.draw();
    };
    for (let i = 0; i < 5; i++) step();
    bestOf(within, 4, step);
    expect(v.view.history.value.back + v.view.history.value.forward).toBe(36);
  });
});
