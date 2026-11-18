/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type Frame, type RangesLay, recordDrawLists } from '@molstar/seqx-core/render';
import { type ExampleView, P53Views } from '@molstar/seqx-examples';
import { nodeFixtureAssets } from '@molstar/seqx-examples/node-assets';
import { GeneFile } from '@molstar/seqx-examples/p53/data';
import { regionSeqText } from '@molstar/seqx-examples/p53/region';
import { DefaultSeqViewerSpec } from '@molstar/seqx-viewer/default-spec';
import { createSeqViewer, type SeqViewer } from '@molstar/seqx-viewer/viewer';
import { afterEach, beforeAll, describe, expect, it, vi } from 'vitest';

// Every drawn mark of a ranges layer has a hit box: on Gene by region and Protein on its gene at the fit view, headless
// on a fake clock, the hit at each mark's centre (in every lane and piece) names an item the mark draws, and hovering
// it outlines that item. At the fit view a lane's marks a device pixel apart merge into one run (`bars`): a run of
// several items is hit, and outlined, by the item under the pointer.

const DelayMs = 300;
const HoverOwner = 'test';

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

/** `delayMs`: each answer's delay, or `staggered` (the story's own: later regions later). */
async function open(
  view: ExampleView,
  delayMs: number | 'staggered' = DelayMs,
  wait = 20 * DelayMs,
): Promise<SeqViewer> {
  vi.useFakeTimers({ toFake: ['setTimeout', 'clearTimeout', 'Date', 'performance'] });
  const text = view.id === 'gene-by-region' ? regionSeqText() : gene;
  const base = DefaultSeqViewerSpec();
  viewer = createSeqViewer({
    ...base,
    ui: undefined,
    registry: [...(base.registry ?? []), ...(view.registry?.(delayMs === 'staggered' ? {} : { delayMs }) ?? [])],
    ports: { assets: () => ({ reader: { readText: async () => text } }) },
  });
  viewer.resize(720);
  void viewer.load(view.build());
  await vi.advanceTimersByTimeAsync(wait);
  viewer.draw();
  return viewer;
}

/** The hover pass's rects. */
const hoverRects = (frame: Frame) =>
  recordDrawLists(frame.lists, ['hover'])
    .split('\n')
    .filter((l) => l.startsWith('rect ')).length;

/** Each drawn mark of every ranges layer hit at its centre: what is wrong, and how many marks were checked. */
function checkMarks(v: SeqViewer) {
  const wrong: string[] = [];
  let marks = 0;
  let merged = 0;
  const layers = v.frame!.layout.sections.flatMap((s) => s.tracks.flatMap((t) => t.layers));
  for (const l of layers) {
    if (l.type.name !== 'ranges') continue;
    const lay = l.lay as RangesLay;
    for (let j = 0; j < lay.m; j++) {
      marks++;
      const x = (lay.mx0[j] + lay.mx1[j]) / 2;
      const y = lay.box.y + (lay.mlane[j] + 0.5) * lay.laneHeight;
      // the items the mark draws: its own, or those of the run's marks
      const drawn = new Set<string | number>();
      if (lay.mitem[j] >= 0) drawn.add(lay.keys[lay.mitem[j]]);
      else {
        merged++;
        for (let k = 0; k < (lay.hn ?? 0); k++)
          if (lay.hlane[k] === lay.mlane[j] && lay.hx1[k] >= lay.mx0[j] && lay.hx0[k] <= lay.mx1[j])
            drawn.add(lay.keys[lay.hitem[k]]);
      }
      const hit = v.hit(x, y);
      const where = `${l.key} mark ${j} (lane ${lay.mlane[j]}, x ${lay.mx0[j]}–${lay.mx1[j]})`;
      if (hit?.layer !== l.key || hit.element === undefined || !drawn.has(hit.element.key)) {
        wrong.push(`${where}: ${hit?.layer} ${JSON.stringify(hit?.element)}`);
        continue;
      }
      // hovered, the item is outlined
      v.hover.set(hit, HoverOwner);
      const f = v.draw() ?? v.frame!;
      if (!hoverRects(f)) wrong.push(`${where}: hovered ${hit.element.key}, nothing outlined`);
      v.hover.set(undefined, HoverOwner);
    }
  }
  return { wrong, marks, merged };
}

describe('every drawn range mark has a hit box', () => {
  it('Gene by region at the fit view: the transcripts and their domains, merged runs among them', async () => {
    const v = await open(P53Views.find((x) => x.id === 'gene-by-region')!);
    const r = checkMarks(v);
    expect(r.marks).toBeGreaterThan(10);
    // the domains' runs of several items (TAD and PRR, NLS and CTD) at this width
    expect(r.merged).toBeGreaterThan(0);
    expect(r.wrong).toEqual([]);
  });

  it('Protein on its gene at the fit view: the domains on the protein and cut into exons on the gene', async () => {
    const v = await open(P53Views.find((x) => x.id === 'protein-on-gene')!);
    const r = checkMarks(v);
    expect(r.marks).toBeGreaterThan(5);
    expect(r.wrong).toEqual([]);
  });

  it('Gene by region while it streams: after each arrival, every mark drawn is hit (the hit index is the frame’s)', async () => {
    const view = P53Views.find((x) => x.id === 'gene-by-region')!;
    // the regions answer one by one (300 ms, then 250 more each 1,024 positions), their domains a delay later
    const v = await open(view, 'staggered', 0);
    let checked = 0;
    for (let t = 0; t < 40; t++) {
      await vi.advanceTimersByTimeAsync(DelayMs / 2);
      v.draw();
      const r = checkMarks(v);
      checked += r.marks;
      expect(r.wrong, `at ${(t + 1) * (DelayMs / 2)} ms`).toEqual([]);
    }
    expect(checked).toBeGreaterThan(100);
  });

  it('Gene by region with no delay (all at once) hits the same marks as streamed', async () => {
    const view = P53Views.find((x) => x.id === 'gene-by-region')!;
    const at = (v: SeqViewer) =>
      v
        .frame!.layout.sections.flatMap((s) => s.tracks.flatMap((t) => t.layers))
        .filter((l) => l.type.name === 'ranges')
        .map((l) => {
          const lay = l.lay as RangesLay;
          return `${l.key} ${Array.from({ length: lay.m }, (_, j) => `${lay.mlane[j]}:${lay.mx0[j]}-${lay.mx1[j]}`).join(' ')}`;
        });
    const streamed = await open(view);
    const marks = at(streamed);
    expect(checkMarks(streamed).wrong).toEqual([]);
    streamed.retire();
    vi.useRealTimers();
    const once = await open(view, 0, 100);
    expect(at(once)).toEqual(marks);
    expect(checkMarks(once).wrong).toEqual([]);
  });
});
