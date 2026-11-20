/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readFileSync } from 'node:fs';
import { FeatureViews } from '@molstar/seqx-examples';
import { PackageFolder } from '@molstar/seqx-examples/node-assets';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { DomHost } from '@molstar/ui/dom';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { projectOntoView } from '@molstar/seqx-core/view';
import {
  composeViewer,
  type QueryReport,
  QueryExamples,
  QueryRunner,
  readableRanges,
  rerunWhileWaiting,
  runQuery,
} from '../compose.js';

// The Query stories, headless on a fake clock, composed as the stories compose them: every example's queries run
// through the query box's own path, each with its printed form, its result's kind and its location. The files are read
// as committed, at once (a read of the disk is not on the fake clock).

const viewOf = (id: string) => FeatureViews.find((v) => v.id === id)!;

let viewer: SeqViewer | undefined;

// a frame tick asked for by an earlier file on timers since dropped would never run: this file's own host function
const raf = DomHost.requestAnimationFrame;
beforeEach(() => {
  DomHost.requestAnimationFrame = (fn) => {
    setTimeout(() => fn(performance.now()), 16);
  };
  vi.useFakeTimers({ toFake: ['setTimeout', 'clearTimeout', 'Date', 'performance'] });
});

afterEach(() => {
  DomHost.requestAnimationFrame = raf;
  viewer?.retire();
  viewer = undefined;
  vi.useRealTimers();
});

async function make(id: string, wait = 200) {
  const view = viewOf(id);
  const read = (uri: string) => readFileSync(`${PackageFolder}${uri}`);
  const v = composeViewer(view, {
    ports: {
      assets: () => ({
        reader: {
          readText: async (uri) => read(uri).toString('utf8'),
          readBytes: async (uri) => new Uint8Array(read(uri)),
        },
      }),
    },
  });
  viewer = v;
  v.resize(720);
  await v.load(view.build());
  await vi.advanceTimersByTimeAsync(wait);
  return v;
}

/** The query alignment's rows band in each section: its box's height and its rows' (every row at its scale). */
function rowsBands(v: SeqViewer) {
  v.draw();
  return v.frame!.layout.sections.flatMap((s) => {
    const lay = s.lays.flat().find((l) => l?.type.name === 'rows')?.lay as
      | { box: { height: number }; rowScale: number; display: { count: number } }
      | undefined;
    return lay ? [{ box: lay.box.height, rows: lay.display.count * lay.rowScale }] : [];
  });
}

describe('Query, headless', () => {
  it('prints what it parses and says what a query is', async () => {
    const v = await make('query');
    const r = runQuery(v, '/data/p53/( range  :to 292   :from 94 )', 'select');
    expect(r).toMatchObject({
      status: 'resolved',
      printed: '/data/p53/(range :from 94 :to 292)',
      kind: 'location',
      size: '199 positions on 1 system',
      systems: [{ system: '/data/p53', ranges: [[93, 292]] }],
    });
    expect(v.addresses.selectionAddress()).toBe('/data/p53/(range :from 94 :to 292)');
    // a later query replaces what the earlier one set
    runQuery(v, '/data/p53/(at :label 175)', 'highlight');
    expect(v.selection.value.value.shown.parts).toEqual([]);
    expect(v.highlight.get('query')?.location.parts[0].ranges).toEqual(new Int32Array([174, 175]));
  });

  it('says a rect’s size as rects, not as positions', async () => {
    const v = await make('query-alignment');
    const r = runQuery(v, '/data/four/(rect :a (row :key "HBA1_BOSGR/26-136") :b (range :from 20 :to 60))', 'select');
    expect(r.size).toBe('1 rect');
  });

  it('alignment: a row form over the seed read from its file names its row, as over the four inline rows', async () => {
    const v = await make('query-alignment');
    const r = runQuery(v, '/data/globins/(row :key "MYG_CYPCA/23-136")', 'select');
    expect(r.status).toBe('resolved');
    expect(r.issues).toEqual([]);
    expect(runQuery(v, '/data/globins/(row :key "NOT_A_ROW")', 'select').issues.map((i) => i.code)).toEqual([
      'seq.query.row.unknown',
    ]);
  });

  it('alignment: a rows band never taller than its rows, the same across a resize of the frame', async () => {
    const v = await make('query-alignment');
    const bands = rowsBands(v);
    expect(bands.length).toBeGreaterThan(0);
    for (const b of bands) expect(b.box).toBeLessThanOrEqual(b.rows);
    for (const w of [1440, 360, 1000, 720]) {
      v.resize(w);
      await vi.advanceTimersByTimeAsync(200);
      v.draw();
      for (const b of rowsBands(v)) expect(b.box).toBeLessThanOrEqual(b.rows);
    }
    // and settled: the canvas as tall as its content, which stays put
    const h = v.frame!.layout.height;
    v.resize(1440);
    await vi.advanceTimersByTimeAsync(200);
    v.resize(720);
    await vi.advanceTimersByTimeAsync(200);
    v.draw();
    expect(v.frame!.layout.height).toBe(h);
  });

  it('alignment: a query waiting for the rows read from the file runs again once they arrive', async () => {
    const v = await make('query-alignment', 0);
    const q = '/data/globins/(row :key "HBA1_BOSGR/26-136")';
    let last: QueryReport | undefined = runQuery(v, q, 'select');
    expect(last.status).toBe('waiting');
    const off = rerunWhileWaiting(
      v,
      () => last,
      () => {
        last = runQuery(v, q, 'select');
      },
    );
    await vi.advanceTimersByTimeAsync(1000);
    off();
    expect(last.status).toBe('resolved');
  });

  it('says why a query does not resolve, with its hint', async () => {
    const v = await make('query');
    expect(runQuery(v, '', 'select').status).toBe('empty');
    const syntax = runQuery(v, 'p53/(range', 'select');
    expect(syntax.status).toBe('failed');
    expect(syntax.printed).toBeUndefined();
    expect(syntax.issues[0].code).toBe('repr.tree.address.syntax');
    const backwards = runQuery(v, '/data/p53/(range :from 292 :to 94)', 'select');
    expect(backwards.status).toBe('failed');
    expect(backwards.issues.map((i) => [i.code, i.hint !== undefined])).toEqual([['seq.query.range.backwards', true]]);
    const none = runQuery(v, '/data/domains/(item :id GONE)', 'select');
    expect(none.status).toBe('missing');
    expect(none.printed).toBe('/data/domains/(item :id GONE)');
    expect(v.selection.value.value.shown.parts).toEqual([]);
  });

  // what each example's queries name, as positions: the headless check of the stories' results
  const P = '/data/p53';
  const Expected: Record<string, readonly { query: string; kind: string; systems?: unknown; rects?: unknown }[]> = {
    category: [
      {
        query: '/data/domains/(item :label DNA-binding)',
        kind: 'elements',
        systems: [{ system: P, ranges: [[93, 292]] }],
      },
      {
        query: '/data/domains/(item :category [binding oligomer])',
        kind: 'elements',
        systems: [
          {
            system: P,
            ranges: [
              [93, 292],
              [323, 355],
            ],
          },
        ],
      },
      {
        query: '/data/domains/(item :id [TAD CTD])',
        kind: 'elements',
        systems: [
          {
            system: P,
            ranges: [
              [0, 61],
              [363, 393],
            ],
          },
        ],
      },
    ],
    ranges: [
      { query: '/data/p53/(range :from 94 :to 292)', kind: 'location', systems: [{ system: P, ranges: [[93, 292]] }] },
      { query: '/data/p53/(at :label 175)', kind: 'location', systems: [{ system: P, ranges: [[174, 175]] }] },
      {
        query: '/data/p53/(range :from 0 :numbering index :to 9)',
        kind: 'location',
        systems: [{ system: P, ranges: [[0, 10]] }],
      },
      {
        query: '/data/tp53/(range :from 101 :to 400)',
        kind: 'location',
        systems: [{ system: '/data/tp53', ranges: [[100, 400]] }],
      },
    ],
    'on-gene': [
      { query: '/data/domains/(item :id DBD)', kind: 'elements', systems: [{ system: P, ranges: [[93, 292]] }] },
      { query: '/data/p53/(range :from 94 :to 292)', kind: 'location', systems: [{ system: P, ranges: [[93, 292]] }] },
      {
        query: '/data/p53/(or :of [(at :label 175) (at :label 248) (at :label 273)])',
        kind: 'location',
        systems: [
          {
            system: P,
            ranges: [
              [174, 175],
              [247, 248],
              [272, 273],
            ],
          },
        ],
      },
    ],
    combined: [
      {
        query: '/data/p53/(and :of [(range :from 94 :to 292) (not :of (range :from 150 :to 200))])',
        kind: 'location',
        systems: [
          {
            system: P,
            ranges: [
              [93, 149],
              [200, 292],
            ],
          },
        ],
      },
      {
        query: '/data/p53/(or :of [(range :from 1 :to 40) (range :from 320 :to 360)])',
        kind: 'location',
        systems: [
          {
            system: P,
            ranges: [
              [0, 40],
              [319, 360],
            ],
          },
        ],
      },
      {
        query: '/data/p53/(and :of [(range :from 94 :to 292) (range :from 250 :to 360)])',
        kind: 'location',
        systems: [{ system: P, ranges: [[249, 292]] }],
      },
    ],
    alignment: [
      {
        query: '/data/globins/(row :key "HBA1_BOSGR/26-136")',
        kind: 'location',
        systems: [{ system: '/data/globins/(rows)', ranges: [[29, 30]] }],
      },
      {
        query: '/data/globins/(range :from 20 :to 60)',
        kind: 'location',
        systems: [{ system: '/data/globins', ranges: [[19, 60]] }],
      },
      {
        query:
          '/data/globins/(rect :a (row :key ["HBA1_BOSGR/26-136" "HBB0_MOUSE/26-142"]) :b (range :from 20 :to 60))',
        kind: 'location',
        rects: [
          {
            a: {
              system: '/data/globins/(rows)',
              ranges: [
                [16, 17],
                [29, 30],
              ],
            },
            b: { system: '/data/globins', ranges: [[19, 60]] },
          },
        ],
      },
      {
        query: '/data/four/(row :key "HBA1_BOSGR/26-136")',
        kind: 'location',
        systems: [{ system: '/data/four/(rows)', ranges: [[0, 1]] }],
      },
    ],
  };

  for (const ex of QueryExamples) {
    it(`${ex.id}: every query of the example resolves, prints as written and names its location`, async () => {
      const v = await make(ex.view);
      const queries = ex.queries.filter((q) => !q.issue);
      for (const q of queries) {
        const r = runQuery(v, q.query, ex.emphasis);
        expect(r.status, `${q.label}: ${r.issues.map((x) => x.message)}`).toBe('resolved');
        expect(r.printed, q.label).toBe(q.query);
        expect(r.issues, q.label).toEqual([]);
        const want = Expected[ex.id].find((w) => w.query === q.query)!;
        expect(want, q.label).toBeDefined();
        expect(r.kind, q.label).toBe(want.kind);
        expect(r.systems, q.label).toEqual(want.systems ?? []);
        expect(r.rects, q.label).toEqual(want.rects ?? []);
        // shown as the example says: selected, or on the query's own channel
        if (ex.emphasis === 'select') expect(v.addresses.selectionAddress(), q.label).toBe(q.query);
        else expect(v.highlight.get('query')?.address, q.label).toBe(q.query);
      }
      expect(queries.length).toBe(Expected[ex.id].length);
    });
  }

  it('combined: a query that mixes result types is an issue, and shows nothing', async () => {
    const v = await make('query');
    const mixed = QueryExamples.find((e) => e.id === 'combined')!.queries.find((q) => q.issue)!;
    const r = runQuery(v, mixed.query, 'highlight');
    expect(r.status).toBe('failed');
    expect(r.issues.map((i) => i.code)).toEqual(['seq.query.form.type']);
    expect(v.highlight.get('query')).toBeUndefined();
  });

  it('on-gene: the protein’s result is carried onto the gene by the codon map, exactly', async () => {
    const v = await make('query');
    const gene = (address: string) => {
      runQuery(v, address, 'select');
      const shown = projectOntoView(v.spec.model.value, v.view.resolved.value!, v.selection.value.value.location);
      return shown.parts.map((p) => [p.system, Array.from(p.ranges)]);
    };
    // residue 175 is codon 174: nucleotides 522 to 524 of the coding sequence, in the second exon (minus strand)
    expect(gene('/data/p53/(at :label 175)')).toEqual([['/data/tp53', [837, 840]]]);
    // residues 94 to 292: the second exon's tail and the third's head, in genomic order
    expect(gene('/data/p53/(range :from 94 :to 292)')).toEqual([['/data/tp53', [486, 900, 960, 1143]]]);
    // the result itself stays on the protein
    expect(v.selection.value.value.location.parts.map((p) => p.system)).toEqual(['/data/p53']);
    // and it is drawn: a rect on each section
    v.view.fit('/views/main/protein');
    v.draw();
    expect(v.frame!.lists.passes.selection.rects.n).toBeGreaterThanOrEqual(2);
    runQuery(v, '', 'select');
    v.draw();
    expect(v.frame!.lists.passes.selection.rects.n).toBe(0);
  });

  it('counts an issue’s offset from the text as typed, leading spaces included', async () => {
    const v = await make('query');
    const bare = runQuery(v, 'p53/(range', 'select').issues[0];
    const typed = runQuery(v, '   p53/(range', 'select').issues[0];
    expect(bare.offset).toBeDefined();
    expect(typed.offset).toBe(bare.offset! + 3);
  });

  it('switching Highlight back to Select moves the result: one is cleared as the other is set', async () => {
    const v = await make('query');
    const q = '/data/p53/(range :from 94 :to 292)';
    runQuery(v, q, 'highlight');
    expect(v.highlight.get('query')?.address).toBe(q);
    expect(v.addresses.selectionAddress()).toBeUndefined();
    runQuery(v, q, 'select');
    expect(v.highlight.get('query')).toBeUndefined();
    expect(v.addresses.selectionAddress()).toBe(q);
    runQuery(v, q, 'highlight');
    expect(v.addresses.selectionAddress()).toBeUndefined();
    expect(v.highlight.get('query')?.address).toBe(q);
  });

  it('leaves the view alone when the query still names what it shows (one change a keystroke)', async () => {
    const v = await make('query');
    runQuery(v, '/data/p53/(range :from 94 :to 292)', 'select');
    const before = v.selection.value;
    runQuery(v, '  /data/p53/( range :from 94  :to 292 ) ', 'select');
    expect(v.selection.value).toBe(before);
  });

  it('a throw is a failed report with its message, and nothing stale', async () => {
    const v = await make('query');
    runQuery(v, '/data/p53/(range :from 94 :to 292)', 'select');
    const boom = vi.spyOn(v.addresses, 'locate').mockImplementation(() => {
      throw new Error('boom');
    });
    const r = runQuery(v, '/data/p53/(at :label 1)', 'select');
    boom.mockRestore();
    expect(r.status).toBe('failed');
    expect(r.issues.map((i) => i.message)).toEqual(['boom']);
    expect(v.addresses.selectionAddress()).toBeUndefined();
  });

  it('reads positions as a reader does: 1-based and inclusive', () => {
    expect(
      readableRanges([
        [93, 292],
        [174, 175],
      ]),
    ).toBe('94–292, 175');
  });
});

describe('the query box’s runner', () => {
  it('typing runs after a rest; Enter runs at once and drops the wait', () => {
    const run = vi.fn();
    const r = new QueryRunner(run, 100);
    r.soon();
    vi.advanceTimersByTime(50);
    expect(run).not.toHaveBeenCalled();
    r.now();
    expect(run).toHaveBeenCalledTimes(1);
    vi.advanceTimersByTime(200);
    expect(run).toHaveBeenCalledTimes(1);
    r.soon();
    r.dispose();
    r.now();
    vi.advanceTimersByTime(200);
    expect(run).toHaveBeenCalledTimes(1);
  });
});
