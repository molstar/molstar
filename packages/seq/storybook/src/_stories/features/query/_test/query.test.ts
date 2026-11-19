/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { FeatureViews } from '@molstar/seqx-examples';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { DomHost } from '@molstar/ui/dom';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { QueryExamples, runQuery, composeViewer } from '../compose.js';

// The Query stories, headless on a fake clock, composed as the stories compose them: every example's queries run
// through the query box's own path, each with its printed form, its result's kind and its location.

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

async function make(id: string) {
  const view = viewOf(id);
  const v = composeViewer(view);
  viewer = v;
  v.resize(720);
  await v.load(view.build());
  await vi.advanceTimersByTimeAsync(200);
  return v;
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
    const r = runQuery(v, '/data/globins/(rect :a (row :key "HBB_HUMAN/1-147") :b (range :from 20 :to 60))', 'select');
    expect(r.size).toBe('1 rect');
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
  const Expected: Record<string, readonly { kind: string; systems?: unknown; rects?: unknown }[]> = {
    category: [
      { kind: 'elements', systems: [{ system: P, ranges: [[93, 292]] }] },
      {
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
      { kind: 'location', systems: [{ system: P, ranges: [[93, 292]] }] },
      { kind: 'location', systems: [{ system: P, ranges: [[174, 175]] }] },
      { kind: 'location', systems: [{ system: P, ranges: [[0, 10]] }] },
      { kind: 'location', systems: [{ system: '/data/tp53', ranges: [[100, 400]] }] },
    ],
    'on-gene': [
      { kind: 'elements', systems: [{ system: P, ranges: [[93, 292]] }] },
      { kind: 'location', systems: [{ system: P, ranges: [[93, 292]] }] },
      {
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
      { kind: 'location', systems: [{ system: P, ranges: [[249, 292]] }] },
    ],
    alignment: [
      { kind: 'location', systems: [{ system: '/data/globins/(rows)', ranges: [[0, 1]] }] },
      { kind: 'location', systems: [{ system: '/data/globins', ranges: [[19, 60]] }] },
      {
        kind: 'location',
        rects: [
          {
            a: { system: '/data/globins/(rows)', ranges: [[0, 2]] },
            b: { system: '/data/globins', ranges: [[19, 60]] },
          },
        ],
      },
    ],
  };

  for (const ex of QueryExamples) {
    it(`${ex.id}: every query of the example resolves, prints as written and names its location`, async () => {
      const v = await make(ex.view);
      const queries = ex.queries.filter((q) => !q.label.includes('(an issue)'));
      for (const [i, q] of queries.entries()) {
        const r = runQuery(v, q.query, ex.emphasis);
        expect(r.status, `${q.label}: ${r.issues.map((x) => x.message)}`).toBe('resolved');
        expect(r.printed, q.label).toBe(q.query);
        expect(r.issues, q.label).toEqual([]);
        const want = Expected[ex.id][i];
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
    const mixed = QueryExamples.find((e) => e.id === 'combined')!.queries.find((q) => q.label.includes('(an issue)'))!;
    const r = runQuery(v, mixed.query, 'highlight');
    expect(r.status).toBe('failed');
    expect(r.issues.map((i) => i.code)).toEqual(['seq.query.form.type']);
    expect(v.highlight.get('query')).toBeUndefined();
  });

  it('on-gene: the protein’s result is drawn on the gene as well, and stays a location on the protein', async () => {
    const v = await make('query');
    v.view.fit('/views/main/protein');
    runQuery(v, '/data/p53/(range :from 94 :to 292)', 'select');
    v.draw();
    const withGene = v.frame!.lists.passes.selection.rects.n;
    runQuery(v, '', 'select');
    v.draw();
    expect(v.frame!.lists.passes.selection.rects.n).toBe(0);
    expect(withGene).toBeGreaterThanOrEqual(2);
    expect(v.selection.value.value.shown.parts.map((p) => p.system)).toEqual([]);
  });
});
