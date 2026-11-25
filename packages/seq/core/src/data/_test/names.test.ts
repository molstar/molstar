/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it, vi } from 'vitest';
import { asyncSeqSource } from '../async-source.js';
import { findNames, type NamedSource } from '../names.js';
import { abortError, type Coverage, type NameHit, type OpenedSource, SeqNodeKeyField } from '../protocol.js';
import { checkCoverage, type SeqData } from '../seq-format.js';
import { openSeq } from '../seq-runtime.js';

const coverage: Coverage = {
  systems: [{ name: 'chr17', length: 1000 }],
  shapes: [{ system: 'chr17', shape: 'ranges' }],
  levels: [1],
  names: ['Name', 'ID'],
};

/** A source finding among its items by `Name`, exact first, then by prefix when asked (what a GFF3 read whole does). */
function source(items: readonly { id: string; name: string; start: number }[], fail?: Error): OpenedSource {
  return {
    coverage,
    read: () => Promise.reject(new Error('not read here')),
    dispose() {},
    findNames(text, options, signal) {
      if (signal.aborted) return Promise.reject(abortError());
      if (fail) return Promise.reject(fail);
      const t = text.toLowerCase();
      const hit = (i: (typeof items)[number]): NameHit => ({
        system: 'chr17',
        start: i.start,
        end: i.start + 10,
        id: i.id,
        label: i.name,
        field: 'Name',
      });
      const exact = items.filter((i) => i.name.toLowerCase() === t).map(hit);
      const prefix = options.prefix
        ? items.filter((i) => i.name.toLowerCase() !== t && i.name.toLowerCase().startsWith(t)).map(hit)
        : [];
      return Promise.resolve([...exact, ...prefix].slice(0, options.max ?? 20));
    },
  };
}

const genes: NamedSource = {
  address: 'genes',
  source: source([
    { id: 'ENSG1', name: 'TP53', start: 100 },
    { id: 'ENSG2', name: 'TP53I3', start: 300 },
    { id: 'ENSG3', name: 'TP63', start: 500 },
  ]),
};
const more: NamedSource = { address: 'more', source: source([{ id: 'X1', name: 'tp53', start: 700 }]) };
const plain: NamedSource = { address: 'plain', source: { coverage, read: genes.source.read, dispose() {} } };

describe('findNames (EC-228)', () => {
  const signal = new AbortController().signal;

  it('merges the sources that serve it, in source order, each hit with its source', async () => {
    const hits = await findNames([plain, genes, more], ' tp53 ', {}, signal);
    expect(hits.map((h) => [h.source, h.id])).toEqual([
      ['genes', 'ENSG1'],
      ['more', 'X1'],
    ]);
    expect(hits[0]).toMatchObject({ system: 'chr17', start: 100, end: 110, label: 'TP53', field: 'Name' });
  });

  it('by prefix when asked, exact first; at most max together', async () => {
    const hits = await findNames([genes, more], 'TP5', { prefix: true }, signal);
    expect(hits.map((h) => h.id)).toEqual(['ENSG1', 'ENSG2', 'X1']);
    expect((await findNames([genes, more], 'TP53', { prefix: true, max: 2 }, signal)).map((h) => h.id)).toEqual([
      'ENSG1',
      'ENSG2',
    ]);
  });

  it('blank text finds nothing; a source that fails is left out; an abort rejects', async () => {
    expect(await findNames([genes], '  ', {}, signal)).toEqual([]);
    const broken: NamedSource = { address: 'broken', source: source([], new Error('read failed')) };
    expect((await findNames([broken, genes], 'TP63', {}, signal)).map((h) => h.id)).toEqual(['ENSG3']);
    const ac = new AbortController();
    ac.abort();
    await expect(findNames([genes], 'TP53', {}, ac.signal)).rejects.toMatchObject({ name: 'AbortError' });
  });

  it('a coverage lists the fields searched as names (strings)', () => {
    expect(checkCoverage(coverage)).toBeUndefined();
    expect(checkCoverage({ ...coverage, names: [1] as unknown as string[] })).toBeDefined();
  });
});

describe('a seq source’s nodes found by key', () => {
  /** Two regions of a system; T1 brought by both, T2 by the second. */
  const data = (names?: readonly string[]): SeqData => ({
    coverage: {
      systems: [{ name: 'chr', length: 200 }],
      shapes: [{ system: 'chr', shape: 'letters' }],
      levels: [1],
      ...(names ? { names } : {}),
    },
    pairs: [0, 100].map((start) => ({
      request: { system: 'chr', start, end: start + 100, level: 1, shapes: ['letters' as const] },
      answer: {
        ok: true as const,
        tables: { letters: { form: 'letters' as const, start, letters: 'A'.repeat(100) } },
        nodes: (start ? ['T1', 'T2'] : ['T1']).map((key) => ({ key, node: { kind: 'ranges' } })),
      },
    })),
  });
  const signal = new AbortController().signal;

  it('finds a key at the first region whose answer brings it, when its coverage lists `key` in `names`', async () => {
    const s = openSeq(data([SeqNodeKeyField]));
    expect(await s.findNames!('T1', {}, signal)).toEqual([
      { system: 'chr', start: 0, end: 100, id: 'T1', label: 'T1', field: 'key' },
    ]);
    expect(await s.findNames!('t2', {}, signal)).toEqual([
      { system: 'chr', start: 100, end: 200, id: 'T2', label: 'T2', field: 'key' },
    ]);
    // a key no answer brings: none, which is the source's word that it has no such node
    expect(await s.findNames!('T9', {}, signal)).toEqual([]);
    expect((await s.findNames!('T', { prefix: true }, signal)).map((h) => h.id)).toEqual(['T1', 'T2']);
  });

  /** One system's pairs, each bringing the keys given, at the level given. */
  const pairs = (rows: readonly [start: number, end: number, keys: string[], level?: number][]): SeqData => ({
    coverage: {
      systems: [{ name: 'chr', length: 10_000 }],
      shapes: [{ system: 'chr', shape: 'letters' }],
      levels: [1, 16],
      names: [SeqNodeKeyField],
    },
    pairs: rows.map(([start, end, keys, level = 1]) => ({
      request: { system: 'chr', start, end, level, shapes: ['letters' as const] },
      answer: {
        ok: true as const,
        tables: { letters: { form: 'letters' as const, start, letters: 'A'.repeat(end - start) } },
        nodes: keys.map((key) => ({ key, node: { kind: 'ranges' } })),
      },
    })),
  });

  it('keeps keys that differ only by case apart: the exact one first, then the others of its case-folded form', async () => {
    const s = openSeq(
      pairs([
        [0, 100, ['T1']],
        [100, 200, ['t1']],
      ]),
    );
    expect((await s.findNames!('t1', {}, signal)).map((h) => [h.id, h.start])).toEqual([
      ['t1', 100],
      ['T1', 0],
    ]);
    expect((await s.findNames!('T1', {}, signal)).map((h) => h.id)).toEqual(['T1', 't1']);
  });

  it('a key brought by far regions is found at the first one, never the span between (bounded work)', async () => {
    const s = openSeq(
      pairs([
        [0, 100, ['T1']],
        [100, 9000, []],
        [9000, 9100, ['T1']],
      ]),
    );
    expect(await s.findNames!('T1', {}, signal)).toEqual([
      { system: 'chr', start: 0, end: 100, id: 'T1', label: 'T1', field: 'key' },
    ]);
  });

  it('a key is found at the finest level that brings it', async () => {
    const s = openSeq(
      pairs([
        [0, 1600, ['T1'], 16],
        [800, 900, ['T1']],
      ]),
    );
    expect((await s.findNames!('T1', {}, signal)).map((h) => [h.start, h.end])).toEqual([[800, 900]]);
  });

  it('serves no find when its coverage does not say so', () => {
    expect(openSeq(data()).findNames).toBeUndefined();
  });

  it('a delayed source finds after its delay, and an abort rejects at once', async () => {
    vi.useFakeTimers();
    try {
      // real timers, faked: the store's clock as a test drives it
      const clock = {
        now: () => Date.now(),
        setTimeout: (fn: () => void, ms: number) => setTimeout(fn, ms),
        clearTimeout: (h: unknown) => clearTimeout(h as number),
      };
      const s = asyncSeqSource(data([SeqNodeKeyField]), { clock, delayMs: 300 });
      let hits: readonly NameHit[] | undefined;
      void s.findNames!('T2', {}, signal).then((h) => {
        hits = h;
      });
      await vi.advanceTimersByTimeAsync(299);
      expect(hits).toBeUndefined();
      await vi.advanceTimersByTimeAsync(1);
      expect(hits?.map((h) => h.id)).toEqual(['T2']);
      const abort = new AbortController();
      const asked = s.findNames!('T1', {}, abort.signal);
      abort.abort();
      await expect(asked).rejects.toMatchObject({ name: 'AbortError' });
    } finally {
      vi.useRealTimers();
    }
  });
});
