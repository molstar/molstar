/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { findNames, type NamedSource } from '../names.js';
import { abortError, type Coverage, type NameHit, type OpenedSource } from '../protocol.js';
import { checkCoverage } from '../seq-format.js';

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
