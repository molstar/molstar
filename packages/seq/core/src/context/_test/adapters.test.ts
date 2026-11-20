/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Data adapters as registry entries (_spec/README.md "Context", Registry entries): one adapter a format, conflicts
// listed; `withAdapter` replaces a default for its formats; the resolver and the overlay bind by format with no order;
// a source's `ref` reaches the adapter (EC-185) and a values node's `field` reads its own channel.

import { channelAddress, type FormatRuntime, SeqRuntime } from '@molstar/seqx-core/data';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { resolveSpec } from '@molstar/seqx-core/resolve';
import { describe, expect, it } from 'vitest';
import { composeSeqContext, runtimesOf, type SeqRegistryEntry, withAdapter } from '../index.js';

const rt = (name: string, ...formats: string[]): FormatRuntime => ({
  name,
  formats,
  open: () => Promise.reject(new Error('unused')),
});
const entry = (runtime: FormatRuntime): SeqRegistryEntry => ({ kind: 'format-runtime', runtime });

describe('one adapter a format', () => {
  it('two runtimes of one format are an error listing every conflict, whatever their order', () => {
    const a = rt('a', 'mmcif', 'bcif');
    const b = rt('b', 'mmcif');
    const c = rt('c', 'bcif', 'pdb');
    expect(() => runtimesOf([a, b, c].map(entry))).toThrow(
      'several format runtimes for one format: mmcif (a, b); bcif (a, c)',
    );
    expect(() => composeSeqContext({ registry: [b, a].map(entry) })).toThrow(/mmcif \(b, a\)/);
    expect(runtimesOf([a, a].map(entry))).toEqual([a]);
  });

  it('withAdapter replaces the default for its formats: a runtime left with others keeps them, other entries kept', () => {
    const records = rt('records', 'tsv', 'csv', 'json');
    const fasta = rt('fasta', 'fasta');
    const other = { kind: 'track-sizing', sizing: { name: 'x', heights: () => [] } } as unknown as SeqRegistryEntry;
    const host = rt('host', 'json', 'fasta');
    const registry = withAdapter([entry(records), other, entry(fasta)], host);
    const runtimes = runtimesOf(registry);
    expect(runtimes.map((r) => [r.name, r.formats])).toEqual([
      ['records', ['tsv', 'csv']],
      ['host', ['json', 'fasta']],
    ]);
    expect(registry.filter((e) => e.kind === 'track-sizing')).toEqual([other]);
    // again with the same adapter: no change, no conflict
    expect(runtimesOf(withAdapter(registry, host)).map((r) => r.name)).toEqual(['records', 'host']);
  });

  it('a runtime listed twice is narrowed once, and an adapter may share the name of the runtime it narrows', () => {
    const records = rt('records', 'tsv', 'json');
    const twice = runtimesOf(withAdapter([entry(records), entry(records)], rt('host', 'json')));
    expect(twice.map((r) => [r.name, r.formats])).toEqual([
      ['records', ['tsv']],
      ['host', ['json']],
    ]);
    const named = runtimesOf(withAdapter([entry(records)], rt('records', 'json')));
    expect(named.map((r) => [r.name, r.formats])).toEqual([
      ['records', ['tsv']],
      ['records', ['json']],
    ]);
  });
});

describe('binding by format', () => {
  function spec(format: string) {
    const b = createSeqSpec();
    const src = b.data.source('s', { format, uri: 'entry.cif', ref: 'model-1' });
    b.data.sequence('a', { alphabet: 'protein', source: src });
    return b.build();
  }

  it('binds the one runtime of a format; several bind none, said, in any order', () => {
    const a = rt('a', 'mmcif');
    const b = rt('b', 'mmcif');
    for (const runtimes of [
      [a, b],
      [b, a],
    ]) {
      const model = resolveSpec(spec('mmcif'), { runtimes });
      expect(model.sources.get('/data/s')?.runtime).toBeUndefined();
      expect(model.issues.map((i) => [i.code, i.message])).toEqual([
        ['seq.resolve.source.runtime', 'Several runtimes read the format of the source'],
      ]);
    }
    expect(resolveSpec(spec('mmcif'), { runtimes: [a, SeqRuntime] }).sources.get('/data/s')?.runtime).toBe(a);
  });

  it('a source’s `ref` rides on its entry, for a host’s adapter (EC-185)', () => {
    const model = resolveSpec(spec('mmcif'), { runtimes: [rt('host', 'mmcif')] });
    expect(model.sources.get('/data/s')?.entry).toMatchObject({ uri: 'entry.cif', ref: 'model-1' });
  });

  it('a values node’s `field` reads its source’s channel; without one, the source', () => {
    const b = createSeqSpec();
    const src = b.data.source('s', { format: 'mmcif', uri: 'entry.cif' });
    const a = b.data.sequence('a', { alphabet: 'protein', letters: 'MEEP' });
    b.data.values('plddt', { target: a, source: src, field: 'plddt' });
    b.data.values('value', { target: a, source: src });
    const model = resolveSpec(b.build(), { runtimes: [rt('host', 'mmcif')] });
    expect(model.annotations.get('/data/plddt')?.access).toMatchObject({ source: channelAddress('/data/s', 'plddt') });
    expect(model.annotations.get('/data/value')?.access).toMatchObject({ source: '/data/s' });
  });
});
