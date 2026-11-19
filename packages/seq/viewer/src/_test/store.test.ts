/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type FormatRuntime, SeqRuntime } from '@molstar/seqx-core/data';
import { RecordsRuntime } from '@molstar/seqx-core/data/records';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { resolveSpec } from '@molstar/seqx-core/resolve';
import { DomHost } from '@molstar/ui/dom';
import { describe, expect, it, vi } from 'vitest';
import { viewerConfig } from '../config.js';
import { DomHostClock, StoreModel } from '../models/store.js';
import { SpecModel } from '../models/spec.js';
import { Gene, Protein, proteinOnGene } from './fixture.js';

const assets = {
  readText: () => Promise.reject(new Error('no network in tests')),
  readBytes: () => Promise.reject(new Error('no network in tests')),
};

describe('StoreModel', () => {
  it('owns a store made with the configuration’s bounds, timed by DomHost', () => {
    const m = new StoreModel(viewerConfig({ store: { inFlight: 2, tileBins: 256 } }));
    expect(m.store.settings).toMatchObject({ maxInFlight: 2, tileBins: 256, maxBytes: 64 * 1024 * 1024 });
    const now = vi.spyOn(DomHost, 'now').mockReturnValue(42);
    expect(DomHostClock.now()).toBe(42);
    now.mockRestore();
    m.retire();
  });

  it('attaches the model’s sources, keeps an unchanged one, opens a changed one again, detaches one gone', async () => {
    const spec = new SpecModel(proteinOnGene());
    const m = new StoreModel();
    await m.sync(spec.model.value, assets);
    const protein = m.store.coverage(Protein);
    expect(protein?.systems[0].length).toBe(4);
    expect(m.store.coverage(Gene)?.systems[0].length).toBe(24);
    // a change elsewhere keeps the protein's source as it was
    spec.apply(proteinOnGene({ geneLength: 21 }));
    await m.sync(spec.model.value, assets);
    expect(m.store.coverage(Protein)).toBe(protein);
    expect(m.store.coverage(Gene)?.systems[0].length).toBe(21);

    const b = createSeqSpec();
    b.data.sequence('p', { alphabet: 'protein', letters: 'MEEP' });
    await m.sync(resolveSpec(b.build()), assets);
    expect(m.store.coverage(Gene)).toBeUndefined();
    expect(m.issues.value).toEqual([]);
    m.retire();
  });

  it('a source that cannot be opened is an issue at its address', async () => {
    const b = createSeqSpec();
    const src = b.data.source('remote', { format: 'seq', uri: 'data/remote.seq' });
    b.data.sequence('g', { alphabet: 'dna', source: src });
    const m = new StoreModel();
    await m.sync(resolveSpec(b.build(), { runtimes: [SeqRuntime] }), assets);
    expect(m.issues.value.map((i) => [i.code, i.address])).toEqual([['seq.source.unreadable', '/data/remote']]);
    m.retire();
  });

  it('an opened source’s own issues (rows a file of records left out) are reported at its address', async () => {
    const b = createSeqSpec();
    const src = b.data.source('table', { format: 'tsv', value: 'start\tend\n1\t2\n5\t3\nx\t4\n', system: 'p' });
    const p = b.data.sequence('p', { alphabet: 'protein', letters: 'MEEPQSD' });
    b.data.ranges('r', { target: p, source: src });
    const m = new StoreModel();
    await m.sync(resolveSpec(b.build(), { runtimes: [SeqRuntime, RecordsRuntime] }), assets);
    expect(m.store.coverage('/data/table')).toBeDefined();
    expect(m.issues.value.map((i) => [i.code, i.address, i.severity])).toEqual([
      ['seq.source.records.order', '/data/table', 'warning'],
      ['seq.source.records.value', '/data/table', 'warning'],
    ]);
    expect(m.issues.value[0].message).toMatch(/1 row \(line 3\)/);
    m.retire();
  });

  it('a failed opening is tried again by the next sync; a source detached takes its issue with it', async () => {
    let opens = 0;
    const flaky: FormatRuntime = {
      name: 'flaky',
      formats: ['flaky'],
      open: async () => {
        if (++opens === 1) throw new Error('not yet');
        return {
          coverage: { systems: [{ name: 'g', length: 10 }], shapes: [], levels: [1] },
          read: () => Promise.reject(new Error('unused')),
          dispose: () => {},
        };
      },
    };
    const spec = (withSource: boolean) => {
      const b = createSeqSpec();
      if (withSource) {
        const src = b.data.source('remote', { format: 'flaky', uri: 'r' });
        b.data.sequence('g', { alphabet: 'dna', source: src });
      } else b.data.sequence('g', { alphabet: 'dna', letters: 'ACGT' });
      return resolveSpec(b.build(), { runtimes: [flaky] });
    };
    const m = new StoreModel();
    const model = spec(true);
    await m.sync(model, assets);
    expect(m.issues.value.map((i) => i.address)).toEqual(['/data/remote']);
    await m.sync(model, assets);
    expect(opens).toBe(2);
    expect(m.issues.value).toEqual([]);
    expect(m.store.coverage('/data/remote')?.systems[0].length).toBe(10);
    // gone: detached
    await m.sync(spec(false), assets);
    expect(m.store.coverage('/data/remote')).toBeUndefined();

    // a failure no longer wanted takes its issue with it
    opens = 0;
    const n = new StoreModel();
    await n.sync(spec(true), assets);
    expect(n.issues.value.length).toBe(1);
    await n.sync(spec(false), assets);
    expect(n.issues.value).toEqual([]);
    m.retire();
    n.retire();
  });
});
