/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A search's letters through the store (_spec/README.md "Search", Letters through the store): a region asked by a want
// of its own and waited for, a pending region never read as empty (EC-169), a failed one its target's issue.

import {
  type Coverage,
  DataStore,
  type OpenedSource,
  type RegionAnswer,
  type RegionRequest,
} from '@molstar/seqx-core/data';
import { describe, expect, it } from 'vitest';
import { storeLetters } from '../search-targets.js';

const Letters = 'ACGTTGCAAC'.repeat(10);
const coverage: Coverage = {
  systems: [{ name: 's', length: Letters.length }],
  shapes: [{ system: 's', shape: 'letters' }],
  levels: [1],
};

/** A source whose answers wait until `release`; `fail` answers a failure. */
function deferredSource(fail = false) {
  const waiting: (() => void)[] = [];
  const source: OpenedSource = {
    coverage,
    read: (r: RegionRequest) =>
      new Promise<RegionAnswer>((resolve) =>
        waiting.push(() =>
          resolve(
            fail
              ? {
                  ok: false,
                  kind: 'not-served',
                  issue: { code: 'seq.source.not_served', severity: 'warning', message: 'no' },
                }
              : {
                  ok: true,
                  tables: { letters: { form: 'letters', start: r.start, letters: Letters.slice(r.start, r.end) } },
                },
          ),
        ),
      ),
  };
  return {
    source,
    release: () => {
      for (const w of waiting.splice(0)) w();
    },
    get waiting() {
      return waiting.length;
    },
  };
}

const turn = async () => {
  for (let i = 0; i < 20; i++) await Promise.resolve();
};

describe('a search’s letters through the store', () => {
  it('EC-169: a pending region is waited for, never read as empty', async () => {
    const store = new DataStore({ settings: { tileBins: 16 } });
    const src = deferredSource();
    store.attach('/data/f', src.source);
    let got: Uint8Array | undefined;
    const read = storeLetters(store, '/data/f', 's')
      .read(5, 45, new AbortController().signal)
      .then((x) => {
        got = x;
      });
    await turn();
    // asked, not answered: nothing read yet, no zeros handed on
    expect(src.waiting).toBeGreaterThan(0);
    expect(got).toBeUndefined();
    src.release();
    await turn();
    src.release();
    await read;
    expect(String.fromCharCode(...got!)).toBe(Letters.slice(5, 45));
  });

  it('a source not open yet is the read’s failure, never letters read as empty', async () => {
    const store = new DataStore({ settings: { tileBins: 16 } });
    await expect(storeLetters(store, '/data/f', 's').read(0, 8, new AbortController().signal)).rejects.toThrow(
      /not open yet/,
    );
  });

  it('a failed region is the read’s failure (its target skipped, said)', async () => {
    const store = new DataStore({ settings: { tileBins: 16 } });
    const src = deferredSource(true);
    store.attach('/data/f', src.source);
    const read = storeLetters(store, '/data/f', 's').read(0, 8, new AbortController().signal);
    await turn();
    src.release();
    await expect(read).rejects.toThrow(/could not be read/);
  });
});
