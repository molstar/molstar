/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it, vi } from 'vitest';

// The spec model over the genome chunk's load (_spec/README.md "SpecModel"): resolved again once when the load lands,
// however many times the spec changed while it loaded, and a failed load is the model's issue, never an unhandled
// rejection. A fresh module registry each case (`vi.resetModules`: this file runs isolated), so the chunk is not loaded.

async function fresh(fail = false) {
  vi.resetModules();
  if (fail)
    vi.doMock('@molstar/seqx-core/resolve/genome-chunk', () => {
      throw new Error('offline');
    });
  const lang = await import('@molstar/seqx-core/lang');
  const { SpecModel } = await import('../models/spec.js');
  const b = lang.createSeqSpec();
  b.data.assembly('g', lang.assembly('toy', [{ name: 'chr1', length: 10 }]));
  return { root: b.build(), SpecModel };
}

describe('the spec model while the genome chunk loads', () => {
  it('resolves again once it lands, a spec swapped meanwhile included', async () => {
    const { root, SpecModel } = await fresh();
    const spec = new SpecModel(root);
    const first = spec.model.value;
    expect(first.loading).toBeDefined();
    let resolved = 0;
    const stop = spec.model.subscribe(() => resolved++);
    // the spec changes while it loads: the same load, never subscribed again
    const then = vi.spyOn(first.loading!, 'then');
    expect(spec.setPreview({ ops: [{ op: 'set', at: ['data', 'g'], key: 'name', value: 'toy2' }] }).issues).toEqual([]);
    expect(spec.model.value).not.toBe(first);
    expect(spec.model.value.loading).toBe(first.loading);
    expect(then).not.toHaveBeenCalled();
    await first.loading;
    await Promise.resolve();
    const after = spec.model.value;
    expect(after.loading).toBeUndefined();
    expect(after.assemblies.get('/data/g')?.systems).toEqual(['/data/g/(sequence :name chr1)']);
    expect(resolved).toBeLessThanOrEqual(3);
    stop?.unsubscribe?.();
  });

  it('a failed load is the model’s issue', async () => {
    const { root, SpecModel } = await fresh(true);
    const spec = new SpecModel(root);
    await spec.model.value.loading;
    await Promise.resolve();
    expect(spec.model.value.issues.map((i) => i.code)).toContain('seq.resolve.genome.load');
    vi.doUnmock('@molstar/seqx-core/resolve/genome-chunk');
  });
});
