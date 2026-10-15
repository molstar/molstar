/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// LM1's letter-map cases, for their heap guard and its quick twin (letter-maps.heap.ts, letter-maps.quick.ts): the
// conformance kit run on `letters`, `rows` and `emphasis` over documents with a complement map, forward and as a
// reverse complement, with the kit's `heap` hook given by the caller (a deep guard's measurement, or a quick case).
// Each call loads the viewer afresh (`vi.resetModules`).

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import { vi } from 'vitest';
import type { AllocOptions } from '../../core/_test/model/alloc.js';
import type { LayerRenderer } from '../src/index.js';

/**
 * The quick twin's counts for each frame the kit draws (there are a dozen a case, each measured with its control): a
 * frame of a few tens of microseconds, windows of 500 steps. A warm-up allowed to end at 1,500 steps stopped on a
 * plateau of 144 bytes a step in one frame (2026-10-15), so it runs at least 4,000.
 */
export const KIT: AllocOptions = { warmSteps: 1_000, steps: 250, windows: 3 };

export const CASES: [string, string][] = [
  ['letters', 'all types'],
  ['letters', '20,000 positions'],
  ['letters', 'gene'],
  ['rows', 'family'],
  ['emphasis', 'emphasized letters'],
];

/** 2,000 nucleotides, every third emphasized, over the letters of their own track. */
function emphasized() {
  const b = createSeqViewSpec();
  const s = b.data.sequence({ ref: 's', label: 'S', alphabet: 'dna', value: 'ACGTTGCAGGCT'.repeat(167) });
  const start = Array.from({ length: 600 }, (_, k) => 3 * k);
  const f = b.data.annotation({
    ref: 'f',
    geometry: 'ranges',
    target: s,
    value: { start, stop: start.map((x) => x + 2) },
  });
  b.view({ ref: 'v', label: 'V' }).track({ label: 'S' }).letters().emphasis({ data: f });
  return b.getDocument();
}

export async function run(
  name: string,
  doc: string,
  reversed: boolean,
  heap: (frame: () => void) => void | Promise<void>,
) {
  vi.resetModules();
  const { KernelRenderers } = await import('../src/index.js');
  const { rendererConformance } = await import('../src/testing/index.js');
  const { asText, directed, mapped } = await import('./conformance-fixtures.js');
  const { opened } = await import('./layout-fixtures.js');
  const { allTypes, gene, scaleDoc } = await import('./prepare-fixtures.js');
  const { case3 } = await import('../../core/_test/format/cases.js');
  const make = {
    'all types': () => allTypes(),
    '20,000 positions': () => scaleDoc({ length: 20_000 }),
    gene: () => gene(),
    family: case3,
    'emphasized letters': emphasized,
  }[doc] as () => ReturnType<typeof allTypes>;
  const all = KernelRenderers as LayerRenderer[];
  const entry = all.find((r) => r.name === name) as LayerRenderer;
  return rendererConformance(entry, {
    open: async ({ geometry, entry: has }) => {
      const d = reversed ? directed(mapped(make()), -1, true) : mapped(make());
      return (
        await opened(geometry === 'wrapped' ? asText(d) : d, {
          renderers: has ? [...all] : all.filter((r) => r !== entry),
        })
      ).ctx;
    },
    heap,
    items: false,
  });
}

// The reverse complement goes through the same map lookups for every letters document and for emphasis, so it is kept
// on one letters document and on rows (dropped on 2026-10-12: letters on 20,000 positions and on gene, emphasis;
// `HEAP_DEEP` runs them).
export const REVERSED = new Set(['letters on all types', 'rows on family']);
