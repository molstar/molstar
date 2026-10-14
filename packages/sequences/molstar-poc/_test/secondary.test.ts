/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `molstar:secondary_structure` over a recorded structure (2OCJ, p53's core domain): helices, strands and turns as
// ranges on the polymer units, read from the model's own assignment; no network.

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import { compile } from '@molstar/seq-core/format';
import { expect, it } from 'vitest';
import { loaders, registries } from './fixtures.js';

it("gives the ranges of the structure's helices, strands and turns on its chains, as the model has them", async () => {
  const b = createSeqViewSpec();
  const m = b.data.complex({
    ref: 'c',
    provider: 'molstar:structure',
    source: { service: 'pdbe', resource: 'structure', id: '2ocj' },
  });
  b.data.computed({ ref: 'ss', provider: 'molstar:secondary_structure', options: { of: m } } as never);
  b.view({ ref: 'v', label: 'V' }).track({ label: '' }).ruler();
  const c = await compile(b.getDocument(), registries(), loaders());
  expect(c.issues.filter((i) => i.severity === 'error')).toEqual([]);
  const a = c.set.annotations.find((x) => x && c.set.strings.get(x.key).endsWith('ss'));
  expect(a?.count).toBeGreaterThan(5);
  const kinds = new Set(Array.from({ length: a?.count ?? 0 }, (_, r) => c.set.strings.get(a?.category?.[r] ?? -1)));
  expect([...kinds].every((k) => ['helix', 'strand', 'turn'].includes(k))).toBe(true);
  expect(kinds.has('strand')).toBe(true);
  // Every range is on a polymer unit and within it.
  for (let r = 0; a && r < a.count; r++) expect(a.stop[r]).toBeGreaterThan(a.start[r]);
});
