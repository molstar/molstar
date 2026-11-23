/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Paul Pillot <paul.pillot@tandemai.com>
 */

import { parseArgs } from 'node:util';
import type { CifBlock } from '@molstar/io/reader/cif/data-model';
import { isCandidate, representative } from './candidates.ts';
import { openCcd } from './ccd.ts';
import { boundComponents, downloadEntries } from './rcsb.ts';

/**
 * Tier C's inputs (plan §3.2): the bound components of every representative entry (data API), then the entries of the
 * candidates (static BinaryCIF), into the cache. `run --tier C` fetches what is missing too; this downloads them ahead.
 */
export async function fetchTierC(ccdIds: string[], read: (id: string) => Promise<CifBlock>, cacheDir: string) {
  const reps = new Map<string, ReturnType<typeof representative>>();
  for (const id of ccdIds) reps.set(id, representative(await read(id)));
  const entries = [...new Set([...reps.values()].flatMap((r) => (r.entry ? [r.entry] : [])))];
  const bound = await boundComponents(entries, cacheDir);
  const candidates = ccdIds.filter((id) => isCandidate(id, reps.get(id)!, bound));
  const toDownload = candidates.map((id) => reps.get(id)!.entry!).filter((e) => bound.get(e) !== null);
  await downloadEntries(toDownload, cacheDir);
  return { reps, bound, candidates };
}

if (import.meta.main) {
  const { values } = parseArgs({ options: { cache: { type: 'string', default: 'cache' } } });
  const ccd = await openCcd(values.cache);
  const t0 = performance.now();
  const { candidates } = await fetchTierC(ccd.ids, (id) => ccd.read(id), values.cache);
  console.log(`${candidates.length} candidates fetched in ${Math.round((performance.now() - t0) / 1000)} s`);
}
