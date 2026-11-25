/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { applyPatch } from '@molstar/core/repr/patch';
import { hashNode, type Node } from '@molstar/core/repr/tree';
import { DataStore, openSeq, type Wanter } from '@molstar/seqx-core/data';
import { addRefPatch, checkSeqLang, SeqLanguage } from '@molstar/seqx-core/lang';
import { resolveAddressIn } from '@molstar/seqx-core/query';
import { buildOverlay, overlayQueryModel, resolveSpec } from '@molstar/seqx-core/resolve';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { build as bookmarksSpec } from '../p53/bookmarks.js';
import { GeneOffset, RegionSource, RegionSystem, regionData } from '../p53/region.js';

// Patches on streamed things, on p53's transcripts by region (the core's Resolving, Status): a bookmark naming a
// transcript, a dynamic node, is a patch like any other whether or not its node is held; evicting the region that
// brought it leaves the bookmark in the spec, waiting (never dangling, never re-bound), and it resolves again, the
// same node, when the region is read again.

/** Positions of the region whose answers bring TP53's transcripts (T1, T2, T3), and a stretch past the gene. */
const Gene = { start: GeneOffset, end: GeneOffset + 1100 };
const Far = { start: 5000, end: 6000 };

let stores: DataStore[] = [];
afterEach(() => {
  for (const s of stores) s.retire();
  stores = [];
  vi.useRealTimers();
});

/** A store over the region's source with a cache that holds about one stretch of it. */
function store(): { data: DataStore; wanter: Wanter } {
  vi.useFakeTimers({ toFake: ['setTimeout', 'clearTimeout', 'Date', 'performance'] });
  const data = new DataStore({ settings: { maxBytes: 6000, tileBins: 256 } });
  stores.push(data);
  expect(data.attach(RegionSource, openSeq(regionData('seq')))).toBe(true);
  return { data, wanter: data.wanter() };
}

async function read(s: { data: DataStore; wanter: Wanter }, at: { start: number; end: number }) {
  s.wanter.want([{ source: RegionSource, system: RegionSystem, ...at, level: 1 }]);
  await vi.advanceTimersByTimeAsync(10);
  expect(s.data.settled.value).toBe(true);
}

function status(root: Node, data: DataStore, address: string) {
  const model = resolveSpec(root);
  return resolveAddressIn(root, overlayQueryModel(model, buildOverlay(model, data)), address).status;
}

describe('a bookmark on a streamed transcript', () => {
  it('is added as one patch while its node is held, and stays when its region is evicted (EC-162, EC-175)', async () => {
    const s = store();
    await read(s, Gene);
    expect(s.data.getNode(RegionSource, 'T2')).toBeDefined();
    const added = applyPatch(
      bookmarksSpec(),
      addRefPatch(bookmarksSpec(), { id: 'iso2', to: `${RegionSource}/(node :key T2)` }),
      {
        language: SeqLanguage,
      },
    );
    expect(added.ok).toBe(true);
    if (!added.ok) return;
    const root = added.root;
    expect(status(root, s.data, '/refs/iso2')).toBe('resolved');

    // the view moves past the gene: the gene's tiles leave the cache and their transcripts go with them
    await read(s, Far);
    expect(s.data.getNode(RegionSource, 'T2')).toBeUndefined();
    // the spec is what it was: arrivals and evictions never touch it, and the bookmark waits, never dangles
    expect(hashNode(root)).toBe(hashNode(added.root));
    expect(status(root, s.data, '/refs/iso2')).toBe('waiting');
    expect(checkSeqLang(root).filter((i) => i.code.startsWith('seq.lang.ref'))).toEqual([]);
  });

  it('is added while its node is not held, and later patches never refuse it as dangling', async () => {
    const s = store();
    await read(s, Far);
    expect(s.data.getNode(RegionSource, 'T2')).toBeUndefined();
    const base = bookmarksSpec();
    const added = applyPatch(base, addRefPatch(base, { id: 't3', to: `${RegionSource}/(node :key T3)` }), {
      language: SeqLanguage,
    });
    expect(added.ok).toBe(true);
    if (!added.ok) return;
    expect(status(added.root, s.data, '/refs/t3')).toBe('waiting');
    // removing another bookmark (`/refs/gene`) is checked for dangling refs: the waiting ones are not
    const removed = applyPatch(
      added.root,
      { base: hashNode(added.root), ops: [{ op: 'remove', at: ['refs', 'gene'] }] },
      { language: SeqLanguage },
    );
    expect(removed.ok ? [] : removed.issues).toEqual([]);
    // a rename of the source rewrites what names its nodes, held or not (refs follow edits)
    const renamed = applyPatch(
      added.root,
      {
        base: hashNode(added.root),
        ops: [{ op: 'move', at: ['data', 'region'], parent: ['data'], place: { after: null }, rename: 'tp53' }],
      },
      { language: SeqLanguage },
    );
    expect(renamed.ok).toBe(true);
    if (!renamed.ok) return;
    const t3 = renamed.root.children?.find((c) => c.id === 'refs')?.children?.find((c) => c.id === 't3');
    expect(t3?.params?.to).toBe('/data/tp53/(node :key T3)');
  });

  it('resolves again when the region is read again: the same node, under the same address (EC-176)', async () => {
    const s = store();
    const root = bookmarksSpec();
    await read(s, Gene);
    const first = s.data.getNode(RegionSource, 'T2');
    expect(status(root, s.data, '/refs/t2')).toBe('resolved');
    await read(s, Far);
    expect(status(root, s.data, '/refs/t2')).toBe('waiting');
    await read(s, Gene);
    expect(status(root, s.data, '/refs/t2')).toBe('resolved');
    expect(s.data.getNode(RegionSource, 'T2')?.node).toEqual(first?.node);
  });
});
