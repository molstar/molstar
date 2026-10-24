/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import type { DocNode } from '@molstar/seq-core/format';
import { filter, firstValueFrom } from 'rxjs';
import { describe, expect, it, vi } from 'vitest';
import type { SequenceViewer } from '../model.js';
import { headlessViewers, hostStructure, ValuesTrack, valuesBehavior } from './fixtures.js';

// The person's edits over the sequence's document (sequence.md §7.1): each structure's, hidden tracks by ref, written
// into each document the model builds and into the tracks a behavior adds (hidden from its first patch), a ref not
// shown kept for later, and a loaded SeqViewSpec's own hidden track shown by an edit that says so.

async function until(check: () => boolean, ms = 3000): Promise<void> {
  const t0 = Date.now();
  while (!check()) {
    if (Date.now() - t0 > ms) throw new Error(`timed out after ${ms} ms`);
    await new Promise((r) => setTimeout(r, 5));
  }
}
const settled = async (h: ReturnType<typeof headlessViewers>) => {
  await firstValueFrom(h.viewer.structure.pipe(filter((s) => s !== undefined)));
  await h.viewer.seq.graph.update();
};
function nodeOf(viewer: SequenceViewer, ref: string): DocNode | undefined {
  const walk = (n: DocNode): DocNode | undefined =>
    n.ref === ref ? n : (n.children ?? []).reduce<DocNode | undefined>((f, c) => f ?? walk(c), undefined);
  const doc = viewer.seq.document.shown.value;
  return doc ? walk(doc.root) : undefined;
}
const hidden = (viewer: SequenceViewer) => viewer.tracks.value.filter((t) => t.hidden).map((t) => t.ref);

describe("the sequence's edits", () => {
  it("hides a behavior's track, which comes back hidden with its structure; another structure keeps its own", async () => {
    const patched: boolean[] = [];
    const h = headlessViewers({ behaviors: [valuesBehavior(patched)] });
    const first = { ...(await hostStructure()), ref: 'a' };
    h.structure.set(first);
    await settled(h);
    await until(() => h.viewer.tracks.value.some((t) => t.ref === ValuesTrack));
    // Named by its behavior's title (its document label is empty).
    expect(h.viewer.tracks.value.find((t) => t.ref === ValuesTrack)?.label).toBe('Values');
    expect(await h.viewer.setTrackHidden(ValuesTrack, true)).toBe(true);
    expect(await h.viewer.setTrackHidden(ValuesTrack, true)).toBe(false);
    expect(await h.viewer.setTrackHidden('nothing', true)).toBe(false);
    expect(hidden(h.viewer)).toEqual([ValuesTrack]);
    expect(h.viewer.edits.value).toEqual({ a: { tracks: { [ValuesTrack]: { hidden: true } } } });

    // Another structure: its own tracks, none hidden (the edits are each structure's, sequence.md §7.1).
    h.structure.set({ ...(await hostStructure('2ocj')), ref: 'b' });
    await settled(h);
    await until(() => patched.length === 2);
    expect(nodeOf(h.viewer, ValuesTrack)?.params?.hidden).toBeUndefined();
    expect(hidden(h.viewer)).toEqual([]);

    // Back: the behavior inserts its track again, hidden as it goes in.
    h.structure.set({ ...(await hostStructure()), ref: 'a' });
    await settled(h);
    await until(() => patched.length === 3);
    expect(nodeOf(h.viewer, ValuesTrack)?.params?.hidden).toBe(true);
    expect(hidden(h.viewer)).toEqual([ValuesTrack]);

    // An edit for a track this document lacks stays for the next that has it.
    await h.viewer.setEdits({ a: { tracks: { [ValuesTrack]: { hidden: true }, later: { hidden: true } } } });
    expect(h.viewer.edits.value.a?.tracks?.later).toEqual({ hidden: true });
    h.viewers.dispose();
  });

  it("an edit keeps what the document shows: a behavior's nodes are not added again", async () => {
    const patched: boolean[] = [];
    const h = headlessViewers({ behaviors: [valuesBehavior(patched)] });
    h.structure.set({ ...(await hostStructure()), ref: 'a' });
    await settled(h);
    await until(() => patched.length === 1);
    expect(await h.viewer.setTrackHidden('ruler', true)).toBe(true);
    expect(await h.viewer.setTrackIn('ruler', false)).toBe(true);
    await new Promise((r) => setTimeout(r, 50));
    expect(patched).toEqual([true]);
    h.viewers.dispose();
  });

  it('hides or shows several tracks in one patch', async () => {
    const h = headlessViewers();
    h.structure.set(await hostStructure());
    await settled(h);
    await until(() => h.viewer.tracks.value.length >= 2);
    const refs = h.viewer.tracks.value.map((t) => t.ref);
    const patch = vi.spyOn(h.viewer.seq, 'applyPatch');
    expect(await h.viewer.setTracksHidden([...refs, 'nothing'], true)).toBe(refs.length);
    expect(patch).toHaveBeenCalledTimes(1);
    expect(hidden(h.viewer)).toEqual(refs);
    expect(await h.viewer.setTracksHidden(refs, true)).toBe(0);
    expect(patch).toHaveBeenCalledTimes(1);
    expect(await h.viewer.setTracksHidden(refs, false)).toBe(refs.length);
    expect(patch).toHaveBeenCalledTimes(2);
    expect(hidden(h.viewer)).toEqual([]);
    expect(h.viewer.edits.value).toEqual({});
    h.viewers.dispose();
  });

  it("a loaded SeqViewSpec's own hidden track: shown by an edit, hidden again by dropping it", async () => {
    const h = headlessViewers();
    const host = { ...(await hostStructure()), ref: 's1' };
    const b = createSeqViewSpec({ title: 'Own hidden' });
    b.data.complex({ ref: 'c', provider: 'molstar:structure', source: { host_file: '1tup.bcif', format: 'bcif' } });
    const section = b.view({ ref: 'v', label: 'V' }).section({ ref: 's', style: 'compact' });
    section.track({ ref: 'shown', label: 'Shown' }).letters({});
    section.track({ ref: 'quiet', label: 'Quiet', hidden: true } as never).letters({});
    h.structure.set(host);
    await settled(h);
    await h.viewer.setOverride({ document: b.getDocument(), structure: 's1' });
    await until(() => h.viewer.structure.value?.section === 's' && h.viewer.tracks.value.length === 2);
    expect(hidden(h.viewer)).toEqual(['quiet']);
    expect(await h.viewer.setTrackHidden('quiet', false)).toBe(true);
    // A change from the document's own value is an edit: shown, kept as one.
    expect(h.viewer.edits.value).toEqual({ s1: { tracks: { quiet: { hidden: false } } } });
    expect(hidden(h.viewer)).toEqual([]);
    expect(await h.viewer.setTrackHidden('quiet', true)).toBe(true);
    expect(h.viewer.edits.value).toEqual({});
    expect(hidden(h.viewer)).toEqual(['quiet']);
    h.viewers.dispose();
  });
});
