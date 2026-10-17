/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { CanvasMenuItem } from '@molstar/prism-ui';
import { ShellModel } from '@molstar/prism-ui';
import type { DocNode } from '@molstar/seq-core/format';
import { firstValueFrom } from 'rxjs';
import { describe, expect, it } from 'vitest';
import { sequenceCommands, sequenceMenu } from '../sequence/commands.js';
import type { SequenceViewer } from '../sequence/model.js';
import { cbs, crambin, expectOk, headlessViewer, until } from './harness.js';
import type { PrismViewer } from '../viewer.js';

// Hiding and showing the sequence's tracks (sequence.md §7–8): `sequence.track.hide` / `.show` as steps Undo takes
// back, the person's edits saved with the sequence's state (version 3) and written into every document built again
// (another structure and back, new colors, a snapshot applied), the commands and the Sequence menu's Tracks. Headless.

const SS = 'ss-track';

async function opened(files = [cbs()]) {
  const v = await headlessViewer({ files });
  const viewer = (await v.views.sequence!.load()).viewers.get();
  await until(() => viewer.tracks.value.length === 3);
  return { v, viewer };
}
const steps = (v: PrismViewer) => v.managers.history.steps.value.length;
/** A track node of the document shown, by ref. */
function nodeOf(viewer: SequenceViewer, ref: string): DocNode | undefined {
  const walk = (n: DocNode): DocNode | undefined =>
    n.ref === ref ? n : (n.children ?? []).reduce<DocNode | undefined>((f, c) => f ?? walk(c), undefined);
  const doc = viewer.seq.document.shown.value;
  return doc ? walk(doc.root) : undefined;
}
const hiddenInDoc = (viewer: SequenceViewer, ref: string) => nodeOf(viewer, ref)?.params?.hidden === true;
const hiddenRefs = (viewer: SequenceViewer) => viewer.tracks.value.filter((t) => t.hidden).map((t) => t.ref);

describe("the sequence's tracks: hide and show", () => {
  it('hides a track as a step (by ref or name), already so, and Undo shows it again', async () => {
    const { v, viewer } = await opened();
    const before = steps(v);
    const r = await v.api.sequence!.hideTrack!({ track: SS });
    expectOk(r);
    expect(r.did).toBe('The Secondary structure track is hidden.');
    expect(steps(v)).toBe(before + 1);
    expect(hiddenRefs(viewer)).toEqual([SS]);
    expect(hiddenInDoc(viewer, SS)).toBe(true);
    expect(viewer.edits.value).toEqual({ tracks: { [SS]: { hidden: true } } });
    expect(await v.api.sequence!.hideTrack!({ track: 'secondary structure' })).toMatchObject({
      ok: true,
      unchanged: true,
    });
    expect(steps(v)).toBe(before + 1);

    expectOk(await v.api.history.undo());
    await until(() => hiddenRefs(viewer).length === 0);
    expect(hiddenInDoc(viewer, SS)).toBe(false);
    expect(viewer.edits.value).toEqual({});
    expectOk(await v.api.history.redo());
    await until(() => hiddenRefs(viewer).length === 1);
  });

  it('shows one track or every hidden one, each one step; names a track it does not know', async () => {
    const { v, viewer } = await opened();
    expectOk(await v.api.sequence!.hideTrack!({ track: 'ruler' }));
    expectOk(await v.api.sequence!.hideTrack!({ track: 'Residues' }));
    expect(hiddenRefs(viewer)).toEqual(['ruler', 'letters']);
    expectOk(await v.api.sequence!.showTrack!({ track: 'ruler' }));
    expect(hiddenRefs(viewer)).toEqual(['letters']);
    expectOk(await v.api.sequence!.hideTrack!({ track: 'ruler' }));
    const before = steps(v);
    const all = await v.api.sequence!.showTrack!({});
    expectOk(all);
    expect(all.did).toBe('The 2 hidden tracks are shown.');
    expect(steps(v)).toBe(before + 1);
    expect(hiddenRefs(viewer)).toEqual([]);
    expect(viewer.edits.value).toEqual({});
    expect(await v.api.sequence!.showTrack!({})).toMatchObject({ ok: true, unchanged: true });
    const missing = await v.api.sequence!.hideTrack!({ track: 'plddt' });
    expect(missing).toMatchObject({ ok: false, code: 'not-found' });
    expect((missing as { hint?: string }).hint).toContain('Secondary structure (ss-track)');
  });

  it('keeps the edits across rebuilds: another structure and back, new letter colors', async () => {
    const { v, viewer } = await opened([cbs(), crambin()]);
    expectOk(await v.api.sequence!.hideTrack!({ track: SS }));
    expectOk(await v.api.sequence!.follow({ structure: '1crn' }));
    await until(() => viewer.structure.value?.title === '1CRN' && viewer.tracks.value.length === 3);
    expect(hiddenRefs(viewer)).toEqual([SS]);
    expect(hiddenInDoc(viewer, SS)).toBe(true);
    expectOk(await v.api.sequence!.follow({ structure: '1cbs' }));
    await until(() => viewer.structure.value?.title === '1CBS' && viewer.tracks.value.length === 3);
    expect(hiddenInDoc(viewer, SS)).toBe(true);
    expectOk(await v.api.sequence!.colors({ scheme: 'taylor' }));
    await until(
      () =>
        (nodeOf(viewer, 'letters')?.children?.[0]?.params?.color as { params?: { name?: string } })?.params?.name ===
        'taylor',
    );
    expect(hiddenInDoc(viewer, SS)).toBe(true);
    expect(hiddenRefs(viewer)).toEqual([SS]);
  });

  it('saves the edits as version 3 of the sequence state, and a snapshot brings them back; version 2 has none', async () => {
    const { v, viewer } = await opened();
    expectOk(await v.api.sequence!.hideTrack!({ track: SS }));
    const snapshot = v.plugin.state.getSnapshot({ camera: false });
    const saved = snapshot.custom?.['prism-viewer.sequence'] as { version: number; edits?: unknown };
    expect(saved).toMatchObject({ version: 3, edits: { tracks: { [SS]: { hidden: true } } } });

    // Another viewer opening the snapshot (a session, a reload): the track is hidden in the document it builds.
    const { v: w, viewer: other } = await opened();
    await w.plugin.state.setSnapshot(snapshot);
    await until(() => other.tracks.value.length === 3 && hiddenRefs(other).length === 1);
    expect(hiddenInDoc(other, SS)).toBe(true);

    // Version 2 (before the edits): read as version 3 with none.
    const { edits: _, ...rest } = saved as Record<string, unknown>;
    await v.plugin.state.setSnapshot({
      ...snapshot,
      custom: { ...snapshot.custom, 'prism-viewer.sequence': { ...rest, version: 2 } },
    });
    await until(() => hiddenRefs(viewer).length === 0);
    expect(viewer.edits.value).toEqual({});
  });

  it('the commands and the Sequence menu hide and show through the API; the Tracks rows follow the tracks', async () => {
    const { v, viewer } = await opened();
    const shell = new ShellModel();
    const commands = sequenceCommands(viewer, shell, v.api);
    const command = (id: string) => commands.find((c) => c.id === id)!;
    const signal = new AbortController().signal;
    expect(command('sequence.track.show').when?.({})).toBe('No sequence track is hidden.');
    const before = steps(v);
    await command('sequence.track.hide').run({ track: SS }, {}, signal);
    expect(steps(v)).toBe(before + 1);
    expect(hiddenRefs(viewer)).toEqual([SS]);
    expect(command('sequence.track.show').when?.({})).toBe(true);

    const menu = sequenceMenu(viewer, shell, () => {}, v.api);
    const tracks = menu.find((m) => m.id === 'track-list')!;
    const rows = async () => (await firstValueFrom(tracks.items!)) as readonly CanvasMenuItem[];
    const checked = async (r: CanvasMenuItem) => firstValueFrom(r.pressed!);
    const list = await rows();
    expect(list.map((r) => r.label)).toEqual(['Ruler', 'Residues', 'Secondary structure', 'Show all']);
    expect(await Promise.all(list.slice(0, 3).map(checked))).toEqual([true, true, false]);
    list[0].run();
    await until(() => hiddenRefs(viewer).length === 2);
    (await rows()).find((r) => r.id === 'show-all')!.run();
    await until(() => hiddenRefs(viewer).length === 0);
    expect((await rows()).map((r) => r.id)).not.toContain('show-all');
    await command('sequence.track.show-all').run({}, {}, signal);
    expect(steps(v)).toBe(before + 3);
  });
});
