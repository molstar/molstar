/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it, vi } from 'vitest';
import { cbs, expectOk, headlessViewer, recordRuns, representationTypes, until } from '../../../_test/harness.js';
import type { PrismViewer } from '../../../viewer.js';
import { PluginStateSnapshotManager } from '@molstar/plugin/state/manager/snapshots';
import { createSnapshotStepper, createSnapshotsModel, createSnapshotTransition } from '../model.js';
import { SnapshotsFeature } from '../index.js';
import type { SnapshotList } from '../operations.js';
import { createSnapshotsSource } from '../source.js';

// The `snapshots` feature (viewer-spec.md §3.1): the `snapshot` area on the plugin's snapshot manager, headless (no
// canvas: the camera is not part of a snapshot here).

const viewer = (files = [cbs()]) => headlessViewer({ extensions: [SnapshotsFeature], files });
const snap = (v: PrismViewer) => v.api.snapshot!;
const steps = (v: PrismViewer) => v.managers.history.steps.value.length;
const list = async (v: PrismViewer) => {
  const r = await snap(v).list();
  expectOk(r);
  return r.value as SnapshotList;
};
const names = async (v: PrismViewer) => (await list(v)).entries.map((e) => e.name);

describe('snapshot.add and list', () => {
  it('adds entries from the current state without a history step', async () => {
    const v = await viewer();
    expect(await list(v)).toMatchObject({ entries: [], playing: false });
    const before = steps(v);
    const a = await snap(v).add({ name: 'Overview', description: 'The **fold**.', key: 'overview' });
    expectOk(a);
    expectOk(await snap(v).add());
    const { entries } = await list(v);
    expect(entries).toHaveLength(2);
    expect(entries[0]).toMatchObject({ position: 1, name: 'Overview', key: 'overview', description: 'The **fold**.' });
    expect(entries[0].id).toBe((a.value as { id: string }).id);
    expect(entries[1]).toMatchObject({ position: 2, name: 'Snapshot 2', current: true });
    expect(entries[0].current).toBe(false);
    // The list is not part of Undo (sequence.md §7.2).
    expect(steps(v)).toBe(before);
  });

  it('records the changes to the list as runs', async () => {
    const v = await viewer();
    const runs = recordRuns(v);
    await snap(v).add({ name: 'A' });
    await until(() => runs.some((r) => r.operation === 'snapshot.add'));
    expect(runs.find((r) => r.operation === 'snapshot.add')?.step).toBeUndefined();
  });

  it('keeps a key unique', async () => {
    const v = await viewer();
    expectOk(await snap(v).add({ name: 'A', key: 'k' }));
    expect(await snap(v).add({ name: 'B', key: 'k' })).toMatchObject({ ok: false, code: 'invalid' });
    expect(await names(v)).toEqual(['A']);
  });

  it('tells the list as lines for the agent', async () => {
    const v = await viewer();
    expect(await snap(v).list()).toMatchObject({ ok: true, did: 'There are no snapshots.' });
    await snap(v).add({ name: 'A' });
    const r = await snap(v).list();
    expect(r.ok && r.did).toMatch(/^1\. A \(current\) \[/);
  });
});

describe('snapshot.update', () => {
  it('renames, describes and keys an entry (by name), keeping what is empty', async () => {
    const v = await viewer();
    await snap(v).add({ name: 'A', description: 'first' });
    expectOk(await snap(v).update({ entry: 'a', name: 'Alpha', key: 'alpha' }));
    expect((await list(v)).entries[0]).toMatchObject({ name: 'Alpha', key: 'alpha', description: 'first' });
    // By key, and cleared on request.
    expectOk(await snap(v).update({ entry: 'alpha', name: 'Alpha', description: '', key: '' }));
    const e = (await list(v)).entries[0];
    expect(e.description).toBeUndefined();
    expect(e.key).toBeUndefined();
  });

  it('is already so when nothing would change, and the current one is the default', async () => {
    const v = await viewer();
    await snap(v).add({ name: 'A' });
    expect(await snap(v).update({ name: 'A' })).toMatchObject({ ok: true, unchanged: true });
    expect(await snap(v).update({ entry: 'A' })).toMatchObject({ ok: true, unchanged: true });
    expectOk(await snap(v).update({ description: 'now' }));
    expect((await list(v)).entries[0].description).toBe('now');
  });

  it('sets how long playback shows it', async () => {
    const v = await viewer();
    await snap(v).add({ name: 'A' });
    expectOk(await snap(v).update({ showFor: 3000 }));
    expect((await list(v)).entries[0].showFor).toBe(3000);
    expect(await snap(v).update({ showFor: 3000 })).toMatchObject({ unchanged: true });
  });

  it('keeps its time when replaced, unless a new one is asked for', async () => {
    const v = await viewer();
    expectOk(await snap(v).add({ name: 'A', showFor: 3000 }));
    // What the tab sends: its save options, with their own time.
    expectOk(await snap(v).replace({ options: { durationInMs: 1500 } }));
    expect((await list(v)).entries[0].showFor).toBe(3000);
    expectOk(await snap(v).replace({ showFor: 5000, options: { durationInMs: 1500 } }));
    expect((await list(v)).entries[0].showFor).toBe(5000);
  });

  it('refuses a key another entry has', async () => {
    const v = await viewer();
    await snap(v).add({ name: 'A', key: 'a' });
    await snap(v).add({ name: 'B' });
    expect(await snap(v).update({ entry: 'B', key: 'a' })).toMatchObject({ ok: false, code: 'invalid' });
  });
});

describe('finding an entry', () => {
  it('fails not-found by name, ambiguous by a shared name, invalid when none is given', async () => {
    const v = await viewer();
    expect(await snap(v).apply({ entry: 'x' })).toMatchObject({ ok: false, code: 'not-found' });
    await snap(v).add({ name: 'Same' });
    await snap(v).add({ name: 'same' });
    expect(await snap(v).apply({ entry: 'nope' })).toMatchObject({ ok: false, code: 'not-found' });
    const r = await snap(v).apply({ entry: 'SAME' });
    expect(r).toMatchObject({ ok: false, code: 'ambiguous' });
    expect(await snap(v).remove({})).toMatchObject({ ok: false, code: 'invalid' });
    // An id always works.
    const id = (await list(v)).entries[1].id;
    expectOk(await snap(v).apply({ entry: id }));
  });
});

describe('snapshot.apply', () => {
  it('makes the entry current only after it is applied: a failed apply leaves the current one', async () => {
    const v = await viewer();
    await snap(v).add({ name: 'A' });
    await snap(v).add({ name: 'B' });
    expectOk(await snap(v).apply({ entry: 'A' }));
    const current = async () => (await list(v)).entries.find((e) => e.current)?.name;
    expect(await current()).toBe('A');

    vi.spyOn(v.plugin.state, 'setSnapshot').mockRejectedValueOnce(new Error('cannot restore'));
    expect(await snap(v).apply({ entry: 'B' })).toMatchObject({ ok: false, code: 'failed' });
    expect(await current()).toBe('A');
    // Stepping is the same: the neighbor is not current unless it was applied.
    vi.spyOn(v.plugin.state, 'setSnapshot').mockRejectedValueOnce(new Error('cannot restore'));
    expect((await snap(v).next()).ok).toBe(false);
    expect(await current()).toBe('A');
    // And it works once the failure is over.
    expectOk(await snap(v).next());
    expect(await current()).toBe('B');
  });

  it('restores the scene it saved with no step: viewing a story, the steps before stay undoable', async () => {
    const v = await viewer();
    const runs = recordRuns(v);
    const original = representationTypes(v);
    expectOk(await snap(v).add({ name: 'Cartoon' }));
    expectOk(await snap(v).add({ name: 'Other' }));
    expectOk(await v.api.styles.apply({ style: 'spacefill' }));
    const changed = representationTypes(v);
    expect(changed).toContain('spacefill');
    const before = steps(v);
    const checkpoints = vi.spyOn(v.plugin.state, 'getSnapshot');

    const r = await snap(v).apply({ entry: 'Cartoon' });
    expectOk(r);
    expect(r.changed.view).toBe(true);
    expect(representationTypes(v)).toEqual(original);
    expect(steps(v)).toBe(before);
    expect(checkpoints).not.toHaveBeenCalled();
    expect(v.managers.history.steps.value.some((x) => !!x.unsafe)).toBe(false);
    // In the record all the same, with no step.
    await until(() => runs.some((x) => x.operation === 'snapshot.apply'));
    expect(runs.find((x) => x.operation === 'snapshot.apply')?.step).toBeUndefined();

    // A step on the snapshot shown, then another snapshot: Undo restores the scene and the snapshot it was shown with.
    expectOk(await v.api.styles.apply({ style: 'spacefill' }));
    expectOk(await snap(v).apply({ entry: 'Other' }));
    expectOk(await v.api.history.undo());
    expect(representationTypes(v)).toEqual(original);
    expect((await list(v)).entries.find((e) => e.current)?.name).toBe('Cartoon');
  });

  it('runs in order with the steps: Undo is refused while a snapshot is being applied', async () => {
    const v = await viewer();
    expectOk(await snap(v).add({ name: 'A' }));
    expectOk(await v.api.styles.apply({ style: 'spacefill' }));
    let open!: () => void;
    const gate = new Promise<void>((r) => (open = r));
    const real = v.plugin.state.setSnapshot.bind(v.plugin.state);
    vi.spyOn(v.plugin.state, 'setSnapshot').mockImplementationOnce(async (s) => {
      await gate;
      return real(s);
    });
    const applying = snap(v).apply({ entry: 'A' });
    await new Promise((r) => setTimeout(r, 0));
    expect(await v.api.history.undo()).toMatchObject({ ok: false });
    open();
    expectOk(await applying);
    expectOk(await v.api.history.undo());
  });

  it("a story's own links, clicks and keys go through the API: no step, the steps before stay undoable", async () => {
    const v = await viewer();
    const runs = recordRuns(v);
    const original = representationTypes(v);
    expectOk(await snap(v).add({ name: 'A', key: 'a' }));
    expectOk(await snap(v).add({ name: 'B', key: 'b' }));
    expectOk(await v.api.styles.apply({ style: 'spacefill' }));
    const before = steps(v);
    const current = async () => (await list(v)).entries.find((e) => e.current)?.name;

    // A description's link (the plugin's markdown extension calls the manager's applyKey).
    expectOk(await snap(v).link({ href: '!apply-snapshot=a', event: 'click' }));
    await until(() => runs.some((r) => r.operation === 'snapshot.apply' && r.origin.kind === 'person'));
    expect(representationTypes(v)).toEqual(original);
    expect(await current()).toBe('A');
    // A clickable primitive's key and the keyboard's next, as the plugin's behaviors call them.
    await v.plugin.managers.snapshot.applyNext(1);
    await until(() => runs.some((r) => r.operation === 'snapshot.next'));
    expect(await current()).toBe('B');
    await v.plugin.managers.snapshot.applyKey('a');
    expect(await current()).toBe('A');

    expect(steps(v)).toBe(before);
    expect(v.managers.history.steps.value.some((x) => !!x.unsafe)).toBe(false);
    expectOk(await v.api.history.undo());
  });

  it('Undo to a scene shown with no snapshot marks none current', async () => {
    const v = await viewer();
    expectOk(await v.api.styles.apply({ style: 'spacefill' }));
    expectOk(await snap(v).add({ name: 'B' }));
    expect((await list(v)).entries.find((e) => e.current)?.name).toBe('B');
    expectOk(await v.api.history.undo());
    expect((await list(v)).entries.some((e) => e.current)).toBe(false);
  });

  it('marks the entry current', async () => {
    const v = await viewer();
    await snap(v).add({ name: 'A' });
    await snap(v).add({ name: 'B' });
    expectOk(await snap(v).apply({ entry: 'A' }));
    expect((await list(v)).entries.map((e) => e.current)).toEqual([true, false]);
  });
});

describe('snapshot.replace', () => {
  it('takes the current state into an entry, keeping its name and key (the id changes)', async () => {
    const v = await viewer();
    const a = await snap(v).add({ name: 'A', key: 'a', description: 'd' });
    expectOk(a);
    expectOk(await v.api.styles.apply({ style: 'spacefill' }));
    const before = steps(v);
    const r = await snap(v).replace({ entry: 'A' });
    expectOk(r);
    expect(steps(v)).toBe(before);
    const e = (await list(v)).entries[0];
    expect(e).toMatchObject({ name: 'A', key: 'a', description: 'd' });
    expect(e.id).not.toBe((a.value as { id: string }).id);
    expect(e.id).toBe((r.value as { id: string }).id);

    // Applying it brings the surface back after the scene is changed again.
    expectOk(await v.api.styles.apply({ style: 'cartoon' }));
    expectOk(await snap(v).apply({ entry: 'a' }));
    expect(representationTypes(v)).toContain('spacefill');
  });

  it('needs an entry', async () => {
    const v = await viewer();
    expect(await snap(v).replace({})).toMatchObject({ ok: false, code: 'not-found' });
    await snap(v).add({ name: 'A' });
    expect(await snap(v).replace({ entry: 'B' })).toMatchObject({ ok: false, code: 'not-found' });
  });
});

describe('snapshot.move', () => {
  it('moves an entry to a position, without applying it', async () => {
    const v = await viewer();
    for (const n of ['A', 'B', 'C', 'D']) await snap(v).add({ name: n });
    expectOk(await snap(v).apply({ entry: 'B' }));
    expectOk(await snap(v).move({ entry: 'A', position: 3 }));
    expect(await names(v)).toEqual(['B', 'C', 'A', 'D']);
    expectOk(await snap(v).move({ entry: 'D', position: 1 }));
    expect(await names(v)).toEqual(['D', 'B', 'C', 'A']);
    // The one applied last stays current.
    expect((await list(v)).entries.find((e) => e.current)?.name).toBe('B');
  });

  it('is already so at the position, and refuses one outside the list', async () => {
    const v = await viewer();
    await snap(v).add({ name: 'A' });
    await snap(v).add({ name: 'B' });
    expect(await snap(v).move({ entry: 'A', position: 1 })).toMatchObject({ ok: true, unchanged: true });
    expect(await snap(v).move({ entry: 'A', position: 3 })).toMatchObject({ ok: false, code: 'invalid' });
    expect(await names(v)).toEqual(['A', 'B']);
  });
});

describe('snapshot.remove and clear', () => {
  it('removes an entry (a destructive operation for the agent), then all', async () => {
    const v = await viewer();
    expect(snap(v).remove.definition.approval).toBe('destructive');
    for (const n of ['A', 'B', 'C']) await snap(v).add({ name: n });
    const before = steps(v);
    expectOk(await snap(v).remove({ entry: 'B' }));
    expect(await names(v)).toEqual(['A', 'C']);
    expect(await snap(v).remove({ entry: 'B' })).toMatchObject({ ok: false, code: 'not-found' });
    expectOk(await snap(v).clear());
    expect(await names(v)).toEqual([]);
    expect(await snap(v).clear()).toMatchObject({ ok: true, unchanged: true });
    expect(steps(v)).toBe(before);
  });
});

describe('snapshot.next and previous', () => {
  it('step through the list, wrapping at the ends as the manager does', async () => {
    const v = await viewer();
    // Nothing to step through: refused before it runs.
    expect(await snap(v).next()).toMatchObject({ ok: false, code: 'blocked' });
    for (const n of ['A', 'B', 'C']) await snap(v).add({ name: n });
    const current = async () => (await list(v)).entries.find((e) => e.current)?.name;
    expect(await current()).toBe('C');
    const before = steps(v);
    const r = await snap(v).next();
    expectOk(r);
    expect(r.value).toMatchObject({ position: 1 });
    expect(await current()).toBe('A');
    expectOk(await snap(v).previous());
    expect(await current()).toBe('C');
    expectOk(await snap(v).previous());
    expect(await current()).toBe('B');
    // Viewing a story: no step.
    expect(steps(v)).toBe(before);
  });

  it('restores each snapshot; Undo takes back the step before, not the stepping', async () => {
    const v = await viewer();
    const original = representationTypes(v);
    await snap(v).add({ name: 'A' });
    await v.api.styles.apply({ style: 'spacefill' });
    await snap(v).add({ name: 'B' });
    const before = steps(v);
    expectOk(await snap(v).previous());
    expect(representationTypes(v)).toEqual(original);
    expectOk(await snap(v).next());
    expect(representationTypes(v)).toContain('spacefill');
    // The last step is the spacefill: Undo restores the scene before it, with A current as it was then.
    expect(steps(v)).toBe(before);
    expectOk(await v.api.history.undo());
    expect(representationTypes(v)).toEqual(original);
    expect((await list(v)).entries.find((e) => e.current)?.name).toBe('A');
  });
});

describe('snapshot.play and stop', () => {
  it('make no step, play needs two snapshots, and stop is already so when not playing', async () => {
    const v = await viewer();
    expect(await snap(v).play()).toMatchObject({ ok: false, code: 'not-found' });
    await snap(v).add({ name: 'A' });
    expect(await snap(v).play()).toMatchObject({ ok: false, code: 'blocked' });
    await snap(v).add({ name: 'B' });
    expect(await snap(v).stop()).toMatchObject({ ok: true, unchanged: true });

    const before = steps(v);
    expectOk(await snap(v).play({ delay: 100 }));
    expect((await list(v)).playing).toBe(true);
    expect(await snap(v).play()).toMatchObject({ ok: true, unchanged: true });
    expectOk(await snap(v).stop());
    expect((await list(v)).playing).toBe(false);
    expect(steps(v)).toBe(before);
  });

  it("restores the manager's delay when a play with its own delay stops", async () => {
    const v = await viewer();
    const manager = v.plugin.managers.snapshot;
    await snap(v).add({ name: 'A' });
    await snap(v).add({ name: 'B' });
    const was = manager.state.nextSnapshotDelayInMs;
    expectOk(await snap(v).play({ delay: 4000 }));
    expect(manager.state.nextSnapshotDelayInMs).toBe(4000);
    expectOk(await snap(v).stop());
    expect(manager.state.nextSnapshotDelayInMs).toBe(was);
  });
});

describe('playing', () => {
  it('is a presentation: no steps, and the steps before it stay undoable', async () => {
    const v = await viewer();
    await snap(v).add({ name: 'A' });
    await v.api.styles.apply({ style: 'spacefill' });
    await snap(v).add({ name: 'B' });
    const before = steps(v);
    expectOk(await snap(v).play());
    await new Promise((r) => setTimeout(r, 50));
    await snap(v).stop();
    expect(steps(v)).toBe(before);
    expect(v.managers.history.steps.value.some((x) => !!x.unsafe)).toBe(false);
  });
});

describe('the views', () => {
  it('follow the list through the state, and the models act through the API', async () => {
    const v = await viewer();
    const area = snap(v);
    const model = createSnapshotsModel(v.api, area, () => {});
    const stepper = createSnapshotStepper(v.api, area);
    expect(model.entries.value).toEqual([]);
    expect(stepper.count.value).toBe(0);

    model.add('One');
    await until(() => model.entries.value.length === 1);
    model.add();
    await until(() => model.entries.value.length === 2);
    expect(model.entries.value.map((e) => e.name)).toEqual(['One', 'Snapshot 2']);
    expect(stepper.count.value).toBe(2);
    expect(stepper.index.value).toBe(1);
    expect(stepper.title.value).toBe('Snapshot 2');
    expect(stepper.titles.value).toEqual(['One', 'Snapshot 2']);

    const first = model.entries.value[0].id;
    model.update(first, { description: 'Hello', key: 'one' });
    await until(() => model.entries.value[0].description === 'Hello');
    expect(model.entries.value[0].key).toBe('one');
    model.update(first, { description: '' });
    await until(() => model.entries.value[0].description === undefined);

    stepper.goTo(0);
    await until(() => stepper.index.value === 0);
    stepper.next();
    await until(() => stepper.index.value === 1);
    stepper.openLink('#one');
    await until(() => stepper.index.value === 0);

    model.moveTo(first, 1);
    await until(() => model.entries.value[1].id === first);
    model.move(first, 1);
    expect(model.entries.value[1].id).toBe(first);
    model.remove(first);
    await until(() => model.entries.value.length === 1);
    model.clear();
    await until(() => model.entries.value.length === 0);
    model.dispose();
  });

  it('save a session through session.save when it exists', async () => {
    const v = await viewer();
    const saved: string[] = [];
    const model = createSnapshotsModel(v.api, snap(v), (_, name) => saved.push(name));
    model.download('state');
    // Without the export-session feature in the spec there is nothing to save.
    await new Promise((r) => setTimeout(r, 20));
    expect(saved).toHaveLength(v.api.describe('session.save') ? 1 : 0);
    if (saved.length) expect(saved[0]).toMatch(/\.molj$/);
  });
});

describe('the Snapshots tab model while hidden', () => {
  /** Records the subscriptions to the list's state, to see when they are open. */
  const spyOnSubscribe = (area: ReturnType<typeof snap>) => {
    const opened: { closed: boolean }[] = [];
    const subscribe = area.state.subscribe.bind(area.state);
    vi.spyOn(area.state, 'subscribe').mockImplementation(((...args: unknown[]) => {
      const sub = (subscribe as (...a: unknown[]) => { closed: boolean })(...args);
      opened.push(sub);
      return sub;
    }) as never);
    return opened;
  };

  it('holds one subscription for all the entry forms, and none while hidden', async () => {
    const v = await viewer();
    const area = snap(v);
    for (const name of ['A', 'B', 'C']) await area.add({ name });
    const model = createSnapshotsModel(v.api, area, () => {});
    const opened = spyOnSubscribe(area);
    try {
      const ids = (await list(v)).entries.map((e) => e.id);
      for (const id of ids) expect(model.entryOptions(id)).toBeDefined();
      // One for the three forms, not one each.
      expect(opened).toHaveLength(1);
      expect(opened[0].closed).toBe(false);

      // Hidden: nothing follows the list. Shown again: one subscription, and the forms are current.
      model.setActive(false);
      expect(opened[0].closed).toBe(true);
      expectOk(await area.update({ entry: ids[0], showFor: 4000 }));
      expect(model.entryOptions(ids[0])?.values.value.durationInMs).not.toBe(4000);
      model.setActive(true);
      expect(opened).toHaveLength(2);
      expect(opened[1].closed).toBe(false);
      expect(model.entryOptions(ids[0])?.values.value.durationInMs).toBe(4000);
    } finally {
      model.dispose();
    }
    expect(opened.every((s) => s.closed)).toBe(true);
  });

  it('drops the forms of removed entries, and the subscription with the last one', async () => {
    const v = await viewer();
    const area = snap(v);
    await area.add({ name: 'A' });
    await area.add({ name: 'B' });
    const model = createSnapshotsModel(v.api, area, () => {});
    const opened = spyOnSubscribe(area);
    try {
      const [a, b] = (await list(v)).entries.map((e) => e.id);
      const formA = model.entryOptions(a);
      model.entryOptions(b);
      expectOk(await area.remove({ entry: a }));
      // The removed entry has no form any more; asking again gives none (it is gone), the other stays the same one.
      expect(model.entryOptions(a)).toBeUndefined();
      expect(model.entryOptions(b)).toBeDefined();
      expect(opened[0].closed).toBe(false);
      expectOk(await area.remove({ entry: b }));
      expect(opened[0].closed).toBe(true);
      expect(formA).toBeDefined();
      // A new entry gets a new subscription, not a leftover.
      await area.add({ name: 'C' });
      const c = (await list(v)).entries[0].id;
      model.entryOptions(c);
      expect(opened).toHaveLength(2);
    } finally {
      model.dispose();
    }
  });

  it('keeps a form with a pending edit while the list changes', async () => {
    const v = await viewer();
    const area = snap(v);
    await area.add({ name: 'A' });
    const model = createSnapshotsModel(v.api, area, () => {});
    try {
      const id = (await list(v)).entries[0].id;
      const form = model.entryOptions(id)!;
      form.set('durationInMs', 2500 as never);
      // Another edit of the list lands before the slider rests: it does not pull the slider back.
      await area.add({ name: 'B' });
      expect(form.values.value.durationInMs).toBe(2500);
      await vi.waitFor(async () => expect((await list(v)).entries[0].showFor).toBe(2500));
    } finally {
      model.dispose();
    }
  });

  it('does not report an edit for a hidden form, and makes the edit a form reports once it rests', async () => {
    const v = await viewer();
    const area = snap(v);
    await area.add({ name: 'A' });
    const model = createSnapshotsModel(v.api, area, () => {});
    try {
      const id = (await list(v)).entries[0].id;
      const form = model.entryOptions(id)!;
      form.set('durationInMs', 2500 as never);
      await vi.waitFor(async () => expect((await list(v)).entries[0].showFor).toBe(2500));
    } finally {
      model.dispose();
    }
  });
});

describe('the snapshots source', () => {
  it('does not read the list again for an event that changed neither entries, current nor playing', async () => {
    const v = await viewer();
    await snap(v).add({ name: 'A' });
    await snap(v).add({ name: 'B' });
    const manager = v.plugin.managers.snapshot;
    const source = createSnapshotsSource(v.plugin);
    try {
      const toArray = vi.spyOn(manager.state.entries, 'toArray');
      // The manager's events on each frame of a transition and step of the scrubber.
      manager.events.changed.next(void 0);
      manager.events.changed.next(void 0);
      expect(toArray).not.toHaveBeenCalled();
      await snap(v).add({ name: 'C' });
      expect(source.state.value.entries.map((e) => e.name)).toEqual(['A', 'B', 'C']);
    } finally {
      source.dispose();
    }
  });
});

describe("a snapshot's transition", () => {
  it('reads its length and where it is, and seeks through it without a history step', async () => {
    const v = await viewer();
    const manager = v.plugin.managers.snapshot;
    // Two frames of 200 ms over the loaded state (Mol*'s `snapshot.transition`).
    const base = v.plugin.state.getSnapshot({ camera: false, canvas3d: false });
    const frames = [0, 1].map(() => ({ durationInMs: 200, data: base.data! }));
    const entry = PluginStateSnapshotManager.Entry({ ...base, transition: { frames } }, { name: 'Moving' });
    manager.add(entry);
    expectOk(await snap(v).apply({ entry: entry.snapshot.id }));
    const model = createSnapshotTransition(v.api, snap(v));
    await until(() => model.duration.value === 400);
    expect(model.step.value).toBe(200);
    expect(model.isPlaying.value).toBe(false);
    const before = steps(v);
    model.seek(300);
    await until(() => model.time.value === 300);
    // Past its end: at its end.
    expectOk(await snap(v).seek({ time: 1000 }));
    expect(model.time.value).toBe(400);
    expect(steps(v)).toBe(before);
  });

  it('has none for a snapshot without one', async () => {
    const v = await viewer();
    expectOk(await snap(v).add());
    expect(createSnapshotTransition(v.api, snap(v)).duration.value).toBe(0);
    expect(await snap(v).seek({ time: 10 })).toMatchObject({ ok: false });
  });
});

describe('stepping through a story quietly', () => {
  it('posts no status message for a person stepping or applying, and says it to screen readers', async () => {
    const v = await viewer();
    expectOk(await snap(v).add({ name: 'One' }));
    expectOk(await snap(v).add({ name: 'Two' }));
    const person = v.api.as({ kind: 'person' });
    const posted: string[] = [];
    const sub = v.managers.notifications.posted.subscribe((m) => posted.push(m.text));
    expectOk(await person.snapshot!.previous());
    expectOk(await person.snapshot!.next());
    const entries = (await list(v)).entries;
    expectOk(await person.snapshot!.apply({ entry: entries[0].id }));
    sub.unsubscribe();
    expect(posted).toEqual([]);
    expect(v.managers.notifications.spoken.value).toMatch(/Applied snapshot 'One'/);
  });
});
