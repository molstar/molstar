/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as fs from 'node:fs';
import * as path from 'node:path';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { cbs, expectOk, headlessViewer } from '../../../_test/harness.js';
import { ExportSessionFeature } from '../index.js';

// The `export-session` feature: `session.save` in the core's `session` area, headless (the download is the page's).

afterEach(() => {
  vi.unstubAllGlobals();
});

describe('export session', () => {
  it('adds session.save and its command, a read of the scene with no history', async () => {
    const v = await headlessViewer({ extensions: [ExportSessionFeature] });
    expect(v.api.session.save).toBeDefined();
    expect(v.api.describe('session.save')?.definition).toMatchObject({ history: false, undoable: false });
    expect(v.managers.commands.get('session.save')).toBeDefined();
    // Nothing loaded, nothing to save.
    expect(v.api.describe('session.save')?.available).toBe('Nothing is loaded.');
    expect(await v.api.session.save!()).toMatchObject({ ok: false, code: 'blocked' });
  });

  it('labels the kinds short enough for the form, the longer words in the description', async () => {
    const v = await headlessViewer({ extensions: [ExportSessionFeature] });
    const kind = v.api.describe('session.save')?.params.kind as unknown as {
      options: [string, string][];
      description?: string;
    };
    expect(kind.options).toEqual([
      ['session', 'Session (.molx), with data'],
      ['state', 'State (.molj)'],
    ]);
    expect(kind.description).toMatch(/keeps the loaded data/);
  });

  it('is only for built-in features', async () => {
    const { builtIn: _, ...outside } = ExportSessionFeature;
    await expect(headlessViewer({ extensions: [{ ...outside, id: 'outside' }] })).rejects.toThrow(/only the viewer/);
  });

  it('saves a session (.molx) and a state (.molj) as non-empty files, with no history step', async () => {
    const v = await headlessViewer({ extensions: [ExportSessionFeature], files: [cbs()] });
    const steps = v.managers.history.steps.value.length;
    const session = await v.api.session.save!({ kind: 'session' });
    expectOk(session);
    expect(session.value!.filename).toBe('1CBS.molx');
    expect(session.value!.blob.size).toBeGreaterThan(0);
    // A .molx is a zip.
    const head = new Uint8Array(await session.value!.blob.slice(0, 2).arrayBuffer());
    expect([...head]).toEqual([0x50, 0x4b]);

    const state = await v.api.session.save!({ kind: 'state', name: 'my scene.molj' });
    expectOk(state);
    expect(state.value!.filename).toBe('my scene.molj');
    expect(state.value!.blob.size).toBeGreaterThan(0);
    expect(JSON.parse(await state.value!.blob.text()).entries.length).toBeGreaterThan(0);
    expect(v.managers.history.steps.value.length).toBe(steps);
  });

  it('saves what the Save options say: no canvas settings or interactivity when they are off', async () => {
    const v = await headlessViewer({ extensions: [ExportSessionFeature], files: [cbs()] });
    const read = async (options?: Record<string, boolean>) => {
      const r = await v.api.session.save!({ kind: 'state', ...(options ? { options: options as never } : {}) });
      expectOk(r);
      return JSON.parse(await r.value!.blob.text()).entries[0].snapshot;
    };
    const all = await read();
    expect(all.data).toBeDefined();
    expect(all.interactivity).toBeDefined();
    const some = await read({ data: true, interactivity: false, canvas3d: false });
    expect(some.data).toBeDefined();
    expect(some.interactivity).toBeUndefined();
    expect(some.canvas3d).toBeUndefined();
  });

  it('opens a saved session: the structures come back, with their data', async () => {
    const v = await headlessViewer({ extensions: [ExportSessionFeature], files: [cbs()] });
    const count = () => v.plugin.managers.structure.hierarchy.current.structures.length;
    expect(count()).toBe(1);
    const saved = await v.api.session.save!({ kind: 'session' });
    expectOk(saved);
    expectOk(await v.api.session.clear());
    expect(count()).toBe(0);
    expectOk(await v.api.session.open({ file: new File([saved.value!.blob], saved.value!.filename) }));
    expect(count()).toBe(1);
  });

  it('opens a saved state when its data can be fetched again (it says where, not what)', async () => {
    const text = fs.readFileSync(
      path.resolve(import.meta.dirname, '../../../../../../../data/examples/1cbs_updated.cif'),
    );
    vi.stubGlobal('fetch', async () => new Response(text));
    const v = await headlessViewer({ extensions: [ExportSessionFeature] });
    const { builders } = v.plugin;
    const data = await builders.data.download({ url: 'https://data.test/1cbs.cif', isBinary: false });
    const trajectory = await builders.structure.parseTrajectory(data, 'mmcif');
    await builders.structure.hierarchy.applyPreset(trajectory, 'default');
    const count = () => v.plugin.managers.structure.hierarchy.current.structures.length;
    expect(count()).toBe(1);
    const saved = await v.api.session.save!({ kind: 'state' });
    expectOk(saved);
    expectOk(await v.api.session.clear());
    expect(count()).toBe(0);
    expectOk(await v.api.session.open({ file: new File([saved.value!.blob], saved.value!.filename) }));
    expect(count()).toBe(1);
  });

  it('opens a state whose one node cannot be made, with a warning; an error logged meanwhile is not a failure', async () => {
    const text = fs.readFileSync(
      path.resolve(import.meta.dirname, '../../../../../../../data/examples/1cbs_updated.cif'),
    );
    let reachable = true;
    const v = await headlessViewer({ extensions: [ExportSessionFeature] });
    vi.stubGlobal('fetch', async (url: string) => {
      if (url.endsWith('gone.cif') && !reachable) {
        v.plugin.log.error('Something else went wrong');
        return new Response('', { status: 404 });
      }
      return new Response(text);
    });
    const { builders } = v.plugin;
    for (const url of ['https://data.test/1cbs.cif', 'https://data.test/gone.cif']) {
      const data = await builders.data.download({ url, isBinary: false });
      const trajectory = await builders.structure.parseTrajectory(data, 'mmcif');
      await builders.structure.hierarchy.applyPreset(trajectory, 'default');
    }
    const saved = await v.api.session.save!({ kind: 'state' });
    expectOk(saved);
    expectOk(await v.api.session.clear());
    reachable = false;
    const opened = await v.api.session.open({ file: new File([saved.value!.blob], saved.value!.filename) });
    expectOk(opened);
    expect(opened.warnings).toEqual([expect.stringMatching(/^1 part of it could not be made: /)]);
    expect(v.plugin.managers.structure.hierarchy.current.structures.filter((s) => s.cell.obj)).toHaveLength(1);
  });

  it("keeps the session's record in the file only when asked, and shows it read-only when the file opens", async () => {
    const v = await headlessViewer({ extensions: [ExportSessionFeature], files: [cbs()] });
    expectOk(await v.api.as({ kind: 'person' }).styles.apply({ style: 'surface' }));
    await v.record.flush();

    // Off by default: no record in the file, none shown.
    const plain = await v.api.session.save!({ kind: 'session' });
    expectOk(plain);
    expectOk(await v.api.session.open({ file: new File([plain.value!.blob], 'plain.molx') }));
    expect(v.keptRecords.opened.value).toBeUndefined();

    const saved = await v.api.session.save!({ kind: 'session', name: 'with record', record: true });
    expectOk(saved);
    expect(saved.did).toMatch(/with its record/);
    const opened = await v.api.session.open({ file: new File([saved.value!.blob], saved.value!.filename) });
    expectOk(opened);
    expect(opened.did).toMatch(/record is in History/);
    const shown = v.keptRecords.opened.value!;
    expect(shown).toMatchObject({ fromFile: true, record: { title: 'with record.molx' } });
    const runs = (await shown.reader.read(50)).items.map((i) => i.entry).filter((e) => e.kind === 'run');
    expect(runs.map((r) => (r.kind === 'run' ? r.operation : ''))).toContain('styles.apply');

    // A state carries it as a field: its entries in their envelopes.
    const state = await v.api.session.save!({ kind: 'state', record: true });
    expectOk(state);
    const json = JSON.parse(await state.value!.blob.text());
    expect(json.entries.length).toBeGreaterThan(0);
    expect(json.record.length).toBeGreaterThan(0);
    expect(json.record[0]).toMatchObject({ codec: 'prism-viewer.record.v1', value: { kind: expect.any(String) } });
  });

  it('saves a state of locally opened files, which cannot be opened without them', async () => {
    const v = await headlessViewer({ extensions: [ExportSessionFeature], files: [cbs()] });
    const saved = await v.api.session.save!({ kind: 'state' });
    expectOk(saved);
    expectOk(await v.api.session.clear());
    expect(await v.api.session.open({ file: new File([saved.value!.blob], saved.value!.filename) })).toMatchObject({
      ok: false,
      code: 'failed',
    });
  });
});

describe('the Snapshots list (gate G5)', () => {
  it('is left as it is by a save: none added to an empty list, and the person’s own not replaced', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const entries = () => v.plugin.managers.snapshot.state.entries.toArray().map((e) => e.snapshot.id);
    expectOk(await v.api.session.save!({ kind: 'state' }));
    expect(entries()).toEqual([]);
    v.plugin.state.setSnapshotParams({ camera: false });
    expectOk(await v.api.run('snapshot.add', { name: 'Mine' }));
    const mine = entries();
    expectOk(await v.api.styles.apply({ style: 'surface' }));
    expectOk(await v.api.session.save!({ kind: 'state' }));
    expect(entries()).toEqual(mine);
  });
});
