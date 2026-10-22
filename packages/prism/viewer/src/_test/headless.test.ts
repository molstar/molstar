/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as fs from 'node:fs';
import * as path from 'node:path';
import { afterEach, describe, expect, it } from 'vitest';
import { DefaultPrismViewerSpec } from '../default-spec.js';
import { createPrismViewerHeadless } from '../headless.js';
import { particlesFile } from '../api/_test/volume-fixtures.js';
import type { PrismViewer } from '../viewer.js';

// The headless viewer (viewer-spec.md §6): the same API as on a page, with no DOM and no canvas.

const dataDir = path.resolve(import.meta.dirname, '../../../../../data/examples');
const crambin = () => new File([fs.readFileSync(path.join(dataDir, '1crn.cif'))], '1crn.cif');
const cbs = () => new File([fs.readFileSync(path.join(dataDir, '1cbs_updated.cif'))], '1cbs.cif');

let viewer: PrismViewer | undefined;
afterEach(async () => {
  await viewer?.dispose();
  viewer = undefined;
});

async function create() {
  viewer = await createPrismViewerHeadless({ spec: DefaultPrismViewerSpec() });
  return viewer;
}

describe('the headless viewer', () => {
  it('is ready, empty, and refuses canvas operations with a reason', async () => {
    const v = await create();
    expect(v.state.ready.value).toBe(true);
    expect(v.state.empty.value).toBe(true);
    const r = await v.api.canvas.background({ color: 'white' });
    expect(r).toMatchObject({ ok: false, code: 'unavailable' });
    if (!r.ok) expect(r.message).toMatch(/canvas/);
  });

  it('loads a file, describes it, and undoes the load', async () => {
    const v = await create();
    const loaded = await v.api.load.files({ files: [crambin()] as never });
    expect(loaded.ok).toBe(true);
    expect(v.state.empty.value).toBe(false);
    expect(v.api.structure.loaded.value).toHaveLength(1);

    const described = await v.api.structure.list();
    expect(described.ok && described.did).toMatch(/327 atoms/);

    expect(v.managers.history.steps.value).toHaveLength(1);
    const undone = await v.api.history.undo();
    expect(undone.ok).toBe(true);
    expect(v.state.empty.value).toBe(true);
    expect(v.api.structure.loaded.value).toHaveLength(0);
  });

  it('is working from the start of a load to its end, before the plugin is (the splash waits on it)', async () => {
    const v = await create();
    expect(v.state.working.value).toBe(false);
    // A particle file: the lazy particles feature loads first, before any plugin work (`busy` is the plugin's)
    const seen: boolean[] = [];
    const sub = v.state.working.subscribe((w) => seen.push(w));
    const loading = v.api.load.files({ files: [particlesFile()] as never });
    await Promise.resolve();
    expect(v.state.working.value).toBe(true);
    expect(v.state.busy.value).toBe(false);
    expect((await loading).ok).toBe(true);
    expect(v.state.working.value).toBe(false);
    sub.unsubscribe();
    // one stretch of work, no gap the splash would take for idle
    expect(seen.filter((w) => w).length).toBeGreaterThan(0);
    expect(seen.join(',')).not.toMatch(/true,false,true/);
  });

  it('checks values: unknown names and wrong types are refused with a hint', async () => {
    const v = await create();
    const unknown = await v.api.run('camera.reset', { speed: 3 });
    expect(unknown).toMatchObject({ ok: false, code: 'invalid' });
    if (!unknown.ok) expect(unknown.hint).toMatch(/animate/);
    const missing = await v.api.run('camera.spin');
    expect(missing).toMatchObject({ ok: false, code: 'not-found' });
  });

  it('runs a group as one step, all or nothing', async () => {
    const v = await create();
    const r = await v.api.group('Load two', async (api) => {
      await api.load.files({ files: [crambin()] as never });
      await api.load.files({ files: [cbs()] as never });
    });
    expect(r.ok).toBe(true);
    expect(v.api.structure.loaded.value).toHaveLength(2);
    expect(v.managers.history.steps.value.map((s) => s.label)).toEqual(['Load two']);

    const failed = await v.api.group('Load and fail', async (api) => {
      await api.session.clear();
      return api.load.files({ files: [] as never });
    });
    expect(failed.ok).toBe(false);
    // The clear inside the failed group was undone with it.
    expect(v.api.structure.loaded.value).toHaveLength(2);

    await v.api.history.undo();
    expect(v.state.empty.value).toBe(true);
  });

  it('a script run is one step; a change outside the API makes earlier steps unsafe', async () => {
    const v = await create();
    await v.run('Load crambin', async (plugin) => {
      const data = await plugin.builders.data.rawData({
        data: fs.readFileSync(path.join(dataDir, '1crn.cif'), 'utf8'),
      });
      const trajectory = await plugin.builders.structure.parseTrajectory(data, 'mmcif');
      await plugin.builders.structure.hierarchy.applyPreset(trajectory, 'default');
    });
    expect(v.managers.history.steps.value.map((s) => s.label)).toEqual(['Load crambin']);

    // Outside the API: straight on the plugin.
    await v.plugin.builders.data.rawData({ data: 'data_x', label: 'outside' });
    expect(v.managers.history.steps.value[0].unsafe).toMatch(/outside/);
    const r = await v.api.history.undo();
    expect(r.ok).toBe(false);
  });

  it('says what exists instead of doing it again, and clears', async () => {
    const v = await create();
    await v.api.load.files({ files: [crambin()] as never });
    const cleared = await v.api.session.clear();
    expect(cleared.ok).toBe(true);
    expect(v.state.empty.value).toBe(true);
    const again = await v.api.session.clear();
    expect(again).toMatchObject({ ok: false, code: 'blocked' });
    await v.api.history.undo();
    expect(v.api.structure.loaded.value).toHaveLength(1);
  });

  it('records runs with their origin', async () => {
    const v = await create();
    const seen: string[] = [];
    v.api.events.finished.subscribe((r) => seen.push(`${r.origin.kind}:${r.operation}:${r.result?.ok}`));
    await v.api.as({ kind: 'person' }).structure.list();
    await v.api.load.files({ files: [crambin()] as never });
    expect(seen).toEqual(['person:structure.list:true', 'script:load.files:true']);
  });

  it('gives sources to read only (api.md §6): marks are set through highlight.mark, formats read when asked', async () => {
    const v = await create();
    expect('set' in v.sources.marks).toBe(false);
    expect('setHover' in v.sources.marks).toBe(false);
    expect(v.sources.marks.marks.value).toEqual([]);
    expect(v.sources.formats.length).toBeGreaterThan(0);
    expect(v.sources.structures).toBe(v.api.structure.loaded);
  });
});
