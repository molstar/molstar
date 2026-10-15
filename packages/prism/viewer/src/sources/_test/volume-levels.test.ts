/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { observableState } from '@molstar/prism-ui/reactive/observable-state';
import { describe, expect, it } from 'vitest';
import { expectOk, headlessViewer, until } from '../../_test/harness.js';
import { mapFile } from '../../api/_test/volume-fixtures.js';
import type { ViewerApi } from '../../api/area.js';
import type { VolumeLevel, VolumeLevelChange, VolumeLevelUnit } from '../../model/volume-levels.js';
import type { PrismViewer } from '../../viewer.js';
import { createVolumeLevels, histogramOf, type VolumeLevelsPort } from '../volume-levels.js';

// A map's levels (model/volume-levels.ts): changes in quick succession are one group (one undo step), applied one at a
// time with the latest winning, shown at once; over a headless viewer, a loaded map's isosurface through the API.

/** A port over plain levels whose changes wait until released, and an API whose groups are recorded. */
function fake() {
  let levels: VolumeLevel[] = [{ id: 'a', label: 'Isosurface', value: 1, color: 0x00ff00, opacity: 1, hidden: false }];
  const applied: VolumeLevelChange[] = [];
  const groups: string[] = [];
  let release: (() => void) | undefined;
  let running = 0;
  const port: VolumeLevelsPort = {
    stats: { min: -2, max: 10, mean: 0, sigma: 1 },
    histogram: () => ({ min: -2, max: 10, counts: [1, 2, 3] }),
    read: () => levels,
    async apply(_api, level, change) {
      running++;
      await new Promise<void>((r) => (release = r));
      running--;
      applied.push(change);
      levels = levels.map((l) => (l.id === level.id ? { ...l, ...change } : l));
    },
  };
  const api = {
    group: async (label: string, work: (g: ViewerApi) => Promise<unknown>) => {
      groups.push(label);
      await work(api);
      return { ok: true, did: label, runs: [] };
    },
  } as unknown as ViewerApi;
  const unit = observableState<VolumeLevelUnit>('sigma');
  const flush = async () => {
    for (let i = 0; i < 10; i++) await Promise.resolve();
  };
  const next = async () => {
    await flush();
    release?.();
    release = undefined;
    await flush();
  };
  return { port, api, unit, applied, groups, next, running: () => running };
}

describe('a map levels model', () => {
  it('applies a drag as one group: one change at a time, the latest winning, shown at once', async () => {
    const f = fake();
    const m = createVolumeLevels(f.api, f.port, f.unit, { quietMs: 10_000 });
    m.change('a', { value: 2 });
    m.change('a', { value: 3 });
    m.change('a', { value: 4 });
    // Where it is being dragged to, before any of it is applied.
    expect(m.levels.value[0].value).toBe(4);
    await f.next();
    m.change('a', { value: 5 });
    m.change('a', { color: 0xff0000 });
    expect(m.levels.value[0]).toMatchObject({ value: 5, color: 0xff0000 });
    await f.next();
    m.settle();
    await f.next();
    // What came while one was applied waited and was merged, the latest values winning.
    expect(f.applied).toEqual([{ value: 4 }, { value: 5, color: 0xff0000 }]);
    expect(f.groups).toEqual(['Volume levels']);
    expect(m.levels.value[0]).toMatchObject({ value: 5, color: 0xff0000 });
  });

  it('starts another step after settle, or after a moment without changes', async () => {
    const f = fake();
    const m = createVolumeLevels(f.api, f.port, f.unit, { quietMs: 5 });
    m.change('a', { value: 2 });
    m.settle();
    await f.next();
    m.change('a', { opacity: 0.5 });
    await f.next();
    await new Promise((r) => setTimeout(r, 20));
    await f.next();
    expect(f.groups).toHaveLength(2);
    expect(f.running()).toBe(0);
    expect(m.levels.value[0]).toMatchObject({ value: 2, opacity: 0.5 });
  });

  it('shares its unit', () => {
    const f = fake();
    const one = createVolumeLevels(f.api, f.port, f.unit);
    const two = createVolumeLevels(f.api, f.port, f.unit);
    one.setUnit('absolute');
    expect(two.unit.value).toBe('absolute');
  });
});

const isoOf = (v: PrismViewer) =>
  (v.plugin.managers.volume.hierarchy.current.volumes[0].representations[0].cell.transform.params as any).type.params
    .isoValue;

describe('the levels of a loaded map', () => {
  it('are on the volume and its isosurface, a drag one undo step, the color uniform', async () => {
    const v = await headlessViewer();
    await v.api.load.files({ files: [mapFile()] });
    await until(() => v.views.scene.roots.value.some((n) => n.kind === 'volume'));
    const volume = v.views.scene.roots.value.find((n) => n.kind === 'volume')!;
    const [levels] = v.views.scene.inspect(volume).levels!;
    const [level] = levels.levels.value;
    expect(level.label).toBe('Isosurface');
    // Mol*'s histogram cached on the grid (1024 bins), folded eight to one: every voxel counted once.
    const grid = v.plugin.managers.volume.hierarchy.current.volumes[0].cell.obj!.data.grid;
    const counts = Array.from(levels.histogram.value!.counts);
    expect(counts).toHaveLength(128);
    expect(counts.reduce((a, b) => a + b, 0)).toBe(grid.cells.data.length);
    // The same object each time it is read (folded once).
    expect(histogramOf(v.plugin.managers.volume.hierarchy.current.volumes[0].cell.obj!.data)).toBe(
      levels.histogram.value,
    );
    expect(Object.keys((grid as any)._histograms)).toEqual(['1024']);
    const { mean, sigma } = levels.stats;
    expect((level.value - mean) / sigma).toBeCloseTo(isoOf(v).relativeValue, 5);

    const steps = v.api.history.steps.value.length;
    for (const s of [1.1, 1.2, 1.3, 1.4]) levels.change(level.id, { value: mean + s * sigma });
    levels.settle();
    await until(() => Math.abs(isoOf(v).relativeValue - 1.4) < 1e-6);
    await until(() => v.api.running.value.length === 0);
    expect(isoOf(v).kind).toBe('relative');
    expect(v.api.history.steps.value.length).toBe(steps + 1);

    const rep = () => v.plugin.managers.volume.hierarchy.current.volumes[0].representations[0];
    const color = (rep().cell.transform.params as any).colorTheme;
    levels.change(level.id, { color: 0x3362b2, opacity: 0.4 });
    levels.settle();
    await until(() => (rep().cell.transform.params as any).colorTheme.params.value === 0x3362b2);
    expect((rep().cell.transform.params as any).colorTheme.name).toBe('uniform');
    expect((rep().cell.transform.params as any).type.params.alpha).toBe(0.4);
    await until(() => v.api.history.steps.value.length === steps + 2);
    // The isosurface's row shows its color; its inspector has the same level.
    const row = volume.children.value[0];
    await until(() => row.swatch?.value === 0x3362b2);
    expect(v.views.scene.inspect(row).levels![0].levels.value.map((l) => l.id)).toEqual([level.id]);

    expectOk(await v.api.history.undo());
    expect((rep().cell.transform.params as any).colorTheme).toEqual(color);
    expectOk(await v.api.history.undo());
    expect(isoOf(v).relativeValue).toBeCloseTo((level.value - mean) / sigma, 5);
  });
});
