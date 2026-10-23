/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { observableState } from '@molstar/prism-ui/reactive/observable-state';
import type { VolumeRef } from '@molstar/plugin/state/manager/volume/hierarchy-state';
import { headlessViewer } from '../../_test/harness.js';
import type { ViewerApi } from '../../api/area.js';
import type { VolumeLevelUnit } from '../../model/volume-levels.js';
import { mapFile } from '../../api/_test/volume-fixtures.js';
import { transferColorAt, transferPeaks } from '../../model/volume-transfer.js';
import { createVolumeTransfer, themeDomain, type VolumeTransferSource } from '../volume-transfer.js';

// A direct volume's transfer function (model/volume-transfer.ts): the peaks of its curve, the colors painted along the
// axis as the renderer does (the theme's domain placed on it, clamped at its ends), and over a real map its points and
// the theme's stops edited through `volume.update`, Colors at peaks, one undo step a gesture.

const stats = { min: -2, max: 8, mean: 0, sigma: 1 };

describe('the transfer function on its axis', () => {
  it('finds the peaks: each run of points above zero, from its first to its last', () => {
    const peaks = transferPeaks([
      { x: 0.1, alpha: 0 },
      { x: 0.2, alpha: 0.05 },
      { x: 0.3, alpha: 0.05 },
      { x: 0.4, alpha: 0 },
      { x: 0.7, alpha: 0.02 },
      { x: 0.8, alpha: 0 },
    ]);
    expect(peaks).toEqual([
      { from: 0.2, to: 0.3 },
      { from: 0.7, to: 0.7 },
    ]);
  });

  it("places the theme's domain on the axis as the renderer does: auto, relative, symmetric, custom", () => {
    const base = { colorList: { kind: 'interpolate', colors: [] }, defaultColor: 0, isRelative: false } as any;
    expect(themeDomain({ ...base, domain: { name: 'auto', params: { symmetric: false } } }, stats)).toEqual([0, 1]);
    expect(themeDomain({ ...base, domain: { name: 'custom', params: [0, 4] } }, stats)).toEqual([0.2, 0.6]);
    expect(themeDomain({ ...base, domain: { name: 'auto', params: { symmetric: true } } }, stats)).toEqual([-0.6, 1]);
    // relative: the domain in sigma units, back to values (perf-audit-2 09/B7); a mean and sigma that are not 0 and 1
    const shifted = { min: -2, max: 8, mean: 1, sigma: 2 };
    const rel = { ...base, isRelative: true };
    expect(themeDomain({ ...rel, domain: { name: 'auto', params: { symmetric: false } } }, shifted)).toEqual([0, 1]);
    expect(themeDomain({ ...rel, domain: { name: 'custom', params: [0, 2] } }, shifted)).toEqual([0.3, 0.7]);
  });

  it('paints the stops interpolated, a set in steps, the ends clamped, within the domain', () => {
    const value = {
      kind: 'value' as const,
      interpolate: true,
      domain: [0.2, 0.6] as const,
      stops: [
        { x: 0.2, color: 0x000000 },
        { x: 0.6, color: 0xff0000 },
      ],
    };
    expect(transferColorAt(value, 0)).toBe(0x000000);
    expect(transferColorAt(value, 0.4)! >> 16).toBeGreaterThan(120);
    expect(transferColorAt(value, 0.9)).toBe(0xff0000);
    expect(transferColorAt({ ...value, interpolate: false }, 0.35)).toBe(0x000000);
    expect(transferColorAt({ ...value, interpolate: false }, 0.45)).toBe(0xff0000);
    expect(transferColorAt({ kind: 'uniform', color: 0x123456 }, 0.5)).toBe(0x123456);
    expect(transferColorAt({ kind: 'theme', label: 'Uniform' }, 0.5)).toBeUndefined();
  });
});

/**
 * A direct volume over a real map (headless has no WebGL to draw one): the loaded volume with a stand-in representation
 * whose params `volume.update` replaces, run by an API that records its groups.
 */
async function directVolume() {
  const v = await headlessViewer();
  await v.api.load.files({ files: [mapFile()] });
  const real = v.plugin.managers.volume.hierarchy.current.volumes[0];
  const data = real.cell.obj!.data;
  const provider = v.plugin.representation.volume.registry.get('direct-volume');
  const defaults = provider.getParams(v.plugin.representation.volume.themes, data) as any;
  const rep = {
    cell: {
      transform: {
        ref: 'dv',
        params: {
          type: { name: 'direct-volume', params: { controlPoints: defaults.controlPoints.defaultValue } },
          colorTheme: {
            name: 'volume-value',
            params: {
              colorList: {
                kind: 'interpolate',
                colors: [
                  [0xffffff, 0],
                  [0xff0000, 0.25],
                  [0xffffff, 0.5],
                  [0x0000ff, 0.75],
                  [0xffffff, 1],
                ],
              },
              domain: { name: 'auto', params: { symmetric: false } },
              isRelative: false,
              defaultColor: 0xcccccc,
            },
          },
        } as any,
      },
      state: {},
    },
  };
  const volume = { ...real, representations: [rep] } as unknown as VolumeRef;
  const groups: string[] = [];
  const updates: any[] = [];
  let source: VolumeTransferSource | undefined;
  const api = {
    group: async (label: string, work: (g: ViewerApi) => Promise<unknown>) => {
      groups.push(label);
      await work(api);
      return { ok: true };
    },
    run: async (id: string, values: any) => {
      expect(id).toBe('volume.update');
      updates.push(values.params);
      rep.cell.transform.params = { ...rep.cell.transform.params, ...values.params };
      source?.refresh();
      return { ok: true };
    },
  } as unknown as ViewerApi;
  source = createVolumeTransfer(v.plugin, api, () => volume, 'dv', observableState<VolumeLevelUnit>('sigma'), {
    quietMs: 10_000,
    throttleMs: 0,
  })!;
  return { tf: source, rep, groups, updates };
}

const flush = () => new Promise((r) => setTimeout(r, 0));

describe('the transfer function of a direct volume', () => {
  it('reads its points and colors, and edits them as one undo step a gesture', async () => {
    const t = await directVolume();
    const points = t.tf.points.value;
    expect(points.length).toBeGreaterThan(1);
    expect(t.tf.colors.value).toMatchObject({ kind: 'value', interpolate: true, domain: [0, 1] });
    expect(t.tf.presets().length).toBeGreaterThan(0);

    for (const a of [0.01, 0.02, 0.03])
      t.tf.change({ points: points.map((p, i) => (i === 1 ? { ...p, alpha: a } : p)) });
    // Shown at once, before it is applied.
    expect(t.tf.points.value[1].alpha).toBe(0.03);
    t.tf.settle();
    await flush();
    expect(t.groups).toEqual(['Transfer function']);
    expect(t.rep.cell.transform.params.type.params.controlPoints[1][1]).toBe(0.03);
  });

  it('emits colors only when they change, not when the points move or a refresh reads the same ones', async () => {
    const t = await directVolume();
    let colorEmits = 0;
    let pointEmits = 0;
    const subs = [t.tf.colors.subscribe(() => colorEmits++), t.tf.points.subscribe(() => pointEmits++)];
    const [c0, p0] = [colorEmits, pointEmits];
    t.tf.refresh();
    expect([colorEmits, pointEmits]).toEqual([c0, p0]);

    const points = t.tf.points.value;
    t.tf.change({ points: points.map((p, i) => (i === 1 ? { ...p, alpha: 0.02 } : p)) });
    expect(pointEmits).toBeGreaterThan(p0);
    expect(colorEmits).toBe(c0);
    t.tf.settle();
    await flush();
    expect(colorEmits).toBe(c0);
    for (const s of subs) s.unsubscribe();
  });

  it('puts colors at the peaks: a stop at each end of each, the theme’s own colors in order, offsets on the axis', async () => {
    const t = await directVolume();
    const peaks = transferPeaks(t.tf.points.value);
    t.tf.colorsAtPeaks();
    await flush();
    const list = t.rep.cell.transform.params.colorTheme.params.colorList.colors as [number, number][];
    expect(list.map((c) => c[1])).toEqual(
      peaks.flatMap((p) => (p.to - p.from > 0.002 ? [p.from, p.to] : [p.from])).map((x) => Math.round(x * 1000) / 1000),
    );
    expect(list[0][0]).toBe(0xff0000);
    expect(t.tf.colors.value.kind === 'value' && t.tf.colors.value.stops.map((s) => s.color)[0]).toBe(0xff0000);
  });

  it('colors the peaks from a palette, and a preset moves the stops onto its peaks with their colors', async () => {
    const t = await directVolume();
    const peakColors = () => {
      const list = t.rep.cell.transform.params.colorTheme.params.colorList.colors as [number, number][];
      return transferPeaks(t.tf.points.value).map((p) => list.find((c) => Math.abs(c[1] - p.from) < 0.002)?.[0]);
    };
    const twoPeaks = [
      { x: 0.1, alpha: 0 },
      { x: 0.2, alpha: 0.05 },
      { x: 0.3, alpha: 0.05 },
      { x: 0.4, alpha: 0 },
      { x: 0.6, alpha: 0 },
      { x: 0.7, alpha: 0.05 },
      { x: 0.8, alpha: 0 },
    ];
    t.tf.change({ points: twoPeaks });
    t.tf.settle();
    await flush();
    t.tf.colorsAtPeaks('viridis');
    await flush();
    // A sequential palette spread over the peaks, low to high.
    expect(peakColors()).toEqual([0x440154, 0xfde725]);

    // A preset (here: the peaks moved): the stops follow, with the colors they had.
    const moved = twoPeaks.map((p) => ({ ...p, x: p.x + 0.05 }));
    t.tf.applyPreset(moved);
    await flush();
    expect(t.rep.cell.transform.params.type.params.controlPoints[1][0]).toBeCloseTo(0.25, 3);
    expect(peakColors()).toEqual([0x440154, 0xfde725]);
  });

  it('on a uniform theme, sets the one color, and Colors at peaks turns it into the volume-value theme', async () => {
    const t = await directVolume();
    t.rep.cell.transform.params = {
      ...t.rep.cell.transform.params,
      colorTheme: { name: 'uniform', params: { value: 0x00ff00 } },
    };
    t.tf.refresh();
    expect(t.tf.colors.value).toEqual({ kind: 'uniform', color: 0x00ff00 });
    t.tf.change({ uniform: 0xcc3311 });
    t.tf.settle();
    await flush();
    expect(t.rep.cell.transform.params.colorTheme).toMatchObject({ name: 'uniform', params: { value: 0xcc3311 } });
    t.tf.colorsAtPeaks();
    await flush();
    expect(t.rep.cell.transform.params.colorTheme.name).toBe('volume-value');
    expect(t.rep.cell.transform.params.colorTheme.params.colorList.colors[0][0]).toBe(0xcc3311);
  });
});
