/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Vec2 } from '@molstar/core/math/linear-algebra';
import { Color } from '@molstar/core/util/color';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { VolumeValueColorThemeParams } from '@molstar/graphics/theme/color/volume-value';
import type { Volume } from '@molstar/model/model/volume';
import type { PluginContext } from '@molstar/plugin/context';
import type { VolumeRef } from '@molstar/plugin/state/manager/volume/hierarchy-state';
import { type ObservableState, observableState } from '@molstar/prism-ui/core';
import type { ViewerApi } from '../api/area.js';
import type { VolumeLevelUnit } from '../model/volume-levels.js';
import {
  type TransferChange,
  type TransferColors,
  type TransferPoint,
  type TransferPreset,
  type TransferStop,
  paletteColors,
  type TransferPalette,
  type TransferPeak,
  transferPeaks,
  type VolumeTransferModel,
} from '../model/volume-transfer.js';
import { gestureSteps } from './gesture-steps.js';
import { histogramOf } from './volume-levels.js';

// A direct volume's transfer function over the API (model/volume-transfer.ts): its control points and its color theme
// through `volume.update`, a gesture one undo step (`gestureSteps`). The theme's stops are placed on the points' axis as
// the renderer places them (volume-value.ts: its domain normalized over the map's range; the shader clamps outside).

export interface VolumeTransferSource extends VolumeTransferModel {
  /** Reads the representation again (the state tree changed). */
  refresh(): void;
}

type DirectParams = {
  type: { name: string; params: { controlPoints?: Vec2[] } & Record<string, unknown> };
  colorTheme?: { name: string; params: Record<string, any> };
} & Record<string, unknown>;

type ValueThemeParams = PD.Values<typeof VolumeValueColorThemeParams>;

/** Colors for peaks when the theme has none to keep: distinct, readable on dark and light. */
const PeakColors = [0x3362b2, 0xd33115, 0x16a34a, 0xf59e0b, 0x9333ea, 0x0891b2, 0xdb2777];

/** The palettes Colors at peaks offers. */
const Palettes: readonly TransferPalette[] = [
  { id: 'distinct', label: 'Distinct', colors: PeakColors },
  { id: 'red-blue', label: 'Red and blue', colors: [0xd33115, 0x3362b2, 0xf59e0b, 0x16a34a] },
  { id: 'pastel', label: 'Pastel', colors: [0x8dd3c7, 0xfb8072, 0x80b1d3, 0xfdb462, 0xb3de69, 0xbebada] },
  {
    id: 'viridis',
    label: 'Viridis',
    colors: [0x440154, 0x414487, 0x2a788e, 0x22a884, 0x7ad151, 0xfde725],
    sequential: true,
  },
  { id: 'magma', label: 'Magma', colors: [0x3b0f70, 0x8c2981, 0xde4968, 0xfe9f6d, 0xfcfdbf], sequential: true },
  { id: 'ice', label: 'Ice', colors: [0xc6dbef, 0x9ecae1, 0x6baed6, 0x3182bd, 0x08519c], sequential: true },
];

/**
 * The theme's own colors for `n` peaks: once each and in order, leaving out whites, greys and blacks (the default's
 * background between its colors); peaks beyond them get distinct ones.
 */
function ownColors(now: TransferColors, n: number): number[] {
  const own = [
    ...new Set(now.kind === 'value' ? now.stops.map((s) => s.color) : now.kind === 'uniform' ? [now.color] : []),
  ].filter(isChromatic);
  const extra = PeakColors.filter((c) => !own.includes(c));
  return Array.from({ length: n }, (_, i) => own[i] ?? extra[(i - own.length) % extra.length]);
}

/** A stop at each end of each peak (one at a single point), in the peak's color. */
function stopsAtPeaks(peaks: readonly TransferPeak[], picks: readonly number[]): TransferStop[] {
  return peaks.flatMap((p, i) => {
    const color = picks[i];
    return p.to - p.from > 0.002
      ? [
          { x: p.from, color },
          { x: p.to, color },
        ]
      : [{ x: p.from, color }];
  });
}

/** Not white, grey or black: its channels differ by more than a little. */
function isChromatic(c: number): boolean {
  const r = (c >> 16) & 255;
  const g = (c >> 8) & 255;
  const b = c & 255;
  return Math.max(r, g, b) - Math.min(r, g, b) > 40;
}

/** Offsets as the palette reads them: three decimals (its texture's resolution follows their precision). */
const round = (x: number) => Math.round(x * 1000) / 1000;

/** Where the volume-value theme's [0, 1] lies on the axis (volume-value.ts, the direct-volume branch). */
export function themeDomain(params: ValueThemeParams, stats: Volume['grid']['stats']): [number, number] {
  const { min, max, mean, sigma } = stats;
  const d: [number, number] = params.domain.name === 'custom' ? [...params.domain.params] : [min, max];
  if (params.domain.name === 'auto' && params.isRelative) {
    d[0] = (d[0] - mean) / sigma;
    d[1] = (d[1] - mean) / sigma;
  }
  if (params.domain.name === 'auto' && params.domain.params.symmetric) {
    const m = Math.max(Math.abs(d[0]), Math.abs(d[1]));
    d[0] = -m;
    d[1] = m;
  }
  // with `isRelative` the domain is in sigma units: back to values, which the shader normalizes by min and max
  if (params.isRelative) {
    d[0] = d[0] * sigma + mean;
    d[1] = d[1] * sigma + mean;
  }
  const span = max - min || 1;
  return [(d[0] - min) / span, (d[1] - min) / span];
}

const sameList = <T>(a: readonly T[], b: readonly T[], same: (x: T, y: T) => boolean) =>
  a === b || (a.length === b.length && a.every((x, i) => same(x, b[i])));

/** Control points by value (each read makes new ones). */
export const samePoints = (a: readonly TransferPoint[], b: readonly TransferPoint[]) =>
  sameList(a, b, (p, q) => p.x === q.x && p.alpha === q.alpha);

/** Colors by value (each read makes new ones). */
export function sameColors(a: TransferColors, b: TransferColors): boolean {
  if (a === b) return true;
  switch (a.kind) {
    case 'uniform':
      return b.kind === 'uniform' && a.color === b.color;
    case 'theme':
      return b.kind === 'theme' && a.label === b.label;
    case 'value':
      return (
        b.kind === 'value' &&
        a.interpolate === b.interpolate &&
        a.domain[0] === b.domain[0] &&
        a.domain[1] === b.domain[1] &&
        sameList(a.stops, b.stops, (p, q) => p.x === q.x && p.color === q.color)
      );
  }
}

/** The transfer function of a direct-volume representation `rep` of `volume`. */
export function createVolumeTransfer(
  plugin: PluginContext,
  api: ViewerApi,
  volume: () => VolumeRef | undefined,
  rep: string,
  unit: ObservableState<VolumeLevelUnit>,
  options: { title?: string; quietMs?: number; throttleMs?: number } = {},
): VolumeTransferSource | undefined {
  const data = volume()?.cell.obj?.data;
  if (!data) return undefined;
  const stats = data.grid.stats;
  const cell = () => volume()?.representations.find((r) => r.cell.transform.ref === rep)?.cell;
  const params = () => cell()?.transform.params as DirectParams | undefined;
  const themeLabel = (name: string) => plugin.representation.volume.themes.colorThemeRegistry.get(name)?.label ?? name;

  const readPoints = (): TransferPoint[] =>
    (params()?.type.params.controlPoints ?? []).map((v) => ({ x: v[0], alpha: v[1] })).sort((a, b) => a.x - b.x);
  const readColors = (): TransferColors => {
    const theme = params()?.colorTheme;
    if (theme?.name === 'uniform') return { kind: 'uniform', color: theme.params.value ?? 0xcccccc };
    if (theme?.name !== 'volume-value') return { kind: 'theme', label: theme ? themeLabel(theme.name) : 'none' };
    const p = theme.params as ValueThemeParams;
    const domain = themeDomain(p, stats);
    const list = p.colorList.colors;
    const stops = list
      .map((c, i): TransferStop => {
        const offset = Array.isArray(c) ? c[1] : list.length > 1 ? i / (list.length - 1) : 0;
        return { x: domain[0] + offset * (domain[1] - domain[0]), color: Array.isArray(c) ? c[0] : (c as number) };
      })
      .sort((a, b) => a.x - b.x);
    return { kind: 'value', stops, interpolate: p.colorList.kind !== 'set', domain };
  };

  const histogram = observableState(histogramOf(data));
  // Compared by value: a drag of the points leaves the colors as they were (and the other way round).
  const points = observableState<readonly TransferPoint[]>([], { compare: samePoints });
  const colors = observableState<TransferColors>(readColors(), { compare: sameColors });
  const publish = () => {
    const o = steps.overlay('tf');
    points.set(o?.points ?? readPoints());
    const c = readColors();
    colors.set(
      o?.uniform !== undefined
        ? { kind: 'uniform', color: o.uniform }
        : o?.stops
          ? c.kind === 'value'
            ? { ...c, stops: o.stops }
            : // Becoming the volume-value theme over the map's range.
              { kind: 'value', stops: o.stops, interpolate: true, domain: [0, 1] }
          : c,
    );
  };

  /** The volume-value theme's params with `stops` (offsets along its domain, or a new theme over the map's range). */
  const valueTheme = (stops: readonly TransferStop[]) => {
    const theme = params()?.colorTheme;
    const current =
      theme?.name === 'volume-value'
        ? (theme.params as ValueThemeParams)
        : {
            ...PD.getDefaultValues(VolumeValueColorThemeParams),
            domain: { name: 'auto' as const, params: { symmetric: false } },
          };
    const [d0, d1] = themeDomain(current, stats);
    const offset = (x: number) => round(d1 > d0 ? Math.min(1, Math.max(0, (x - d0) / (d1 - d0))) : 0);
    return {
      name: 'volume-value',
      params: {
        ...current,
        colorList: {
          kind: current.colorList?.kind ?? 'interpolate',
          colors: [...stops].sort((a, b) => a.x - b.x).map((s) => [Color(s.color), offset(s.x)]),
        },
      },
    };
  };

  const steps = gestureSteps<TransferChange>(api, {
    label: options.title ? `${options.title} transfer function` : 'Transfer function',
    apply: async (g, _key, change) => {
      const p = params();
      if (!p) return;
      const type =
        change.points === undefined
          ? p.type
          : {
              ...p.type,
              params: {
                ...p.type.params,
                controlPoints: [...change.points]
                  .sort((a, b) => a.x - b.x)
                  .map((q) => [round(q.x), Math.round(q.alpha * 1e5) / 1e5] as unknown as Vec2),
              },
            };
      const colorTheme =
        change.uniform !== undefined
          ? {
              name: 'uniform',
              params: {
                ...(p.colorTheme?.name === 'uniform' ? p.colorTheme.params : {}),
                value: Color(change.uniform),
              },
            }
          : change.stops
            ? valueTheme(change.stops)
            : p.colorTheme;
      await g.run('volume.update', { node: rep, params: { type, ...(colorTheme ? { colorTheme } : {}) } });
    },
    changed: () => publish(),
    ...(options.quietMs !== undefined ? { quietMs: options.quietMs } : {}),
    ...(options.throttleMs !== undefined ? { throttleMs: options.throttleMs } : {}),
  });

  const provider = plugin.representation.volume.registry.get('direct-volume');
  const logMin =
    (provider?.getParams(plugin.representation.volume.themes, data) as { controlPoints?: PD.LineGraph } | undefined)
      ?.controlPoints?.yAxis?.logMin ?? 1e-3;

  publish();
  return {
    ...(options.title ? { title: options.title } : {}),
    stats,
    histogram: histogram.asReadonly,
    points: points.asReadonly,
    colors: colors.asReadonly,
    logMin,
    presets(): TransferPreset[] {
      const d = volume()?.cell.obj?.data;
      const def =
        d && (provider?.getParams(plugin.representation.volume.themes, d) as { controlPoints?: PD.LineGraph });
      return (def?.controlPoints?.getPresets?.() ?? []).map(([pts, label]) => ({
        label,
        points: pts.map((v) => ({ x: v[0], alpha: v[1] })),
      }));
    },
    unit: unit.asReadonly,
    setUnit: (u) => unit.set(u),
    change: (change) => steps.change('tf', change),
    settle: () => steps.settle(),
    palettes: Palettes,
    colorsAtPeaks(palette) {
      const peaks = transferPeaks(points.value);
      if (!peaks.length) return;
      const chosen = palette ? Palettes.find((p) => p.id === palette) : undefined;
      const picks = chosen ? paletteColors(chosen, peaks.length) : ownColors(colors.value, peaks.length);
      steps.change('tf', { stops: stopsAtPeaks(peaks, picks) });
      steps.settle();
    },
    applyPreset(next) {
      const now = colors.value;
      // The stops follow the curve's peaks on the volume-value theme; another theme's colors are not on the axis.
      const peaks = transferPeaks(next);
      const stops =
        now.kind === 'value' && peaks.length ? stopsAtPeaks(peaks, ownColors(now, peaks.length)) : undefined;
      steps.change('tf', { points: next, ...(stops ? { stops } : {}) });
      steps.settle();
    },
    refresh() {
      const d = volume()?.cell.obj?.data;
      if (d) histogram.set(histogramOf(d));
      publish();
    },
  };
}
