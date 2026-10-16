/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { calculateHistogram } from '@molstar/core/math/histogram';
import { Color } from '@molstar/core/util/color';
import type { VolumeStreaming } from '@molstar/plugin/behavior/dynamic/volume-streaming/behavior';
import type { PluginContext } from '@molstar/plugin/context';
import type { VolumeRef } from '@molstar/plugin/state/manager/volume/hierarchy-state';
import { type ObservableState, observableState } from '@molstar/prism-ui/core';
import { Grid, Volume } from '@molstar/model/model/volume';
import type { ViewerApi } from '../api/area.js';
import { type StreamingState, streamingChannels } from '../api/volume.js';
import { gestureSteps } from './gesture-steps.js';
import type {
  VolumeLevel,
  VolumeLevelChange,
  VolumeLevelHistogram,
  VolumeLevelStats,
  VolumeLevelsModel,
  VolumeLevelUnit,
} from '../model/volume-levels.js';

// A map's levels over the API (model/volume-levels.ts): a loaded map's isosurfaces through `volume.isovalue` and
// `volume.update`, streamed density's channels through `volume.streaming-view`. Changes in quick succession are one
// undo step (`gestureSteps`).

/** What a map's levels read and how a change is made. */
export interface VolumeLevelsPort {
  readonly title?: string;
  readonly stats: VolumeLevelStats;
  histogram(): VolumeLevelHistogram | undefined;
  read(): readonly VolumeLevel[];
  /** Applies one level's change, as part of the step in progress (`api` is the group's). */
  apply(api: ViewerApi, level: VolumeLevel, change: VolumeLevelChange): Promise<unknown>;
}

export interface VolumeLevelsSource extends VolumeLevelsModel {
  /** Reads the levels again (the state tree changed). */
  refresh(): void;
}

/** Levels by value (each read makes new ones). */
export const sameLevels = (a: readonly VolumeLevel[], b: readonly VolumeLevel[]) =>
  a === b ||
  (a.length === b.length &&
    a.every((l, i) => {
      const m = b[i];
      return (
        l.id === m.id &&
        l.label === m.label &&
        l.value === m.value &&
        l.color === m.color &&
        l.opacity === m.opacity &&
        l.hidden === m.hidden
      );
    }));

export function createVolumeLevels(
  api: ViewerApi,
  port: VolumeLevelsPort,
  unit: ObservableState<VolumeLevelUnit>,
  options: { quietMs?: number; throttleMs?: number } = {},
): VolumeLevelsSource {
  const histogram = observableState(port.histogram());
  const levels = observableState<readonly VolumeLevel[]>([], { compare: sameLevels });
  // What is shown: the tree's levels with the changes of steps still in progress over them (a marker being dragged
  // stays under the pointer while earlier positions land).
  const publish = () =>
    levels.set(
      port.read().map((l) => {
        const c = steps.overlay(l.id);
        return c ? { ...l, ...c } : l;
      }),
    );
  const steps = gestureSteps<VolumeLevelChange>(api, {
    label: port.title ? `${port.title} levels` : 'Volume levels',
    apply: async (g, id, change) => {
      const level = port.read().find((l) => l.id === id);
      if (level) await port.apply(g, level, change);
    },
    changed: () => publish(),
    ...(options.quietMs !== undefined ? { quietMs: options.quietMs } : {}),
    ...(options.throttleMs !== undefined ? { throttleMs: options.throttleMs } : {}),
  });

  publish();
  return {
    ...(port.title ? { title: port.title } : {}),
    stats: port.stats,
    histogram: histogram.asReadonly,
    levels: levels.asReadonly,
    unit: unit.asReadonly,
    setUnit: (u) => unit.set(u),
    change: (id, change) => steps.change(id, change),
    settle: () => steps.settle(),
    refresh() {
      histogram.set(port.histogram());
      publish();
    },
  };
}

/** Bins of a map's histogram: enough to see its shape in a narrow panel; Mol*'s cached 1024 bins, eight to one. */
const Bins = 128;
const folded = new WeakMap<object, VolumeLevelHistogram>();
const streamed = new WeakMap<object, VolumeLevelHistogram>();

/**
 * A loaded map's values over its range: Mol*'s histogram cached on the grid (`Grid.getHistogram`, the one its robust
 * statistics and the direct-volume transfer function use), folded once.
 */
export function histogramOf(volume: Volume): VolumeLevelHistogram {
  const fine = Grid.getHistogram(volume.grid, Grid.RobustStatsBinCount);
  let h = folded.get(fine);
  if (!h) {
    const step = Math.max(1, Math.round(fine.counts.length / Bins));
    const counts = new Int32Array(Math.ceil(fine.counts.length / step));
    for (let i = 0; i < fine.counts.length; i++) counts[Math.floor(i / step)] += fine.counts[i];
    h = { min: fine.min, max: fine.max, counts };
    folded.set(fine, h);
  }
  return h;
}

/**
 * Streamed density's values over the whole map's range (the server's, not the box's), once per box that arrives.
 */
function streamedHistogramOf(volume: Volume, min: number, max: number): VolumeLevelHistogram {
  const grid = volume.grid;
  let h = streamed.get(grid);
  if (!h) {
    const { counts } = calculateHistogram(grid.cells.data, Bins, { min, max });
    h = { min, max, counts };
    streamed.set(grid, h);
  }
  return h;
}

type IsoParams = {
  type: { name: string; params: { isoValue?: Volume.IsoValue; alpha?: number } & Record<string, unknown> };
  colorTheme?: { name: string; params: { value?: Color } & Record<string, unknown> };
} & Record<string, unknown>;

/** A loaded map: its isosurfaces (the representations with an iso value), or one of them. */
export function mapLevels(plugin: PluginContext, volume: () => VolumeRef | undefined, only?: string): VolumeLevelsPort {
  const data = volume()?.cell.obj?.data;
  const stats = data ? data.grid.stats : { min: 0, max: 0, mean: 0, sigma: 0 };
  const label = (name: string) => plugin.representation.volume.registry.get(name)?.label ?? name;
  return {
    stats,
    histogram: () => {
      const d = volume()?.cell.obj?.data;
      return d ? histogramOf(d) : undefined;
    },
    read: () =>
      (volume()?.representations ?? []).flatMap((r) => {
        const ref = r.cell.transform.ref;
        const p = r.cell.transform.params as IsoParams | undefined;
        const iso = p?.type.params.isoValue;
        if (!p || !iso || (only && ref !== only)) return [];
        return [
          {
            id: ref,
            label: label(p.type.name),
            value: Volume.IsoValue.toAbsolute(iso, stats).absoluteValue,
            color: p.colorTheme?.name === 'uniform' ? p.colorTheme.params.value : undefined,
            opacity: p.type.params.alpha ?? 1,
            hidden: !!r.cell.state.isHidden,
          },
        ];
      }),
    async apply(api, level, change) {
      const r = volume()?.representations.find((x) => x.cell.transform.ref === level.id);
      const p = r?.cell.transform.params as IsoParams | undefined;
      if (!r || !p) return;
      if (change.value !== undefined) {
        // In the unit the representation keeps (a relative level stays relative as the map changes).
        const relative = p.type.params.isoValue?.kind !== 'absolute';
        await api.run('volume.isovalue', {
          value: relative ? (change.value - stats.mean) / (stats.sigma || 1) : change.value,
          kind: relative ? 'relative' : 'absolute',
          volume: volume()!.cell.transform.ref,
          representation: level.id,
        });
      }
      if (change.color === undefined && change.opacity === undefined) return;
      // Read again: the iso value above changed the params.
      const now = r.cell.transform.params as IsoParams;
      const type =
        change.opacity === undefined
          ? now.type
          : { ...now.type, params: { ...now.type.params, alpha: change.opacity } };
      const colorTheme =
        change.color === undefined
          ? now.colorTheme
          : {
              name: 'uniform',
              params: {
                ...(now.colorTheme?.name === 'uniform' ? now.colorTheme.params : {}),
                value: Color(change.color),
              },
            };
      await api.run('volume.update', {
        representation: level.id,
        params: { type, ...(colorTheme ? { colorTheme } : {}) },
      });
    },
  };
}

const ChannelLabels: Readonly<Record<VolumeStreaming.ChannelType, string>> = {
  em: 'EM',
  '2fo-fc': '2Fo-Fc',
  'fo-fc(+ve)': 'Fo-Fc(+ve)',
  'fo-fc(-ve)': 'Fo-Fc(-ve)',
};

/**
 * Streamed density (v5's Volume Streaming channels): one map for EM, two for X-ray (2Fo-Fc; Fo-Fc with its two signs).
 * σ is the whole map's (the server's), the histogram the part streamed.
 */
export function streamingLevels(
  plugin: PluginContext,
  structure: string,
  state: () => StreamingState | undefined,
): VolumeLevelsPort[] {
  const on = state();
  if (!on) return [];
  const has = streamingChannels(on.entry);
  const maps: { title: string; channels: VolumeStreaming.ChannelType[]; info: number }[] =
    on.entry.kind === 'em'
      ? [{ title: 'EM', channels: ['em'], info: 0 }]
      : [
          { title: '2Fo-Fc', channels: ['2fo-fc'], info: 0 },
          { title: 'Fo-Fc', channels: ['fo-fc(+ve)', 'fo-fc(-ve)'], info: 1 },
        ];
  return maps.map(({ title, channels, info }) => {
    const v = on.entry.header.sampling[0]?.valuesInfo[info];
    const stats: VolumeLevelStats = v
      ? { min: v.min, max: v.max, mean: v.mean, sigma: v.sigma }
      : { min: -1, max: 1, mean: 0, sigma: 1 };
    const params = () =>
      state()?.params.entry.params.channels as Partial<Record<string, VolumeStreaming.ChannelParams>> | undefined;
    return {
      title,
      stats,
      histogram: () => {
        const s = state();
        const behavior =
          s && (plugin.state.data.cells.get(s.behavior)?.obj?.data as VolumeStreaming.Behavior | undefined);
        const data = behavior?.channels[channels[0]]?.data;
        return data && data.grid.cells.data.length ? streamedHistogramOf(data, stats.min, stats.max) : undefined;
      },
      read: () => {
        const p = params();
        const view = state()?.params.entry.params.view.name;
        return channels.flatMap((c) => {
          const ch = p?.[c];
          if (!ch || !has.includes(c)) return [];
          return [
            {
              id: c,
              label: ChannelLabels[c],
              value: Volume.IsoValue.toAbsolute(ch.isoValue, stats).absoluteValue,
              color: ch.color,
              opacity: ch.opacity,
              hidden: view === 'off',
            },
          ];
        });
      },
      apply: (api, level, change) =>
        api.run('volume.streaming-view', {
          structure,
          channels: {
            [level.id]: {
              ...(change.value !== undefined ? { isoValue: (change.value - stats.mean) / (stats.sigma || 1) } : {}),
              ...(change.color !== undefined ? { color: change.color } : {}),
              ...(change.opacity !== undefined ? { opacity: change.opacity } : {}),
            },
          },
        }),
    };
  });
}
