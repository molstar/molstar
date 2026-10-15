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
import { Volume } from '@molstar/model/model/volume';
import type { ViewerApi } from '../api/area.js';
import { type StreamingState, streamingChannels } from '../api/volume.js';
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
// `api.group` (one undo step), the latest winning while one is applied.

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

/** A step in progress: its changes not applied yet, the ones shown until it ends, and its loop's wake-up. */
interface Step {
  readonly pending: Map<string, VolumeLevelChange>;
  readonly shown: Map<string, VolumeLevelChange>;
  wake?: () => void;
  closed: boolean;
}

/** After this long without a change (and no gesture held), the changes so far are one step. */
const QuietMs = 600;

export function createVolumeLevels(
  api: ViewerApi,
  port: VolumeLevelsPort,
  unit: ObservableState<VolumeLevelUnit>,
  options: { quietMs?: number } = {},
): VolumeLevelsSource {
  const quietMs = options.quietMs ?? QuietMs;
  const steps: Step[] = [];
  let open: Step | undefined;
  let quiet: ReturnType<typeof setTimeout> | undefined;
  const histogram = observableState(port.histogram());
  const levels = observableState<readonly VolumeLevel[]>([]);

  // What is shown: the tree's levels with the changes of steps still in progress over them (a marker being dragged
  // stays under the pointer while earlier positions land).
  const publish = () => {
    let list = port.read();
    for (const step of steps) {
      if (!step.shown.size) continue;
      list = list.map((l) => {
        const c = step.shown.get(l.id);
        return c ? { ...l, ...c } : l;
      });
    }
    levels.set(list);
  };

  const run = (step: Step) =>
    api
      .group(
        port.title ? `${port.title} levels` : 'Volume levels',
        async (g) => {
          for (;;) {
            const next = step.pending.entries().next();
            if (!next.done) {
              const [id, change] = next.value;
              step.pending.delete(id);
              const level = port.read().find((l) => l.id === id);
              if (level) await port.apply(g, level, change);
              continue;
            }
            if (step.closed) return;
            await new Promise<void>((resolve) => (step.wake = resolve));
            step.wake = undefined;
          }
        },
        { onFailure: 'continue' },
      )
      .catch(() => undefined)
      .finally(() => {
        steps.splice(steps.indexOf(step), 1);
        publish();
      });

  const close = () => {
    clearTimeout(quiet);
    const step = open;
    open = undefined;
    if (!step) return;
    step.closed = true;
    step.wake?.();
  };

  publish();
  return {
    ...(port.title ? { title: port.title } : {}),
    stats: port.stats,
    histogram: histogram.asReadonly,
    levels: levels.asReadonly,
    unit: unit.asReadonly,
    setUnit: (u) => unit.set(u),
    change(id, change) {
      let step = open;
      if (!step) {
        step = open = { pending: new Map(), shown: new Map(), closed: false };
        steps.push(step);
        void run(step);
      }
      step.pending.set(id, { ...step.pending.get(id), ...change });
      step.shown.set(id, { ...step.shown.get(id), ...change });
      step.wake?.();
      publish();
      clearTimeout(quiet);
      quiet = setTimeout(close, quietMs);
    },
    settle: close,
    refresh() {
      histogram.set(port.histogram());
      publish();
    },
  };
}

/** Bins of a map's histogram: enough to see its shape in a narrow panel. */
const Bins = 96;
const histograms = new WeakMap<object, VolumeLevelHistogram>();

/** A grid's values over [min, max], once per grid (a map's data does not change). */
export function histogramOf(volume: Volume, min: number, max: number): VolumeLevelHistogram {
  const grid = volume.grid;
  let h = histograms.get(grid);
  if (!h || h.min !== min || h.max !== max) {
    const { counts } = calculateHistogram(grid.cells.data, Bins, { min, max });
    h = { min, max, counts };
    histograms.set(grid, h);
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
      return d ? histogramOf(d, d.grid.stats.min, d.grid.stats.max) : undefined;
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
        return data && data.grid.cells.data.length ? histogramOf(data, stats.min, stats.max) : undefined;
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
