/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Color } from '@molstar/core/util/color';
import { stripTags } from '@molstar/core/util/string';
import type { PluginContext } from '@molstar/plugin/context';
import { StructureMeasurementParams } from '@molstar/plugin/state/manager/structure/measurement';
import { type IconName, observableState, type ReadonlyState } from '@molstar/prism-ui/core';
import { ParamsFormModel } from '@molstar/prism-ui/forms/params-form';
import { merge, skip } from 'rxjs';
import type { ViewerApi } from '../../api/area.js';
import type { MeasurementKind, MeasurementPoint, MeasurementsModel } from '../../model/measurements.js';
import type { SceneNode } from '../../model/scene.js';
import { sameItems } from '../../sources/state.js';
import { watchedState } from '../../sources/watched.js';
import { debounce } from '../../ui/debounce.js';

// The measurements as the views read them (state.md §1, `managers.structure.measurement` ↔ `MeasurementsModel`): the
// measurement manager's cells described in words, and the Measure section's model over them and over the selection's
// additions history (the ordered picks). Reads the plugin, as `sources/` does; every change of measurements goes
// through `measure.*`. The picks themselves are view state, not steps (state.md §3): ordering and dropping them is done
// here on the selection manager, and nothing records it.

export {
  findMeasurement,
  MeasureKinds,
  type MeasuredCell,
  type MeasurementInfo,
  type MeasureKind,
  readMeasurements,
} from '../../sources/measurements.js';
import {
  type MeasurementInfo,
  type MeasureKind,
  MeasureKinds,
  measurementLocis,
  readMeasurements,
} from '../../sources/measurements.js';

const Kinds: readonly MeasurementKind[] = (Object.keys(MeasureKinds) as MeasureKind[]).map((id) => ({
  id,
  ...MeasureKinds[id],
}));

/** The model, and what releases its subscriptions. */
export interface MeasurementsSource {
  readonly model: MeasurementsModel;
  dispose(): void;
}

class MeasurementNode implements SceneNode {
  readonly kind = 'measurement' as const;
  readonly #label = observableState('');
  readonly #description = observableState<string | undefined>(undefined);
  readonly #hidden = observableState(false);
  readonly #expanded = observableState(false);
  readonly #children = observableState<readonly SceneNode[]>([]);
  icon: IconName = 'ruler';
  removable = true;

  constructor(
    readonly id: string,
    private readonly actions: { visibility(hide: boolean): void; remove(): void; highlight(on: boolean): void },
  ) {}

  get label(): ReadonlyState<string> {
    return this.#label.asReadonly;
  }
  get description(): ReadonlyState<string | undefined> {
    return this.#description.asReadonly;
  }
  get isHidden(): ReadonlyState<boolean> {
    return this.#hidden.asReadonly;
  }
  get isExpanded(): ReadonlyState<boolean> {
    return this.#expanded.asReadonly;
  }
  get children(): ReadonlyState<readonly SceneNode[]> {
    return this.#children.asReadonly;
  }

  set(m: MeasurementInfo): void {
    this.icon = MeasureKinds[m.kind].icon;
    this.#label.set(m.label);
    this.#description.set(m.value || undefined);
    this.#hidden.set(m.hidden);
  }
  toggleVisibility(): void {
    this.actions.visibility(!this.#hidden.value);
  }
  toggleExpanded(): void {
    this.#expanded.set(!this.#expanded.value);
  }
  highlight(on: boolean): void {
    this.actions.highlight(on);
  }
  focus(): void {}
  remove(): void {
    this.actions.remove();
  }
}

/** The values of the options form: Mol*'s own (`StructureMeasurementParams`). */
const OptionsParams = {
  distanceUnitLabel: StructureMeasurementParams.distanceUnitLabel,
  textColor: StructureMeasurementParams.textColor,
};

/**
 * Mol*'s label of a pick lists what it is in, most specific first ("REA 200 | B [auth A] | Instance ASM-1 | Model 1 |
 * 1CBS"); a row shows the first and the last part, the pick and its structure.
 */
export function shortLabel(label: string): string {
  const parts = label
    .split('|')
    .map((p) => p.trim())
    .filter(Boolean);
  return parts.length > 1 ? `${parts[0]} · ${parts[parts.length - 1]}` : label;
}

const sameEntries = (a: { id: string; label: string }, b: { id: string; label: string }) =>
  a.id === b.id && a.label === b.label;

/**
 * The model of the Measure section over the plugin. `api` is bound to the person: adding, removing and showing or
 * hiding a measurement are its operations; the points are the selection's additions history (Mol*'s, newest first).
 *
 * Nothing here listens to the plugin until a view binds to it (prism AGENTS.md, hidden UI): `points`, `added` and the
 * options' values are watched states, `added` re-read once per animation frame however many events a transaction (or a
 * playing trajectory) fires.
 */
export function createMeasurementsSource(plugin: PluginContext, api: ViewerApi): MeasurementsSource {
  const selection = plugin.managers.structure.selection;
  const manager = plugin.managers.structure.measurement;

  const points = watchedState(
    (): readonly MeasurementPoint[] =>
      selection.additionsHistory.map((e) => {
        const label = stripTags(e.label).trim();
        return { id: e.id, label, short: shortLabel(label) };
      }),
    () => selection.events.additionsHistoryUpdated,
    { compare: (a, b) => sameItems(a, b, sameEntries), coalesce: 'none' },
  );
  const entry = (id: string) => selection.additionsHistory.find((e) => e.id === id);

  const nodes = new Map<string, MeasurementNode>();
  const added = watchedState(
    (): readonly SceneNode[] => {
      const list = readMeasurements(plugin);
      const live = new Set(list.map((m) => m.ref));
      for (const ref of [...nodes.keys()]) if (!live.has(ref)) nodes.delete(ref);
      return list.map((m) => {
        let node = nodes.get(m.ref);
        if (!node) {
          node = new MeasurementNode(m.ref, {
            visibility: (hide) => void api.measure?.visibility({ measurement: m.ref, visible: !hide }),
            remove: () => void api.measure?.remove({ measurement: m.ref }),
            highlight: (on) => {
              const loci = measurementLocis(plugin, m.ref)[0];
              if (on && loci) void api.highlight.set({ target: { loci } });
              else void api.highlight.clear();
            },
          });
          nodes.set(m.ref, node);
        }
        node.set(m);
        return node;
      });
    },
    // The manager's state replays on subscribe: the first read is the watcher's own.
    () =>
      merge(
        manager.behaviors.state.pipe(skip(1)),
        plugin.state.data.events.cell.stateUpdated,
        plugin.state.data.events.changed,
      ),
    { compare: (a, b) => sameItems(a, b) },
  );

  // The options form follows the manager (Undo, a session, a command) while it shows; an edit waits for its debounce
  // and its run, and the form keeps what the person is editing until they have landed.
  const optionValues = watchedState(
    () => manager.state.options,
    () => manager.behaviors.state.pipe(skip(1)),
    { compare: (a, b) => a.distanceUnitLabel === b.distanceUnitLabel && a.textColor === b.textColor, coalesce: 'none' },
  );
  let sending = 0;
  const send = debounce((values: PD.Values<typeof OptionsParams>) => {
    sending++;
    void Promise.resolve(
      api.measure?.options({
        distanceUnitLabel: values.distanceUnitLabel,
        textColor: Color.toHexStyle(values.textColor),
      }),
    ).finally(() => {
      sending--;
      if (!sending && !send.pending) options.load(manager.state.options as PD.Values<typeof OptionsParams>);
    });
  });
  const options = new ParamsFormModel(OptionsParams, manager.state.options, { onChange: (values) => send(values) });

  const model: MeasurementsModel = {
    points,
    movePoint: (id, by) => {
      const e = entry(id);
      if (e) selection.modifyHistory(e, by < 0 ? 'up' : 'down', 4);
    },
    movePointTo: (id, index) => {
      // Dragged to a position: the history moves the entry there, in place (Mol*'s list is the state).
      const e = entry(id);
      if (!e) return;
      const history = selection.additionsHistory;
      history.splice(history.indexOf(e), 1);
      history.splice(Math.max(0, Math.min(history.length, index)), 0, e);
      selection.events.additionsHistoryUpdated.next(undefined);
    },
    removePoint: (id) => {
      const e = entry(id);
      if (e) selection.modifyHistory(e, 'remove');
    },
    focusPoint: (id) => {
      const e = entry(id);
      if (e) void api.camera.frame({ target: { loci: e.loci } });
    },
    highlightPoint: (id) => {
      const e = id ? entry(id) : undefined;
      if (e) void api.highlight.set({ target: { loci: e.loci } });
      else void api.highlight.clear();
    },
    kinds: Kinds,
    add: (kind) => {
      const op = api.measure?.[kind];
      if (op) void op({});
    },
    options,
    optionValues,
    // An edit the person made is not undone by an older value arriving: it waits for its run to land.
    loadOptions: (values) => {
      if (!send.pending && !sending) options.load(values as PD.Values<typeof OptionsParams>);
    },
    added,
    styleOf: () => undefined,
  };
  return { model, dispose: () => send.dispose() };
}
