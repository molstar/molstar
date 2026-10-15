/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StateSelection } from '@molstar/core/state';
import type { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { PluginContext } from '@molstar/plugin/context';
import { WiggleStructureRepresentation3DFromBundle } from '@molstar/plugin/state/transforms/structure/effects/wiggle';
import type { ReadonlyState } from '@molstar/prism-ui/core';
import { merge } from 'rxjs';
import { sameItems } from '../../sources/state.js';
import { watchedState } from '../../sources/watched.js';

// The animations as the views read them (state.md §4, row 2): the plugin's animation manager (the animations
// registered, which plays) and the structures' procedural animation (the wiggle). Sources by role (contexts.md §7):
// this one may read the plugin, so the feature's views never do. The animation area owns it (`api.animation.state`);
// it listens to the plugin only while something watches it (prism AGENTS.md, hidden UI).

/** A registered animation. */
export interface AnimationInfo {
  /** Its registered name (`built-in.animate-camera-spin`). */
  readonly name: string;
  readonly title: string;
  readonly description?: string;
  /** True, or why it cannot play now (a trajectory to animate, snapshots to play). */
  readonly available: true | string;
}

/** Wiggle: Mol*'s procedural animation of a structure's atoms. */
export type WiggleMode = 'dynamics' | 'uncertainty' | 'off';

export interface WiggleState {
  readonly mode: WiggleMode;
  /** The amplitude of the dynamics (0 when off or by uncertainty). */
  readonly amplitude: number;
}

export interface AnimationState {
  readonly animations: readonly AnimationInfo[];
  /** The animation playing now, by name. */
  readonly playing?: string;
  readonly wiggle: WiggleState;
}

export interface AnimationSource {
  readonly state: ReadonlyState<AnimationState>;
  dispose(): void;
}

export function animationInfos(plugin: PluginContext): AnimationInfo[] {
  return plugin.managers.animation.animations.map((a) => {
    const can = a.canApply?.(plugin);
    return {
      name: a.name,
      title: a.display.name,
      ...(a.display.description ? { description: a.display.description } : {}),
      available: !can || can.canApply ? true : (can.reason ?? 'It cannot play now.'),
    };
  });
}

const words = (key: string) => key.replace(/([a-z0-9])([A-Z])/g, '$1 $2').replace(/^./, (c) => c.toUpperCase());

/** "Duration In Ms" is "Duration (ms)": the unit goes in brackets, in lower case. */
function unitLabel(key: string): string | undefined {
  const m = /^(.*?)In(Ms|S)$/.exec(key);
  return m ? `${words(m[1] || key)} (${m[2] === 'Ms' ? 'ms' : 's'})` : undefined;
}

/** The params with their duration labels in words; the rest as they are (a param with a label keeps it). */
export function relabel(params: PD.Params): PD.Params {
  const out: PD.Params = {};
  for (const [key, p] of Object.entries(params)) {
    const label = p.label ? undefined : unitLabel(key);
    if (p.type === 'group') out[key] = { ...p, params: relabel(p.params), ...(label ? { label } : {}) } as PD.Any;
    else out[key] = label ? ({ ...p, label } as PD.Any) : p;
  }
  return out;
}

/** An animation's params (with defaults), by name; undefined when none is registered under it. */
export function paramsOf(plugin: PluginContext, name: string): PD.Params | undefined {
  const params = plugin.managers.animation.animations.find((a) => a.name === name)?.params(plugin) as
    | PD.Params
    | undefined;
  return params && relabel(params);
}

/** The structures' wiggle cells: the per-group layers of the uncertainty wiggle. */
export function wiggleCells(plugin: PluginContext) {
  return plugin.state.data.select(StateSelection.Generators.ofTransformer(WiggleStructureRepresentation3DFromBundle));
}

/** The strength of the uncertainty wiggle: its largest layer (the most uncertain atoms move that far). */
export function uncertaintyScale(plugin: PluginContext): number {
  let scale = 0;
  for (const cell of wiggleCells(plugin)) {
    const layers = (cell.transform.params as { layers?: readonly { value: number }[] } | undefined)?.layers ?? [];
    for (const l of layers) if (l.value > scale) scale = l.value;
  }
  return scale;
}

export function wiggleOf(plugin: PluginContext): WiggleState {
  const amplitude = plugin.managers.structure.component.state.options.animation.wiggleAmplitude;
  if (amplitude > 0) return { mode: 'dynamics', amplitude };
  return { mode: wiggleCells(plugin).length ? 'uncertainty' : 'off', amplitude: 0 };
}

const sameInfo = (a: AnimationInfo, b: AnimationInfo) =>
  a.name === b.name && a.title === b.title && a.description === b.description && a.available === b.available;

/**
 * The animation state, followed only while it is watched (a view bound to it, or a caller subscribed): the plugin is not
 * listened to otherwise, and `value` reads it fresh. Bursts of events are read once per animation frame.
 */
export function createAnimationSource(plugin: PluginContext): AnimationSource {
  const manager = plugin.managers.animation;
  const read = (): AnimationState => ({
    animations: animationInfos(plugin),
    ...(manager.isAnimating && manager.current ? { playing: manager.current.anim.name } : {}),
    wiggle: wiggleOf(plugin),
  });
  const state = watchedState(
    read,
    () =>
      merge(
        manager.events.updated,
        plugin.state.data.events.changed,
        plugin.managers.structure.component.events.optionsUpdated,
      ),
    {
      compare: (a, b) =>
        a.playing === b.playing &&
        a.wiggle.mode === b.wiggle.mode &&
        a.wiggle.amplitude === b.wiggle.amplitude &&
        sameItems(a.animations, b.animations, sameInfo),
    },
  );
  return { state, dispose: () => {} };
}
