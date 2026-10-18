/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StateTransform, StateTransformer } from '@molstar/core/state';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Canvas3DContext, Canvas3DParams, type Canvas3DProps } from '@molstar/graphics/canvas3d/canvas3d';
import { PluginBehavior } from '@molstar/plugin/behavior/behavior';
import type { PluginContext } from '@molstar/plugin/context';
import { InteractivityManager } from '@molstar/plugin/state/manager/interactivity';
import type { IconName, ReadonlyState } from '@molstar/prism-ui/core';
import { watchedState } from './watched.js';

// The Settings tab's groups (parity.md §5, api.md §2 `settings`): Canvas3D's params, the canvas context's,
// the interactivity manager's and the behaviors' params, mapped into groups by intent so no param is lost. The mapping
// is one table (`SettingsMapping`): what it does not name lands in Advanced, so a param added to Canvas3D shows up
// there. This file reads the plugin; the API (`api/settings.ts`) writes it.

export type SettingsGroupId =
  | 'appearance'
  | 'camera'
  | 'interaction'
  | 'performance'
  | 'viewer'
  | 'behaviors'
  | 'advanced';

/** One group of settings: its params, and the current values. */
export interface SettingsGroupState {
  readonly id: SettingsGroupId;
  readonly title: string;
  readonly icon: IconName;
  readonly description: string;
  /** Labels are long here (behaviors): the view gives them the row. */
  readonly wideLabels?: boolean;
  readonly params: PD.Params;
  readonly values: Record<string, unknown>;
  /** Where each param lives (the writer applies by it): Canvas3D's props, the canvas context's, the interactivity's. */
  readonly canvas: readonly string[];
  readonly context: readonly string[];
  readonly interactivity: readonly string[];
  /** A behavior's params by param name (the behaviors group): the behavior's state ref. */
  readonly behaviors?: Readonly<Record<string, string>>;
  /** The viewer's own (the viewer group). */
  readonly viewer?: boolean;
}

interface GroupDef {
  readonly id: SettingsGroupId;
  readonly title: string;
  readonly icon: IconName;
  readonly description: string;
  readonly wideLabels?: boolean;
  /** Canvas3D param paths: `key.child` for a child of a group param, `key` for the whole param (all its children). */
  readonly canvas?: readonly string[];
  readonly context?: boolean;
  readonly interactivity?: boolean;
  readonly behaviors?: boolean;
  readonly viewer?: boolean;
  /** Every Canvas3D param no other group names. */
  readonly rest?: boolean;
}

/** The mapping, once: group → which Canvas3D param paths (and which other sources). Order is the tab's. */
export const SettingsMapping: readonly GroupDef[] = [
  {
    id: 'appearance',
    title: 'Appearance',
    icon: 'palette',
    description: 'Background, colors, lighting and effects.',
    canvas: [
      'renderer.backgroundColor',
      'renderer.highlightColor',
      'renderer.selectColor',
      'renderer.dimColor',
      'renderer.highlightStrength',
      'renderer.selectStrength',
      'renderer.dimStrength',
      'transparentBackground',
      'checkeredTransparentBackground',
      'postprocessing',
      'cameraFog',
      'marking',
    ],
  },
  {
    id: 'camera',
    title: 'Camera',
    icon: 'camera',
    description: 'Projection, field of view, clipping and reset.',
    canvas: [
      'camera',
      'cameraClipping',
      'cameraResetDurationMs',
      'cameraResetEasing',
      'cameraResetTrajectory',
      'sceneRadiusFactor',
    ],
  },
  {
    id: 'interaction',
    title: 'Interaction',
    icon: 'mouse-pointer-2',
    description: 'Mouse and touch controls, picking.',
    canvas: ['trackball', 'interaction', 'pickPadding', 'userInteractionReleaseMs'],
    interactivity: true,
  },
  {
    id: 'performance',
    title: 'Performance',
    icon: 'gauge',
    description: 'Resolution, transparency and sampling.',
    canvas: ['multiSample', 'hiZ', 'dpoitIterations'],
    context: true,
  },
  {
    id: 'viewer',
    title: 'Viewer',
    icon: 'app-window',
    description: "The viewer's own: Undo.",
    viewer: true,
  },
  {
    id: 'behaviors',
    title: 'Behaviors',
    icon: 'toggle-right',
    description: "The plugin's behaviors and their params.",
    wideLabels: true,
    behaviors: true,
  },
  {
    id: 'advanced',
    title: 'Advanced',
    icon: 'wrench',
    description: 'Every other canvas param.',
    rest: true,
  },
];

export const SettingsGroupIds: readonly SettingsGroupId[] = SettingsMapping.map((g) => g.id);

/** Canvas3D's param paths, one level down: `key.child` for a group param's children, `key` for the others. */
export function canvasLeaves(): string[] {
  const out: string[] = [];
  for (const [key, p] of Object.entries(Canvas3DParams) as [string, PD.Any][]) {
    if (p.type === 'group') for (const child of Object.keys(p.params)) out.push(`${key}.${child}`);
    else out.push(key);
  }
  return out;
}

const claims = (entry: string, leaf: string) => leaf === entry || leaf.startsWith(`${entry}.`);

/** The Canvas3D param paths of each group, rest included (a path in exactly one group). */
export function canvasPathsByGroup(): Record<SettingsGroupId, string[]> {
  const leaves = canvasLeaves();
  const named = SettingsMapping.flatMap((g) => g.canvas ?? []);
  const out = {} as Record<SettingsGroupId, string[]>;
  for (const g of SettingsMapping) {
    out[g.id] = g.rest
      ? leaves.filter((l) => !named.some((n) => claims(n, l)))
      : leaves.filter((l) => g.canvas?.some((n) => claims(n, l)));
  }
  return out;
}

/**
 * Labels of params the Canvas card names too, by Canvas3D path: the Settings tab says them the same ("Fog", "Color",
 * "Transparent"), not by Mol*'s own ("Camera Fog", "Background Color", "Transparent Background").
 */
export const CanvasLabels: Readonly<Record<string, string>> = {
  cameraFog: 'Fog',
  'renderer.backgroundColor': 'Color',
  transparentBackground: 'Transparent',
};

const words = (key: string) => key.replace(/([a-z0-9])([A-Z])/g, '$1 $2').replace(/^./, (c) => c.toUpperCase());

const renamed = (p: PD.Any, label: string | undefined): PD.Any => (label ? { ...p, label } : p);

/** The params of a set of Canvas3D leaves: a whole param when all its children are in, else a group of the ones that are. */
function canvasParamsOf(leaves: readonly string[]): PD.Params {
  const byKey = new Map<string, string[]>();
  for (const leaf of leaves) {
    const [key, child] = leaf.split('.');
    const list = byKey.get(key);
    if (list) list.push(child);
    else byKey.set(key, child === undefined ? [] : [child]);
  }
  const params: PD.Params = {};
  for (const [key, original] of Object.entries(Canvas3DParams) as [string, PD.Any][]) {
    const children = byKey.get(key);
    if (!children) continue;
    if (original.type !== 'group' || children.length === Object.keys(original.params).length) {
      params[key] = renamed(original, CanvasLabels[key]);
    } else {
      const picked: PD.Params = {};
      for (const c of Object.keys(original.params)) {
        if (children.includes(c)) picked[c] = renamed(original.params[c], CanvasLabels[`${key}.${c}`]);
      }
      params[key] = PD.Group(picked, { label: original.label ?? words(key) });
    }
  }
  return params;
}

// Param definitions only, made from Mol*'s static `Canvas3DParams`: no viewer's state (values live in each viewer's plugin),
// so one cache for the page is harmless with several viewers (README decision 29).
const canvasGroupParams = new Map<SettingsGroupId, PD.Params>();
function paramsOfCanvasGroup(id: SettingsGroupId): PD.Params {
  let params = canvasGroupParams.get(id);
  if (!params) {
    params = canvasParamsOf(canvasPathsByGroup()[id]);
    canvasGroupParams.set(id, params);
  }
  return params;
}

/** A group's values from the props they live in: a group param's own children only. */
function pickValues(params: PD.Params, from: Record<string, any>, into: Record<string, unknown>) {
  for (const [key, p] of Object.entries(params)) {
    const v = from[key];
    if (v === undefined) continue;
    if (p.type === 'group' && v && typeof v === 'object') {
      const picked: Record<string, unknown> = {};
      for (const child of Object.keys(p.params)) if (child in v) picked[child] = v[child];
      into[key] = picked;
    } else into[key] = v;
  }
}

/** The behaviors with params the person can set, as a group: params by behavior, values, and refs. */
interface BehaviorEntry {
  readonly key: string;
  readonly ref: string;
  readonly label: string;
  readonly params: PD.Params;
  readonly values: Record<string, unknown>;
}

function behaviorEntries(plugin: PluginContext): BehaviorEntry[] {
  const state = plugin.state.behaviors;
  const root = state.cells.get(StateTransform.RootRef)?.obj;
  const out: BehaviorEntry[] = [];
  const used = new Set<string>();
  for (const cell of state.cells.values()) {
    const t = cell.transform;
    if (t.ref === StateTransform.RootRef || !cell.obj || !PluginBehavior.Behavior.is(cell.obj)) continue;
    let all: PD.Params;
    try {
      all = StateTransformer.getParamDefinition(t.transformer, root, plugin) as PD.Params;
    } catch {
      continue;
    }
    const params: PD.Params = {};
    const values: Record<string, unknown> = {};
    for (const [name, p] of Object.entries(all)) {
      if (p.isHidden) continue;
      params[name] = p;
      values[name] = (t.params as Record<string, unknown> | undefined)?.[name] ?? p.defaultValue;
    }
    if (!Object.keys(params).length) continue;
    // The behavior's transformer name, as agents say it (`representation-highlight-loci`); the ref when it repeats.
    const id = t.transformer.id;
    const name = id.slice(id.indexOf('.') + 1);
    const key = used.has(name) ? t.ref : name;
    used.add(key);
    out.push({ key, ref: t.ref, label: t.transformer.definition.display.name, params, values });
  }
  return out;
}

const behaviorGroups = new WeakMap<PluginContext, Map<string, PD.Group<any>>>();
/** A behavior's params as a sub-group, made once per behavior (the form compares params by identity). */
function behaviorGroup(plugin: PluginContext, e: BehaviorEntry): PD.Group<any> {
  let groups = behaviorGroups.get(plugin);
  if (!groups) behaviorGroups.set(plugin, (groups = new Map()));
  const id = `${e.key}:${Object.keys(e.params).join(',')}`;
  let g = groups.get(id);
  if (!g) {
    g = PD.Group(e.params, { label: e.label, isExpanded: false });
    groups.set(id, g);
  }
  return g;
}

/** The viewer's own params: whether every step takes a checkpoint (history.setCheckpointing). */
export const ViewerSettingsParams = {
  checkpointing: PD.Boolean(true, {
    label: 'Undo each step',
    description:
      'Take a checkpoint before each step so Undo can restore it. Turned off for a large scene, where it is slow: turn it back on here (earlier steps stay unsafe).',
  }),
};

/**
 * The groups with their params and current values, those the host can show: the canvas ones need a canvas, Performance
 * its context. `checkpointing` is the history's (the viewer group's value).
 */
export function readSettings(plugin: PluginContext, checkpointing: boolean): SettingsGroupState[] {
  const canvas = plugin.canvas3d?.props as Canvas3DProps | undefined;
  const context = plugin.canvas3dContext?.props as Readonly<Canvas3DContext.Props> | undefined;
  const out: SettingsGroupState[] = [];
  for (const def of SettingsMapping) {
    const base = { id: def.id, title: def.title, icon: def.icon, description: def.description } as const;
    const wide = def.wideLabels ? { wideLabels: true as const } : {};
    if (def.viewer) {
      out.push({
        ...base,
        params: ViewerSettingsParams,
        values: { checkpointing },
        canvas: [],
        context: [],
        interactivity: [],
        viewer: true,
      });
      continue;
    }
    if (def.behaviors) {
      const entries = behaviorEntries(plugin);
      if (!entries.length) continue;
      const params: PD.Params = {};
      const values: Record<string, unknown> = {};
      const refs: Record<string, string> = {};
      for (const e of entries) {
        params[e.key] = behaviorGroup(plugin, e);
        values[e.key] = e.values;
        refs[e.key] = e.ref;
      }
      out.push({ ...base, ...wide, params, values, canvas: [], context: [], interactivity: [], behaviors: refs });
      continue;
    }
    // The canvas groups need a canvas (Performance its context too): without them they are not listed.
    if (!canvas || (def.context && !context)) continue;
    const params: PD.Params = {};
    const values: Record<string, unknown> = {};
    const canvasKeys: string[] = [];
    const contextKeys: string[] = [];
    const interactivityKeys: string[] = [];
    const own = paramsOfCanvasGroup(def.id);
    Object.assign(params, own);
    canvasKeys.push(...Object.keys(own));
    pickValues(own, canvas, values);
    if (def.context && context) {
      Object.assign(params, Canvas3DContext.Params);
      contextKeys.push(...Object.keys(Canvas3DContext.Params));
      pickValues(Canvas3DContext.Params, context, values);
    }
    if (def.interactivity) {
      Object.assign(params, InteractivityManager.Params);
      interactivityKeys.push(...Object.keys(InteractivityManager.Params));
      pickValues(InteractivityManager.Params, plugin.managers.interactivity.props, values);
    }
    if (!Object.keys(params).length) continue;
    out.push({ ...base, params, values, canvas: canvasKeys, context: contextKeys, interactivity: interactivityKeys });
  }
  return out;
}

/** Whether two reads are the same, group by group (the source emits on change only). */
export function sameSettings(a: readonly SettingsGroupState[], b: readonly SettingsGroupState[]): boolean {
  if (a === b) return true;
  if (a.length !== b.length) return false;
  for (let i = 0; i < a.length; i++) {
    const x = a[i];
    const y = b[i];
    if (x.id !== y.id) return false;
    const keys = Object.keys(x.params);
    if (keys.length !== Object.keys(y.params).length || keys.some((k) => x.params[k] !== y.params[k])) return false;
    if (!PD.areEqual(x.params, x.values, y.values)) return false;
  }
  return true;
}

export interface SettingsSources {
  /**
   * The groups the host can show, with their current values: it follows the canvas (settings, a new canvas), the
   * interactivity and the behaviors. The viewer group's value is the history's: the reader sets it
   * (`withCheckpointing`). Followed only while subscribed (the Settings tab hidden holds no subscription: a slider drag
   * or a camera focus elsewhere does not read the settings), read once per animation frame.
   */
  readonly groups: ReadonlyState<readonly SettingsGroupState[]>;
}

export function createSettingsSources(plugin: PluginContext): SettingsSources {
  const groups = watchedState(
    () => readSettings(plugin, true),
    () => [
      plugin.events.canvas3d.settingsUpdated,
      plugin.behaviors.canvas3d.initialized,
      plugin.managers.interactivity.events.propsUpdated,
      plugin.state.behaviors.events.changed,
    ],
    { compare: sameSettings },
  );
  return { groups };
}

/** The groups with the viewer group's value set (it is the history's, which the sources do not hold). */
export function withCheckpointing(groups: readonly SettingsGroupState[], on: boolean): SettingsGroupState[] {
  return groups.map((g) => (g.viewer ? { ...g, values: { checkpointing: on } } : g));
}
