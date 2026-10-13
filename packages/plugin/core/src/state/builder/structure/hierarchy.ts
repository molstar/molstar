/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { PluginContext } from '@molstar/plugin/context';
import { StateObjectRef } from '@molstar/core/state';
import { Task } from '@molstar/core/task';
import { isProductionMode } from '@molstar/core/util/debug';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { PluginStateObject } from '../../objects.js';
import type {
  BuiltInTrajectoryHierarchyPresetAlias,
  BuiltInTrajectoryHierarchyPresetId,
  PresetTrajectoryHierarchy,
} from './hierarchy-presets/catalog.js';
import type { TrajectoryHierarchyPresetProvider } from './hierarchy-presets/types.js';
import { PluginConfig } from '@molstar/plugin/config';
import { PresetRegistry } from '../preset-registry.js';

// TODO factor out code shared with StructureRepresentationBuilder?

export type TrajectoryHierarchyPresetProviderRef =
  | BuiltInTrajectoryHierarchyPresetId
  | BuiltInTrajectoryHierarchyPresetAlias
  | TrajectoryHierarchyPresetProvider
  | string;

type BuiltInPresetOf<P extends { id: string; alias?: string }, K extends string> = P extends P
  ? K extends P['id'] | NonNullable<P['alias']>
    ? P
    : never
  : never;
type BuiltInPreset<K extends string> = BuiltInPresetOf<PresetTrajectoryHierarchy[keyof PresetTrajectoryHierarchy], K>;

export class TrajectoryHierarchyBuilder {
  private registry = new PresetRegistry<TrajectoryHierarchyPresetProvider>('Trajectory hierarchy preset registry');

  resolveProvider(ref: TrajectoryHierarchyPresetProviderRef) {
    return typeof ref === 'string' ? this.registry.resolve(ref) : ref;
  }

  hasPreset(t: PluginStateObject.Molecule.Trajectory) {
    for (const p of this.registry.providers) {
      if (!p.isApplicable || p.isApplicable(t, this.plugin)) return true;
    }
    return false;
  }

  get providers(): ReadonlyArray<TrajectoryHierarchyPresetProvider> {
    return this.registry.providers;
  }

  getPresets(t?: PluginStateObject.Molecule.Trajectory) {
    if (!t) return this.providers;
    const ret = [];
    for (const p of this.registry.providers) {
      if (p.isApplicable && !p.isApplicable(t, this.plugin)) continue;
      ret.push(p);
    }
    return ret;
  }

  getPresetSelect(t?: PluginStateObject.Molecule.Trajectory): PD.Select<string> {
    const options: [string, string][] = [];
    for (const p of this.registry.providers) {
      if (t && p.isApplicable && !p.isApplicable(t, this.plugin)) continue;
      options.push([p.id, p.display.name]);
    }
    const configured = this.registry.resolve(
      this.plugin.config.get(PluginConfig.Structure.DefaultHierarchyPreset) ?? '',
    )?.id;
    const defaultId = options.some((o) => o[0] === configured) ? configured! : (options[0]?.[0] ?? '');
    return PD.Select(defaultId, options);
  }

  getPresetsWithOptions(t: PluginStateObject.Molecule.Trajectory) {
    const options: [string, string][] = [];
    const map: { [K in string]: PD.Any } = Object.create(null);
    for (const p of this.registry.providers) {
      if (p.isApplicable && !p.isApplicable(t, this.plugin)) continue;

      options.push([p.id, p.display.name]);
      map[p.id] = p.params ? PD.Group(p.params(t, this.plugin)) : PD.EmptyGroup();
    }
    if (options.length === 0) return PD.MappedStatic('', { '': PD.EmptyGroup() });
    return PD.MappedStatic(options[0][0], map, { options });
  }

  /** The message `registerPreset` would throw for this provider, or `undefined` if it would succeed. */
  findConflict(provider: TrajectoryHierarchyPresetProvider): string | undefined {
    return this.registry.findConflict(provider);
  }

  /** Whether a preset is registered under this id or alias. */
  has(idOrAlias: string) {
    return this.registry.has(idOrAlias);
  }

  registerPreset(provider: TrajectoryHierarchyPresetProvider) {
    this.registry.register(provider);
  }

  unregisterPreset(provider: TrajectoryHierarchyPresetProvider | string) {
    this.registry.unregister(provider);
  }

  applyPreset<K extends BuiltInTrajectoryHierarchyPresetId | BuiltInTrajectoryHierarchyPresetAlias>(
    parent: StateObjectRef<PluginStateObject.Molecule.Trajectory>,
    preset: K,
    params?: Partial<TrajectoryHierarchyPresetProvider.Params<BuiltInPreset<K>>>,
  ): Promise<TrajectoryHierarchyPresetProvider.State<BuiltInPreset<K>>> | undefined;
  applyPreset<P = any, S = {}>(
    parent: StateObjectRef<PluginStateObject.Molecule.Trajectory>,
    provider: TrajectoryHierarchyPresetProvider<P, S>,
    params?: P,
  ): Promise<S> | undefined;
  applyPreset(
    parent: StateObjectRef<PluginStateObject.Molecule.Trajectory>,
    presetIdOrAlias: string,
    params?: any,
  ): Promise<any> | undefined;
  applyPreset(
    parent: StateObjectRef,
    providerRef: string | TrajectoryHierarchyPresetProvider,
    params?: any,
  ): Promise<any> | undefined {
    const provider = this.resolveProvider(providerRef);
    if (!provider) throw new Error(`Preset '${providerRef}' is not registered in this plugin`);

    const state = this.plugin.state.data;
    const cell = StateObjectRef.resolveAndCheck(state, parent);
    if (!cell) {
      if (!isProductionMode) console.warn(`Applying hierarchy preset provider to bad cell.`);
      return;
    }

    const prms =
      params || (provider.params ? PD.getDefaultValues(provider.params(cell.obj, this.plugin) as PD.Params) : {});

    const task = Task.create(`${provider.display.name}`, () => provider.apply(cell, prms, this.plugin) as Promise<any>);
    return this.plugin.runTask(task);
  }

  constructor(public plugin: PluginContext) {}
}
