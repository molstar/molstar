/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { PluginContext } from '@molstar/plugin/context';
import { StateBuilder, StateObjectRef, StateObjectSelector, StateTransform } from '@molstar/core/state';
import { Task } from '@molstar/core/task';
import { isProductionMode } from '@molstar/core/util/debug';
import { objectForEach } from '@molstar/core/util/object';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import {
  createStructureRepresentationParams,
  type StructureColorThemeRef,
  type StructureRepresentationProps,
  type StructureRepresentationRef,
  type StructureSizeThemeRef,
} from '../../helpers/structure-representation-params.js';
import type { PluginStateObject } from '../../objects.js';
import { StructureRepresentation3D } from '@molstar/plugin/state/transforms/structure/representation';
import {
  type BuiltInStructureRepresentationPresetAlias,
  type BuiltInStructureRepresentationPresetId,
  PresetStructureRepresentations,
} from './representation-presets/catalog.js';
import type { StructureRepresentationPresetProvider } from './representation-presets/types.js';
import { PresetRegistry } from '../preset-registry.js';
import { PluginConfig } from '@molstar/plugin/config';

// TODO factor out code shared with TrajectoryHierarchyBuilder?

export type StructureRepresentationPresetProviderRef =
  | BuiltInStructureRepresentationPresetId
  | BuiltInStructureRepresentationPresetAlias
  | StructureRepresentationPresetProvider
  | string;

type BuiltInPresetOf<P extends { id: string; alias?: string }, K extends string> = P extends P
  ? K extends P['id'] | NonNullable<P['alias']>
    ? P
    : never
  : never;
type BuiltInPreset<K extends string> = BuiltInPresetOf<
  PresetStructureRepresentations[keyof PresetStructureRepresentations],
  K
>;

export class StructureRepresentationBuilder {
  private registry = new PresetRegistry<StructureRepresentationPresetProvider>(
    'Structure representation preset registry',
  );
  private get dataState() {
    return this.plugin.state.data;
  }

  resolveProvider(ref: StructureRepresentationPresetProviderRef) {
    return typeof ref === 'string' ? this.registry.resolve(ref) : ref;
  }

  hasPreset(s: PluginStateObject.Molecule.Structure) {
    for (const p of this.registry.providers) {
      if (!p.isApplicable || p.isApplicable(s, this.plugin)) return true;
    }
    return false;
  }

  get providers(): ReadonlyArray<StructureRepresentationPresetProvider> {
    return this.registry.providers;
  }

  getPresets(s?: PluginStateObject.Molecule.Structure) {
    if (!s) return this.providers;
    const ret = [];
    for (const p of this.registry.providers) {
      if (p.isApplicable && !p.isApplicable(s, this.plugin)) continue;
      ret.push(p);
    }
    return ret;
  }

  getPresetSelect(s?: PluginStateObject.Molecule.Structure): PD.Select<string> {
    const options: [string, string, string | undefined][] = [];
    for (const p of this.registry.providers) {
      if (s && p.isApplicable && !p.isApplicable(s, this.plugin)) continue;
      options.push([p.id, p.display.name, p.display.group]);
    }
    const configured = this.registry.resolve(
      this.plugin.config.get(PluginConfig.Structure.DefaultRepresentationPreset) ?? '',
    )?.id;
    const defaultId = options.some((o) => o[0] === configured) ? configured! : (options[0]?.[0] ?? '');
    return PD.Select(defaultId, options);
  }

  getPresetsWithOptions(s: PluginStateObject.Molecule.Structure) {
    const options: [string, string][] = [];
    const map: { [K in string]: PD.Any } = Object.create(null);
    for (const p of this.registry.providers) {
      if (p.isApplicable && !p.isApplicable(s, this.plugin)) continue;

      options.push([p.id, p.display.name]);
      map[p.id] = p.params ? PD.Group(p.params(s, this.plugin)) : PD.EmptyGroup();
    }
    if (options.length === 0) return PD.MappedStatic('', { '': PD.EmptyGroup() });
    return PD.MappedStatic(options[0][0], map, { options });
  }

  /** The message `registerPreset` would throw for this provider, or `undefined` if it would succeed. */
  findConflict(provider: StructureRepresentationPresetProvider): string | undefined {
    return this.registry.findConflict(provider);
  }

  /** Whether a preset is registered under this id or alias. */
  has(idOrAlias: string) {
    return this.registry.has(idOrAlias);
  }

  registerPreset(provider: StructureRepresentationPresetProvider) {
    this.registry.register(provider);
  }

  unregisterPreset(provider: StructureRepresentationPresetProvider | string) {
    this.registry.unregister(provider);
  }

  applyPreset<K extends BuiltInStructureRepresentationPresetId | BuiltInStructureRepresentationPresetAlias>(
    parent: StateObjectRef<PluginStateObject.Molecule.Structure>,
    preset: K,
    params?: StructureRepresentationPresetProvider.Params<BuiltInPreset<K>>,
  ): Promise<StructureRepresentationPresetProvider.State<BuiltInPreset<K>>> | undefined;
  applyPreset<P = any, S extends {} = {}>(
    parent: StateObjectRef<PluginStateObject.Molecule.Structure>,
    provider: StructureRepresentationPresetProvider<P, S>,
    params?: P,
  ): Promise<S> | undefined;
  applyPreset(
    parent: StateObjectRef<PluginStateObject.Molecule.Structure>,
    presetIdOrAlias: string,
    params?: any,
  ): Promise<any> | undefined;
  applyPreset(
    parent: StateObjectRef,
    providerRef: string | StructureRepresentationPresetProvider,
    params?: any,
  ): Promise<any> | undefined {
    const provider = this.resolveProvider(providerRef);
    if (!provider) throw new Error(`Preset '${providerRef}' is not registered in this plugin`);

    const state = this.plugin.state.data;
    const cell = StateObjectRef.resolveAndCheck(state, parent);
    if (!cell) {
      if (!isProductionMode) console.warn(`Applying structure repr. provider to bad cell.`);
      return;
    }

    const pd = (provider.params?.(cell.obj, this.plugin) as PD.Params) || {};
    let prms = params || (provider.params ? PD.getDefaultValues(pd) : {});

    const defaults = this.plugin.config.get(PluginConfig.Structure.DefaultRepresentationPresetParams);
    prms = PD.merge(pd, defaults, prms);

    const task = Task.create(`${provider.display.name}`, () => provider.apply(cell, prms, this.plugin) as Promise<any>);
    return this.plugin.runTask(task);
  }

  /** `type`, `color`, and `size` are each a provider or a name (see `StructureRepresentationProps`). */
  async addRepresentation<
    R extends StructureRepresentationRef = StructureRepresentationRef,
    C extends StructureColorThemeRef = StructureColorThemeRef,
    S extends StructureSizeThemeRef = StructureSizeThemeRef,
  >(
    structure: StateObjectRef<PluginStateObject.Molecule.Structure>,
    props: StructureRepresentationProps<R, C, S>,
    options?: Partial<StructureRepresentationBuilder.AddRepresentationOptions>,
  ): Promise<StateObjectSelector<PluginStateObject.Molecule.Structure.Representation3D>>;
  async addRepresentation(
    structure: StateObjectRef<PluginStateObject.Molecule.Structure>,
    props: StructureRepresentationProps<any, any, any>,
    options?: Partial<StructureRepresentationBuilder.AddRepresentationOptions>,
  ) {
    const repr = this.dataState.build();
    const selector = this.buildRepresentation(repr, structure, props, options);
    if (!selector) return;

    await repr.commit();
    return selector;
  }

  /** `type`, `color`, and `size` are each a provider or a name (see `StructureRepresentationProps`). */
  buildRepresentation<
    R extends StructureRepresentationRef = StructureRepresentationRef,
    C extends StructureColorThemeRef = StructureColorThemeRef,
    S extends StructureSizeThemeRef = StructureSizeThemeRef,
  >(
    builder: StateBuilder.Root,
    structure: StateObjectRef<PluginStateObject.Molecule.Structure> | undefined,
    props: StructureRepresentationProps<R, C, S>,
    options?: Partial<StructureRepresentationBuilder.AddRepresentationOptions>,
  ): StateObjectSelector<PluginStateObject.Molecule.Structure.Representation3D>;
  buildRepresentation(
    builder: StateBuilder.Root,
    structure: StateObjectRef<PluginStateObject.Molecule.Structure> | undefined,
    props: StructureRepresentationProps<any, any, any>,
    options?: Partial<StructureRepresentationBuilder.AddRepresentationOptions>,
  ) {
    if (!structure) return;
    const data = StateObjectRef.resolveAndCheck(this.dataState, structure)?.obj?.data;
    if (!data) return;

    const params = createStructureRepresentationParams(this.plugin, data, props);
    return options?.tag
      ? builder
          .to(structure)
          .applyOrUpdateTagged(options.tag, StructureRepresentation3D, params, { state: options?.initialState })
          .selector
      : builder.to(structure).apply(StructureRepresentation3D, params, { state: options?.initialState }).selector;
  }

  constructor(public plugin: PluginContext) {
    objectForEach(PresetStructureRepresentations, (r) => this.registerPreset(r));
  }
}

export namespace StructureRepresentationBuilder {
  export interface AddRepresentationOptions {
    initialState?: Partial<StateTransform.State>;
    tag?: string;
  }
}
