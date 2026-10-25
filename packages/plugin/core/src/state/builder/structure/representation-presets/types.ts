/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Adam Midlik <midlik@gmail.com>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import type { PresetProvider } from '../../preset-provider.js';
import type { PluginStateObject } from '../../../objects.js';
import { deepEqual } from '@molstar/core/util';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { type VisualQuality, VisualQualityOptions } from '@molstar/graphics/geo/geometry/base';
import type { ColorTheme } from '@molstar/graphics/theme/color';
import type { Structure } from '@molstar/model/model/structure';
import type { PluginContext } from '@molstar/plugin/context';
import type { StateObjectRef, StateObjectSelector } from '@molstar/core/state';
import type { StaticStructureComponentType } from '../../../helpers/structure-component.js';
import type { StructureSelectionQuery } from '@molstar/plugin/state/queries/structure/query';
import type { StructureFocusRepresentationProps } from '@molstar/plugin/behavior/dynamic/selection/structure-focus-representation';
import { StructureFocusRepresentationId } from '@molstar/plugin/behavior/dynamic/selection/structure-focus-representation/id';
import { createStructureColorThemeParams } from '../../../helpers/structure-representation-params.js';
import { ChainIdColorThemeProvider } from '@molstar/graphics/theme/color/chain-id';
import { OperatorNameColorThemeProvider } from '@molstar/graphics/theme/color/operator-name';

export interface StructureRepresentationPresetProvider<
  P = any,
  S extends _Result = _Result,
  Id extends string = string,
  Alias extends string = string,
> extends PresetProvider<PluginStateObject.Molecule.Structure, P, S, Id, Alias> {}
export function StructureRepresentationPresetProvider<
  P,
  S extends _Result,
  const Id extends string,
  const Alias extends string = string,
>(repr: StructureRepresentationPresetProvider<P, S, Id, Alias>) {
  return repr;
}
export namespace StructureRepresentationPresetProvider {
  export type Params<P extends StructureRepresentationPresetProvider> =
    P extends StructureRepresentationPresetProvider<infer T> ? T : never;
  export type State<P extends StructureRepresentationPresetProvider> =
    P extends StructureRepresentationPresetProvider<infer _, infer S> ? S : never;
  export type Result = {
    components?: { [name: string]: StateObjectSelector | undefined };
    representations?: { [name: string]: StateObjectSelector | undefined };
  };

  export const CommonParams = {
    ignoreHydrogens: PD.Optional(PD.Boolean(false)),
    ignoreHydrogensVariant: PD.Optional(PD.Select('all', PD.arrayToOptions(['all', 'non-polar'] as const))),
    ignoreLight: PD.Optional(PD.Boolean(false)),
    quality: PD.Optional(PD.Select<VisualQuality>('auto', VisualQualityOptions)),
    theme: PD.Optional(
      PD.Group({
        globalName: PD.Optional(PD.Text<ColorTheme.BuiltIn>('')),
        globalColorParams: PD.Optional(PD.Value<any>({}, { isHidden: true })),
        carbonColor: PD.Optional(
          PD.Select('chain-id', PD.arrayToOptions(['chain-id', 'operator-name', 'element-symbol'] as const)),
        ),
        symmetryColor: PD.Optional(PD.Text<ColorTheme.BuiltIn>('')),
        symmetryColorParams: PD.Optional(PD.Value<any>({}, { isHidden: true })),
        focus: PD.Optional(
          PD.Group({
            name: PD.Optional(PD.Text<ColorTheme.BuiltIn>('')),
            params: PD.Optional(PD.Value<ColorTheme.BuiltInParams<ColorTheme.BuiltIn>>({} as any)),
          }),
        ),
      }),
    ),
  };
  export type CommonParams = PD.ValuesFor<typeof CommonParams>;

  function getCarbonColorParams(name: 'chain-id' | 'operator-name' | 'element-symbol') {
    return name === 'chain-id'
      ? { name, params: ChainIdColorThemeProvider.defaultValues }
      : name === 'operator-name'
        ? { name, params: OperatorNameColorThemeProvider.defaultValues }
        : { name, params: {} };
  }

  function isSymmetry(structure: Structure) {
    return structure.units.some((u) => !u.conformation.operator.assembly && u.conformation.operator.spgrOp >= 0);
  }

  export function reprBuilder(plugin: PluginContext, params: CommonParams, structure?: Structure) {
    const update = plugin.state.data.build();
    const builder = plugin.builders.structure.representation;
    const options = plugin.managers.structure.component.state.options;
    const h = options.hydrogens;
    // the component options, as the component manager gives its representations: a representation follows them
    // (StructureComponentManager.setOptions) unless the preset sets its own value below
    const typeParams = {
      quality: options.visualQuality,
      ignoreHydrogens: h !== 'all',
      ignoreHydrogensVariant: (h === 'only-polar' ? 'non-polar' : 'all') as 'all' | 'non-polar',
      ignoreLight: options.ignoreLight,
      material: options.materialStyle,
      clip: options.clipObjects,
      interior: options.interior,
      animation: options.animation,
    };
    if (params.quality && params.quality !== 'auto') typeParams.quality = params.quality;
    if (params.ignoreHydrogens !== void 0) typeParams.ignoreHydrogens = !!params.ignoreHydrogens;
    if (params.ignoreHydrogensVariant !== void 0) typeParams.ignoreHydrogensVariant = params.ignoreHydrogensVariant;
    if (params.ignoreLight !== void 0) typeParams.ignoreLight = !!params.ignoreLight;
    const color: ColorTheme.BuiltIn | undefined = params.theme?.globalName ? params.theme?.globalName : void 0;
    const ballAndStickColor: ColorTheme.BuiltInParams<'element-symbol'> =
      params.theme?.carbonColor !== undefined
        ? { carbonColor: getCarbonColorParams(params.theme?.carbonColor), ...params.theme?.globalColorParams }
        : { ...params.theme?.globalColorParams };
    const symmetryColor: ColorTheme.BuiltIn | undefined =
      structure && params.theme?.symmetryColor ? (isSymmetry(structure) ? params.theme?.symmetryColor : color) : color;
    const symmetryColorParams = params.theme?.symmetryColorParams
      ? { ...params.theme?.globalColorParams, ...params.theme?.symmetryColorParams }
      : { ...params.theme?.globalColorParams };
    const globalColorParams = params.theme?.globalColorParams ? { ...params.theme?.globalColorParams } : undefined;
    const surfaceTypeParams = {
      ...typeParams,
      solidInterior: options.solidSurface,
    };

    return {
      update,
      builder,
      color,
      symmetryColor,
      symmetryColorParams,
      globalColorParams,
      typeParams,
      surfaceTypeParams,
      ballAndStickColor,
    };
  }

  /**
   * Sets the color theme of the focus representation behavior's target and surroundings representations. Does nothing
   * when the behavior is absent, and skips a representation whose type or the requested theme is not registered. The
   * theme defaults to the default color theme of each representation type.
   */
  export function updateFocusRepr<T extends ColorTheme.BuiltIn>(
    plugin: PluginContext,
    structure: Structure,
    themeName: T | undefined,
    themeParams: ColorTheme.BuiltInParams<T> | undefined,
  ) {
    if (!plugin.state.hasBehavior(StructureFocusRepresentationId)) return;

    const current = plugin.state.behaviors.cells.get(StructureFocusRepresentationId)?.params?.values as
      | StructureFocusRepresentationProps
      | undefined;
    if (!current) return;

    const { registry, themes } = plugin.representation.structure;
    if (themeName && !themes.colorThemeRegistry.has(themeName)) return;

    const colorTheme = (type: string) => {
      if (!registry.has(type)) return;
      const name = themeName || registry.get(type).defaultColorTheme.name;
      if (!themes.colorThemeRegistry.has(name)) return;
      return createStructureColorThemeParams(plugin, structure, type, name, themeParams);
    };
    const surroundings = colorTheme(current.surroundingsParams.type.name);
    const target = colorTheme(current.targetParams.type.name);
    if (!surroundings && !target) return;
    // presets call this once per structure, mostly with these themes already set: the update would only apply the
    // same params again (and walk the whole state for the focus representations). Compared with the params of the
    // last update started (the tree), and never while an earlier call's update waits: the last call wins, as before.
    const latest = plugin.state.behaviors.tree.transforms.get(StructureFocusRepresentationId)?.params as
      | StructureFocusRepresentationProps
      | undefined;
    const isSet = (theme: typeof target, old: unknown) => !theme || deepEqual(theme, old);
    if (
      latest &&
      !pendingFocusUpdates.has(plugin) &&
      isSet(surroundings, latest.surroundingsParams.colorTheme) &&
      isSet(target, latest.targetParams.colorTheme)
    ) {
      return;
    }

    const update = plugin.state.updateBehavior<StructureFocusRepresentationProps>(
      StructureFocusRepresentationId,
      (p) => {
        if (surroundings) p.surroundingsParams.colorTheme = surroundings;
        if (target) p.targetParams.colorTheme = target;
      },
    );
    pendingFocusUpdates.set(plugin, update);
    const settled = () => {
      if (pendingFocusUpdates.get(plugin) === update) pendingFocusUpdates.delete(plugin);
    };
    update.then(settled, settled);
    return update;
  }
}

/** the last focus theme update of `updateFocusRepr` per plugin, while it has not settled */
const pendingFocusUpdates = new WeakMap<PluginContext, Promise<void>>();

type _Result = StructureRepresentationPresetProvider.Result;

export const BuiltInPresetGroupName = 'Basic';

export function presetStaticComponent(
  plugin: PluginContext,
  structure: StateObjectRef<PluginStateObject.Molecule.Structure>,
  type: StaticStructureComponentType,
  params?: { label?: string; tags?: string[] },
) {
  return plugin.builders.structure.tryCreateComponentStatic(structure, type, params);
}

/**
 * Creates a component from a selection query. The component's key is `selection-<tag>`; import the query from its
 * module (for example `protein` from `@molstar/plugin/state/queries/structure/type`) instead of the query catalog.
 */
export function presetSelectionComponent(
  plugin: PluginContext,
  structure: StateObjectRef<PluginStateObject.Molecule.Structure>,
  query: StructureSelectionQuery,
  tag: string,
  params?: { label?: string; tags?: string[] },
) {
  return plugin.builders.structure.tryCreateComponentFromSelection(structure, query, `selection-${tag}`, params);
}
