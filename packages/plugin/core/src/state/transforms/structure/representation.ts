/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginStateTransform, PluginStateObject as SO } from '@molstar/plugin/state/objects';
import type { PluginContext } from '@molstar/plugin/context';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { arrayEqual, deepEqual } from '@molstar/core/util';
import { Structure } from '@molstar/model/model/structure';
import { Task } from '@molstar/core/task';
import { Theme } from '@molstar/graphics/theme/theme';
import { StateTransformer } from '@molstar/core/state';
import { Color } from '@molstar/core/util/color';
import {
  assertRepresentationScope,
  emptyMappedParams,
  hasColorThemes,
  hasSizeThemes,
} from '../../helpers/representation-registry.js';
import { getRepresentationSignatures, registerTransformSignatures } from '@molstar/plugin/state/precompile';

export { StructureRepresentation3D };
type StructureRepresentation3D = typeof StructureRepresentation3D;

type ThemeCache = { themeProperties?: unknown[] };

/**
 * The property data the themes of a representation of `s` may read: the structure's, its root's and its models', as
 * the containers holding the values with their versions (a property computed again or with other props puts a new
 * container). Equal arrays: nothing a theme reads was computed again.
 */
function propertyState(s: Structure): unknown[] {
  const out: unknown[] = [];
  const add = (bag: { [name: string]: any }) => {
    for (const name of Object.keys(bag)) out.push(name, bag[name], bag[name]?.data?.version);
  };
  add(s.currentPropertyData);
  if (s.parent) add(s.inheritedPropertyData);
  for (const m of s.models) {
    add(m._staticPropertyData);
    add(m._dynamicPropertyData);
  }
  return out;
}
const StructureRepresentation3D = PluginStateTransform.BuiltIn({
  name: 'structure-representation-3d',
  display: '3D Representation',
  from: SO.Molecule.Structure,
  to: SO.Molecule.Structure.Representation3D,
  params: (a, ctx: PluginContext) => {
    const scope = ctx.representation.structure;
    const { registry, themes: themeCtx } = scope;
    // Param definitions never throw; without a registered representation or theme the params are empty mapped params
    const type = registry.default?.provider;
    const hasColor = !!type && hasColorThemes(scope);
    const hasSize = !!type && hasSizeThemes(scope);

    if (!a) {
      const colorThemeInfo = {
        help: (value: { name: string; params: {} }) => {
          const { name, params } = value;
          const p = themeCtx.colorThemeRegistry.get(name);
          const ct = p.factory({}, params);
          return { description: ct.description, legend: ct.legend };
        },
      };

      return {
        type: registry.default
          ? PD.Mapped<any>(registry.default.name, registry.types, (name) =>
              PD.Group<any>(registry.get(name).getParams(themeCtx, Structure.Empty)),
            )
          : emptyMappedParams(),
        colorTheme: hasColor
          ? PD.Mapped<any>(
              type.defaultColorTheme.name,
              themeCtx.colorThemeRegistry.types,
              (name) => PD.Group<any>(themeCtx.colorThemeRegistry.get(name).getParams({ structure: Structure.Empty })),
              colorThemeInfo,
            )
          : emptyMappedParams(),
        sizeTheme: hasSize
          ? PD.Mapped<any>(type.defaultSizeTheme.name, themeCtx.sizeThemeRegistry.types, (name) =>
              PD.Group<any>(themeCtx.sizeThemeRegistry.get(name).getParams({ structure: Structure.Empty })),
            )
          : emptyMappedParams(),
      };
    }

    const dataCtx = { structure: a.data };
    const colorThemeInfo = {
      help: (value: { name: string; params: {} }) => {
        const { name, params } = value;
        const p = themeCtx.colorThemeRegistry.get(name);
        const ct = p.factory(dataCtx, params);
        return { description: ct.description, legend: ct.legend };
      },
    };

    return {
      type: registry.default
        ? PD.Mapped<any>(registry.default.name, registry.getApplicableTypes(a.data), (name) =>
            PD.Group<any>(registry.get(name).getParams(themeCtx, a.data)),
          )
        : emptyMappedParams(),
      colorTheme: hasColor
        ? PD.Mapped<any>(
            type.defaultColorTheme.name,
            themeCtx.colorThemeRegistry.getApplicableTypes(dataCtx),
            (name) => PD.Group<any>(themeCtx.colorThemeRegistry.get(name).getParams(dataCtx)),
            colorThemeInfo,
          )
        : emptyMappedParams(),
      sizeTheme: hasSize
        ? PD.Mapped<any>(type.defaultSizeTheme.name, themeCtx.sizeThemeRegistry.getApplicableTypes(dataCtx), (name) =>
            PD.Group<any>(themeCtx.sizeThemeRegistry.get(name).getParams(dataCtx)),
          )
        : emptyMappedParams(),
    };
  },
})({
  canAutoUpdate({ a, oldParams, newParams }) {
    // TODO: other criteria as well?
    return (
      a.data.elementCount < 10000 ||
      (oldParams.type.name === newParams.type.name && newParams.type.params.quality !== 'custom')
    );
  },
  apply({ a, params, cache }, plugin: PluginContext) {
    assertRepresentationScope('structure', plugin.representation.structure);
    return Task.create('Structure Representation', async (ctx) => {
      const propertyCtx = { runtime: ctx, assetManager: plugin.managers.asset, errorContext: plugin.errorContext };
      const provider = plugin.representation.structure.registry.get(params.type.name);
      const data = provider.getData?.(a.data, params.type.params) || a.data;
      if (provider.ensureCustomProperties) await provider.ensureCustomProperties.attach(propertyCtx, data);
      const repr = provider.factory(
        { graphics: plugin.canvas3d?.graphics, ...plugin.representation.structure.themes },
        provider.getParams,
      );
      await Theme.ensureDependencies(propertyCtx, plugin.representation.structure.themes, { structure: data }, params);
      repr.setTheme(Theme.create(plugin.representation.structure.themes, { structure: data }, params));

      const props = params.type.params || {};
      await repr.createOrUpdate(props, data).runInContext(ctx);
      // after the visuals: a theme or a visual may read a property only then (its container made on first read)
      (cache as ThemeCache).themeProperties = propertyState(data);
      return new SO.Molecule.Structure.Representation3D({ repr, sourceData: a.data }, { label: provider.label });
    });
  },
  update({ a, b, oldParams, newParams, cache }, plugin: PluginContext) {
    assertRepresentationScope('structure', plugin.representation.structure);
    return Task.create('Structure Representation', async (ctx) => {
      if (newParams.type.name !== oldParams.type.name) return StateTransformer.UpdateResult.Recreate;

      const provider = plugin.representation.structure.registry.get(newParams.type.name);
      if (provider.mustRecreate?.(oldParams.type.params, newParams.type.params))
        return StateTransformer.UpdateResult.Recreate;

      const data = provider.getData?.(a.data, newParams.type.params) || a.data;
      const propertyCtx = { runtime: ctx, assetManager: plugin.managers.asset, errorContext: plugin.errorContext };
      if (provider.ensureCustomProperties) await provider.ensureCustomProperties.attach(propertyCtx, data);

      // An animation frame that changes the representation's own params (opacity, emissive, clipping) on the same
      // structure keeps the themes, and their dependencies' references. A theme holds the version of a property it
      // reads (secondary structure, ASA, a validation report): with any property of the structure or its models
      // computed again since the themes were made, they are made again, and so they are when the representation is
      // evaluated with its params as they were (an ancestor or a dependency, e.g. a volume a theme refers to, changed).
      const themeCache = cache as ThemeCache;
      const sameThemes =
        a.data === b.data.sourceData &&
        !deepEqual(oldParams.type, newParams.type) &&
        deepEqual(oldParams.colorTheme, newParams.colorTheme) &&
        deepEqual(oldParams.sizeTheme, newParams.sizeTheme) &&
        !!themeCache.themeProperties &&
        arrayEqual(themeCache.themeProperties, propertyState(data));
      if (!sameThemes) {
        Theme.releaseDependencies(plugin.representation.structure.themes, { structure: b.data.sourceData }, oldParams);
        await Theme.ensureDependencies(
          propertyCtx,
          plugin.representation.structure.themes,
          { structure: data },
          newParams,
        );
        b.data.repr.setTheme(Theme.create(plugin.representation.structure.themes, { structure: data }, newParams));
      }

      const props = { ...b.data.repr.props, ...newParams.type.params };
      await b.data.repr.createOrUpdate(props, data).runInContext(ctx);
      if (!sameThemes) themeCache.themeProperties = propertyState(data);
      b.data.sourceData = a.data;
      return StateTransformer.UpdateResult.Updated;
    });
  },
  dispose({ b, params }, plugin: PluginContext) {
    if (!b || !params) return;

    const structure = b.data.sourceData;
    const provider = plugin.representation.structure.registry.get(params.type.name);
    if (provider.ensureCustomProperties) provider.ensureCustomProperties.detach(structure);
    Theme.releaseDependencies(plugin.representation.structure.themes, { structure }, params);
  },
  interpolate(src, tar, t) {
    if (src.colorTheme.name !== 'uniform' || tar.colorTheme.name !== 'uniform') {
      return t <= 0.5 ? src : tar;
    }
    const from = src.colorTheme.params.value as Color,
      to = tar.colorTheme.params.value as Color;
    const value = Color.interpolate(from, to, t);
    return {
      type: t <= 0.5 ? src.type : tar.type,
      colorTheme: { name: 'uniform', params: { value } },
      sizeTheme: t <= 0.5 ? src.sizeTheme : tar.sizeTheme,
    };
  },
});

registerTransformSignatures(StructureRepresentation3D, (params, ctx) => {
  const { registry, themes } = ctx.plugin.representation.structure;
  return getRepresentationSignatures(registry.get(params.type.name), themes, params, ctx);
});
