/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginStateTransform, PluginStateObject as SO } from '@molstar/plugin/state/objects';
import type { PluginContext } from '@molstar/plugin/context';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
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

export { StructureRepresentation3D };
type StructureRepresentation3D = typeof StructureRepresentation3D;
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
        { webgl: plugin.canvas3d?.webgl, ...plugin.representation.structure.themes },
        provider.getParams,
      );
      await Theme.ensureDependencies(propertyCtx, plugin.representation.structure.themes, { structure: data }, params);
      repr.setTheme(Theme.create(plugin.representation.structure.themes, { structure: data }, params));

      const props = params.type.params || {};
      await repr.createOrUpdate(props, data).runInContext(ctx);
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

      // TODO: if themes had a .needsUpdate method the following block could
      //       be optimized and only executed conditionally
      Theme.releaseDependencies(plugin.representation.structure.themes, { structure: b.data.sourceData }, oldParams);
      await Theme.ensureDependencies(
        propertyCtx,
        plugin.representation.structure.themes,
        { structure: data },
        newParams,
      );
      b.data.repr.setTheme(Theme.create(plugin.representation.structure.themes, { structure: data }, newParams));

      const props = { ...b.data.repr.props, ...newParams.type.params };
      await b.data.repr.createOrUpdate(props, data).runInContext(ctx);
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
