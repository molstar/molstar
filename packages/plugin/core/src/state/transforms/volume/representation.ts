/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginStateTransform, PluginStateObject as SO } from '@molstar/plugin/state/objects';
import type { PluginContext } from '@molstar/plugin/context';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Volume } from '@molstar/model/model/volume';
import { Task } from '@molstar/core/task';
import { Theme } from '@molstar/graphics/theme/theme';
import { VolumeRepresentation3DHelpers } from './representation-helpers.js';
import { StateTransformer } from '@molstar/core/state';

export { VolumeRepresentation3D };
type VolumeRepresentation3D = typeof VolumeRepresentation3D;
const VolumeRepresentation3D = PluginStateTransform.BuiltIn({
  name: 'volume-representation-3d',
  display: '3D Representation',
  from: SO.Volume.Data,
  to: SO.Volume.Representation3D,
  params: (a, ctx: PluginContext) => {
    const { registry, themes: themeCtx } = ctx.representation.volume;
    const type = registry.get(registry.default?.name ?? '');

    if (!a) {
      return {
        type: PD.Mapped<any>(registry.default?.name ?? '', registry.types, (name) =>
          PD.Group<any>(registry.get(name).getParams(themeCtx, Volume.One)),
        ),
        colorTheme: PD.Mapped<any>(type.defaultColorTheme.name, themeCtx.colorThemeRegistry.types, (name) =>
          PD.Group<any>(themeCtx.colorThemeRegistry.get(name).getParams({ volume: Volume.One })),
        ),
        sizeTheme: PD.Mapped<any>(type.defaultSizeTheme.name, themeCtx.sizeThemeRegistry.types, (name) =>
          PD.Group<any>(themeCtx.sizeThemeRegistry.get(name).getParams({ volume: Volume.One })),
        ),
      };
    }

    const dataCtx = { volume: a.data };
    return {
      type: PD.Mapped<any>(registry.default?.name ?? '', registry.types, (name) =>
        PD.Group<any>(registry.get(name).getParams(themeCtx, a.data)),
      ),
      colorTheme: PD.Mapped<any>(
        type.defaultColorTheme.name,
        themeCtx.colorThemeRegistry.getApplicableTypes(dataCtx),
        (name) => PD.Group<any>(themeCtx.colorThemeRegistry.get(name).getParams(dataCtx)),
      ),
      sizeTheme: PD.Mapped<any>(
        type.defaultSizeTheme.name,
        themeCtx.sizeThemeRegistry.getApplicableTypes(dataCtx),
        (name) => PD.Group<any>(themeCtx.sizeThemeRegistry.get(name).getParams(dataCtx)),
      ),
    };
  },
})({
  canAutoUpdate({ oldParams, newParams }) {
    return oldParams.type.name === newParams.type.name;
  },
  apply({ a, params }, plugin: PluginContext) {
    return Task.create('Volume Representation', async (ctx) => {
      const propertyCtx = { runtime: ctx, assetManager: plugin.managers.asset, errorContext: plugin.errorContext };
      const provider = plugin.representation.volume.registry.get(params.type.name);
      if (provider.ensureCustomProperties) await provider.ensureCustomProperties.attach(propertyCtx, a.data);
      const repr = provider.factory(
        { webgl: plugin.canvas3d?.webgl, ...plugin.representation.volume.themes },
        provider.getParams,
      );
      await Theme.ensureDependencies(propertyCtx, plugin.representation.volume.themes, { volume: a.data }, params);
      repr.setTheme(
        Theme.create(
          plugin.representation.volume.themes,
          { volume: a.data, locationKinds: provider.locationKinds },
          params,
        ),
      );

      const props = params.type.params || {};
      await repr.createOrUpdate(props, a.data).runInContext(ctx);
      return new SO.Volume.Representation3D(
        { repr, sourceData: a.data },
        { label: provider.label, description: VolumeRepresentation3DHelpers.getDescription(props) },
      );
    });
  },
  update({ a, b, oldParams, newParams }, plugin: PluginContext) {
    return Task.create('Volume Representation', async (ctx) => {
      if (newParams.type.name !== oldParams.type.name) return StateTransformer.UpdateResult.Recreate;

      const provider = plugin.representation.volume.registry.get(newParams.type.name);
      if (provider.mustRecreate?.(oldParams.type.params, newParams.type.params))
        return StateTransformer.UpdateResult.Recreate;

      const propertyCtx = { runtime: ctx, assetManager: plugin.managers.asset, errorContext: plugin.errorContext };
      if (provider.ensureCustomProperties) await provider.ensureCustomProperties.attach(propertyCtx, a.data);

      // TODO: if themes had a .needsUpdate method the following block could
      //       be optimized and only executed conditionally
      Theme.releaseDependencies(plugin.representation.volume.themes, { volume: b.data.sourceData }, oldParams);
      await Theme.ensureDependencies(propertyCtx, plugin.representation.volume.themes, { volume: a.data }, newParams);
      b.data.repr.setTheme(
        Theme.create(
          plugin.representation.volume.themes,
          { volume: a.data, locationKinds: provider.locationKinds },
          newParams,
        ),
      );

      const props = { ...b.data.repr.props, ...newParams.type.params };
      await b.data.repr.createOrUpdate(props, a.data).runInContext(ctx);
      b.data.sourceData = a.data;
      b.description = VolumeRepresentation3DHelpers.getDescription(props);
      return StateTransformer.UpdateResult.Updated;
    });
  },
  dispose({ b, params }, plugin: PluginContext) {
    if (!b || !params) return;

    const volume = b.data.sourceData;
    const provider = plugin.representation.volume.registry.get(params.type.name);
    if (provider.ensureCustomProperties) provider.ensureCustomProperties.detach(volume);
    Theme.releaseDependencies(plugin.representation.volume.themes, { volume }, params);
  },
});
