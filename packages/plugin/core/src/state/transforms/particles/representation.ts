/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginStateTransform, PluginStateObject as SO } from '@molstar/plugin/state/objects';
import type { PluginContext } from '@molstar/plugin/context';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Task } from '@molstar/core/task';
import { Theme } from '@molstar/graphics/theme/theme';
import { StateTransformer } from '@molstar/core/state';

export { ParticlesRepresentation3D };
type ParticlesRepresentation3D = typeof ParticlesRepresentation3D;
const ParticlesRepresentation3D = PluginStateTransform.BuiltIn({
  name: 'particles-representation-3d',
  display: '3D Representation',
  from: SO.Particle.List,
  to: SO.Particle.Representation3D,
  params: (a, ctx: PluginContext) => {
    const { registry, themes: themeCtx } = ctx.representation.particles;
    const type = registry.get(registry.default.name);

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
        type: PD.Mapped<any>(registry.default.name, registry.types, (name) =>
          PD.Group<any>(registry.get(name).getParams(themeCtx, undefined as any)),
        ),
        colorTheme: PD.Mapped<any>(
          type.defaultColorTheme.name,
          themeCtx.colorThemeRegistry.types,
          (name) => PD.Group<any>(themeCtx.colorThemeRegistry.get(name).getParams({})),
          colorThemeInfo,
        ),
        sizeTheme: PD.Mapped<any>(type.defaultSizeTheme.name, themeCtx.sizeThemeRegistry.types, (name) =>
          PD.Group<any>(themeCtx.sizeThemeRegistry.get(name).getParams({})),
        ),
      };
    }

    const dataCtx = { particles: a.data };
    const colorThemeInfo = {
      help: (value: { name: string; params: {} }) => {
        const { name, params } = value;
        const p = themeCtx.colorThemeRegistry.get(name);
        const ct = p.factory(dataCtx, params);
        return { description: ct.description, legend: ct.legend };
      },
    };

    return {
      type: PD.Mapped<any>(registry.default.name, registry.getApplicableTypes(a.data), (name) =>
        PD.Group<any>(registry.get(name).getParams(themeCtx, a.data)),
      ),
      colorTheme: PD.Mapped<any>(
        type.defaultColorTheme.name,
        themeCtx.colorThemeRegistry.getApplicableTypes(dataCtx),
        (name) => PD.Group<any>(themeCtx.colorThemeRegistry.get(name).getParams(dataCtx)),
        colorThemeInfo,
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
    return Task.create('Particles Representation', async (ctx) => {
      const themes = plugin.representation.particles.themes;
      const provider = plugin.representation.particles.registry.get(params.type.name);
      const repr = provider.factory({ webgl: plugin.canvas3d?.webgl, ...themes }, provider.getParams);
      repr.setTheme(Theme.create(themes, { particles: a.data }, params));
      const props = params.type.params || {};
      await repr.createOrUpdate(props, a.data).runInContext(ctx);
      return new SO.Particle.Representation3D({ repr, sourceData: a.data }, { label: provider.label });
    });
  },
  update({ a, b, oldParams, newParams }, plugin: PluginContext) {
    return Task.create('Particles Representation', async (ctx) => {
      if (newParams.type.name !== oldParams.type.name) return StateTransformer.UpdateResult.Recreate;

      const provider = plugin.representation.particles.registry.get(newParams.type.name);
      if (provider.mustRecreate?.(oldParams.type.params, newParams.type.params))
        return StateTransformer.UpdateResult.Recreate;

      const themes = plugin.representation.particles.themes;
      b.data.repr.setTheme(Theme.create(themes, { particles: a.data }, newParams));
      const props = { ...b.data.repr.props, ...newParams.type.params };
      await b.data.repr.createOrUpdate(props, a.data).runInContext(ctx);
      b.data.sourceData = a.data;
      return StateTransformer.UpdateResult.Updated;
    });
  },
});
