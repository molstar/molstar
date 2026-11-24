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
import { ParticleTileDrawings } from '@molstar/model/model/particles/particle-tile-drawings';
import { ParticleNesting } from '@molstar/model/model/particles/particle-nesting';
import type { ParticleList } from '@molstar/model/model/particles/particle-list';

import {
  assertRepresentationScope,
  emptyMappedParams,
  hasColorThemes,
  hasSizeThemes,
} from '../../helpers/representation-registry.js';
import {
  getRepresentationSignatures,
  prepareRepresentation,
  registerTransformSignatures,
} from '@molstar/plugin/state/precompile';

/**
 * What a list is themed by: a list with templates by its leaves by level (`ParticleNesting.getThemeList`: a slot list of
 * its top list's and its templates' leaves, no copy expanded; its representation draws it so, or flattened where it
 * cannot, its loci the list's by path), any other as it is
 */
function getThemeData(list: ParticleList) {
  return ParticleNesting.hasTemplates(list) ? ParticleNesting.getThemeList(list) : list;
}

export { ParticlesRepresentation3D };
type ParticlesRepresentation3D = typeof ParticlesRepresentation3D;
const ParticlesRepresentation3D = PluginStateTransform.BuiltIn({
  name: 'particles-representation-3d',
  display: '3D Representation',
  from: SO.Particle.List,
  to: SO.Particle.Representation3D,
  params: (a, ctx: PluginContext) => {
    const scope = ctx.representation.particles;
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
              PD.Group<any>(registry.get(name).getParams(themeCtx, undefined as any)),
            )
          : emptyMappedParams(),
        colorTheme: hasColor
          ? PD.Mapped<any>(
              type.defaultColorTheme.name,
              themeCtx.colorThemeRegistry.types,
              (name) => PD.Group<any>(themeCtx.colorThemeRegistry.get(name).getParams({})),
              colorThemeInfo,
            )
          : emptyMappedParams(),
        sizeTheme: hasSize
          ? PD.Mapped<any>(type.defaultSizeTheme.name, themeCtx.sizeThemeRegistry.types, (name) =>
              PD.Group<any>(themeCtx.sizeThemeRegistry.get(name).getParams({})),
            )
          : emptyMappedParams(),
      };
    }

    const data = getThemeData(a.data);
    const dataCtx = { particles: data };
    const colorThemeInfo = {
      help: (value: { name: string; params: {} }) => {
        const { name, params } = value;
        const p = themeCtx.colorThemeRegistry.get(name);
        const ct = p.factory(dataCtx, params);
        return { description: ct.description, legend: ct.legend };
      },
    };

    // a type of the theme list's (its params of it: the target's, by the leaves' targets) or the list's own, which
    // the types drawing a list whole take (flattened) and a slot list's refuse: not of the theme list's alone (a
    // spacefill's params else mapped to the default type's, a target representation more)
    const applies = (name: string) => {
      const provider = registry.get(name);
      return provider.isApplicable(data) ? data : provider.isApplicable(a.data) ? a.data : undefined;
    };
    const types = registry.types.filter(([name]) => applies(name));
    return {
      type: registry.default
        ? PD.Mapped<any>(registry.default.name, types, (name) =>
            PD.Group<any>(registry.get(name).getParams(themeCtx, applies(name) ?? data)),
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
  canAutoUpdate({ oldParams, newParams }) {
    return oldParams.type.name === newParams.type.name;
  },
  apply({ a, params }, plugin: PluginContext) {
    assertRepresentationScope('particles', plugin.representation.particles);
    return Task.create('Particles Representation', async (ctx) => {
      const themes = plugin.representation.particles.themes;
      const provider = plugin.representation.particles.registry.get(params.type.name);
      const repr = provider.factory({ graphics: plugin.canvas3d?.graphics, ...themes }, provider.getParams);
      const data = getThemeData(a.data);
      const theme = Theme.create(themes, { particles: data }, params);
      repr.setTheme(theme);
      const props = params.type.params || {};
      // compiling what it draws while it builds
      prepareRepresentation(plugin, repr, { ...repr.props, ...props }, theme, data);
      await repr.createOrUpdate(props, a.data).runInContext(ctx);
      // how it draws a tile view (its clip objects, its opacity): what the view picks its tiles by
      ParticleTileDrawings.set(repr, a.data, { ...repr.props, ...props });
      return new SO.Particle.Representation3D({ repr, sourceData: a.data }, { label: provider.label });
    });
  },
  update({ a, b, oldParams, newParams }, plugin: PluginContext) {
    assertRepresentationScope('particles', plugin.representation.particles);
    return Task.create('Particles Representation', async (ctx) => {
      if (newParams.type.name !== oldParams.type.name) return StateTransformer.UpdateResult.Recreate;

      const provider = plugin.representation.particles.registry.get(newParams.type.name);
      if (provider.mustRecreate?.(oldParams.type.params, newParams.type.params))
        return StateTransformer.UpdateResult.Recreate;

      const themes = plugin.representation.particles.themes;
      const data = getThemeData(a.data);
      b.data.repr.setTheme(Theme.create(themes, { particles: data }, newParams));
      // its data alone changed (a tile view's step), its params the same: the representation's props as they are, not
      // the params', which would undo a clipping gizmo's preview (a drag sets the clip on the representation alone)
      const dataOnly = JSON.stringify(oldParams.type.params) === JSON.stringify(newParams.type.params);
      const props = dataOnly ? { ...b.data.repr.props } : { ...b.data.repr.props, ...newParams.type.params };
      prepareRepresentation(plugin, b.data.repr, props, b.data.repr.theme, data);
      await b.data.repr.createOrUpdate(dataOnly ? {} : props, a.data).runInContext(ctx);
      b.data.sourceData = a.data;
      ParticleTileDrawings.set(b.data.repr, a.data, props);
      return StateTransformer.UpdateResult.Updated;
    });
  },
  dispose({ b }) {
    if (b) ParticleTileDrawings.set(b.data.repr, undefined, undefined);
  },
});

registerTransformSignatures(ParticlesRepresentation3D, (params, ctx) => {
  const { registry, themes } = ctx.plugin.representation.particles;
  return getRepresentationSignatures(registry.get(params.type.name), themes, params, ctx);
});
