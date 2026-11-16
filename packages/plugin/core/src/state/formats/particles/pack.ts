/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Task } from '@molstar/core/task';
import { Sniffable } from '@molstar/core/util/sniffable';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import {
  createParticleListFromPack,
  type ParticlePackDocument,
  ParticlePackMagic,
  readParticlePack,
} from '@molstar/model/formats/particles/pack';
import { PluginStateTransform, PluginStateObject as SO } from '@molstar/plugin/state/objects';
import { DataFormatProvider } from '@molstar/plugin/state/formats/provider';
import type { PluginContext } from '@molstar/plugin/context';
import type { PluginRegistryEntry } from '@molstar/plugin/spec';
import { Asset } from '@molstar/core/util/assets';
import { isParticlePackStreamable, openParticlePackStream } from '@molstar/model/formats/particles/pack-stream';
import { getParticlePackSource } from '@molstar/plugin/state/helpers/particle-pack-source';
import { ParticleTileSetFromPack, ParticleTilesView } from '@molstar/plugin/state/transforms/particles/tiles';
import { ParticlesRepresentation3D } from '@molstar/plugin/state/transforms/particles/representation';
import { ParticlesFormatCategory } from './category.js';
import { type ParticleFormatData, complexVisuals, ComplexParticleVisuals } from './provider.js';

export { ParticleListFromPack };
type ParticleListFromPack = typeof ParticleListFromPack;
/**
 * The particle list of a particle pack (`formats/particles/pack.ts`): its arrays views of the data's bytes, which stay
 * the list's storage (not copied)
 */
const ParticleListFromPack = PluginStateTransform.BuiltIn({
  name: 'particle-list-from-pack',
  display: { name: 'Particle List from Pack', description: 'Create ParticleList from a particle pack (.mspk).' },
  from: [SO.Data.Binary, SO.Data.Response],
  to: SO.Particle.List,
  params: { label: PD.Optional(PD.Text('')) },
})({
  apply({ a, params }) {
    return Task.create('Read Particle Pack', async () => {
      const list = createParticleListFromPack(readParticlePack(SO.Data.readBytes(a)));
      const label = params.label || list.label || a.label || 'Particles';
      return new SO.Particle.List(list, { label, description: 'Particle Pack' });
    });
  },
});

export const ParticlePackProvider = DataFormatProvider({
  name: 'particle_pack',
  label: 'Particle Pack',
  description: 'Particle Pack (a particle list and its targets, read as views)',
  category: ParticlesFormatCategory,
  binaryExtensions: ['mspk'],
  sniffsResponse: true,
  isApplicable: (_info, data) => Sniffable.of(data).sniffBytes(ParticlePackMagic, 0),
  parse: async (plugin, data, params?: { label?: string }) => {
    const list = plugin.state.data.build().to(data).apply(ParticleListFromPack, { label: params?.label });
    await list.commit({ revertOnError: true });
    return { format: data, list: list.selector };
  },
  visuals: (plugin: PluginContext, data: ParticleFormatData) => complexVisuals(plugin, data),
  stream: (plugin: PluginContext, asset: Asset, options: { readonly visuals: boolean }) =>
    openParticlePackStreamed(plugin, asset, options),
});

/**
 * A tile-ordered pack opened by range (`ParticleTileSetFromPack`): its tiles read as a view following the camera wants
 * them (its default visuals: that view and its targets' representation, colored as a whole pack's), never whole; false
 * where it is not tile-ordered, not a pack as it is (gzipped), or not read by range (a file not held, a server ignoring
 * ranges): read whole then
 */
async function openParticlePackStreamed(plugin: PluginContext, asset: Asset, options: { readonly visuals: boolean }) {
  const read = await getParticlePackSource(asset, plugin.managers.asset, (m) => plugin.log.warn(`Particle pack: ${m}`));

  if (!read) return false;
  // not a pack as it is (gzipped: `.mspk.gz`; a server's error page): read whole then, which inflates it or says why not
  let document: ParticlePackDocument;
  try {
    ({ document } = await openParticlePackStream(read));
  } catch {
    return false;
  }
  if (!isParticlePackStreamable(document)) return false;

  const source = Asset.isUrl(asset)
    ? { name: 'url' as const, params: asset }
    : { name: 'file' as const, params: asset as Asset.File };
  const tiles = await plugin
    .build()
    .toRoot()
    .apply(ParticleTileSetFromPack, { source })
    .commit({ revertOnError: true });
  if (!options.visuals || !tiles.data) return true;
  const { catalog } = tiles.data.set;
  const colorTheme =
    catalog.entityInfo && catalog.compartmentInfo
      ? 'particle-hierarchy'
      : catalog.entityInfo
        ? 'particle-entity'
        : catalog.compartmentInfo
          ? 'particle-compartment'
          : 'particle-index';
  await plugin
    .build()
    .to(tiles.ref)
    .apply(ParticleTilesView, {})
    .apply(ParticlesRepresentation3D, {
      type: { name: 'target', params: {} },
      colorTheme: { name: colorTheme, params: {} },
    })
    .commit();
  return true;
}

/** The particle pack data format with its action and the representations its visuals use. */
export const ParticlePack: PluginRegistryEntry = {
  formats: [ParticlePackProvider],
  actions: [ParticleListFromPack],
  ...ComplexParticleVisuals,
};
