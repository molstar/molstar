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
import { Particle, type ParticleList } from '@molstar/model/model/particles/particle-list';
import {
  isParticlePackListStreamable,
  isParticlePackStreamable,
  openParticleListFromPackStream,
  openParticlePackStream,
} from '@molstar/model/formats/particles/pack-stream';
import {
  getHeldParticlePackSource,
  getParticlePackSource,
  type ParticlePackHold,
} from '@molstar/plugin/state/helpers/particle-pack-source';
import { canCullOnGpu } from '@molstar/graphics/gl/renderable/merged';
import { ParticleTileSetFromPack, ParticleTilesView } from '@molstar/plugin/state/transforms/particles/tiles';
import { ParticlesRepresentation3D } from '@molstar/plugin/state/transforms/particles/representation';
import { ParticlesFormatCategory } from './category.js';
import { type ParticleFormatData, complexVisuals, ComplexParticleVisuals } from './provider.js';

export { ParticleListFromPack, ParticleListFromPackStream };
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

type ParticleListFromPackStream = typeof ParticleListFromPackStream;
/**
 * The nested list of a pack read by range (`openParticleListFromPackStream`): its top list and its templates' runs,
 * tile directories and coarse levels read, their particles read by page as the GPU's cull wants them (WebGPU's
 * nested draw: what reads all of a template's particles refuses it, read whole then). Read whole where the GPU's cull
 * is not there (`read: 'whole'`, or none at apply); the backend changed, by range or whole again as it can be drawn
 * (`syncPackStreams`)
 */
const ParticleListFromPackStream = PluginStateTransform.BuiltIn({
  name: 'particle-list-from-pack-stream',
  display: {
    name: 'Particle List from Pack (by range)',
    description:
      "A nested particle pack read by range: its templates' particles read as the view draws them, the pack never whole.",
  },
  from: SO.Root,
  to: SO.Particle.List,
  params: {
    source: PD.MappedStatic('url', {
      url: PD.Url('', { description: "The pack's URL: its server answering range requests." }),
      file: PD.File({ accept: '.mspk' }),
    }),
    label: PD.Optional(PD.Text('')),
    read: PD.Select<'range' | 'whole'>(
      'range',
      [
        ['range', 'By range'],
        ['whole', 'Whole'],
      ],
      { isHidden: true, description: "Read whole: the backend's (no GPU cull to draw it by range)." },
    ),
  },
})({
  apply({ params, cache }, plugin: PluginContext) {
    syncPackStreams(plugin);
    return Task.create('Read Particle Pack by range', async () => {
      const read = await getHeldParticlePackSource(
        params.source,
        cache as ParticlePackHold,
        plugin.managers.asset,
        (m) => plugin.log.warn(`Particle pack: ${m}`),
      );
      // no GPU cull to draw it by range (a session restored on WebGL, a canvas without one): read whole, as a pack
      // opened without range reads is
      if (params.read === 'whole' || !canCullOnGpu(plugin.canvas3d?.graphics)) {
        plugin.log.warn('Particle pack: read whole, its templates not drawn by range without WebGPU');
        const list = createParticleListFromPack(readParticlePack(await read.read(0, read.byteLength)));
        const label = params.label || list.label || 'Particles';
        return new SO.Particle.List(list, { label, description: 'Particle Pack' });
      }
      const list = await openParticleListFromPackStream(await openParticlePackStream(read));
      const label = params.label || list.label || 'Particles';
      return new SO.Particle.List(list, { label, description: 'Particle Pack · read by range' });
    });
  },
  dispose({ cache }, plugin: PluginContext) {
    const held = (cache as ParticlePackHold | undefined)?.held;
    if (held) plugin.managers.asset.delete(held);
  },
});

/** whether `list`'s templates are not all read (a pack opened by range: `Particle.DeferredParticles`) */
function hasUnreadTemplates(list: ParticleList) {
  for (const t of list.targetMapping?.values() ?? [])
    if (t.kind === 'particles' && Particle.getDeferred(t.particles)) return true;
  return false;
}

/**
 * How a pack read by range (`ParticleListFromPackStream`) is read again after a backend change: whole where its
 * templates are not read and the GPU's cull is gone, by range again where it was read whole for that and the cull is
 * back; undefined: as it is
 */
export function getPackStreamRead(culls: boolean, read: 'range' | 'whole', unread: boolean) {
  if (!culls && unread) return 'whole';
  if (culls && read === 'whole') return 'range';
  return undefined;
}

/** the plugins whose backend changes `syncPackStreams` follows */
const syncedPlugins = new WeakSet<PluginContext>();

/**
 * Once a plugin: its backend changed, each pack read by range (`ParticleListFromPackStream`) read again as it can be
 * drawn: whole where its templates are not read and the GPU's cull is gone (their particles what another backend
 * draws: it drew nothing, not refused), by range again where it was read whole for that and the cull is back
 */
function syncPackStreams(plugin: PluginContext) {
  if (syncedPlugins.has(plugin)) return;
  syncedPlugins.add(plugin);
  plugin.events.canvas3d.backendChanged.subscribe(() => {
    const culls = canCullOnGpu(plugin.canvas3d?.graphics);
    const state = plugin.state.data;
    const update = state.build();
    let changed = false;
    for (const cell of state.cells.values()) {
      if (cell.transform.transformer !== ParticleListFromPackStream) continue;
      const params = cell.transform.params as { readonly read?: 'range' | 'whole' };
      const list = cell.obj?.data as ParticleList | undefined;
      const read = getPackStreamRead(culls, params.read ?? 'range', !!list && hasUnreadTemplates(list));
      if (!read) continue;
      update.to(cell.transform.ref).update({ ...params, read });
      changed = true;
    }
    if (!changed) return;
    plugin.runTask(state.updateTree(update)).catch((e) => plugin.log.error(`Particle pack: read again failed: ${e}`));
  });
}

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
  const source = Asset.isUrl(asset)
    ? { name: 'url' as const, params: asset }
    : { name: 'file' as const, params: asset as Asset.File };
  // a nested pack read by range where the GPU's cull draws it (its templates' pages as the view wants them), else whole
  if (!isParticlePackStreamable(document)) {
    if (!isParticlePackListStreamable(document) || !canCullOnGpu(plugin.canvas3d?.graphics)) return false;
    const list = await plugin
      .build()
      .toRoot()
      .apply(ParticleListFromPackStream, { source })
      .commit({ revertOnError: true });
    if (options.visuals && list.data) await complexVisuals(plugin, { format: list.ref, list: list.ref });
    return true;
  }

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
