/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Box3D } from '@molstar/core/math/geometry/primitives/box3d';
import { Vec3 } from '@molstar/core/math/linear-algebra/3d/vec3';
import { Task } from '@molstar/core/task';
import { UserTiming } from '@molstar/core/task/util/user-timing';
import { StateTransformer } from '@molstar/core/state';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { createParticleTileLoader, type ParticleTileLoader } from '@molstar/model/model/particles/particle-tile-loader';
import { ParticleTileOpaque } from '@molstar/model/model/particles/particle-tile-opaque';
import type { ParticleList } from '@molstar/model/model/particles/particle-list';
import { ParticleSlots } from '@molstar/model/model/particles/particle-slots';
import { ParticleTileView } from '@molstar/model/model/particles/particle-tile-view';
import {
  getParticleTileSetRadius,
  ParticleTileVisibility,
} from '@molstar/model/model/particles/particle-tile-visibility';
import { ParticleTileDrawings } from '@molstar/model/model/particles/particle-tile-drawings';
import { ParticleTiles, type ParticleTileSet } from '@molstar/model/model/particles/particle-tiles';
import type { PluginContext } from '@molstar/plugin/context';
import { getParticleTileCamera, type ParticleTileCameraSnapshot } from '@molstar/plugin/state/helpers/particle-tiles';
import { getHeldParticlePackSource, type ParticlePackHold } from '@molstar/plugin/state/helpers/particle-pack-source';

import {
  createParticleTileSetFromPack,
  isParticlePackStreamable,
  openParticlePackStream,
} from '@molstar/model/formats/particles/pack-stream';
import { PluginStateTransform, PluginStateObject as SO } from '@molstar/plugin/state/objects';

// Particle tiles in the state: a list or trajectory subdivided into tiles loaded on
// demand, and views of some of them as particle lists (the representations as they are). Examples to test the model
// with: a format's tiles would be its own transform giving a tile set.

/** a tile set's cells along each axis to choose from: powers of two (a tile-ordered pack's tiles coarser by them read as they are) */
const CellOptions: [number, string][] = [1, 2, 4, 8, 16, 32, 64].map((n) => [n, `${n}`]);

export { ParticleTileSetFromList };
type ParticleTileSetFromList = typeof ParticleTileSetFromList;
const ParticleTileSetFromList = PluginStateTransform.BuiltIn({
  name: 'particle-tile-set-from-list',
  display: {
    name: 'Particle Tiles',
    description:
      'Subdivide a particle list or trajectory into tiles on a grid, their particles and targets loaded on demand.',
  },
  from: [SO.Particle.List, SO.Particle.Trajectory],
  to: SO.Particle.TileSet,
  params: {
    cells: PD.Select<'auto' | number>('auto', [['auto', 'Auto'], ...CellOptions], {
      description:
        "The grid cells along each axis. Auto: a tile-ordered pack's written grid (its tiles read as they are; another orders every particle again), else about 8k particles a tile.",
    }),
    delay: PD.Numeric(
      0,
      { min: 0, max: 5000, step: 10 },
      { description: 'Serve a tile or target after this many milliseconds, as a fetch would (0: at once).' },
    ),
    concurrency: PD.Numeric(4, { min: 1, max: 32, step: 1 }, { description: 'The tiles loading at once.' }),
    opacity: PD.Boolean(false, {
      description:
        "Compute the tiles' opacities where the data has none (a tile-ordered particle pack has them), so that tiles hide those behind them in a view by the camera.",
    }),
  },
})({
  apply({ a, params }, plugin: PluginContext) {
    return Task.create('Particle Tiles', async (ctx) => {
      // auto: what the data is ordered for (a tile-ordered pack's grid), else about a store page a tile
      const first = SO.Particle.List.is(a) ? a.data : (a as SO.Particle.Trajectory).data.representative;
      const cells = params.cells === 'auto' ? ParticleTiles.getAutoCells(first) : params.cells;
      const dimensions: [number, number, number] = [cells, cells, cells];
      const options = { delay: params.delay || undefined };
      let set = SO.Particle.List.is(a)
        ? ParticleTiles.fromList(a.data, dimensions, options)
        : ParticleTiles.fromTrajectory((a as SO.Particle.Trajectory).data, dimensions, options);
      // a tile-ordered pack's tiles have theirs (written with it): computed where none has
      const known = [...set.tiles.values()].some((t) => t.opacity);
      if (params.opacity && !known) set = await ParticleTileOpaque.computeOpacities(set).runInContext(ctx);
      const loader = createParticleTileLoader(set, { concurrency: params.concurrency, onError: tileError(plugin) });
      return new SO.Particle.TileSet(
        { set, loader },
        { label: a.label, description: `${cells}³ grid · ${set.tiles.size} tiles · ${set.count} particles` },
      );
    });
  },
  dispose({ b }) {
    b?.data.loader.dispose();
  },
});

export { ParticleTileSetFromPack };
/** A target's sites kept (`ParticleTileSetFromPackOptions.siteStride`) */
const SiteStrideOptions: [number, string][] = [
  [1, 'All'],
  [10, 'Every 10th'],
  [100, 'Every 100th'],
];

type ParticleTileSetFromPack = typeof ParticleTileSetFromPack;
const ParticleTileSetFromPack = PluginStateTransform.BuiltIn({
  name: 'particle-tile-set-from-pack',
  display: {
    name: 'Particle Pack Tiles',
    description:
      'A tile-ordered particle pack read by range as a tile set: its tiles and targets read as views want them, the pack never whole.',
  },
  from: SO.Root,
  to: SO.Particle.TileSet,
  params: {
    source: PD.MappedStatic('url', {
      url: PD.Url('', { description: "The pack's URL: its server answering range requests." }),
      file: PD.File({ accept: '.mspk' }),
    }),
    label: PD.Optional(PD.Text('')),
    concurrency: PD.Numeric(4, { min: 1, max: 32, step: 1 }, { description: 'The tiles loading at once.' }),
    levels: PD.Boolean(true, {
      description:
        "The pack's coarse levels where written: far tiles drawn as a sample of their particles (every 8th, 64th, ...), larger, so a whole cell fits few particles; picked by a view by the camera.",
    }),
    siteStride: PD.Select(1, SiteStrideOptions, {
      description:
        "A target's sites kept: all, or every 10th or 100th (drawn larger, as a coarser level of detail is), a tenth or a hundredth of the memory its targets take.",
    }),
  },
})({
  apply({ params, cache }, plugin: PluginContext) {
    return Task.create('Particle Pack Tiles', async () => {
      // a file held by the asset manager, not read: a session saved keeps it, one opened reads it from there
      const read = await getHeldParticlePackSource(params.source, cache as PackCache, plugin.managers.asset, (m) =>
        plugin.log.warn(`Particle pack: ${m}`),
      );
      const stream = await openParticlePackStream(read);
      if (!isParticlePackStreamable(stream.document)) throw new Error('particle pack not tile-ordered: read it whole');
      const set = await createParticleTileSetFromPack(stream, {
        label: params.label || undefined,
        siteStride: params.siteStride,
        levels: params.levels,
      });
      const loader = createParticleTileLoader(set, { concurrency: params.concurrency, onError: tileError(plugin) });
      const n = set.grid.dimensions[0];
      const base = set.levels?.base ?? set;
      return new SO.Particle.TileSet(
        { set, loader },
        {
          label: params.label || set.label || 'Particles',
          description: `${n}³ grid · ${base.tiles.size} tiles${set.levels ? ` · ${set.levels.count} levels` : ''} · ${base.count} particles · read by range${params.siteStride > 1 ? ` · every ${params.siteStride}th site` : ''}`,
        },
      );
    });
  },
  dispose({ b, cache }, plugin: PluginContext) {
    b?.data.loader.dispose();
    const held = (cache as PackCache | undefined)?.held;
    if (held) plugin.managers.asset.delete(held);
  },
});

/**
 * a tile set's loader's failures logged (a dropped request, a short answer): the tile loaded again where wanted again,
 * a few times (`ParticleTileLoadAttempts`)
 */
function tileError(plugin: PluginContext) {
  return (id: number, frame: number, error: unknown) =>
    plugin.log.warn(`Particle tiles: tile ${id}${frame > 0 ? ` at frame ${frame}` : ''} not read: ${error}`);
}

/** A pack tile set's cache (its cell's, kept as it is made again): its pack's file's hold */
type PackCache = ParticlePackHold;

/** A canvas camera as a param: its snapshot and its viewport's size (none: the canvas' own at the time) */
const CameraParam = PD.Value<
  | {
      snapshot: ParticleTileCameraSnapshot;
      width: number;
      height: number;
    }
  | undefined
>(undefined, { isHidden: true });

const TilesParam = PD.MappedStatic(
  'camera',
  {
    all: PD.EmptyGroup(),
    ids: PD.Group({ ids: PD.Text('', { description: 'Tile ids, e.g. "0, 3, 8-12".' }) }, { isFlat: true }),
    box: PD.Group(
      {
        min: PD.Vec3(Vec3.create(0, 0, 0), {}, { description: "The box's corner of least coordinates, Å." }),
        max: PD.Vec3(Vec3.create(100, 100, 100), {}, { description: "The box's corner of most coordinates, Å." }),
      },
      { isFlat: true },
    ),
    camera: PD.Group(
      {
        minPixels: PD.Numeric(
          1,
          { min: 0, max: 64, step: 0.5 },
          { description: 'The smallest size on screen worth drawing.' },
        ),
        occlusion: PD.Boolean(true, {
          description: 'Leave out tiles hidden behind others (where opacities are known).',
        }),
        opticalDepth: PD.Numeric(
          3,
          { min: 0.5, max: 12, step: 0.5 },
          {
            description:
              'The occlusion: the optical depth that hides what is behind (the light let through exp(-depth): 1 37%, 3 5%, 7 0.1%). Lower leaves out more tiles behind thin content, faster but at times a gap; higher only behind dense content.',
          },
        ),
        follow: PD.Boolean(true, {
          description:
            'Pick the tiles again as the camera moves (with the particle tiles behavior), unless a camera is given: the view stepping toward them as they load.',
        }),
        budget: PD.Numeric(
          1_000_000,
          { min: 0, max: 10_000_000, step: 10_000 },
          {
            description:
              'Following, the particles added at most per update (at least a tile; 0: all at once): the view adding the tiles in steps, frames drawn between.',
          },
        ),
        detailPixels: PD.Numeric(
          24,
          { min: 0, max: 64, step: 0.5 },
          {
            description:
              "With levels (a tile set's coarse levels): the radius on screen, in pixels, of the set's largest molecule at a tile from which its particles are all drawn; smaller, a coarse level's sample of them drawn larger (0: all). Alike for tiles alike far, whatever they hold.",
          },
        ),
        maxParticles: PD.Numeric(
          0,
          { min: 0, max: 50_000_000, step: 100_000 },
          {
            description:
              'The particles shown at most (0: no cap): the tiles past it, by priority (the least occluded, the largest on screen, the nearest first), left out. What bounds the memory a view takes (on a phone), whatever the data.',
          },
        ),
        stepMs: PD.Numeric(
          0,
          { min: 0, max: 100, step: 1 },
          {
            description:
              "Following, the time an update aims at (ms; 0: by the budget alone): its particles by how long the last updates took, at most the budget. Off by default: smaller steps are more of them, each with an update's fixed cost (a step a frame at 4 or 8 ms: three times the frames past a display's frame).",
          },
        ),
        camera: CameraParam,
      },
      { isFlat: true },
    ),
  },
  {
    options: [
      ['all', 'All'],
      ['ids', 'By Id'],
      ['box', 'In a Box'],
      ['camera', 'Seen by the Camera'],
    ],
  },
);

/** "0, 3, 8-12" as ids */
function parseIds(text: string): number[] {
  const out: number[] = [];
  for (const part of text.split(',')) {
    const [a, b] = part.split('-').map((s) => Number.parseInt(s.trim(), 10));
    if (Number.isNaN(a)) continue;
    // a single id, or a range
    const end = b === undefined || Number.isNaN(b) ? a : b;
    for (let i = Math.min(a, end); i <= Math.max(a, end); ++i) out.push(i);
  }
  return out;
}

export { pickTiles as pickParticleTiles };

/**
 * the occlusion's rays a side (`ParticleTileVisibilityOptions.resolution`): settled, and while what draws the tiles
 * changes (a cut dragged: a pick a frame or so, a quarter of the rays, about a quarter of the time) or the camera
 * moves (a pick every `intervalMs`)
 */
const SettledResolution = 64;
const MovingResolution = 32;

/** per canvas camera when it last changed (`setParticleTilesCameraChanged`, the follow behavior's: every change) */
const cameraChangedAt = new WeakMap<object, number>();

/**
 * `camera` (a canvas') changed `at`: its picks the quicker ones while it changed within `ParticleTileDrawings.SettleMs`,
 * a full pick once it settled (the follow behavior's pass after it)
 */
export function setParticleTilesCameraChanged(camera: object, at = performance.now()) {
  cameraChangedAt.set(camera, at);
}

function isCameraMoving(camera: object) {
  const at = cameraChangedAt.get(camera);
  return at !== undefined && performance.now() - at < ParticleTileDrawings.SettleMs;
}
/** a block a clip surface cuts split in eight once while the cut moves (settled: the visibility's default, 3) */
const MovingClipSplitDepth = 1;

/** the tiles `params` picks of `set`, by priority where by the camera, else ascending */
function pickTiles(
  set: ParticleTileSet,
  params: PD.Values<{ tiles: typeof TilesParam }>['tiles'],
  plugin: PluginContext,
  /**
   * the set's latest pick of these params (`latestPicks`), whatever camera it was of, where there is one: a pass not
   * the camera's (a tile loaded, a step's next) steps toward it, the camera's passes pick again
   */
  latest = false,
) {
  switch (params.name) {
    case 'all':
      return [...tilesOf(set).keys()].sort((x, y) => x - y);
    case 'ids':
      return parseIds(params.params.ids).filter((id) => tilesOf(set).has(id));
    case 'box': {
      const box = Box3D.create(Vec3.clone(params.params.min), Vec3.clone(params.params.max));
      return [...tilesOf(set).values()]
        .filter((t) => Box3D.overlaps(t.box, box))
        .map((t) => t.id)
        .sort((x, y) => x - y);
    }
    case 'camera': {
      const given = params.params.camera;
      const canvas = plugin.canvas3d;
      const snapshot = given?.snapshot ?? canvas?.camera.getSnapshot();
      // no camera (headless, none given): all
      if (!snapshot) return [...tilesOf(set).keys()].sort((x, y) => x - y);
      const width = given ? given.width : canvas!.camera.viewport.width;
      const height = given ? given.height : canvas!.camera.viewport.height;
      const { minPixels, occlusion } = params.params;
      const opticalDepth = params.params.opticalDepth ?? 3;
      // the set's last picks kept by their cameras: a following view's steps (and the behavior's check) of one camera
      // pick once (the visibility a pass over the tiles, their occlusion), each of a few views of the set its own
      const detailPixels = params.params.detailPixels ?? 24;
      // what draws its tiles (their representations' clip objects and opacity): picked again where it changes
      const drawing = ParticleTileDrawings.get(set);
      // a cut being moved (its drawings changing): the quicker pick, a coarser occlusion and a block a surface cuts
      // split once, picked again in full once it settles (the follow behavior's pass after `SettleMs`)
      const moving = ParticleTileDrawings.isChanging(set);
      // the canvas' camera moving (none given): the coarser occlusion too, picked again in full once it settles
      const cameraMoving = !given && !!canvas && isCameraMoving(canvas.camera);
      const resolution = moving || cameraMoving ? MovingResolution : SettledResolution;
      const paramsKey = `${minPixels},${occlusion},${opticalDepth},${params.params.maxParticles ?? 0},${detailPixels},${drawing.version},${resolution}`;
      if (latest && !given) {
        const last = latestPicks.get(set)?.get(paramsKey);
        if (last) return last;
      }
      const key = `${getCameraKey(snapshot, width, height, minPixels, occlusion, opticalDepth)},${paramsKey}`;
      const kept = getPick(set, key);
      if (kept) {
        setLatestPick(set, paramsKey, kept);
        return kept;
      }
      const v = UserTiming.span('tiles: visibility', () =>
        ParticleTileVisibility.get(set.levels?.base ?? set, getParticleTileCamera(snapshot, width, height), {
          minPixels,
          occlusion,
          opticalDepth,
          drawings: drawing.list.length > 0 ? drawing.list : undefined,
          resolution,
          clipSplitDepth: moving ? MovingClipSplitDepth : undefined,
        }),
      );
      // the tiles by priority within the cap of particles
      const max = params.params.maxParticles ?? 0;
      if (set.levels) {
        const ids = UserTiming.span('tiles: cut', () => pickNodes(set, v.tiles, detailPixels, max));
        keepPick(set, key, ids);
        setLatestPick(set, paramsKey, ids);
        return ids;
      }
      let ids = v.tiles.map((t) => t.id);
      if (max > 0) {
        let count = 0;
        let n = 0;
        while (n < ids.length && count + set.tiles.get(ids[n])!.count <= max) count += set.tiles.get(ids[n++])!.count;
        ids = ids.slice(0, Math.max(1, n));
      }
      keepPick(set, key, ids);
      setLatestPick(set, paramsKey, ids);
      return ids;
    }
  }
}

/** the tiles of `set` a camera sees and picks by: its levels' base's, else its own */
function tilesOf(set: ParticleTileSet) {
  return (set.levels?.base ?? set).tiles;
}

/**
 * How far the cap lowers the detail at most (`pickNodes`): a node's molecules drawn at most this many times
 * `detailPixels` in radius on screen, the least important tiles left out past it. Why: coarser near the camera draws
 * molecules filling the screen many layers deep (a phone's GPU stopped by the overdraw), where leaving out tiles behind
 * keeps what is near as it is
 */
const MaxDetailRaise = 4;

/**
 * The coarsest level of a visible tile whose sample, drawn larger, keeps the set's largest molecule's radius on screen
 * there (`size`: `getParticleTileSetRadius` by the tile's scale) at `detailPixels`: level `L` keeps `8^-L` of the particles, at least the fraction `(size / detailPixels)³`;
 * 0 (its own particles) from `detailPixels` on, or nearer than its reach
 */
function getWantedLevel(size: number, detailPixels: number, levelCount: number) {
  if (!(detailPixels > 0) || !Number.isFinite(size)) return 0;
  const fraction = Math.pow(Math.max(0, size) / detailPixels, 3);
  let level = levelCount;
  while (level > 0 && Math.pow(8, -level) < fraction) --level;
  return level;
}

/**
 * The nodes of `set`'s levels (`ParticleTileLevels`) and its tiles that show its visible tiles `visible` (by
 * priority): a cut of the octree, each visible tile shown by its coarsest node all of whose visible tiles want it
 * coarse enough (`getWantedLevel`), the nodes by the priority of their first visible tile. Within `max` particles:
 * the detail lowered (`detailPixels` raised, by bisection, at most `MaxDetailRaise` times) until the cut fits, then the
 * least important left out (at least one)
 */
function pickNodes(
  set: ParticleTileSet,
  visible: readonly { readonly id: number; readonly size: number; readonly scale: number }[],
  detailPixels: number,
  max: number,
) {
  const levels = set.levels!;
  // a tile's detail by the set's largest molecule at the tile (alike far, alike coarse), its own size where none known
  const radius = getParticleTileSetRadius(set);
  const sizes = visible.map((t) => (radius === undefined ? t.size : radius * t.scale));
  const count = levels.count;
  const dimensions = levels.base.grid.dimensions[0];
  // per visible tile its nodes, by level (0 itself; -1 none)
  const stride = count + 1;
  const nodes = new Int32Array(visible.length * stride);
  for (let i = 0; i < visible.length; ++i) {
    const id = visible[i].id;
    nodes[i * stride] = id;
    for (let level = 1; level <= count; ++level) {
      nodes[i * stride + level] = levels.getNode(level, ParticleTiles.getLevelCell(dimensions, id, level)) ?? -1;
    }
  }
  const cut = (detail: number) => {
    // per node the least level its visible tiles want
    const least = new Map<number, number>();
    for (let i = 0; i < visible.length; ++i) {
      const wanted = getWantedLevel(sizes[i], detail, count);
      for (let level = 1; level <= count; ++level) {
        const node = nodes[i * stride + level];
        if (node < 0) continue;
        const l = least.get(node);
        if (l === undefined || wanted < l) least.set(node, wanted);
      }
    }
    // a visible tile's shown by its coarsest node wanted at least that coarse (the nodes above it are not: fewer
    // visible tiles want a node the coarser it is), the nodes by their first visible tile
    const chosen = new Map<number, number>();
    let total = 0;
    for (let i = 0; i < visible.length; ++i) {
      let shown = nodes[i * stride];
      for (let level = count; level >= 1; --level) {
        const node = nodes[i * stride + level];
        if (node >= 0 && least.get(node)! >= level) {
          shown = node;
          break;
        }
      }
      if (chosen.has(shown)) continue;
      chosen.set(shown, i);
      total += set.tiles.get(shown)?.count ?? 0;
    }
    return { ids: [...chosen.keys()], total };
  };
  let picked = cut(detailPixels);
  if (max > 0 && picked.total > max) {
    // the least detail that fits: raised by bisection (geometric), the coarsest where none does
    let lo = detailPixels > 0 ? detailPixels : 1e-3;
    let hi = Math.max(lo, 1) * MaxDetailRaise;
    const coarsest = cut(hi);
    if (coarsest.total > max) {
      picked = coarsest;
    } else {
      picked = coarsest;
      for (let k = 0; k < 24; ++k) {
        const mid = Math.sqrt(lo * hi);
        const c = cut(mid);
        if (c.total <= max) {
          hi = mid;
          picked = c;
        } else {
          lo = mid;
        }
      }
    }
    if (picked.total > max) {
      const ids = picked.ids;
      let total = picked.total;
      let end = ids.length;
      while (end > 1 && total > max) total -= set.tiles.get(ids[--end])?.count ?? 0;
      return ids.slice(0, end);
    }
  }
  return picked.ids;
}

/**
 * per tile set its last picks by a camera (`pickTiles`), by key, the least recent first. Why several: views of one set
 * by other cameras (a box and a following view, two canvases) would each pick again at every step otherwise
 */
const lastPicks = new WeakMap<ParticleTileSet, Map<string, number[]>>();
/** the picks a set keeps */
const KeptPicks = 4;

/** the pick of `set` by `key` where kept, made the most recent (the same array: `wantTiles` compares it by identity) */
function getPick(set: ParticleTileSet, key: string) {
  const picks = lastPicks.get(set);
  const ids = picks?.get(key);
  if (!picks || !ids) return undefined;
  picks.delete(key);
  picks.set(key, ids);
  return ids;
}

/**
 * per tile set its latest pick by params other than the camera's (`pickTiles`' `latest`). Why: a pass woken by a tile
 * loaded or a step's next, while the camera moves, would pick again for the camera where it is now: a pick a frame
 * (the visibility a pass over the tiles), the camera's own passes throttled for nothing
 */
const latestPicks = new WeakMap<ParticleTileSet, Map<string, number[]>>();

function setLatestPick(set: ParticleTileSet, paramsKey: string, ids: number[]) {
  let picks = latestPicks.get(set);
  if (!picks) {
    picks = new Map();
    latestPicks.set(set, picks);
  }
  picks.set(paramsKey, ids);
}

function keepPick(set: ParticleTileSet, key: string, ids: number[]) {
  let picks = lastPicks.get(set);
  if (!picks) {
    picks = new Map();
    lastPicks.set(set, picks);
  }
  picks.set(key, ids);
  while (picks.size > KeptPicks) picks.delete(picks.keys().next().value!);
}

/**
 * what a camera's pick is of: its pose and lens, its near planes too (`getParticleTileCamera`'s: `minNear`, and the
 * clipping's, of its focus' radius)
 */
function getCameraKey(
  s: ParticleTileCameraSnapshot,
  width: number,
  height: number,
  minPixels: number,
  occlusion: boolean,
  opticalDepth: number,
) {
  return [
    s.mode,
    s.fov,
    s.minNear,
    s.radius,
    s.radiusMax,
    s.clipFar,
    s.minFar,
    ...s.position,
    ...s.target,
    ...s.up,
    width,
    height,
    minPixels,
    occlusion,
    opticalDepth,
  ].join(',');
}

type TilesValues = PD.Values<{ tiles: typeof TilesParam }>['tiles'];

/**
 * per tile set the particles of its following views' next step by time (`ParticleTilesFollowCamera`: the last steps'
 * particles a ms times `stepMs`), at most their budget
 */
const stepParticles = new WeakMap<ParticleTileSet, number>();

/** the particles of `set`'s following views' next step (by the time the last took), none: their budget alone */

export function setParticleTilesStep(set: ParticleTileSet, particles: number | undefined) {
  if (particles === undefined) stepParticles.delete(set);
  else stepParticles.set(set, particles);
}

/**
 * the particle budget of a step where a view of `tiles` follows the camera in steps (`follow`, no camera given), else 0:
 * its budget, less where its steps go by time (`stepMs`, `setParticleTilesStep`)
 */
function getStepBudget(tiles: TilesValues, set: ParticleTileSet) {
  if (tiles.name !== 'camera' || !tiles.params.follow || tiles.params.camera) return 0;
  const budget = Math.max(0, tiles.params.budget ?? 0);
  const byTime = (tiles.params.stepMs ?? 0) > 0 ? stepParticles.get(set) : undefined;
  if (byTime === undefined) return budget;
  return budget > 0 ? Math.max(1, Math.min(budget, byTime)) : Math.max(1, byTime);
}

/**
 * A following view's next tiles (`getStepBudget`): those of `shown` still picked, then the picked ones ready (by
 * priority) up to its budget of particles (at least one tile, or one region switched), not waiting for the rest;
 * `shown` as it is where that is none (a pick all new, none of it ready). Why: an empty view a blank frame, and what
 * draws it made again for it and again after. A region switched (`getSwitchGroups`: a level switch, a node for its
 * tiles or the tiles for their node) steps whole: its picked tiles in once all have arrived (ready or failed), the shown
 * ones it replaces out in the same step, kept until then. Why: no region empty while it loads, none drawn twice (a
 * node's sample over its tiles: more particles to draw, and the two looks at once)
 */
function stepTiles(
  set: ParticleTileSet,
  loader: ParticleTileLoader,
  picked: readonly number[],
  frame: number,
  shown: ArrayLike<number>,
  budget: number,
) {
  const wanted = new Set(picked);
  const next: number[] = [];
  const kept = new Set<number>();
  for (let i = 0; i < shown.length; ++i) {
    if (!wanted.has(shown[i])) continue;
    next.push(shown[i]);
    kept.add(shown[i]);
  }
  const groups = getSwitchGroups(set, picked, shown, kept);
  const ready = (id: number) => loader.getState(id, frame) === 'ready';
  const arrived = (id: number) => {
    const state = loader.getState(id, frame);
    return state === 'ready' || state === 'failed';
  };
  const visited = new Set<SwitchGroup>();
  const taken = new Set<SwitchGroup>();
  let particles = 0;
  for (const id of picked) {
    if (kept.has(id)) continue;
    const group = groups?.get(id);
    if (!group) {
      if (!ready(id)) continue;
      const count = set.tiles.get(id)!.count;
      if (particles > 0 && particles + count > budget) break;
      next.push(id);
      kept.add(id);
      particles += count;
      continue;
    }
    // a region switched: at its first picked tile (by priority), all of them or none
    if (visited.has(group)) continue;
    visited.add(group);
    if (!group.picked.every(arrived)) continue;
    const ids = group.picked.filter(ready);
    let count = 0;
    for (const t of ids) count += set.tiles.get(t)!.count;
    if (particles > 0 && particles + count > budget) break;
    for (const t of ids) {
      next.push(t);
      kept.add(t);
    }
    taken.add(group);
    particles += count;
  }
  if (next.length === 0) return Array.from(shown);
  // the regions not switched yet keep what they show
  if (groups) {
    for (const group of new Set(groups.values())) {
      if (!taken.has(group)) for (const t of group.shown) next.push(t);
    }
  }
  return next;
}

export { stepTiles as stepParticleTiles };

/** A region switched (`getSwitchGroups`): the shown tiles leaving it, the picked tiles over it */
interface SwitchGroup {
  readonly shown: number[];
  readonly picked: number[];
}

/**
 * The regions a step switches (with levels): the tiles of `shown` no longer picked and the picked tiles not `kept`
 * over them, joined where their regions share a base cell (`getTileRegions`; a node's its block of them), per picked
 * tile its group; none without levels (its tiles' regions apart) or where none leaves. A shown tile no picked tile is
 * over (the camera turned away) in none: dropped. By the base grid's cells, not tile against tile: thousands of each
 * at a switch
 */
function getSwitchGroups(
  set: ParticleTileSet,
  picked: readonly number[],
  shown: ArrayLike<number>,
  kept: ReadonlySet<number>,
) {
  if (!set.levels) return undefined;
  const regions = getTileRegions(set);
  const wanted = new Set(picked);
  const leaving: number[] = [];
  for (let i = 0; i < shown.length; ++i) {
    if (!wanted.has(shown[i]) && regions.has(shown[i])) leaving.push(shown[i]);
  }
  if (leaving.length === 0) return undefined;
  const d = set.levels.base.grid.dimensions[0];
  // per base cell the leaving tile over it (its index), -1 none
  const owner = new Int32Array(d * d * d).fill(-1);
  for (let i = 0; i < leaving.length; ++i) {
    someRegionCell(d, regions.get(leaving[i])!, (cell) => {
      owner[cell] = i;
      return false;
    });
  }
  // the leaving tiles joined by the picked tiles over them (union-find)
  const parent = Int32Array.from({ length: leaving.length }, (_, i) => i);
  const root = (i: number) => {
    while (parent[i] !== i) i = parent[i] = parent[parent[i]];
    return i;
  };
  const over = new Map<number, number>();
  for (const id of picked) {
    if (kept.has(id)) continue;
    const region = regions.get(id);
    if (region === undefined) continue;
    let first = -1;
    someRegionCell(d, region, (cell) => {
      const o = owner[cell];
      if (o < 0) return false;
      if (first < 0) first = o;
      else parent[root(o)] = root(first);
      return false;
    });
    if (first >= 0) over.set(id, first);
  }
  if (over.size === 0) return undefined;
  const byRoot = new Map<number, SwitchGroup>();
  const groups = new Map<number, SwitchGroup>();
  for (const [id, o] of over) {
    const r = root(o);
    let group = byRoot.get(r);
    if (!group) byRoot.set(r, (group = { shown: [], picked: [] }));
    group.picked.push(id);
    groups.set(id, group);
  }
  for (let i = 0; i < leaving.length; ++i) byRoot.get(root(i))?.shown.push(leaving[i]);
  return groups;
}

/** per set with levels its tiles' regions (`getTileRegions`) */
const tileRegions = new WeakMap<ParticleTileSet, Map<number, number>>();

/**
 * The region of each tile of `set` (with levels): its level and its cell in that level's grid, as one number
 * (`cell * 32 + level`); a base tile level 0 at its own cell, a node by its level's grid (`ParticleTileLevels.getNode`)
 */
function getTileRegions(set: ParticleTileSet) {
  let regions = tileRegions.get(set);
  if (regions) return regions;
  regions = new Map();
  const levels = set.levels!;
  for (const id of levels.base.tiles.keys()) regions.set(id, id * 32);
  const dimensions = levels.base.grid.dimensions[0];
  for (let level = 1; level <= levels.count; ++level) {
    const d = dimensions >> level;
    for (let cell = 0, n = d * d * d; cell < n; ++cell) {
      const node = levels.getNode(level, cell);
      if (node !== undefined) regions.set(node, cell * 32 + level);
    }
  }
  tileRegions.set(set, regions);
  return regions;
}

/**
 * Whether `test` holds for a base cell (of a grid of `dimensions` cells an axis) within `region` (`getTileRegions`: a
 * level's cell, `2^level` base cells an axis), stopping at the first that does
 */
function someRegionCell(dimensions: number, region: number, test: (cell: number) => boolean) {
  const level = region % 32;
  const cell = (region - level) / 32;
  const d = dimensions >> level;
  const x0 = (cell % d) << level;
  const y0 = (Math.floor(cell / d) % d) << level;
  const z0 = Math.floor(cell / (d * d)) << level;
  const n = 1 << level;
  for (let z = z0; z < z0 + n; ++z) {
    for (let y = y0; y < y0 + n; ++y) {
      for (let x = x0; x < x0 + n; ++x) {
        if (test(x + y * dimensions + z * dimensions * dimensions)) return true;
      }
    }
  }
  return false;
}

/**
 * resolves where one of `ids` (`any`), else all, is ready or failed, or the loader is idle. Why not idle alone: the
 * loader wants the tiles of every view of its set (`wantTiles`), a view not waiting for the others'
 */
function whenSettled(loader: ParticleTileLoader, ids: readonly number[], frame: number, any: boolean) {
  const isSettled = (id: number) => {
    const state = loader.getState(id, frame);
    return state === 'ready' || state === 'failed';
  };
  const settled = () => (any ? ids.some(isSettled) : ids.every(isSettled));
  if (ids.length === 0 || settled()) return Promise.resolve();
  return new Promise<void>((resolve) => {
    let done = false;
    let unsubscribe = () => {};
    const finish = () => {
      if (done) return;
      done = true;
      unsubscribe();
      resolve();
    };
    unsubscribe = loader.subscribe(() => {
      if (settled()) finish();
    });
    loader.whenIdle().then(finish);
  });
}

/** A view's cache (its cell's, `StateTransformer.ApplyParams.cache`): the key of its wants on its set's loader */
interface ViewCache {
  loader?: ParticleTileLoader;
}

/** what a view of a loader's set wants (its picked tiles at a frame, `pickTiles`' array) and shows */
interface ViewWant {
  ids: readonly number[];
  frame: number;
  shown: ArrayLike<number>;
}

/**
 * per loader its views' wants, by view. Why: the loader is the tile set's, shared by its views; a view setting what it
 * wants alone cancels another's loads, and letting go of what it neither wants nor shows (`trimLoader`) drops another's
 * tiles
 */
const viewWants = new WeakMap<ParticleTileLoader, Map<ViewCache, ViewWant>>();

/** `loader` wanting what its views want: their picks by priority, the views' in turn (none waiting on another's) */
function setLoaderWanted(loader: ParticleTileLoader, views: ReadonlyMap<ViewCache, ViewWant>) {
  const wants = [...views.values()];
  const requests: { id: number; frame: number }[] = [];
  const longest = wants.reduce((n, w) => Math.max(n, w.ids.length), 0);
  for (let i = 0; i < longest; ++i) {
    for (const w of wants) if (i < w.ids.length) requests.push({ id: w.ids[i], frame: w.frame });
  }
  loader.setWanted(requests);
}

/** view `view` of `loader`'s set wanting `ids` at `frame` (by priority), where not the pick it wants already */
function wantTiles(loader: ParticleTileLoader, view: ViewCache, ids: readonly number[], frame: number) {
  // a view moved to another set: its wants on the last one's loader dropped
  if (view.loader && view.loader !== loader) forgetView(view);
  view.loader = loader;
  let views = viewWants.get(loader);
  if (!views) {
    views = new Map();
    viewWants.set(loader, views);
  }
  const last = views.get(view);
  if (last?.ids === ids && last.frame === frame) return;
  views.set(view, { ids, frame, shown: last?.shown ?? [] });
  setLoaderWanted(loader, views);
}

/** view `view`'s wants dropped from its loader (the view disposed, or of another set), the other views' kept */
function forgetView(view: ViewCache) {
  const loader = view.loader;
  if (!loader) return;
  view.loader = undefined;
  const views = viewWants.get(loader);
  if (!views?.delete(view)) return;
  setLoaderWanted(loader, views);
}

/**
 * Whether a following view (`getStepBudget`) is behind what the camera picks now: a step (`stepTiles`, any budget)
 * changing what it shows (what `ParticleTilesFollowCamera` steps it toward). Not behind while what a step would add is
 * not there yet: a picked tile not ready, a region switched not all arrived, a view none of whose tiles is picked and
 * none of the picked ready (the loader's change waking the behavior); not stepped to the same tiles every frame
 * meanwhile (forever where they fail)
 */
export function isParticleTilesViewBehind(
  data: SO.Particle.TileSet['data'],
  params: { tiles: TilesValues; frame: number },
  list: ParticleList,
  plugin: PluginContext,
  /** the view's cell's cache (`StateObjectCell.cache`): its wants on the loader the other views' kept */
  cache: unknown,
  /** the set's latest pick, not the camera's now (`pickTiles`): a pass the camera's change did not wake */
  latest = false,
) {
  const { set, loader } = data;
  const view = ParticleTileView.of(list);
  if (!view || view.set !== set) return true;
  const frame = Math.max(0, Math.min(params.frame, set.frameCount - 1));
  const picked = pickTiles(set, params.tiles, plugin, latest);
  // the loader wanting what the camera picks now, whether the view steps or not: else a pick that keeps the tiles shown
  // (the first, before the camera framed the set) and brings others never loads them
  if (cache) wantTiles(loader, cache as ViewCache, picked, frame);
  // the step itself, unbounded: behind where it changes the tiles (a bounded step changes them too, at least a tile or
  // a region)
  const next = stepTiles(set, loader, picked, frame, view.tiles, Infinity);
  const shown = new Set(view.tiles);
  if (next.length !== shown.size) return true;
  for (const id of next) if (!shown.has(id)) return true;
  return false;
}

export { ParticleTilesView };
type ParticleTilesView = typeof ParticleTilesView;
const ParticleTilesView = PluginStateTransform.BuiltIn({
  name: 'particle-tiles-view',
  display: {
    name: 'Particle Tiles View',
    description:
      'Some tiles of a tile set at a frame as one particle list: all, by id, in a box or seen by the camera.',
  },
  from: SO.Particle.TileSet,
  to: SO.Particle.List,
  params: (a) => ({
    tiles: TilesParam,
    frame: PD.Numeric(
      0,
      { min: 0, max: Math.max(0, (a?.data.set.frameCount ?? 1) - 1), step: 1 },
      { description: 'Zero-based frame index', immediateUpdate: true },
    ),
    slots: PD.Boolean(true, {
      description:
        "A slot list: a tile's particles a page at its own slots, a tile added or let go updates what draws it (the target representation on WebGPU), not all of it.",
    }),
  }),
})({
  apply({ a, params, cache }, plugin: PluginContext) {
    return Task.create('Particle Tiles View', async () => {
      const view = await createView(a.data, params, plugin, cache as ViewCache);
      return new SO.Particle.List(view.list, describe(view));
    });
  },
  update({ a, b, newParams, cache }, plugin: PluginContext) {
    return Task.create('Particle Tiles View', async () => {
      // the previous view: a frame step of the same tiles the same particles posed otherwise
      const view = await createView(a.data, newParams, plugin, cache as ViewCache, ParticleTileView.of(b.data), b.data);
      b.data = view.list;
      const { label, description } = describe(view);
      b.label = label;
      b.description = description;
      return StateTransformer.UpdateResult.Updated;
    });
  },
  dispose({ cache }) {
    if (cache) forgetView(cache as ViewCache);
  },
});

async function createView(
  data: SO.Particle.TileSet['data'],
  params: { tiles: PD.Values<{ tiles: typeof TilesParam }>['tiles']; frame: number; slots?: boolean },
  plugin: PluginContext,
  cache: ViewCache,
  previous?: ParticleTileView,
  previousList?: ParticleList,
) {
  const { set, loader } = data;
  const frame = Math.max(0, Math.min(params.frame, set.frameCount - 1));
  // a following view's step toward the set's latest pick: the behavior's pass before it picked for the camera where it
  // woke for it (`isParticleTilesViewBehind`), not again for a camera moved since
  const following = params.tiles.name === 'camera' && params.tiles.params.follow && !params.tiles.params.camera;
  const ids = pickTiles(set, params.tiles, plugin, following);
  wantTiles(loader, cache, ids, frame);
  let ready: number[];
  const budget = getStepBudget(params.tiles, set);
  if (budget > 0) {
    // a step toward the tiles picked: those ready up to the budget, the others as they load (the next steps)
    const shown = previous?.set === set && previous.frame === frame ? previous.tiles : [];
    if (shown.length === 0) await whenSettled(loader, ids, frame, true);
    ready = UserTiming.span('tiles: step', () => stepTiles(set, loader, ids, frame, shown, budget));
  } else {
    await whenSettled(loader, ids, frame, false);
    // the failed left out: a view of what is there
    ready = ids.filter((id) => loader.getState(id, frame) === 'ready');
  }
  // a slot list the next version of the last (its changes by page) where that is one of the same set
  const view = UserTiming.span('tiles: view', () =>
    params.slots
      ? ParticleTileView.createSlots(set, loader, ready, frame, previousList)
      : ParticleTileView.create(set, loader, ready, frame, previous?.set === set ? previous : undefined),
  );
  const want = viewWants.get(loader)?.get(cache);
  if (want) want.shown = view.tiles;
  const max = params.tiles.name === 'camera' ? (params.tiles.params.maxParticles ?? 0) : 0;
  if (max > 0) UserTiming.span('tiles: trim', () => trimLoader(loader, frame, LoadedParticlesPerCap * max));
  return view;
}

/** the particles a loader keeps at most, by the cap (`trimLoader`): what the views show, what they want, and some more */
const LoadedParticlesPerCap = 2;

/**
 * `loader`'s ready tiles at `frame` that no view of its set wants or shows let go where its particles are past
 * `budget`. Why: a view following the camera loads the tiles it passes, kept otherwise (the loader keeps what is ready).
 * Its targets kept: one read again a new object, the representation's members all made again
 */
function trimLoader(loader: ParticleTileLoader, frame: number, budget: number) {
  if (loader.getLoaded().particles <= budget) return;
  const keep = new Set<number>();
  for (const w of viewWants.get(loader)?.values() ?? []) {
    // a view at another frame wants its tiles there, none of this frame's
    if (w.frame !== frame) continue;
    for (const id of w.ids) keep.add(id);
    for (let i = 0; i < w.shown.length; ++i) keep.add(w.shown[i]);
  }
  const release = loader.getReady(frame).filter((r) => !keep.has(r.id));
  if (release.length === 0) return;
  loader.release(release);
}

function describe(view: ParticleTileView) {
  const { set } = view;
  const frames = set.frameCount > 1 ? ` · frame ${view.frame + 1} of ${set.frameCount}` : '';
  return {
    label: set.label || 'Particle Tiles',
    description: `${view.tiles.length} of ${set.tiles.size} tiles · ${ParticleSlots.getLiveCount(view.list)} particles${frames}`,
  };
}
