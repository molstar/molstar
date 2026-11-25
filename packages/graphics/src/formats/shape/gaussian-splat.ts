/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { type RuntimeContext, Task } from '@molstar/core/task';
import { Color } from '@molstar/core/util/color';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Mat4 } from '@molstar/core/math/linear-algebra/3d/mat4';
import { Vec3 } from '@molstar/core/math/linear-algebra/3d/vec3';
import type { ReadonlyUint8Array } from '@molstar/core/util/type-helpers';
import type { PlyFile } from '@molstar/io/reader/ply/schema';
import { parseSplat } from '@molstar/io/reader/splat/parser';
import { parseSpz } from '@molstar/io/reader/spz/parser';
import { decodeImageBitmap } from '@molstar/io/common/image';
import { decodeSogImage, parseSog, parseSogArchive } from '@molstar/io/reader/sog/parser';
import type { SogImage, SogImageBytes, SogMeta } from '@molstar/io/reader/sog/schema';
import { Sphere3D } from '@molstar/core/math/geometry/primitives/sphere3d';
import { isSplatSortSupported } from '@molstar/graphics/gl/compute/splat-sort';
import type { Shape } from '@molstar/model/model/shape';
import {
  GaussianSplat,
  type GaussianSplatData,
  type GaussianSplatReadOptions,
} from '@molstar/model/formats/shape/gaussian-splat/data';
import { gaussianSplatsFromPly } from '@molstar/model/formats/shape/gaussian-splat/ply';
import { gaussianSplatsFromSplat } from '@molstar/model/formats/shape/gaussian-splat/splat';
import { gaussianSplatsFromSpz } from '@molstar/model/formats/shape/gaussian-splat/spz';
import {
  createSogMeansReader,
  gaussianSplatsFromSog,
  getSogShLayout,
} from '@molstar/model/formats/shape/gaussian-splat/sog';
import * as ShapeUtils from '@molstar/graphics/geo/shape/shape';
import type { ShapeProvider } from '@molstar/graphics/geo/shape/provider';
import { GaussianSplats } from '@molstar/graphics/geo/geometry/gaussian-splats/gaussian-splats';
import { encodeSplatSh, NoSplatSh, setSplat } from '@molstar/graphics/geo/geometry/gaussian-splats/splat-data';
import { int32ToFloatBits, PositionGroupStride } from '@molstar/graphics/geo/geometry/center-group';
import type { GraphicsContext } from '@molstar/graphics/gl/context';

/** The params of a gaussian splat file's shape: how it is read (the turn, the view-dependent color), the geometry's */
export const GaussianSplatFileParams = {
  ...GaussianSplats.Params,
  shDegree: PD.Numeric(
    0,
    { min: 0, max: 3, step: 1 },
    {
      description:
        "The degree of view-dependent color read (spherical harmonics, at most the file's): 0 the base color only; 3 " +
        'adds 60 bytes a splat.',
    },
  ),
  upright: PD.Boolean(true, {
    description:
      'Turns the scene 180° about x: 3DGS scenes are y-down (the camera convention of their training), upside down as stored.',
  }),
};
export type GaussianSplatFileParams = typeof GaussianSplatFileParams;

/**
 * The splats and their base colors (rgb bytes, a group each; none where the shader decodes them: SOG read in place), in
 * the file's frame (`GaussianSplatData.upright`)
 */
export interface GaussianSplatShapeData {
  splats: GaussianSplats;
  colors?: Uint8Array;
  upright: boolean;
}

/** The turn of 180° about x, between the 3DGS PLY's y-down frame and right-up-back */
const UprightTurn: Mat4 = Mat4.fromRotation(Mat4(), Math.PI, Vec3.unitX);

/**
 * The shape's transforms of splats stored `stored` upright for the file's `transforms` (none: one instance), turned 180°
 * about x first where `upright` differs (the turn about x as the data's own: the view-dependent color, evaluated in the
 * splats' frame, turns with it)
 */
export function getSplatTransforms(
  transforms: Mat4[] | undefined,
  stored: boolean,
  upright: boolean,
): Mat4[] | undefined {
  if (stored === upright) return transforms;
  if (!transforms) return [UprightTurn];
  return transforms.map((t) => Mat4.mul(Mat4(), t, UprightTurn));
}

function toByte(v: number) {
  return Math.round(Math.min(1, Math.max(0, v)) * 255);
}

const UpdateChunk = 100_000;

/**
 * The geometry of `data`, each splat decoded once into the layouts the shaders read (no other copy), the geometry of
 * `splats` reused where given
 */
export async function buildGaussianSplats(
  ctx: RuntimeContext,
  data: GaussianSplatData,
  splats?: GaussianSplats,
): Promise<GaussianSplatShapeData> {
  const { count } = data;
  const { positionGroup, splatData } = GaussianSplats.createArrays(count);
  const colors = new Uint8Array(count * 3);
  const s = GaussianSplat();
  for (let i = 0; i < count; ++i) {
    data.getSplat(i, s);
    positionGroup[i * 4] = s.x;
    positionGroup[i * 4 + 1] = s.y;
    positionGroup[i * 4 + 2] = s.z;
    positionGroup[i * 4 + 3] = int32ToFloatBits(i);
    setSplat(splatData, i, s.logScaleX, s.logScaleY, s.logScaleZ, s.qx, s.qy, s.qz, s.qw, s.opacity);
    colors[i * 3] = toByte(s.r);
    colors[i * 3 + 1] = toByte(s.g);
    colors[i * 3 + 2] = toByte(s.b);
    if (i % UpdateChunk === 0 && ctx.shouldUpdate) {
      await ctx.update({ message: 'Building gaussian splats', current: i, max: count });
    }
  }
  const sh = data.sh ? encodeSplatSh(count, data.sh.degree, data.sh.get) : NoSplatSh;
  return { splats: GaussianSplats.create(positionGroup, splatData, count, splats, sh), colors, upright: data.upright };
}

/** The shape of splats read from `source` (`name` its kind), each splat a group of its base color */
export function createSplatShape(
  name: string,
  source: unknown,
  transforms: Mat4[] | undefined,
  data: GaussianSplatShapeData,
): Shape<GaussianSplats> {
  const { splats, colors } = data;
  return ShapeUtils.create(
    name,
    source,
    splats,
    // without colors the shader's (`direct`, see `ShapeUtils.getTheme`)
    colors
      ? (group: number) => Color.fromRgb(colors[group * 3], colors[group * 3 + 1], colors[group * 3 + 2])
      : () => Color(0xffffff),
    () => 1,
    (group: number) => `Splat ${group}`,
    transforms,
    splats.splatCount,
  );
}

/** A gaussian splat file: its parsed data (PLY) or bytes (`.splat`, `.spz`, `.sog`) */
export type GaussianSplatSource =
  | { readonly kind: 'ply'; readonly file: PlyFile }
  | { readonly kind: 'splat' | 'spz' | 'sog'; readonly bytes: ReadonlyUint8Array };

export type GaussianSplatShapeSource = GaussianSplatSource & { readonly transforms?: Mat4[] };

/**
 * The splats of `source` read with `options` (`.splat` and `.spz` parsed first: their arrays as views; `.sog` its images
 * decoded by the browser)
 */
export async function readGaussianSplats(
  ctx: RuntimeContext,
  source: GaussianSplatSource,
  options: GaussianSplatReadOptions,
): Promise<GaussianSplatData> {
  if (source.kind === 'ply') return gaussianSplatsFromPly(source.file, options);
  if (source.kind === 'splat') {
    const parsed = await parseSplat(source.bytes).runInContext(ctx);
    if (parsed.isError) throw new Error(parsed.message);
    return gaussianSplatsFromSplat(parsed.result, options);
  }
  if (source.kind === 'sog') {
    const parsed = await parseSog(source.bytes).runInContext(ctx);
    if (parsed.isError) throw new Error(parsed.message);
    return gaussianSplatsFromSog(parsed.result, options);
  }
  const parsed = await parseSpz(source.bytes).runInContext(ctx);
  if (parsed.isError) throw new Error(parsed.message);
  return gaussianSplatsFromSpz(parsed.result, options);
}

/** Whether SOG's images can be uploaded as they are (decoded by the browser into `ImageBitmap`s) */
function canReadSogOnGpu() {
  return typeof createImageBitmap !== 'undefined';
}

/** A lossless WebP of one pixel, rgba (200, 100, 50, 128): its rgb as stored where nothing premultiplies it */
const PremultiplyProbe = new Uint8Array([
  82, 73, 70, 70, 30, 0, 0, 0, 87, 69, 66, 80, 86, 80, 56, 76, 17, 0, 0, 0, 47, 0, 0, 0, 16, 7, 80, 178, 34, 87, 166,
  128, 129, 136, 232, 127, 0, 0,
]);

const UploadProbes = new WeakMap<GraphicsContext, Promise<boolean>>();

/**
 * Whether images decoded by the browser upload to `graphics` as stored: some browsers premultiply them despite
 * `premultiplyAlpha: 'none'` (WebKit under WebGPU, the SOG spec warns), which corrupts `sh0`'s codebook indices below
 * alpha 255. The probe image through the same path, read back; once per context
 */
function canUploadImagesAsStored(graphics: GraphicsContext) {
  let probe = UploadProbes.get(graphics);
  if (!probe) {
    probe = (async () => {
      const bitmap = await decodeBitmap({ name: 'probe', bytes: PremultiplyProbe });
      const texture = graphics.resources.texture('image-uint8', 'rgba', 'ubyte', 'nearest');
      try {
        texture.load(bitmap);
        const { array } = await graphics.readTextureData(texture);
        return [200, 100, 50, 128].every((v, i) => Math.abs(array[i] - v) <= 1);
      } finally {
        bitmap.close();
        texture.destroy();
      }
    })().catch(() => false);
    UploadProbes.set(graphics, probe);
  }
  return probe;
}

let warnedPremultiplied = false;

/** The browser's decode of a SOG image's lossless WebP as stored: not premultiplied, no color conversion (`sh0`'s alpha) */
function decodeBitmap(image: SogImageBytes) {
  return decodeImageBitmap(image.bytes, 'image/webp', image.name);
}

/** The centers of `count` splats decoded from SOG's means (`meansLow`, `meansHigh`: RGBA bytes) */
function decodeSogCenters(meta: SogMeta, count: number, meansLow: Uint8Array, meansHigh: Uint8Array) {
  const positionGroup = new Float32Array(count * PositionGroupStride);
  const position = createSogMeansReader(meta, meansLow, meansHigh);
  for (let i = 0; i < count; ++i) {
    positionGroup[i * 4] = position(i, 0);
    positionGroup[i * 4 + 1] = position(i, 1);
    positionGroup[i * 4 + 2] = position(i, 2);
    positionGroup[i * 4 + 3] = int32ToFloatBits(i);
  }
  return positionGroup;
}

/** The bounds of SOG's centers from the means' ranges (`meta.json`, in the log domain) */
function getSogCentersBounds(meta: SogMeta) {
  const out = (n: number) => Math.sign(n) * (Math.exp(Math.abs(n)) - 1);
  const { mins, maxs } = meta.means;
  const min = Vec3.create(out(mins[0]), out(mins[1]), out(mins[2]));
  const max = Vec3.create(out(maxs[0]), out(maxs[1]), out(maxs[2]));
  return Sphere3D.create(Vec3.scale(Vec3(), Vec3.add(Vec3(), min, max), 0.5), Vec3.distance(min, max) / 2);
}

/**
 * The splats of the SOG file `bytes` read in place on `graphics` (the `sog` encoding): the images uploaded as the browser
 * decodes them, none read back but the means for the CPU's centers (the sorts, the bounds): before the first frame
 * where the CPU sorts (WebGL), after it where the GPU does (it decodes its own centers: `decodeSogCenters`); the
 * geometry of `splats` reused where given. None where `uploadsAsStored` (checked once decoded, not uploaded before)
 * resolves false: the CPU's path then
 */
export async function buildSogSplats(
  ctx: RuntimeContext,
  bytes: ReadonlyUint8Array,
  graphics: GraphicsContext,
  shDegree: number,
  splats?: GaussianSplats,
  uploadsAsStored?: Promise<boolean>,
): Promise<GaussianSplatShapeData | undefined> {
  const parsed = await parseSogArchive(bytes).runInContext(ctx);
  if (parsed.isError) throw new Error(parsed.message);
  const archive = parsed.result;
  const { meta, count } = archive;
  const shN = shDegree > 0 ? archive.shN : undefined;
  // the CPU's centers later where the GPU sorts: the means uploaded as the others
  const later = !!graphics.compute && isSplatSortSupported(graphics.compute);
  const [meansLow, meansHigh, scales, quats, sh0, labels, centroids] = await Promise.all([
    later ? decodeBitmap(archive.meansLow) : decodeSogImage(archive, archive.meansLow),
    later ? decodeBitmap(archive.meansHigh) : decodeSogImage(archive, archive.meansHigh),
    decodeBitmap(archive.scales),
    decodeBitmap(archive.quats),
    decodeBitmap(archive.sh0),
    shN ? decodeBitmap(shN.labels) : undefined,
    shN ? decodeBitmap(shN.centroids) : undefined,
  ]);
  if (uploadsAsStored && !(await uploadsAsStored)) {
    for (const image of [meansLow, meansHigh, scales, quats, sh0, labels, centroids]) {
      if (image instanceof ImageBitmap) image.close();
    }
    return;
  }
  if (scales.width * scales.height < count)
    throw new Error(`SOG images of ${scales.width}x${scales.height} for ${count} splats`);

  const upload = (image: ImageBitmap | SogImage) => {
    const texture = graphics.resources.texture('image-uint8', 'rgba', 'ubyte', 'nearest');
    if (image instanceof ImageBitmap) {
      texture.load(image);
      // the GPU's copy kept, the browser's decode not (a restore decodes the file's bytes again)
      image.close();
    } else {
      texture.load({ array: image.data as Uint8Array, width: image.width, height: image.height });
    }
    return texture;
  };
  const codebook = new Float32Array(256 * 4);
  const scaleCodes = meta.scales.codebook;
  const colorCodes = meta.sh0.codebook;
  const shCodes = meta.shN?.codebook ?? [];
  let maxLogScale = -Infinity;
  for (let i = 0; i < 256; ++i) {
    codebook[i * 4] = scaleCodes[i] ?? 0;
    codebook[i * 4 + 1] = colorCodes[i] ?? 0;
    codebook[i * 4 + 2] = shCodes[i] ?? 0;
    if (i < scaleCodes.length) maxLogScale = Math.max(maxLogScale, scaleCodes[i]);
  }
  const layout = shN && centroids ? getSogShLayout(meta, centroids.width, centroids.height) : undefined;
  const { mins, maxs } = meta.means;
  const meansLowTexture = upload(meansLow);
  const geometry = GaussianSplats.createSog(
    {
      graphics,
      count,
      positionGroup:
        meansLow instanceof ImageBitmap
          ? undefined
          : decodeSogCenters(meta, count, meansLow.data, (meansHigh as SogImage).data),
      centersBounds: getSogCentersBounds(meta),
      textures: {
        meansLow: meansLowTexture,
        meansHigh: upload(meansHigh),
        scales: upload(scales),
        quats: upload(quats),
        sh0: upload(sh0),
        shLabels: labels && upload(labels),
        shCentroids: centroids && upload(centroids),
      },
      meansMin: Vec3.create(mins[0], mins[1], mins[2]),
      meansMax: Vec3.create(maxs[0], maxs[1], maxs[2]),
      codebook,
      version: archive.version,
      sh: layout && layout.bands > 0 ? { degree: Math.min(shDegree, layout.bands), ...layout } : undefined,
      maxScale: Number.isFinite(maxLogScale) ? Math.exp(maxLogScale) : 0,
    },
    splats,
  );
  if (later) {
    // after the first frame: the GPU sorts and draws without them
    setTimeout(async () => {
      try {
        const [lo, hi] = await Promise.all([
          decodeSogImage(archive, archive.meansLow),
          decodeSogImage(archive, archive.meansHigh),
        ]);
        // the same images still (not replaced, not released)
        if (geometry.sog.meansLow.ref.value !== meansLowTexture) return;
        GaussianSplats.setSogCenters(geometry, decodeSogCenters(meta, count, lo.data, hi.data));
      } catch (e) {
        console.error('Reading the centers of SOG failed', e);
      }
    }, 0);
  }
  return { splats: geometry, upright: false };
}

function makeShapeGetter() {
  let _source: GaussianSplatShapeSource | undefined;
  let _upright: boolean | undefined;
  let _shDegree: number | undefined;
  let _graphics: GraphicsContext | undefined;
  let _built: GaussianSplatShapeData;
  let _shape: Shape<GaussianSplats>;
  return async (
    ctx: RuntimeContext,
    source: GaussianSplatShapeSource,
    props: PD.Values<GaussianSplatFileParams>,
    shape?: Shape<GaussianSplats>,
    graphics?: GraphicsContext,
  ) => {
    // the splats depend on the file and the view-dependent color's degree only (the other params are the geometry's,
    // applied by the representation), SOG read in place on the context too (its textures); the turn is the shape's
    // transform, a new shape of the same geometry
    const gpu = source.kind === 'sog' && !!graphics && canReadSogOnGpu();
    // SOG read in place is the representation's own (its textures released with it, `dispose`): another
    // representation of the provider (no shape, or another's) builds its own; callers reading the shape alone (no
    // context: e.g. its center) take the cached one
    const owner = !!shape && shape.geometry === _built?.splats;
    const rebuild = _source !== source || _shDegree !== props.shDegree || (gpu && (_graphics !== graphics || !owner));
    if (rebuild) {
      // the upload probed while the file is read (`canUploadImagesAsStored`)
      const sog =
        gpu && source.kind === 'sog'
          ? await buildSogSplats(
              ctx,
              source.bytes,
              graphics!,
              props.shDegree,
              shape?.geometry,
              canUploadImagesAsStored(graphics!),
            )
          : undefined;
      if (gpu && !sog && !warnedPremultiplied) {
        warnedPremultiplied = true;
        console.warn('SOG read on the CPU: this browser premultiplies the images it uploads');
      }
      if (sog) {
        _built = sog;
      } else {
        const data = await readGaussianSplats(ctx, source, { shDegree: props.shDegree });
        _built = await buildGaussianSplats(ctx, data, shape?.geometry);
      }
      _graphics = graphics;
    }
    if (rebuild || _upright !== props.upright) {
      _shape = createSplatShape(
        `${source.kind}-gaussian-splats`,
        source.kind === 'ply' ? source.file : source.bytes,
        getSplatTransforms(source.transforms, _built.upright, props.upright),
        _built,
      );
      _source = source;
      _upright = props.upright;
      _shDegree = props.shDegree;
    }
    return _shape;
  };
}

/** The shape provider of a gaussian splat file (`source`) */
export function shapeFromGaussianSplats(source: GaussianSplatShapeSource) {
  return Task.create<ShapeProvider<GaussianSplatShapeSource, GaussianSplats, GaussianSplatFileParams>>(
    'Shape Provider',
    async () => ({
      label: 'Gaussian Splats',
      data: source,
      params: GaussianSplatFileParams,
      getShape: makeShapeGetter(),
      geometryUtils: GaussianSplats.Utils,
      dispose: GaussianSplats.dispose,
    }),
  );
}
