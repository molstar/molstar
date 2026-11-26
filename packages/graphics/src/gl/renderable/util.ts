/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Sphere3D } from '@molstar/core/math/geometry';
import { Vec3, Mat4 } from '@molstar/core/math/linear-algebra';
import { BoundaryHelper } from '@molstar/core/math/geometry/boundary-helper';
import { InstanceTransforms } from '@molstar/core/math/geometry/instance-transforms';
import type { InstanceCenters } from '@molstar/core/math/geometry/instance-grid';
import type { TextureFilter } from '../texture.js';
import { arrayMinMax } from '@molstar/core/util/array';
import { SegmentedArray } from '../segmented-array.js';
import { getObjectTransforms, type TransformValues } from '../transform-parts.js';
import { ValueCell } from '@molstar/core/util/value-cell';
import type { ComputeContext, StorageBuffer } from '../compute.js';

// avoiding namespace lookup improved performance in Chrome (Aug 2020)
const v3fromArray = Vec3.fromArray;

export function calculateTextureInfo(n: number, itemSize: number) {
  n = Math.max(n, 2); // observed issues with 1 pixel textures
  const sqN = Math.sqrt(n);
  let width = Math.ceil(sqN);
  width = width + ((itemSize - (width % itemSize)) % itemSize);
  const height = width > 0 ? Math.ceil(n / width) : 0;
  return { width, height, length: width * height * itemSize };
}

export interface TextureImage<T extends Uint8Array | Float32Array | Int32Array> {
  readonly array: T;
  readonly width: number;
  readonly height: number;
  readonly flipY?: boolean;
  readonly filter?: TextureFilter;
  /**
   * Incremented by whoever writes `array` in place (the font atlas appending a glyph): a shared texture
   * (`TextureSpecOptions.shared`) uploads the image again only when it changed
   */
  version?: number;
}

export interface TextureVolume<T extends Uint8Array | Float32Array> {
  readonly array: T;
  readonly width: number;
  readonly height: number;
  readonly depth: number;
}

export function createTextureImage<T extends Uint8Array | Float32Array>(
  n: number,
  itemSize: number,
  arrayCtor: new (length: number) => T,
  array?: T,
): TextureImage<T> {
  const { length, width, height } = calculateTextureInfo(n, itemSize);
  array = array && array.length >= length ? array : new arrayCtor(length);
  return { array, width, height };
}

/**
 * The value of a `DataSpec`: `count` elements by index in `array`, which may be longer (reused by producers, rounded
 * up to whole words); no 2D layout, backends lay it out (WebGL data textures, WebGPU storage pages)
 */
export interface DataArray<T extends Uint8Array | Float32Array | Int16Array> {
  readonly array: T;
  readonly count: number;
  /**
   * Its elements written on the GPU (a smoothed layer of a merged surface sampled per vertex, `computeLayerSampling`):
   * `array` then zeros (`getSharedZeros`), the render item writing the elements by this where it draws them
   */
  readonly gpu?: GpuDataWriter;
  /**
   * an array of its own its producer may reuse for the next data of the same cell (`createDataArray`): with `gpu`, the
   * input it was sampled from (a smoothed layer's by group and instance), `array` the shared zeros. Why: else every set
   * of the layer allocates it again (groups times instances, a segment of 60 copies about 0.5 MB)
   */
  readonly spare?: T;
}

/** What writes a data array's elements on the GPU (`DataArray.gpu`) */
export interface GpuDataWriter {
  /** the context it writes with: of another one, its data array foreign there (`hasForeignResources`) */
  readonly compute: ComputeContext;
  /** bytes per element where written (4 for `rgba8`, 1 for `r8`) */
  readonly elementBytes: number;
  /** Writes the elements into `buffer` (a storage page) from `byteOffset`, the bytes around them kept */
  write(buffer: StorageBuffer, byteOffset: number): void;
  /** The elements as last written; zeros where never written */
  read(): Promise<Uint8Array>;
  /** Releases its GPU resources: the data array replaced or dropped */
  dispose(): void;
}

/** the buffers of `getSharedZeros`, by element size */
const SharedZeroBuffers = new Map<number, ArrayBuffer>();
const SharedZeroBufferSet = new WeakSet<ArrayBuffer>();

/**
 * `length` zeros of `ctor`, a view of a buffer shared by all such arrays (the arrays of data written on the GPU,
 * `DataArray.gpu`): never written, and not reused by `createDataArray`
 */
export function getSharedZeros<T extends Uint8Array | Float32Array | Int16Array>(
  ctor: new (buffer: ArrayBuffer, byteOffset: number, length: number) => T,
  length: number,
): T {
  const bpe = (ctor as unknown as { BYTES_PER_ELEMENT: number }).BYTES_PER_ELEMENT;
  const words = Math.ceil((length * bpe) / 4) * 4;
  let buffer = SharedZeroBuffers.get(bpe);
  if (!buffer || buffer.byteLength < words) {
    // a new one where longer needed: the views of the old one keep it
    buffer = new ArrayBuffer(Math.max(words, (buffer?.byteLength ?? 0) * 2));
    SharedZeroBuffers.set(bpe, buffer);
    SharedZeroBufferSet.add(buffer);
  }
  return new ctor(buffer, 0, words / bpe);
}

export function isSharedZeros(array: ArrayBufferView) {
  return SharedZeroBufferSet.has(array.buffer as ArrayBuffer);
}

/** Releases what writes `data` on the GPU (`DataArray.gpu`), where it is replaced */
export function releaseGpuData(data: DataArray<any> | undefined) {
  data?.gpu?.dispose();
}

/**
 * `values` with the data arrays written on the GPU (`DataArray.gpu`) read back, for CPU readers of them (the geometry
 * export): a view of `values` with those cells replaced, `values` itself unchanged
 */
export async function readGpuDataArrays<V extends object>(values: V): Promise<V> {
  let view: { [k: string]: unknown } | undefined;
  for (const k in values) {
    const cell = (values as { [k: string]: unknown })[k] as ValueCell<DataArray<any>> | undefined;
    const data = cell && typeof cell === 'object' && 'ref' in cell ? cell.ref.value : undefined;
    if (!data?.gpu) continue;
    const bytes = await data.gpu.read();
    view ??= Object.create(values) as { [k: string]: unknown };
    view[k] = ValueCell.create({ array: bytes, count: data.count });
  }
  return (view ?? values) as V;
}

/**
 * `count` elements of `itemSize` components, in `array` where long enough (its contents kept), else a new one of whole
 * words
 */
export function createDataArray<T extends Uint8Array | Float32Array | Int16Array>(
  count: number,
  itemSize: number,
  arrayCtor: new (length: number) => T,
  array?: T,
): DataArray<T> {
  const length = count * itemSize;
  // shared zeros (`getSharedZeros`) never: they are no producer's to write
  if (!array || array.length < length || isSharedZeros(array)) {
    const bpe = (arrayCtor as unknown as { BYTES_PER_ELEMENT: number }).BYTES_PER_ELEMENT;
    array = new arrayCtor(Math.ceil((length * bpe) / 4) * (4 / bpe));
  }
  return { array, count };
}

/**
 * `data` laid out as a 2D image of `itemSize` components per texel (`calculateTextureInfo`), for consumers that load it
 * into a texture of their own (e.g. compute passes); a copy, padded
 */
export function getDataTextureImage<T extends Uint8Array | Float32Array>(
  data: DataArray<T>,
  itemSize: number,
): TextureImage<T> {
  const { width, height, length } = calculateTextureInfo(data.count, itemSize);
  const source = data.array as T | SegmentedArray;
  const segmented = SegmentedArray.is(source);
  const ctor = (segmented ? source.ctor : source.constructor) as new (length: number) => T;
  const array = new ctor(length);
  const n = Math.min(length, data.count * itemSize);
  array.set(segmented ? (SegmentedArray.read(source, 0, n) as T) : source.subarray(0, n));
  return { array, width, height };
}

const DefaultPrintImageOptions = {
  scale: 1,
  pixelated: false,
  id: 'molstar.debug.image',
  normalize: false,
  useCanvas: false,
  flipY: false,
};
export type PrintImageOptions = typeof DefaultPrintImageOptions;

export function printTextureImage(textureImage: TextureImage<any>, options: Partial<PrintImageOptions> = {}) {
  const { array, width, height } = textureImage;
  const itemSize = array.length / (width * height);
  const data = new Uint8ClampedArray(width * height * 4);
  const [min, max] = arrayMinMax(array);
  if (itemSize === 1) {
    data.fill(255);
    for (let y = 0; y < height; ++y) {
      for (let x = 0; x < width; ++x) {
        const i = y * width + x;
        if (options.normalize) {
          data[i * 4 + 0] = ((array[i] - min) / (max - min)) * 255;
        } else {
          data[i * 4 + 0] = array[i] * 255;
        }
      }
    }
  } else if (itemSize === 4) {
    if (options.normalize) {
      for (let i = 0, il = width * height * 4; i < il; i += 4) {
        data[i] = ((array[i] - min) / (max - min)) * 255;
        data[i + 1] = ((array[i + 1] - min) / (max - min)) * 255;
        data[i + 2] = ((array[i + 2] - min) / (max - min)) * 255;
        data[i + 3] = 255;
      }
    } else {
      data.set(array);
    }
  } else if (itemSize === 3) {
    for (let i = 0, il = width * height; i < il; ++i) {
      data[i * 4] = array[i * 3];
      data[i * 4 + 1] = array[i * 3 + 1];
      data[i * 4 + 2] = array[i * 3 + 2];
      data[i * 4 + 3] = 255;
    }
  } else {
    console.warn(`itemSize '${itemSize}' not supported`);
  }
  return printImageData(new ImageData(data, width, height), options);
}

let tmpCanvas: HTMLCanvasElement;
let tmpCanvasCtx: CanvasRenderingContext2D;
let tmpContainer: HTMLDivElement;

export function printImageData(imageData: ImageData, options: Partial<PrintImageOptions> = {}) {
  const o = { ...DefaultPrintImageOptions, ...options };

  if (!tmpContainer) {
    tmpContainer = document.createElement('div');
    tmpContainer.style.position = 'absolute';
    tmpContainer.style.bottom = '0px';
    tmpContainer.style.right = '0px';
    tmpContainer.style.border = 'solid orange';
    tmpContainer.style.pointerEvents = 'none';
    if (options.flipY) tmpContainer.style.transform = 'scaleY(-1)';
    document.body.appendChild(tmpContainer);
  }

  if (o.useCanvas) {
    const existingCanvas = document.getElementById(o.id) as HTMLCanvasElement;
    const outCanvas = existingCanvas || document.createElement('canvas');
    outCanvas.width = imageData.width;
    outCanvas.height = imageData.height;
    const outCtx = outCanvas.getContext('2d');
    if (!outCtx) throw new Error('Could not create canvas 2d context');
    outCtx.putImageData(imageData, 0, 0);
    outCanvas.id = o.id;
    outCanvas.style.width = imageData.width * o.scale + 'px';
    outCanvas.style.height = imageData.height * o.scale + 'px';
    if (o.pixelated) {
      outCanvas.style.imageRendering = 'pixelated';
    }
    outCanvas.style.position = 'relative';
    outCanvas.style.border = 'solid grey';
    outCanvas.style.pointerEvents = 'none';
    if (!existingCanvas) tmpContainer.appendChild(outCanvas);
  } else {
    const canvas = tmpCanvas || document.createElement('canvas');
    tmpCanvas = canvas;
    canvas.width = imageData.width;
    canvas.height = imageData.height;
    const ctx = tmpCanvasCtx || canvas.getContext('2d');
    tmpCanvasCtx = ctx;
    if (!ctx) throw new Error('Could not create canvas 2d context');
    ctx.putImageData(imageData, 0, 0);

    canvas.toBlob((imgBlob) => {
      const objectURL = URL.createObjectURL(imgBlob!);
      const existingImg = document.getElementById(o.id) as HTMLImageElement;
      const img = existingImg || document.createElement('img');
      img.id = o.id;
      img.src = objectURL;
      img.style.width = imageData.width * o.scale + 'px';
      img.style.height = imageData.height * o.scale + 'px';
      if (o.pixelated) {
        // not supported in Firefox and IE
        img.style.imageRendering = 'pixelated';
      }
      img.style.position = 'relative';
      img.style.border = 'solid grey';
      img.style.pointerEvents = 'none';
      if (!existingImg) tmpContainer.appendChild(img);
    }, 'image/png');
  }
}

//

const v = Vec3();
const boundaryHelperCoarse = new BoundaryHelper('14');
const boundaryHelperFine = new BoundaryHelper('98');

function getHelper(count: number) {
  return count > 100_000 ? boundaryHelperCoarse : boundaryHelperFine;
}

/** `stride`: floats per position (3, or 4 where a group follows each, as `tPositionGroup`) */
export function calculateInvariantBoundingSphere(
  position: Float32Array,
  positionCount: number,
  stepFactor: number,
  stride = 3,
): Sphere3D {
  const step = stepFactor * stride;
  const boundaryHelper = getHelper(positionCount);

  boundaryHelper.reset();
  for (let i = 0, _i = positionCount * stride; i < _i; i += step) {
    v3fromArray(v, position, i);
    boundaryHelper.includePosition(v);
  }
  boundaryHelper.finishedIncludeStep();
  for (let i = 0, _i = positionCount * stride; i < _i; i += step) {
    v3fromArray(v, position, i);
    boundaryHelper.radiusPosition(v);
  }

  const sphere = boundaryHelper.getSphere();

  // no more positions than the fine helper's extrema (49 directions, two each): the positions themselves, never more
  // extrema than positions
  if (positionCount <= 98) {
    const extrema: Vec3[] = [];
    for (let i = 0, _i = positionCount * stride; i < _i; i += step) {
      extrema.push(v3fromArray(Vec3(), position, i));
    }
    Sphere3D.setExtrema(sphere, extrema);
  }

  return sphere;
}

const _mat4 = Mat4();

/**
 * The sphere around the box of the transformed centers with the largest distance from its center, plus the invariant
 * radius: encloses, two passes over the transforms, many times cheaper than the helper's. For many instances as tight
 * or tighter: the helper's center is the extrema's centroid, pulled off where they pile up unevenly (a crystal lattice
 * 8% wider, 3PQR; a capsid 1%, 1M4X)
 */
function calculateBoxTransformBoundingSphere(invariantBoundingSphere: Sphere3D, t: InstanceTransforms): Sphere3D {
  const { center, radius } = invariantBoundingSphere;
  InstanceTransforms.writePoints(t, center, undefined, undefined, undefined, boxMin, boxMax);
  const c = Vec3.create((boxMin[0] + boxMax[0]) / 2, (boxMin[1] + boxMax[1]) / 2, (boxMin[2] + boxMax[2]) / 2);
  return Sphere3D.create(c, Math.sqrt(InstanceTransforms.getMaxDistanceSq(t, center, c)) + radius);
}
const boxMin = Vec3();
const boxMax = Vec3();

/** `calculateBoxTransformBoundingSphere`'s center and largest distance per transforms' cell, by their versions */
const boxSpheres = new WeakMap<
  object,
  { transforms: number; poses: number; count: number; invariant: Vec3; center: Vec3; distance: number }
>();

/**
 * Seeds `calculateTransformDataBoundingSphere` with the instances' centers' sphere of `data` as it is (its grid's
 * pass over them, `InstanceCenters`): the first bounding sphere then without a pass of its own (many instances only)
 */
export function setTransformDataCenters(data: TransformValues, centers: InstanceCenters) {
  const count = data.instanceCount.ref.value;
  if (count < BoxTransformBoundingSphereMinCount && !data.instanceSphereFast?.ref.value) return;
  boxSpheres.set(data.aTransform, {
    transforms: data.aTransform.ref.version,
    poses: data.instancePoses?.ref.version ?? -1,
    count,
    invariant: Vec3.clone(centers.invariantCenter),
    center: Vec3.clone(centers.center),
    distance: centers.distance,
  });
}

/**
 * `calculateInstancesBoundingSphere` of the instances of `data` (its poses or matrices, `getObjectTransforms`), the
 * many instances' box sphere's center and largest distance kept per transforms' version and `invariantBoundingSphere`'s
 * center: a radius changed alone (a size's padding) adds to the distance, no pass over the instances again
 */
export function calculateTransformDataBoundingSphere(
  invariantBoundingSphere: Sphere3D,
  data: TransformValues,
): Sphere3D {
  const t = getObjectTransforms(data);
  // a few instances' the boundary helper's, unless the object asks for the box sphere (`instanceSphereFast`)
  if (t.count < BoxTransformBoundingSphereMinCount && !data.instanceSphereFast?.ref.value)
    return getFewInstancesBoundingSphere(invariantBoundingSphere, data, t);
  const { center, radius } = invariantBoundingSphere;
  const transforms = data.aTransform.ref.version;
  const poses = data.instancePoses?.ref.version ?? -1;
  const e = boxSpheres.get(data.aTransform);
  if (
    e &&
    e.transforms === transforms &&
    e.poses === poses &&
    e.count === t.count &&
    Vec3.exactEquals(e.invariant, center)
  ) {
    return Sphere3D.create(Vec3.clone(e.center), e.distance + radius);
  }
  InstanceTransforms.writePoints(t, center, undefined, undefined, undefined, boxMin, boxMax);
  const c = Vec3.create((boxMin[0] + boxMax[0]) / 2, (boxMin[1] + boxMax[1]) / 2, (boxMin[2] + boxMax[2]) / 2);
  const distance = Math.sqrt(InstanceTransforms.getMaxDistanceSq(t, center, c));
  boxSpheres.set(data.aTransform, {
    transforms,
    poses,
    count: t.count,
    invariant: Vec3.clone(center),
    center: Vec3.clone(c),
    distance,
  });
  return Sphere3D.create(c, distance + radius);
}

/** `calculateInstancesBoundingSphere`'s sphere of few instances, per transforms' cell by its versions and the invariant sphere */
const fewSpheres = new WeakMap<
  object,
  { transforms: number; poses: number; count: number; invariant: Sphere3D; sphere: Sphere3D }
>();

/**
 * The sphere of a few instances (`calculateInstancesBoundingSphere`: the boundary helper's, tighter than the box sphere,
 * not additive in the radius) kept per transforms' version and the whole invariant sphere: an object's values and its
 * bounds' update computing it for the same transforms (a pack's targets: thousands of small ones) read it once
 */
/** the same center, radius and extrema (a few instances' sphere reads them) */
function isSameSphere(a: Sphere3D, b: Sphere3D) {
  if (!Sphere3D.exactEquals(a, b)) return false;
  const ea = a.extrema;
  const eb = b.extrema;
  if (!ea || !eb) return ea === eb;
  return ea.length === eb.length && ea.every((v, i) => Vec3.exactEquals(v, eb[i]));
}

function getFewInstancesBoundingSphere(invariant: Sphere3D, data: TransformValues, t: InstanceTransforms): Sphere3D {
  const transforms = data.aTransform.ref.version;
  const poses = data.instancePoses?.ref.version ?? -1;
  const e = fewSpheres.get(data.aTransform);
  if (
    e &&
    e.transforms === transforms &&
    e.poses === poses &&
    e.count === t.count &&
    isSameSphere(e.invariant, invariant)
  ) {
    return Sphere3D.clone(e.sphere);
  }
  const sphere = calculateInstancesBoundingSphere(invariant, t);
  fewSpheres.set(data.aTransform, {
    transforms,
    poses,
    count: t.count,
    invariant: Sphere3D.clone(invariant),
    sphere: Sphere3D.clone(sphere),
  });
  return sphere;
}

/** Sets whether the instances' bounding sphere of `data` is the box sphere at any count (`instanceSphereFast`) */
export function setInstanceSphereFast(data: { instanceSphereFast?: ValueCell<boolean> }, fast: boolean) {
  if (data.instanceSphereFast) ValueCell.updateIfChanged(data.instanceSphereFast, fast);
  else data.instanceSphereFast = ValueCell.create(fast);
}

/** The values `getDrawnBoundingSphere` reads */
export interface DrawnBoundsValues {
  readonly boundingSphere: ValueCell<Sphere3D>;
  /** the identity where missing */
  readonly uObjectTransform?: ValueCell<Mat4>;
}

const drawnBounds = new WeakMap<DrawnBoundsValues, { bounds: number; transform: number; sphere: Sphere3D }>();

/**
 * The bounds of an object as drawn: `boundingSphere` (its instances as placed) with its transient transform
 * (`uObjectTransform`) applied, `boundingSphere` itself at rest. What culls and bounds the scene while drawing reads
 * this; what fits the camera to the data reads `boundingSphere`
 */
export function getDrawnBoundingSphere(values: DrawnBoundsValues): Sphere3D {
  const cell = values.uObjectTransform;
  if (!cell || Mat4.isIdentity(cell.ref.value)) return values.boundingSphere.ref.value;
  const transform = cell.ref.value;
  let e = drawnBounds.get(values);
  const bounds = values.boundingSphere.ref.version;
  const version = cell.ref.version;
  if (e && e.bounds === bounds && e.transform === version) return e.sphere;
  if (!e) {
    e = { bounds, transform: version, sphere: Sphere3D() };
    drawnBounds.set(values, e);
  }
  e.bounds = bounds;
  e.transform = version;
  // center and radius only: the extrema transformed are allocations per call, and every union of them (the scene's
  // bounds, a merged object's) costs per point; a transient view's bounds need not be tight
  const s = values.boundingSphere.ref.value;
  Vec3.transformMat4(e.sphere.center, s.center, transform);
  e.sphere.radius = s.radius * Mat4.getMaxScaleOnAxis(transform);
  return e.sphere;
}

/**
 * Sphere `i` of `spheres` (center and radius, four floats each: an instance grid's cells and batches) into `out`,
 * transformed by `transform` whose largest scale is `scale` (`Mat4.getMaxScaleOnAxis`)
 */
export function transformSphereAt(out: Sphere3D, spheres: Float32Array, i: number, transform: Mat4, scale: number) {
  const o = i * 4;
  const x = spheres[o],
    y = spheres[o + 1],
    z = spheres[o + 2];
  out.center[0] = transform[0] * x + transform[4] * y + transform[8] * z + transform[12];
  out.center[1] = transform[1] * x + transform[5] * y + transform[9] * z + transform[13];
  out.center[2] = transform[2] * x + transform[6] * y + transform[10] * z + transform[14];
  out.radius = spheres[o + 3] * scale;
  return out;
}

/**
 * From how many instances the sphere is the box's (`calculateBoxTransformBoundingSphere`): few instances spread wide (a
 * dozen particles over a cell's surface) put the box's center far off and the sphere 1.5x too wide (MG_1184's scene);
 * the helper's is cheap for them
 */
const BoxTransformBoundingSphereMinCount = 1000;

export function calculateTransformBoundingSphere(
  invariantBoundingSphere: Sphere3D,
  transform: Float32Array,
  transformCount: number,
  transformOffset: number,
): Sphere3D {
  return calculateInstancesBoundingSphere(
    invariantBoundingSphere,
    InstanceTransforms.ofMatrices(transform, transformCount, transformOffset),
  );
}

/** The sphere around the instances `t` (matrices or poses, `InstanceTransforms`) of `invariantBoundingSphere` */
export function calculateInstancesBoundingSphere(invariantBoundingSphere: Sphere3D, t: InstanceTransforms): Sphere3D {
  const transformCount = t.count;
  if (transformCount === 1) {
    InstanceTransforms.getMatrix(t, 0, _mat4);
    const s = Sphere3D.clone(invariantBoundingSphere);
    return Mat4.isIdentity(_mat4) ? s : Sphere3D.transform(s, s, _mat4);
  }
  if (transformCount >= BoxTransformBoundingSphereMinCount) {
    return calculateBoxTransformBoundingSphere(invariantBoundingSphere, t);
  }

  const boundaryHelper = getHelper(transformCount);
  boundaryHelper.reset();

  const { center, radius, extrema } = invariantBoundingSphere;
  // the balls around the extrema (`Sphere3D.extremaRadius`): instance transforms are rigid
  const extremaRadius = Sphere3D.getExtremaRadius(invariantBoundingSphere);

  // only use extrema if there are not too many transforms
  if (extrema && transformCount <= 14) {
    for (let i = 0, _i = transformCount; i < _i; ++i) {
      for (const e of extrema) {
        InstanceTransforms.getPoint(t, i, e, v);
        boundaryHelper.includePositionRadius(v, extremaRadius);
      }
    }
    boundaryHelper.finishedIncludeStep();
    for (let i = 0, _i = transformCount; i < _i; ++i) {
      for (const e of extrema) {
        InstanceTransforms.getPoint(t, i, e, v);
        boundaryHelper.radiusPositionRadius(v, extremaRadius);
      }
    }
    // with balls, few points: the transformed points themselves (the helper's extrema would be points on their balls,
    // as many as its directions)
    if (extremaRadius > 0 && extrema.length * transformCount <= 98) {
      const points: Vec3[] = [];
      for (let i = 0, _i = transformCount; i < _i; ++i) {
        for (const e of extrema) {
          const p = Vec3();
          InstanceTransforms.getPoint(t, i, e, p);
          points.push(p);
        }
      }
      return Sphere3D.setExtrema(boundaryHelper.getSphere(), points, extremaRadius);
    }
  } else {
    for (let i = 0, _i = transformCount; i < _i; ++i) {
      InstanceTransforms.getPoint(t, i, center, v);
      boundaryHelper.includePositionRadius(v, radius);
    }
    boundaryHelper.finishedIncludeStep();
    for (let i = 0, _i = transformCount; i < _i; ++i) {
      InstanceTransforms.getPoint(t, i, center, v);
      boundaryHelper.radiusPositionRadius(v, radius);
    }
  }

  return boundaryHelper.getSphere();
}

export function calculateBoundingSphere(
  position: Float32Array,
  positionCount: number,
  transform: Float32Array,
  transformCount: number,
  padding = 0,
  stepFactor = 1,
): { boundingSphere: Sphere3D; invariantBoundingSphere: Sphere3D } {
  const invariantBoundingSphere = calculateInvariantBoundingSphere(position, positionCount, stepFactor);
  const boundingSphere = calculateTransformBoundingSphere(invariantBoundingSphere, transform, transformCount, 0);
  Sphere3D.expand(boundingSphere, boundingSphere, padding);
  Sphere3D.expand(invariantBoundingSphere, invariantBoundingSphere, padding);
  return { boundingSphere, invariantBoundingSphere };
}
