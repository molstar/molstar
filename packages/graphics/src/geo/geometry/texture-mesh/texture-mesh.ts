/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Cai Huiyu <szmun.caihy@gmail.com>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import { ValueCell } from '@molstar/core/util';
import { Box3D, Sphere3D } from '@molstar/core/math/geometry';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { LocationIterator, PositionLocation } from '../../util/location-iterator.js';
import type { TransformData } from '../transform-data.js';
import { createColors } from '../color-data.js';
import { createMarkers } from '../marker-data.js';
import type { GeometryUtils } from '../geometry.js';
import type { Theme } from '@molstar/graphics/theme/theme';
import type { Color } from '@molstar/core/util/color';
import { BaseGeometry, resolveInstanceGranularity } from '../base.js';
import { createEmptyOverpaint } from '../overpaint-data.js';
import { createEmptyTransparency } from '../transparency-data.js';
import type { TextureMeshValues } from '@molstar/graphics/gl/renderable/texture-mesh';
import { calculateTransformBoundingSphere } from '@molstar/graphics/gl/renderable/util';
import { createNullTexture, isNullTexture, type Texture } from '@molstar/graphics/gl/texture';
import { type Mat4, Vec2, Vec3, Vec4 } from '@molstar/core/math/linear-algebra';
import { createEmptyClipping } from '../clipping-data.js';
import { NullLocation } from '@molstar/model/model/location';
import { createEmptySubstance } from '../substance-data.js';
import type { RenderableState } from '@molstar/graphics/gl/renderable';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import type { GeometryRequest } from '../request.js';
import type { ComputeContext, DrawIndirect, StorageBuffer } from '@molstar/graphics/gl/compute';
import {
  type ChunkedSurface,
  destroyChunkedSurface,
  releaseRetiredTables,
  onVertexArenaGrow,
  readChunkedVertices,
} from '@molstar/graphics/gl/compute/marching-cubes/chunks';
import { createEmptyEmissive } from '../emissive-data.js';
import { createEmptyWiggle } from '../wiggle-data.js';
import { createInteriorValues, getInteriorParam, updateInteriorValues } from '../interior.js';
import { getAnimationParam, createAnimationValues, updateAnimationValues } from '../animation.js';

export interface TextureMesh {
  readonly kind: 'texture-mesh';

  /** Number of vertices in the texture-mesh, an upper bound with `drawIndirect` */
  vertexCount: number;
  /** Number of groups in the texture-mesh */
  groupCount: number;

  readonly geoTextureDim: ValueCell<Vec2>;
  readonly vertexTexture: ValueCell<Texture>;
  readonly groupTexture: ValueCell<Texture>;
  readonly normalTexture: ValueCell<Texture>;
  readonly varyingGroup: ValueCell<boolean>;
  readonly doubleBuffer: TextureMesh.DoubleBuffer;
  /** `vertexCount`, shared with the render values (`drawCount`, `uVertexCount`) */
  readonly vertexCountCell: ValueCell<number>;
  /** The vertex count in a GPU buffer (WebGPU), the vertices beyond it up to `vertexCount` are not drawn */
  readonly drawIndirect: ValueCell<DrawIndirect | undefined>;
  /**
   * Whether the vertices are in pages of the vertex arena (WebGPU, see `computeChunkedIsosurface`) instead of the
   * textures: the page table of the surface and the arena's buffers
   */
  readonly isoChunks: ValueCell<boolean>;
  readonly isoTable: ValueCell<StorageBuffer | undefined>;
  readonly isoVertices: ValueCell<StorageBuffer | undefined>;

  readonly boundingSphere: Sphere3D;

  readonly meta: {
    /** The vertices and normals read back for the position iterator (`TextureMesh.Utils.prepare`) */
    vertexData?: { vertices: Float32Array; normals: Float32Array };
    /** Called after the vertices were extracted again (data derived from them, e.g. smoothing), by kind */
    reextracted?: Map<string, () => void>;
    /** The exact vertex count of indirectly drawn meshes, once read back (`vertexCount` is an upper bound) */
    readVertexCount?: number;
    /** Bytes of the pages and page table of a mesh in the vertex arena, once read back */
    chunkByteCount?: number;
    /** Settles once the vertices are final (indirectly drawn meshes: after the count check and any re-extraction) */
    extracted?: Promise<void>;
    /**
     * A box holding every untransformed vertex (e.g. the grid the surface was extracted from), which lets layer
     * smoothing leave out the instances a layer does not reach (`getTextureMeshSplatInstances`); unknown: all drawn
     */
    vertexBox?: Box3D;
    [k: string]: unknown;
  };
}

export namespace TextureMesh {
  export class DoubleBuffer {
    private index = 0;
    private textures: ({ vertex: Texture; group: Texture; normal: Texture } | undefined)[] = [];
    /** the surface in the vertex arena of a compute context (WebGPU), instead of the textures */
    private chunks: { surface: ChunkedSurface; compute: ComputeContext; unsubscribe: () => void } | undefined;

    /** The surface in the vertex arena of `graphics` to reuse; one of another context is dropped, not destroyed */
    getChunks(graphics: GraphicsContext) {
      if (this.chunks && this.chunks.compute !== graphics.compute) {
        this.chunks.unsubscribe();
        this.chunks = undefined;
      }
      return this.chunks?.surface;
    }

    /** `surface` in an arena of `compute`, `onGrow` called with the current surface when an arena grew (its buffers replaced) */
    setChunks(surface: ChunkedSurface, compute: ComputeContext, onGrow: (surface: ChunkedSurface) => void) {
      if (this.chunks?.compute !== compute) {
        this.chunks?.unsubscribe();
        const unsubscribe = onVertexArenaGrow(compute, () => {
          if (this.chunks) onGrow(this.chunks.surface);
        });
        this.chunks = { surface, compute, unsubscribe };
      } else {
        this.chunks.surface = surface;
      }
    }

    /**
     * The textures to reuse (WebGL); textures of another graphics context than `graphics` (e.g. before a backend
     * switch, or a restore that replaced the context) are not reused, nor destroyed (their context does it)
     */
    get(graphics: GraphicsContext) {
      if (this.textures.some((t) => t && !graphics.resources.ownsTexture(t.vertex))) {
        this.textures = [];
        this.index = 0;
      }
      return this.textures[this.index];
    }

    /** Forgets the textures (the mesh filled into the vertex arena instead): of the previous context, not destroyed */
    dropTextures() {
      this.textures = [];
      this.index = 0;
    }

    set(vertex: Texture, group: Texture, normal: Texture) {
      this.textures[this.index] = Object.assign(this.textures[this.index] || {}, {
        vertex,
        group,
        normal,
      });
      this.index = (this.index + 1) % 2;
    }

    destroy() {
      for (const buffer of this.textures) {
        buffer!.vertex.destroy();
        buffer!.group.destroy();
        buffer!.normal.destroy();
      }
      if (this.chunks) {
        // the pages back to the arena; the surface's draw arguments are destroyed with it
        destroyChunkedSurface(this.chunks.compute, this.chunks.surface);
        this.chunks.unsubscribe();
        this.chunks = undefined;
      }
    }
  }

  export function create(
    vertexCount: number,
    groupCount: number,
    vertexTexture: Texture,
    groupTexture: Texture,
    normalTexture: Texture,
    boundingSphere: Sphere3D,
    textureMesh?: TextureMesh,
    drawIndirect?: DrawIndirect,
  ): TextureMesh {
    const width = vertexTexture.getWidth();
    const height = vertexTexture.getHeight();
    if (textureMesh) {
      textureMesh.vertexCount = vertexCount;
      textureMesh.groupCount = groupCount;
      textureMesh.meta.vertexData = undefined;
      textureMesh.meta.extracted = undefined;
      textureMesh.meta.vertexBox = undefined;
      ValueCell.updateIfChanged(textureMesh.vertexCountCell, vertexCount);
      ValueCell.updateIfChanged(textureMesh.drawIndirect, drawIndirect);
      ValueCell.update(textureMesh.geoTextureDim, Vec2.set(textureMesh.geoTextureDim.ref.value, width, height));
      ValueCell.update(textureMesh.vertexTexture, vertexTexture);
      ValueCell.update(textureMesh.groupTexture, groupTexture);
      ValueCell.update(textureMesh.normalTexture, normalTexture);
      textureMesh.doubleBuffer.set(vertexTexture, groupTexture, normalTexture);
      // last filled into the vertex arena (before a switch to WebGL, or emptied): the textures now; a surface kept in
      // the arena is reused (`getChunks`) or destroyed with the mesh
      if (textureMesh.isoChunks.ref.value) {
        ValueCell.update(textureMesh.isoChunks, false);
        ValueCell.update(textureMesh.isoTable, undefined);
        ValueCell.update(textureMesh.isoVertices, undefined);
      }
      Sphere3D.copy(textureMesh.boundingSphere, boundingSphere);
      return textureMesh;
    } else {
      return {
        kind: 'texture-mesh',
        vertexCount,
        groupCount,
        geoTextureDim: ValueCell.create(Vec2.create(width, height)),
        vertexTexture: ValueCell.create(vertexTexture),
        groupTexture: ValueCell.create(groupTexture),
        normalTexture: ValueCell.create(normalTexture),
        varyingGroup: ValueCell.create(false),
        doubleBuffer: new DoubleBuffer(),
        vertexCountCell: ValueCell.create(vertexCount),
        drawIndirect: ValueCell.create(drawIndirect),
        isoChunks: ValueCell.create(false),
        isoTable: ValueCell.create(undefined),
        isoVertices: ValueCell.create(undefined),
        boundingSphere: Sphere3D.clone(boundingSphere),
        meta: {},
      };
    }
  }

  /**
   * The texture mesh of a surface in the vertex arena of `compute` (see `computeChunkedIsosurface`), drawn indirectly;
   * `vertexCount` is an upper bound (for draws without the GPU count, e.g. multi-draw lists)
   */
  export function createChunked(
    vertexCount: number,
    groupCount: number,
    surface: ChunkedSurface,
    compute: ComputeContext,
    boundingSphere: Sphere3D,
    textureMesh?: TextureMesh,
  ): TextureMesh {
    const drawIndirect = { buffer: surface.drawArgs, byteOffset: 0 };
    let mesh: TextureMesh;
    if (textureMesh) {
      mesh = textureMesh;
      mesh.vertexCount = vertexCount;
      mesh.groupCount = groupCount;
      mesh.meta.vertexData = undefined;
      mesh.meta.extracted = undefined;
      mesh.meta.vertexBox = undefined;
      ValueCell.updateIfChanged(mesh.vertexCountCell, vertexCount);
      ValueCell.updateIfChanged(mesh.drawIndirect, drawIndirect);
      Sphere3D.copy(mesh.boundingSphere, boundingSphere);
      // last filled into textures (before a switch from WebGL): of the previous context, not destroyed (it does it)
      if (!isNullTexture(mesh.vertexTexture.ref.value)) {
        const nullTexture = createNullTexture();
        ValueCell.update(mesh.vertexTexture, nullTexture);
        ValueCell.update(mesh.groupTexture, nullTexture);
        ValueCell.update(mesh.normalTexture, nullTexture);
        mesh.doubleBuffer.dropTextures();
      }
    } else {
      const nullTexture = createNullTexture();
      mesh = create(
        vertexCount,
        groupCount,
        nullTexture,
        nullTexture,
        nullTexture,
        boundingSphere,
        undefined,
        drawIndirect,
      );
    }
    bindChunkedTable(mesh, surface);
    ValueCell.updateIfChanged(mesh.isoChunks, true);
    // an arena's buffers replaced when it grew
    mesh.doubleBuffer.setChunks(surface, compute, (s) =>
      ValueCell.updateIfChanged(mesh.isoVertices, s.arena.vertices.ref.value),
    );
    return mesh;
  }

  /**
   * Binds the current page table of `surface` (grown for a larger count, see `computeChunkedIsosurface`) and the
   * vertices of its arena (another after a move), and releases the tables it grew out of
   */
  export function bindChunkedTable(mesh: TextureMesh, surface: ChunkedSurface) {
    ValueCell.updateIfChanged(mesh.isoTable, surface.table);
    ValueCell.updateIfChanged(mesh.isoVertices, surface.arena.vertices.ref.value);
    releaseRetiredTables(surface);
  }

  /**
   * The box of the points of a grid of `gridDim` points placed by `transform`, which holds every vertex of an isosurface
   * extracted from it (marching cubes puts them on the edges between grid points), for `meta.vertexBox`
   */
  export function getGridVertexBox(gridDim: Vec3, transform: Mat4) {
    const box = Box3D.create(Vec3.create(0, 0, 0), Vec3.create(gridDim[0] - 1, gridDim[1] - 1, gridDim[2] - 1));
    return Box3D.transform(box, box, transform);
  }

  /** Sets the handler of `kind` called after the vertices were extracted again */
  export function setReextractionHandler(meta: TextureMesh['meta'], kind: string, handler: () => void) {
    if (!meta.reextracted) meta.reextracted = new Map();
    meta.reextracted.set(kind, handler);
  }

  /**
   * The vertices read back where the position iterator is used (see `createPositionIterator`), once they are final (an
   * indirectly drawn mesh may be extracted again with the exact count)
   */
  function prepare(textureMesh: TextureMesh, request: GeometryRequest, graphics: GraphicsContext | undefined) {
    if (!request.positions || !graphics || textureMesh.meta.vertexData) return;
    const { extracted } = textureMesh.meta;
    const read = () => readVertexData(textureMesh, graphics);
    return extracted ? extracted.then(read) : read();
  }

  /** Reads the vertices and normals back into `meta.vertexData`, for the position iterator */
  async function readVertexData(textureMesh: TextureMesh, graphics: GraphicsContext) {
    const table = textureMesh.isoTable.ref.value;
    const arena = textureMesh.isoVertices.ref.value;
    if (textureMesh.isoChunks.ref.value && table && arena && graphics.compute) {
      const count = textureMesh.meta.readVertexCount ?? textureMesh.vertexCount;
      const { vertices, normals } = await readChunkedVertices(graphics.compute, { table, vertices: arena }, count);
      textureMesh.meta.vertexData = { vertices, normals };
      return;
    }
    const [vertices, normals] = await Promise.all([
      graphics.readTextureData(textureMesh.vertexTexture.ref.value),
      graphics.readTextureData(textureMesh.normalTexture.ref.value),
    ]);
    textureMesh.meta.vertexData = {
      vertices: vertices.array as Float32Array,
      normals: normals.array as Float32Array,
    };
  }

  /**
   * The vertex count of a mesh of the vertex arena, read back after its extraction, into its draw arguments: drawn
   * directly with it from then on, until the next extraction draws indirectly again (`createChunked`). An indirect draw
   * costs far more GPU time than a direct one in Chrome; with hundreds of surfaces, several passes a frame, they took
   * most of it (3j3q, 1356 surfaces: 105–540 ms a frame against 25 ms drawn directly). `vertexCount` stays the upper
   * bound (the stride of vertex-instance data).
   */
  export function setDrawnVertexCount(textureMesh: TextureMesh, count: number) {
    const args = textureMesh.drawIndirect.ref.value;
    if (args && args.count !== count) ValueCell.update(textureMesh.drawIndirect, { ...args, count });
  }

  /** Updates the vertex count (an upper bound with `drawIndirect`) after the textures were resized */
  export function setVertexCount(textureMesh: TextureMesh, vertexCount: number) {
    const texture = textureMesh.vertexTexture.ref.value;
    textureMesh.vertexCount = vertexCount;
    textureMesh.meta.vertexData = undefined;
    ValueCell.updateIfChanged(textureMesh.vertexCountCell, vertexCount);
    ValueCell.update(
      textureMesh.geoTextureDim,
      Vec2.set(textureMesh.geoTextureDim.ref.value, texture.getWidth(), texture.getHeight()),
    );
  }

  export function createEmpty(textureMesh?: TextureMesh): TextureMesh {
    const vt = textureMesh ? textureMesh.vertexTexture.ref.value : createNullTexture();
    const gt = textureMesh ? textureMesh.groupTexture.ref.value : createNullTexture();
    const nt = textureMesh ? textureMesh.normalTexture.ref.value : createNullTexture();
    const bs = textureMesh ? textureMesh.boundingSphere : Sphere3D();
    return create(0, 0, vt, gt, nt, bs, textureMesh);
  }

  export const Params = {
    ...BaseGeometry.Params,
    doubleSided: PD.Boolean(false, BaseGeometry.CustomQualityParamInfo),
    flipSided: PD.Boolean(false, BaseGeometry.ShadingCategory),
    flatShaded: PD.Boolean(false, BaseGeometry.ShadingCategory),
    ignoreLight: PD.Boolean(false, BaseGeometry.ShadingCategory),
    celShaded: PD.Boolean(false, BaseGeometry.ShadingCategory),
    xrayShaded: PD.Select<boolean | 'inverted'>(
      false,
      [
        [false, 'Off'],
        [true, 'On'],
        ['inverted', 'Inverted'],
      ],
      BaseGeometry.ShadingCategory,
    ),
    transparentBackfaces: PD.Select(
      'off',
      PD.arrayToOptions(['off', 'on', 'opaque'] as const),
      BaseGeometry.ShadingCategory,
    ),
    solidInterior: PD.Boolean(false, {
      ...BaseGeometry.ShadingCategory,
      description: 'Render a solid cap where the camera near plane or a clip object cuts a closed surface',
    }),
    bumpFrequency: PD.Numeric(0, { min: 0, max: 10, step: 0.1 }, BaseGeometry.ShadingCategory),
    bumpAmplitude: PD.Numeric(1, { min: 0, max: 5, step: 0.1 }, BaseGeometry.ShadingCategory),
    interior: getInteriorParam(),
    animation: getAnimationParam(),
  };
  export type Params = typeof Params;

  export const Utils: GeometryUtils<TextureMesh, Params> = {
    Params,
    createEmpty,
    createValues,
    createValuesSimple,
    updateValues,
    updateBoundingSphere,
    createRenderableState,
    updateRenderableState,
    createPositionIterator,
    prepare,
  };

  /**
   * The position iterator reads the vertices from `meta.vertexData`, read back asynchronously before (`prepare`, for a
   * request with `positions`); without them its locations are `NullLocation`
   */
  function createPositionIterator(textureMesh: TextureMesh, transform: TransformData): LocationIterator {
    // looked up per call: the read back may follow the creation of the iterator (e.g. a new color theme)
    const getData = () => textureMesh.meta.vertexData;

    const groupCount = textureMesh.vertexCount;
    const instanceCount = transform.instanceCount.ref.value;
    const location = PositionLocation();
    const p = location.position;
    const n = location.normal;
    const m = transform.aTransform.ref.value;
    const getLocation = (groupIndex: number, instanceIndex: number) => {
      const data = getData();
      if (!data) return NullLocation;
      const { vertices, normals } = data;
      if (instanceIndex < 0) {
        Vec3.fromArray(p, vertices, groupIndex * 4);
        Vec3.fromArray(n, normals, groupIndex * 4);
      } else {
        Vec3.transformMat4Offset(p, vertices, m, 0, groupIndex * 4, instanceIndex * 16);
        Vec3.transformDirectionOffset(n, normals, m, 0, groupIndex * 4, instanceIndex * 16);
      }
      return location;
    };
    return LocationIterator(groupCount, instanceCount, 1, getLocation);
  }

  function createValues(
    textureMesh: TextureMesh,
    transform: TransformData,
    locationIt: LocationIterator,
    theme: Theme,
    props: PD.Values<Params>,
  ): TextureMeshValues {
    const { instanceCount, groupCount } = locationIt;
    const positionIt = Utils.createPositionIterator(textureMesh, transform);

    const color = createColors(locationIt, positionIt, theme.color);
    const marker = resolveInstanceGranularity(props.instanceGranularity, groupCount, instanceCount)
      ? createMarkers(instanceCount, 'instance')
      : createMarkers(instanceCount * groupCount, 'groupInstance');
    const overpaint = createEmptyOverpaint();
    const transparency = createEmptyTransparency();
    const emissive = createEmptyEmissive();
    const substance = createEmptySubstance();
    const clipping = createEmptyClipping();
    const wiggle = createEmptyWiggle();

    const counts = {
      drawCount: textureMesh.vertexCount,
      vertexCount: textureMesh.vertexCount,
      groupCount,
      instanceCount,
    };

    const invariantBoundingSphere = Sphere3D.clone(textureMesh.boundingSphere);
    const boundingSphere = calculateTransformBoundingSphere(
      invariantBoundingSphere,
      transform.aTransform.ref.value,
      instanceCount,
      0,
    );

    return {
      dGeometryType: ValueCell.create('textureMesh'),

      uGeoTexDim: textureMesh.geoTextureDim,
      tPosition: textureMesh.vertexTexture,
      tGroup: textureMesh.groupTexture,
      tNormal: textureMesh.normalTexture,
      dVaryingGroup: textureMesh.varyingGroup,

      boundingSphere: ValueCell.create(boundingSphere),
      invariantBoundingSphere: ValueCell.create(invariantBoundingSphere),
      uInvariantBoundingSphere: ValueCell.create(Vec4.ofSphere(invariantBoundingSphere)),

      ...color,
      ...marker,
      ...overpaint,
      ...transparency,
      ...emissive,
      ...substance,
      ...clipping,
      ...wiggle,
      ...transform,

      ...BaseGeometry.createValues(props, counts),
      drawCount: textureMesh.vertexCountCell,
      uVertexCount: textureMesh.vertexCountCell,
      drawIndirect: textureMesh.drawIndirect,
      dIsoChunks: textureMesh.isoChunks,
      sIsoTable: textureMesh.isoTable,
      sIsoVertices: textureMesh.isoVertices,
      uDoubleSided: ValueCell.create(props.doubleSided),
      dFlatShaded: ValueCell.create(props.flatShaded),
      dFlipSided: ValueCell.create(props.flipSided),
      dIgnoreLight: ValueCell.create(props.ignoreLight),
      dCelShaded: ValueCell.create(props.celShaded),
      dXrayShaded: ValueCell.create(
        props.xrayShaded === 'inverted' ? 'inverted' : props.xrayShaded === true ? 'on' : 'off',
      ),
      dTransparentBackfaces: ValueCell.create(props.transparentBackfaces),
      dSolidInterior: ValueCell.create(props.solidInterior),
      uBumpFrequency: ValueCell.create(props.bumpFrequency),
      uBumpAmplitude: ValueCell.create(props.bumpAmplitude),

      meta: ValueCell.create(textureMesh.meta),

      ...createInteriorValues(props.interior),
      ...createAnimationValues(props.animation),
    };
  }

  function createValuesSimple(
    textureMesh: TextureMesh,
    props: Partial<PD.Values<Params>>,
    colorValue: Color,
    sizeValue: number,
    transform?: TransformData,
  ) {
    const s = BaseGeometry.createSimple(colorValue, sizeValue, transform);
    const p = { ...PD.getDefaultValues(Params), ...props };
    return createValues(textureMesh, s.transform, s.locationIterator, s.theme, p);
  }

  function updateValues(values: TextureMeshValues, props: PD.Values<Params>) {
    BaseGeometry.updateValues(values, props);
    ValueCell.updateIfChanged(values.uDoubleSided, props.doubleSided);
    ValueCell.updateIfChanged(values.dFlatShaded, props.flatShaded);
    ValueCell.updateIfChanged(values.dFlipSided, props.flipSided);
    ValueCell.updateIfChanged(values.dIgnoreLight, props.ignoreLight);
    ValueCell.updateIfChanged(values.dCelShaded, props.celShaded);
    ValueCell.updateIfChanged(
      values.dXrayShaded,
      props.xrayShaded === 'inverted' ? 'inverted' : props.xrayShaded === true ? 'on' : 'off',
    );
    ValueCell.updateIfChanged(values.dTransparentBackfaces, props.transparentBackfaces);
    ValueCell.updateIfChanged(values.dSolidInterior, props.solidInterior);
    ValueCell.updateIfChanged(values.uBumpFrequency, props.bumpFrequency);
    ValueCell.updateIfChanged(values.uBumpAmplitude, props.bumpAmplitude);
    updateInteriorValues(values, props.interior);
    updateAnimationValues(values, props.animation);
  }

  function updateBoundingSphere(values: TextureMeshValues, textureMesh: TextureMesh) {
    const invariantBoundingSphere = Sphere3D.clone(textureMesh.boundingSphere);
    const boundingSphere = calculateTransformBoundingSphere(
      invariantBoundingSphere,
      values.aTransform.ref.value,
      values.instanceCount.ref.value,
      0,
    );

    if (!Sphere3D.equals(boundingSphere, values.boundingSphere.ref.value)) {
      ValueCell.update(values.boundingSphere, boundingSphere);
    }
    if (!Sphere3D.equals(invariantBoundingSphere, values.invariantBoundingSphere.ref.value)) {
      ValueCell.update(values.invariantBoundingSphere, invariantBoundingSphere);
      ValueCell.update(
        values.uInvariantBoundingSphere,
        Vec4.fromSphere(values.uInvariantBoundingSphere.ref.value, invariantBoundingSphere),
      );
    }
  }

  function createRenderableState(props: PD.Values<Params>): RenderableState {
    const state = BaseGeometry.createRenderableState(props);
    updateRenderableState(state, props);
    return state;
  }

  function updateRenderableState(state: RenderableState, props: PD.Values<Params>) {
    BaseGeometry.updateRenderableState(state, props);
    state.opaque = state.opaque && !props.xrayShaded;
    state.writeDepth = state.opaque;
  }
}
