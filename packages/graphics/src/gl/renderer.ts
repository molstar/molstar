/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import type { OcclusionPyramid } from './compute/cell-cull.js';
import { type SplatOrderGpu, updateSplatOrder } from './renderable/gaussian-splats-order.js';
import { createSplatSorter, isSplatSortSupported, SplatSortKernels } from './compute/splat-sort.js';
import type { GaussianSplatsValues } from './renderable/gaussian-splats.js';
import { Viewport } from '@molstar/core/util/viewport';
import type { RenderCamera } from './camera.js';
import type { Scene } from './scene.js';
import type { GraphicsContext, PipelineRequest, TargetDescriptor } from './context.js';
import { Mat4, Vec3, Vec4, Vec2 } from '@molstar/core/math/linear-algebra';
import { type GraphicsRenderable, type Frame, createFrame, LodView } from './renderable.js';
import { getDrawnBoundingSphere } from './renderable/util.js';
import {
  DrawFact,
  DrawRow,
  DrawState,
  type DrawTable,
  getDrawFacts,
  getDrawState,
  hasCameraBounds,
} from './draw-table.js';
import { Color } from '@molstar/core/util/color';
import { ValueCell, deepEqual } from '@molstar/core/util';
import { GlobalUniformSchema, type GlobalUniformValues, MaxLightCount } from './renderable/schema.js';
import type { GraphicsRenderVariant, ProgramItem } from './render-item.js';
import type { ObjectSignature } from './precompile.js';
import { getRenderableProgramSource } from './render-object.js';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { degToRad } from '@molstar/core/math/misc';
import type { Texture, Textures } from './texture.js';
import { arrayMapUpsert } from '@molstar/core/util/array';
import { Frustum3D } from '@molstar/core/math/geometry/primitives/frustum3d';
import { Plane3D } from '@molstar/core/math/geometry/primitives/plane3d';
import { Sphere3D } from '@molstar/core/math/geometry';
import { Clip } from '@molstar/core/util/clip';
import { RenderState } from './render-state.js';
import {
  Flag,
  getObjectFacts,
  getObjectState,
  getPassDraws,
  getSolidInteriorFillStates,
  getSolidInteriorMarkState,
  isPassDrawn,
  OpaqueState,
  type Pass as PassDrawsPass,
  PassDraws,
} from './renderer/pass-draws.js';
import { createRendererDrawLists, type PassPart, WholePass } from './renderer/draw-lists.js';
import { createSolidInteriorCaps } from './renderer/solid-interior-caps.js';
import { createGpuCulling } from './renderer/gpu-culling.js';
import { createRenderPassAttachments } from './render-pass-attachments.js';

export interface RendererStats {
  programCount: number;
  shaderCount: number;

  attributeCount: number;
  elementsCount: number;
  framebufferCount: number;
  renderbufferCount: number;
  textureCount: number;
  vertexArrayCount: number;

  drawCount: number;
  instanceCount: number;
  instancedDrawCount: number;
}

export enum PickType {
  None = 0,
  Object = 1,
  Instance = 2,
  Group = 3,
}

export enum MarkingType {
  None = 0,
  Depth = 1,
  Mask = 2,
}

interface Renderer {
  readonly stats: RendererStats;
  readonly props: Readonly<RendererProps>;
  readonly light: Readonly<Light>;
  readonly ambientColor: Vec3;

  clear: (toBackgroundColor: boolean, ignoreTransparentBackground?: boolean, forceToTransparency?: boolean) => void;
  clearDepth: (packed?: boolean) => void;
  /**
   * `side`: a frame of a pass's own within the frame (culled anew, e.g. the tracing's thickness), not counted as one:
   * the frame updated before it and after it are the same one
   */
  update: (camera: RenderCamera, scene: Scene, frame: Frame, side?: boolean) => void;
  /**
   * Culls on the GPU (merged render objects, `Renderable.prepareCull`) with the planes of the last `update`: outside of
   * passes (compute), before the passes of the frame that draw with it; others draw the CPU culling's lists
   */
  /** `clipGhost`: the ghost pass draws this frame (`renderClipGhost`), its lists culled on the GPU too */
  prepareCull: (scene: Scene.Group, twoPhase?: boolean, clipGhost?: boolean) => void;
  /**
   * The gaussian splats of `scene` back to front for `camera` (their order, `tSplatOrder` or on the GPU `sSplatOrder`),
   * sorted again where the view turned: outside of passes (compute), before the passes of the view that draw them
   */
  prepareGaussianSplats: (scene: Scene, camera: RenderCamera) => void;
  /**
   * Culls on the GPU for a view of its own (a frame of its own: the thickness pass, a ray's pick) with the planes of the
   * last `update`, no pass open: the objects culled on the GPU (`Renderable.prepareCullView`); the others draw the CPU
   * culling's lists there. Not tested against a pyramid
   */
  prepareCullView: (scene: Scene.Group) => void;
  /** The objects whose second phase is to run after the last `prepareCull` with `twoPhase` (none: no second part) */
  readonly cullSecondPending: number;
  /**
   * After the opaque pass's first part, no pass open: the second phase of GPU culling against `pyramid`, of the depth so
   * far, and of the sphere raster; none: left out. Then `renderOpaqueSecond`
   */
  prepareCullSecond: (scene: Scene.Group, pyramid: OcclusionPyramid | undefined) => void;
  /** The opaque pass's second part: what the second phase found, on the first part's target */
  renderOpaqueSecond: (group: Scene.Group, camera: RenderCamera) => void;
  /**
   * The tracing pass's second part (illumination's input), as `renderOpaqueSecond`; the passes after it draw both
   * phases (the back depth)
   */
  renderTracingSecond: (group: Scene.Group, camera: RenderCamera) => void;
  /**
   * Draws the small opaque spheres of `scene` by the sphere raster (`sphereRaster`, WebGPU with `atomicMinMax64`) for
   * the next `renderOpaque`, which resolves them; after `update` and `prepareCull`, no pass open (compute)
   */
  prepareSphereRaster: (scene: Scene.Group, camera: RenderCamera) => void;
  /** The objects the last `prepareCull` culled on the GPU */
  readonly gpuCulledCount: number;
  setTime: (time: number) => void;

  /** `occluder`: for an overlay (`dOverlay`), the scene's opaque depth its occlusion tests against */
  renderPick: (
    group: Scene.Group,
    camera: RenderCamera,
    variant: 'pick' | 'depth',
    pickType: PickType,
    occluder?: Renderer.Occluder | null,
  ) => void;
  renderDepth: (group: Scene.Group, camera: RenderCamera) => void;
  renderDepthOpaque: (group: Scene.Group, camera: RenderCamera) => void;
  renderDepthOpaqueBack: (group: Scene.Group, camera: RenderCamera) => void;
  renderDepthTransparent: (group: Scene.Group, camera: RenderCamera, depthTexture: Texture) => void;
  /**
   * The ghost pass: what the clip objects remove, as depth (the `depth` variant, `uClipGhost`), of
   * the objects with a ghost (`DrawFact.ClipGhost`) the clip objects reach; `scale`: of a target at that fraction of
   * the resolution (the viewport set so), sizes in pixels scaled with it
   */
  renderClipGhost: (group: Scene.Group, camera: RenderCamera, scale?: number) => void;
  renderMarkingDepth: (group: Scene.Group, camera: RenderCamera) => void;
  renderMarkingMask: (group: Scene.Group, camera: RenderCamera, depthTexture: Texture | null) => void;
  renderEmissiveOpaque: (
    group: Scene.Group,
    camera: RenderCamera,
    depthTexture: Texture,
    occludeWithOpaqueDepth: boolean,
  ) => void;
  renderEmissiveTransparent: (group: Scene.Group, camera: RenderCamera, depthTexture: Texture) => void;
  renderTracing: (group: Scene.Group, camera: RenderCamera) => void;
  renderBlended: (group: Scene, camera: RenderCamera) => void;
  /**
   * An overlay (`dOverlay`), as `renderBlended`, over the scene's opaque depth (`occluder`): the transparent pass also
   * draws the occluded parts of opaque objects (`uOccludedAlpha` < 1)
   */
  renderOverlay: (group: Scene.Group, camera: RenderCamera, occluder: Renderer.Occluder | null) => void;
  renderOpaque: (group: Scene.Group, camera: RenderCamera) => void;
  renderBlendedTransparent: (group: Scene.Group, camera: RenderCamera) => void;
  /**
   * Gaussian splats, blended over what is drawn (after the other transparency's composite), tested against the
   * attached opaque depth, not writing it
   */
  renderGaussianSplats: (group: Scene.Group, camera: RenderCamera) => void;
  renderVolume: (group: Scene.Group, camera: RenderCamera, depthTexture: Texture) => void;
  renderWboitTransparent: (group: Scene.Group, camera: RenderCamera, depthTexture: Texture) => void;
  renderDpoitTransparent: (
    group: Scene.Group,
    camera: RenderCamera,
    depthTexture: Texture,
    dpoitTextures: { depth: Texture; frontColor: Texture },
  ) => void;

  /**
   * The pipelines `pass` draws the objects of `group` with into a target like `target`, for
   * `resources.preparePipelines`: every object the pass would draw regardless of visibility and culling, and with `all`
   * also those that the objects' current state leaves out (e.g. unmarked objects for the marking passes)
   */
  collect: (
    pass: Renderer.Pass,
    group: Scene.Group,
    target: TargetDescriptor,
    options?: { all?: boolean },
  ) => PipelineRequest[];
  /**
   * The pipelines `pass` draws objects of the signatures with, through their program items
   * (`createSignatureProgramItem`): what drawing those objects would request, without them, for each way they may draw
   * (`ObjectSignature.opaque`, ...); `all` as for `collect`
   */
  collectSignatures: (
    pass: Renderer.Pass,
    items: readonly { readonly signature: ObjectSignature; readonly item: ProgramItem }[],
    target: TargetDescriptor,
    options?: { all?: boolean },
  ) => PipelineRequest[];

  setProps: (props: Partial<RendererProps>) => void;
  setViewport: (x: number, y: number, width: number, height: number) => void;
  setTransparentBackground: (value: boolean) => void;
  setDrawingBufferSize: (width: number, height: number) => void;
  setPixelRatio: (value: number) => void;
  setOcclusionTest: (f: ((s: Sphere3D) => boolean) | null) => void;
  /** The test of `setOcclusionTest`, e.g. to draw what the scene's depth does not hide (an overlay) and restore it */
  readonly occlusionTest: ((s: Sphere3D) => boolean) | null;
  /**
   * The depth pyramid GPU culling tests occlusion against (`prepareCull`, the previous frame's), null for none (e.g.
   * where passes need what is hidden, illumination)
   */
  setOcclusionPyramid: (get: (() => OcclusionPyramid | undefined) | null) => void;

  dispose: () => void;
}

export const RendererParams = {
  backgroundColor: PD.Color(Color(0x000000), { description: 'Background color of the 3D canvas' }),

  pickingAlphaThreshold: PD.Numeric(
    0.5,
    { min: 0.0, max: 1.0, step: 0.01 },
    { description: 'The minimum opacity value needed for an object to be pickable.' },
  ),

  colorMarker: PD.Boolean(true, {
    description:
      'Tint marked and dim unmarked objects, by the marking pass when enabled, otherwise in the material color.',
  }),
  highlightColor: PD.Color(Color.fromNormalizedRgb(1.0, 0.4, 0.6)),
  selectColor: PD.Color(Color.fromNormalizedRgb(0.2, 1.0, 0.1)),
  dimColor: PD.Color(Color.fromNormalizedRgb(1.0, 1.0, 1.0)),
  highlightStrength: PD.Numeric(0.3, { min: 0.0, max: 1.0, step: 0.1 }),
  selectStrength: PD.Numeric(0.3, { min: 0.0, max: 1.0, step: 0.1 }),
  dimStrength: PD.Numeric(0.0, { min: 0.0, max: 1.0, step: 0.1 }),
  markerPriority: PD.Select(1, [
    [1, 'Highlight'],
    [2, 'Select'],
  ]),

  xrayEdgeFalloff: PD.Numeric(1, { min: 0.0, max: 3.0, step: 0.1 }),
  celSteps: PD.Numeric(5, { min: 2, max: 16, step: 1 }),
  exposure: PD.Numeric(1, { min: 0.0, max: 3.0, step: 0.01 }),

  light: PD.ObjectList(
    {
      inclination: PD.Numeric(150, { min: 0, max: 180, step: 1 }),
      azimuth: PD.Numeric(320, { min: 0, max: 360, step: 1 }),
      color: PD.Color(Color.fromNormalizedRgb(1.0, 1.0, 1.0)),
      intensity: PD.Numeric(0.6, { min: 0.0, max: 5.0, step: 0.01 }),
    },
    (o) => Color.toHexString(o.color),
    {
      defaultValue: [
        {
          inclination: 150,
          azimuth: 320,
          color: Color.fromNormalizedRgb(1.0, 1.0, 1.0),
          intensity: 0.6,
        },
      ],
    },
  ),
  ambientColor: PD.Color(Color.fromNormalizedRgb(1.0, 1.0, 1.0)),
  ambientIntensity: PD.Numeric(0.4, { min: 0.0, max: 2.0, step: 0.01 }),

  enableAnimation: PD.Boolean(true, { description: 'Enable time-based animations.' }),
  drawLists: PD.Boolean(true, {
    description:
      'Replay the passes from recordings while they draw the same objects (WebGPU render bundles); objects changing every frame are drawn directly.',
  }),
  sphereRaster: PD.Select('auto', PD.arrayToOptions(['off', 'auto', 'on', 'only'] as const), {
    description:
      'Draw small opaque spheres by compute kernels instead of impostor quads (WebGPU with 64-bit atomics, perspective views). "auto" draws those of objects culled on the GPU (merged, many instances: their instances split), where measured faster; "on" also those of unmerged objects (split per sphere, their quads drawn all the same: no faster); merged objects need GPU culling either way (`gpuCulling`); "only" (for testing) draws all spheres of the eligible objects by the kernels, whatever their size, and none as quads (spheres through the near plane are missing).',
  }),
  sphereRasterRadius: PD.Numeric(
    24,
    { min: 1, max: 64, step: 1 },
    {
      description:
        'The largest sphere radius on screen, in CSS pixels, the sphere raster draws; larger ones are quads.',
    },
  ),
  sphereRasterObjects: PD.Numeric(
    8,
    { min: 0, max: 256, step: 1 },
    {
      description:
        'With `sphereRaster` "auto": of the objects culled on the GPU for the raster only (many spheres, few instances), the most culled there and drawn by the raster, the largest visible ones; the others culled on the CPU and drawn as quads. Each costs its own dispatches a frame, whatever its size.',
    },
  ),
};
export type RendererProps = PD.Values<typeof RendererParams>;

export type Light = {
  count: number;
  direction: number[];
  color: number[];
};

/**
 * The ghost pass's levels of detail (`Renderer.renderClipGhost`), chosen as if a CSS pixel were this many times
 * larger: coarser, a removed particle as its few spheres
 */
const ClipGhostLevelScale = 4;

const tmpDir = Vec3();
const tmpColor = Vec3();
function getLight(props: RendererProps['light'], light?: Light): Light {
  const count = Math.min(props.length, MaxLightCount);
  const { direction, color } = light || {
    direction: new Array(count * 3).fill(0),
    color: new Array(count * 3).fill(0),
  };
  for (let i = 0; i < count; ++i) {
    const p = props[i];
    Vec3.directionFromSpherical(tmpDir, degToRad(p.inclination), degToRad(p.azimuth), 1);
    Vec3.toArray(tmpDir, direction, i * 3);
    Vec3.scale(tmpColor, Color.toVec3Normalized(tmpColor, p.color), p.intensity);
    Vec3.toArray(tmpColor, color, i * 3);
  }
  return { count, direction, color };
}

export function getTransformedLightDirection(light: Light, t: Mat4): Light['direction'] {
  const tld = new Array(light.count * 3);
  for (let i = 0, il = light.count; i < il; ++i) {
    Vec3.fromArray(tmpDir, light.direction, i * 3);
    Vec3.transformDirection(tmpDir, tmpDir, t);
    Vec3.toArray(tmpDir, tld, i * 3);
  }
  return tld;
}

namespace Renderer {
  /**
   * The scene's depth an overlay is drawn over: the opaque depth's texture (`tDepth`), the front transparent depth with
   * its alpha (`tDepthAlpha`, where the frame wrote it) and the planes they were drawn with
   */
  export type Occluder = {
    readonly depth: Texture;
    readonly depthAlpha?: Texture;
    readonly near: number;
    readonly far: number;
  };

  export type Pass = PassDrawsPass;

  const enum Mask {
    All = 0,
    Opaque = 1,
    Transparent = 2,
  }

  /** 1x1 depth texture cleared to the far plane */
  function createFarDepthTexture(ctx: GraphicsContext) {
    const texture = ctx.resources.texture('image-depth', 'depth', 'ushort', 'nearest');
    texture.define(1, 1);
    const target = createRenderPassAttachments(ctx.resources, [], texture);
    ctx.beginRenderPass(target, { clear: { depth: 0 } });
    ctx.endRenderPass();
    target.destroy();
    return texture;
  }

  export function create(ctx: GraphicsContext, props: Partial<RendererProps> = {}): Renderer {
    const { stats } = ctx;
    const p = PD.merge(RendererParams, PD.getDefaultValues(RendererParams), props);
    const light = getLight(p.light);

    const viewport = Viewport();
    const initialDrawingBufferSize = ctx.getDrawingBufferSize();
    const drawingBufferSize = Vec2.create(initialDrawingBufferSize.width, initialDrawingBufferSize.height);
    const bgColor = Color.toVec3Normalized(Vec3(), p.backgroundColor);

    let transparentBackground = false;
    let isOccluded: ((s: Sphere3D) => boolean) | null = null;

    // placeholders at the far plane (depth is reversed: 0, packed as zeros)
    const emptyDepthAlphaTexture = ctx.resources.texture('image-uint8', 'rgba', 'ubyte', 'nearest');
    emptyDepthAlphaTexture.define(1, 1);
    emptyDepthAlphaTexture.load({ array: new Uint8Array([0, 0, 0, 0]), width: 1, height: 1 });
    const emptyDepthTexture = ctx.capabilities.depthTexture ? createFarDepthTexture(ctx) : emptyDepthAlphaTexture;
    // every global texture from the start: the `dpoit` mode's programs sample the dpoit textures in all their color
    // draws, also those before the first dpoit pass sets them (the opaque pass), which would read whatever their units
    // hold, e.g. a texture attached to the target ('invalid operation' on WebGL)
    const sharedTexturesList: Textures = [
      ['tDepth', emptyDepthTexture],
      ['tDepthAlpha', emptyDepthAlphaTexture],
      ['tDpoitDepth', emptyDepthAlphaTexture],
      ['tDpoitFrontColor', emptyDepthAlphaTexture],
    ];

    /** The values each cell was last set to (`updateArray`) */
    const lastArrays = new Map<ValueCell<ArrayLike<number>>, Float64Array>();
    /**
     * Sets `cell` to the matrix or vector `a` (computed in place, into the arrays the cells hold) where it differs from
     * the values last set: unchanged, as at rest, no new box (`ValueCell.update` makes one) and no new version (the
     * uniform blocks' writers then write nothing)
     */
    const updateArray = <T extends ArrayLike<number>>(cell: ValueCell<T>, a: T) => {
      const n = a.length;
      let last = lastArrays.get(cell);
      if (last !== undefined && cell.ref.value === a && last.length === n) {
        let same = true;
        for (let i = 0; i < n && same; ++i) same = last[i] === a[i];
        if (same) return;
      }
      if (last === undefined || last.length !== n) lastArrays.set(cell, (last = new Float64Array(n)));
      for (let i = 0; i < n; ++i) last[i] = a[i];
      ValueCell.update(cell, a);
    };

    const model = Mat4();
    const view = Mat4();
    const invView = Mat4();
    const modelView = Mat4();
    const invModelView = Mat4();
    const invProjection = Mat4();
    const modelViewProjection = Mat4();
    const invModelViewProjection = Mat4();
    const invHeadRotation = Mat4();
    const modelViewEye = Mat4();
    const invModelViewEye = Mat4();

    const cameraDir = Vec3();
    const cameraPosition = Vec3();
    const cameraTarget = Vec3();
    const cameraPlane = Plane3D();
    const viewOffset = Vec2();
    const frustum = Frustum3D();
    /**
     * `cameraPlane` and `frustum` in the unscaled world (their constants over the camera scale), the instances' culling's
     * (`Renderable.cull`, the GPU's): of a point `p` the scaled planes' distances of `s p` over `s`, so culled and the
     * levels chosen as unscaled (XR's world scale); the same as them unscaled
     */
    const cullPlane = Plane3D();
    const cullFrustum = Frustum3D();
    const lodView = LodView();
    /** the culling view of the last update (camera plane, frustum, level scales), see `LodView.version` */
    const lastView = new Float64Array(30);
    let viewVersions = 0;
    /** whether the view was compared since the last update (`updateViewVersion`) */
    let viewChecked = false;
    /** compares the culling view with that of the last cull: a new `lodView.version` where it differs */
    const updateViewVersion = () => {
      let changed = lodView.version === 0;
      const set = (i: number, v: number) => {
        if (lastView[i] !== v) {
          lastView[i] = v;
          changed = true;
        }
      };
      const setPlane = (p: Plane3D, o: number) => {
        set(o, p.normal[0]);
        set(o + 1, p.normal[1]);
        set(o + 2, p.normal[2]);
        set(o + 3, p.constant);
      };
      setPlane(cullPlane, 0);
      for (let i = 0; i < 6; ++i) setPlane(cullFrustum[i as Frustum3D.PlaneIndex], (i + 1) * 4);
      set(28, lodView.depthScale);
      set(29, lodView.pixelSize);
      if (changed) lodView.version = ++viewVersions;
    };

    let modelScale = 1;
    /** of the object being drawn, with `drawFacts` (`isVisible`, `setRowBoundingSphere`) */
    const boundingSphere = Sphere3D();
    /** the facts of the object being drawn (`DrawFact`): its render states and culling read them */
    let drawFacts = 0;
    let currentFrame: Frame = createFrame();
    /**
     * the ghost pass's culling (`renderClipGhost`), into lists of their own (the renderables' list sets): a view of its
     * own, its levels coarser (`ClipGhostLevelScale`), a version of its own
     */
    const ghostLodView = LodView();
    let ghostLodOf = -1;
    const updateGhostCull = () => {
      if (ghostLodOf !== lodView.version || lodView.version === 0) {
        ghostLodView.depthScale = lodView.depthScale * ClipGhostLevelScale;
        ghostLodView.pixelSize = lodView.pixelSize * ClipGhostLevelScale;
        ghostLodView.version = lodView.version === 0 ? 0 : ++viewVersions;
        ghostLodOf = lodView.version;
      }
    };
    /** frames updated (`update` with another frame, not a side one): what the visible sets' changes are counted in */
    let frameCount = 0;
    /** the last frame counted (`frameCount`): a pass's side frame and back to it is one frame */
    let countedFrame: Frame | undefined;

    const ambientColor = Vec3();
    Vec3.scale(ambientColor, Color.toArrayNormalized(p.ambientColor, ambientColor, 0), p.ambientIntensity);

    const globalUniforms: GlobalUniformValues = {
      uDrawId: ValueCell.create(0),

      uModel: ValueCell.create(Mat4.identity()),
      uView: ValueCell.create(view),
      uInvView: ValueCell.create(invView),
      uModelView: ValueCell.create(modelView),
      uInvModelView: ValueCell.create(invModelView),
      uInvProjection: ValueCell.create(invProjection),
      uProjection: ValueCell.create(Mat4()),
      uModelViewProjection: ValueCell.create(modelViewProjection),
      uInvModelViewProjection: ValueCell.create(invModelViewProjection),
      uHasHeadRotation: ValueCell.create(false),
      uInvHeadRotation: ValueCell.create(invHeadRotation),
      uHasEyeCamera: ValueCell.create(false),
      uModelViewEye: ValueCell.create(modelViewEye),
      uInvModelViewEye: ValueCell.create(invModelViewEye),

      uIsOrtho: ValueCell.create(1),
      uViewOffset: ValueCell.create(viewOffset),
      uModelScale: ValueCell.create(1),

      uPixelRatio: ValueCell.create(ctx.pixelRatio),
      uViewport: ValueCell.create(Viewport.toVec4(Vec4(), viewport)),
      uDrawingBufferSize: ValueCell.create(drawingBufferSize),

      uCameraPosition: ValueCell.create(cameraPosition),
      uCameraDir: ValueCell.create(cameraDir),
      uCameraPlane: ValueCell.create(Plane3D.toArray(cameraPlane, Vec4(), 0)),
      uNear: ValueCell.create(1),
      uFar: ValueCell.create(10000),
      uFog: ValueCell.create(true),
      uFogNear: ValueCell.create(1),
      uFogFar: ValueCell.create(10000),
      uFogColor: ValueCell.create(bgColor),

      uRenderMask: ValueCell.create(0),
      uMarkingDepthTest: ValueCell.create(false),
      uDepthBack: ValueCell.create(false),
      uClipGhost: ValueCell.create(false),
      uPickType: ValueCell.create(PickType.None),
      uMarkingType: ValueCell.create(MarkingType.None),
      uSolidInteriorPass: ValueCell.create(0),
      uSolidInteriorClip: ValueCell.create(-1),

      uTransparentBackground: ValueCell.create(false),

      uLightDirection: ValueCell.create(light.direction),
      uLightColor: ValueCell.create(light.color),
      uAmbientColor: ValueCell.create(ambientColor),

      uPickingAlphaThreshold: ValueCell.create(p.pickingAlphaThreshold),

      uHighlightColor: ValueCell.create(Color.toVec3Normalized(Vec3(), p.highlightColor)),
      uSelectColor: ValueCell.create(Color.toVec3Normalized(Vec3(), p.selectColor)),
      uDimColor: ValueCell.create(Color.toVec3Normalized(Vec3(), p.dimColor)),
      uHighlightStrength: ValueCell.create(p.highlightStrength),
      uSelectStrength: ValueCell.create(p.selectStrength),
      uDimStrength: ValueCell.create(p.dimStrength),
      uMarkerPriority: ValueCell.create(p.markerPriority),
      uMarkerAverage: ValueCell.create(0),

      uXrayEdgeFalloff: ValueCell.create(p.xrayEdgeFalloff),
      uCssHeight: ValueCell.create(1),
      uOccluderNear: ValueCell.create(1),
      uOccluderFar: ValueCell.create(10000),
      uCelSteps: ValueCell.create(p.celSteps),
      uExposure: ValueCell.create(p.exposure),

      uTime: ValueCell.create(0),
      uEnableAnimation: ValueCell.create(p.enableAnimation),
    };
    const globals = ctx.createGlobalBindings(GlobalUniformSchema, Object.entries(globalUniforms), sharedTexturesList);

    /**
     * Whether clip objects remove all of `r` within `sphere` (unscaled): with clip objects and no clipping layer
     * (`DrawFact.ClipCull`), its bounding sphere clipped whole (`Clip.classifySphere`)
     */
    const clipSphere = Sphere3D();
    const isClippedWhole = (r: GraphicsRenderable, sphere: Sphere3D) => {
      const values = r.values as { dClipObjects?: ValueCell<boolean>; dClipping?: ValueCell<boolean> };
      if (!values.dClipObjects?.ref.value || values.dClipping?.ref.value) return false;
      const objects = (r.values as { clipObjects?: ValueCell<Clip.Objects> }).clipObjects?.ref.value;
      return !!objects && Clip.classifySphere(objects, sphere) === Clip.SphereClass.All;
    };

    /** drawing the ghost pass (`renderClipGhost`): what the clip objects remove, the visibility inverted */
    let clipGhost = false;
    /**
     * whether `r` has a ghost (`DrawFact.ClipGhost`) the clip objects reach: its bounding sphere clipped whole or cut
     * (`Clip.classifySphere` not none); with a clipping layer (flags exempting groups) not classified, so drawn
     */
    const isClipGhosted = (r: GraphicsRenderable, sphere: Sphere3D) => {
      const values = r.values as {
        dClipObjects?: ValueCell<boolean>;
        dClipping?: ValueCell<boolean>;
        clipGhost?: ValueCell<boolean>;
      };
      if (!values.dClipObjects?.ref.value || values.clipGhost?.ref.value === false) return false;
      if (values.dClipping?.ref.value) return true;
      const objects = (r.values as { clipObjects?: ValueCell<Clip.Objects> }).clipObjects?.ref.value;
      return !!objects && Clip.classifySphere(objects, sphere) !== Clip.SphereClass.None;
    };

    /** also sets `boundingSphere` to the scaled bounding sphere of `r`, and `drawFacts` */
    const isVisible = (r: GraphicsRenderable, variant: GraphicsRenderVariant) => {
      if (r.state.disposed || !r.state.visible || (!r.state.pickable && variant === 'pick')) {
        return false;
      }

      if (!r.values.drawCount.ref.value) {
        return false;
      }
      // the ghost pass: objects with a ghost (`DrawFact.ClipGhost`)
      if (clipGhost && (getDrawFacts(r.values) & DrawFact.ClipGhost) === 0) return false;

      // set before any `return true`: drawing reads them (culling, the render state), not the last object's
      drawFacts = getDrawFacts(r.values);
      // as drawn: moved by a transient transform
      const drawnBoundingSphere = getDrawnBoundingSphere(r.values);
      Sphere3D.scaleNX(boundingSphere, drawnBoundingSphere, modelScale);

      // its size in the world depends on the camera, its bounds don't hold
      if (hasCameraBounds(r.values)) return true;

      if (!Frustum3D.intersectsSphere3D(frustum, boundingSphere)) {
        return false;
      }
      // clipped whole (the unscaled sphere: clip objects test model space divided by the model scale); the ghost
      // pass draws what the clip objects remove
      if (clipGhost ? !isClipGhosted(r, drawnBoundingSphere) : isClippedWhole(r, drawnBoundingSphere)) return false;

      const [minDistance, maxDistance] = r.values.uLod.ref.value;
      if (minDistance !== 0 || maxDistance !== 0) {
        const { center, radius } = boundingSphere;
        const d = Plane3D.distanceToPoint(cameraPlane, center);
        if (d + radius < minDistance * modelScale) return false;
        if (d - radius > maxDistance * modelScale) return false;
      }

      if (modelScale === 1 && isOccluded !== null && isOccluded(boundingSphere)) {
        return false;
      }

      return true;
    };

    /**
     * sets `boundingSphere` to the scaled bounding sphere of a table's row, and `drawFacts` (as `isVisible` does for an
     * object): from the row, without reading the object's values
     */
    const setRowBoundingSphere = (table: DrawTable, slot: number) => {
      const o = slot * DrawRow.Stride;
      const { floats } = table;
      const s = modelScale;
      Vec3.set(
        boundingSphere.center,
        floats[o + DrawRow.CenterX] * s,
        floats[o + DrawRow.CenterY] * s,
        floats[o + DrawRow.CenterZ] * s,
      );
      boundingSphere.radius = floats[o + DrawRow.Radius] * s;
      drawFacts = table.bits[o + DrawRow.Facts];
    };

    /**
     * The planes deciding visibility, the frustum's six and the camera plane (normal, constant), as of the last `update`;
     * with the model scale and the occlusion test, the view what was derived from it (the tables' visibility, the draw
     * lists' snapshots) compares by value: passes of several cameras alternate in a frame (e.g. the camera helper's)
     */
    const viewPlanes = new Float64Array(28);
    const setViewPlane = (o: number, plane: Plane3D) => {
      const { normal, constant } = plane;
      viewPlanes[o] = normal[0];
      viewPlanes[o + 1] = normal[1];
      viewPlanes[o + 2] = normal[2];
      viewPlanes[o + 3] = constant;
    };
    const setViewPlanes = () => {
      for (let i = 0; i < 6; ++i) setViewPlane(i * 4, frustum[i as Frustum3D.PlaneIndex]);
      setViewPlane(24, cameraPlane);
    };
    const isSameView = (planes: Float64Array) => {
      for (let i = 0; i < 28; ++i) if (planes[i] !== viewPlanes[i]) return false;
      return true;
    };

    /**
     * The visibility of a table's rows for a view, by slot (1 visible): state (not disposed, visible), draw count,
     * frustum, lod, occlusion; as `isVisible` without the pass (`isPassDrawn`). Its key: the table's and the view's
     * versions, the model scale, the occlusion test and with it the frame (Hi-Z changes between frames)
     */
    type RowVisibility = {
      tableVersion: number;
      readonly view: Float64Array;
      modelScale: number;
      isOccluded: ((s: Sphere3D) => boolean) | null;
      frame: Frame | undefined;
      visible: Uint8Array;
      /** as `visible` without the frustum and occlusion tests: drawn, the GPU clips what is outside */
      present: Uint8Array;
      visibleCount: number;
      presentCount: number;
      /** the `frameCount` at which `visible` last changed */
      changedAt: number;
    };
    const rowVisibilities = new WeakMap<DrawTable, RowVisibility>();
    /** the ghost pass's (`renderClipGhost`), the clip test inverted */
    const ghostRowVisibilities = new WeakMap<DrawTable, RowVisibility>();
    const occlusionSphere = Sphere3D();
    /** drawing a pass from a draw list (`RendererDrawLists.render`): its rows unculled while the view moves */
    let drawingFromList = false;
    /** frames since the visible set last changed that count as moving (a rotation changes it nearly every frame) */
    const MovingFrames = 8;
    /**
     * the share of the present rows visible from which a moving view draws them unculled: zoomed in to half of 1000
     * boxes visible, unculled replays cost less than recording per frame; to 8%, about the same
     */
    const UnculledShare = 0.15;
    /**
     * The rows to draw by the view (`getRowVisibilities`): while the visible set changed within `MovingFrames` frames
     * and at least `UnculledShare` of the present rows are visible, a pass replayed from a draw list (`drawingFromList`)
     * draws the present rows: unculled, its recording stays valid as objects enter and leave the view (re-recorded per
     * frame otherwise), the GPU drawing what is off-screen at under 1 µs a draw; at rest culled again (recorded once)
     */
    const getRowVisibility = (table: DrawTable) => {
      const v = getRowVisibilities(table);
      if (
        drawingFromList &&
        frameCount - v.changedAt <= MovingFrames &&
        v.visibleCount >= UnculledShare * v.presentCount
      ) {
        return v.present;
      }
      return v.visible;
    };
    const getRowVisibilities = (table: DrawTable) => {
      table.refresh();
      const visibilities = clipGhost ? ghostRowVisibilities : rowVisibilities;
      let v = visibilities.get(table);
      if (!v) {
        v = {
          tableVersion: -1,
          view: new Float64Array(28).fill(NaN),
          modelScale: NaN,
          isOccluded: null,
          frame: undefined,
          visible: new Uint8Array(0),
          present: new Uint8Array(0),
          visibleCount: 0,
          presentCount: 0,
          changedAt: -Infinity,
        };
        visibilities.set(table, v);
      }
      if (
        v.tableVersion === table.version &&
        isSameView(v.view) &&
        v.modelScale === modelScale &&
        v.isOccluded === isOccluded &&
        (isOccluded === null || v.frame === currentFrame)
      ) {
        return v;
      }
      v.tableVersion = table.version;
      v.view.set(viewPlanes);
      v.modelScale = modelScale;
      v.isOccluded = isOccluded;
      v.frame = currentFrame;
      if (v.visible.length < table.capacity) {
        v.visible = new Uint8Array(table.capacity);
        v.present = new Uint8Array(table.capacity);
      }
      const { visible, present } = v;
      let visibleCount = 0;
      let presentCount = 0;
      let changed = false;
      const { order, count, floats, bits } = table;
      const s = modelScale;
      const occlude = s === 1 ? isOccluded : null;
      const p = viewPlanes;
      // inlined plane tests: the loop runs over every row per view
      for (let k = 0; k < count; ++k) {
        const slot = order[k];
        const o = slot * DrawRow.Stride;
        const was = visible[slot];
        visible[slot] = 0;
        present[slot] = 0;
        const drawable =
          (bits[o + DrawRow.State] & (DrawState.Disposed | DrawState.Visible)) === DrawState.Visible &&
          (bits[o + DrawRow.Facts] & DrawFact.DrawCount) !== 0 &&
          (!clipGhost || (bits[o + DrawRow.Facts] & DrawFact.ClipGhost) !== 0);
        if (drawable && (bits[o + DrawRow.Facts] & DrawFact.CameraBounds) !== 0) {
          // its size in the world depends on the camera, its bounds don't hold
          present[slot] = 1;
          presentCount += 1;
          visible[slot] = 1;
          visibleCount += 1;
        } else if (drawable) {
          const x = floats[o + DrawRow.CenterX] * s;
          const y = floats[o + DrawRow.CenterY] * s;
          const z = floats[o + DrawRow.CenterZ] * s;
          const r = floats[o + DrawRow.Radius] * s;
          let drawn = true;
          if (clipGhost) {
            // the ghost pass: objects with a ghost the clip objects reach
            if ((bits[o + DrawRow.Facts] & DrawFact.ClipGhost) === 0) {
              drawn = false;
            } else {
              Vec3.set(
                clipSphere.center,
                floats[o + DrawRow.CenterX],
                floats[o + DrawRow.CenterY],
                floats[o + DrawRow.CenterZ],
              );
              clipSphere.radius = floats[o + DrawRow.Radius];
              if (!isClipGhosted(table.objects[slot]!, clipSphere)) drawn = false;
            }
          } else if ((bits[o + DrawRow.Facts] & DrawFact.ClipCull) !== 0) {
            Vec3.set(
              clipSphere.center,
              floats[o + DrawRow.CenterX],
              floats[o + DrawRow.CenterY],
              floats[o + DrawRow.CenterZ],
            );
            clipSphere.radius = floats[o + DrawRow.Radius];
            if (isClippedWhole(table.objects[slot]!, clipSphere)) drawn = false;
          }
          const minDistance = floats[o + DrawRow.LodMin];
          const maxDistance = floats[o + DrawRow.LodMax];
          if (drawn && (minDistance !== 0 || maxDistance !== 0)) {
            const d = p[24] * x + p[25] * y + p[26] * z + p[27];
            if (d + r < minDistance * s || d - r > maxDistance * s) drawn = false;
          }
          if (drawn) {
            present[slot] = 1;
            presentCount += 1;
            let inside = true;
            for (let i = 0; i < 24; i += 4) {
              if (p[i] * x + p[i + 1] * y + p[i + 2] * z + p[i + 3] < -r) {
                inside = false;
                break;
              }
            }
            if (inside && occlude !== null) {
              Vec3.set(occlusionSphere.center, x, y, z);
              occlusionSphere.radius = r;
              if (occlude(occlusionSphere)) inside = false;
            }
            if (inside) {
              visible[slot] = 1;
              visibleCount += 1;
            }
          }
        }
        if (was !== visible[slot]) changed = true;
      }
      v.visibleCount = visibleCount;
      v.presentCount = presentCount;
      if (changed) v.changedAt = frameCount;
      return v;
    };

    /** the instances and levels of detail of `r` to draw; assumes `isVisible` returned true (`boundingSphere`) */
    const cullObject = (r: GraphicsRenderable) => {
      if (clipGhost) updateGhostCull();
      const view = clipGhost ? ghostLodView : lodView;
      const lodLevels = (drawFacts & DrawFact.LodLevels) !== 0;
      const hasInstanceGrid = (drawFacts & DrawFact.InstanceGrid) !== 0;
      const hasMultipleInstances = (drawFacts & DrawFact.MultipleInstances) !== 0;
      if (hasInstanceGrid && (hasMultipleInstances || lodLevels)) {
        if (!viewChecked) {
          updateViewVersion();
          viewChecked = true;
        }
        // the marking mask drawing the marked instances alone needs no culling
        if (drawingMarked && r.prepareMarked?.(cullPlane, cullFrustum, lodView)) return;
        // in the unscaled world (`cullFrustum`); the occlusion pyramid's spheres scaled: not tested scaled
        const occluded = modelScale === 1 ? isOccluded : null;
        r.cull(cullPlane, cullFrustum, occluded, ctx.stats, currentFrame, view, clipGhost);
      } else {
        r.uncull(view, clipGhost);
      }
    };

    /**
     * assumes `isVisible` returned true, which also set `boundingSphere`; `segments`: of an object drawing several, those
     * to draw (`Renderable.renderSegments`), all where empty
     */
    const drawObject = (
      r: GraphicsRenderable,
      variant: GraphicsRenderVariant,
      flag: Flag,
      passState: RenderState,
      segments?: readonly number[],
    ) => {
      // culled all the same: the views of a frame culled already draw (and raster) by its lists (`getRasterInput`)
      cullObject(r);
      // the raster draws all its objects' spheres (`sphereRaster` 'only')
      if (culling.isRasterDrawn(r)) return;
      if (drawingSecond && !r.renderSecond) return;
      const state = getObjectState(drawFacts, variant, flag, passState);
      if (!state) return;
      ctx.setRenderState(state);

      if (drawingSecond) r.renderSecond!(variant, globals);
      else if (segments && segments.length > 0) r.renderSegments!(variant, globals, segments);
      else if (drawingMarked && r.renderMarked) r.renderMarked(variant, globals);
      else r.render(variant, globals);
    };

    const renderObject = (
      r: GraphicsRenderable,
      variant: GraphicsRenderVariant,
      flag: Flag,
      passState: RenderState,
    ) => {
      if (isVisible(r, variant)) drawObject(r, variant, flag, passState);
    };

    const {
      supported: solidInteriorCapSupported,
      canRender: canRenderSolidInteriorCaps,
      has: hasSolidInteriorCap,
      get: getSolidInteriorCaps,
      segments: capSegments,
      draw: drawSolidInteriorCaps,
      drawAll: renderSolidInteriorCaps,
    } = createSolidInteriorCaps({
      ctx,
      globals,
      boundingSphere,
      cameraPlane,
      cameraPosition,
      get modelScale() {
        return modelScale;
      },
      get near() {
        return globalUniforms.uNear.ref.value;
      },
      get cullSecondPending() {
        return culling.secondPending;
      },
      drawObject,
      isVisible,
    });

    /** drawing the opaque pass's second part (`drawObject`) */
    let drawingSecond = false;
    /** drawing the marking mask: objects draw their marked instances (`Renderable.renderMarked`) */
    let drawingMarked = false;
    const culling = createGpuCulling({
      ctx,
      props: p,
      globals,
      globalUniforms,
      viewport,
      cameraPlane: cullPlane,
      frustum: cullFrustum,
      lodView,
      get frame() {
        return currentFrame;
      },
    });
    const update = (camera: RenderCamera, scene: Scene, frame: Frame, side = false) => {
      if (!side && frame !== countedFrame) {
        frameCount += 1;
        countedFrame = frame;
      }
      currentFrame = frame;
      culling.setView(camera);

      updateArray(globalUniforms.uView, camera.view);
      updateArray(globalUniforms.uInvView, Mat4.invert(invView, camera.view));
      updateArray(globalUniforms.uProjection, camera.projection);
      updateArray(globalUniforms.uInvProjection, Mat4.invert(invProjection, camera.projection));

      ValueCell.updateIfChanged(globalUniforms.uIsOrtho, camera.state.mode === 'orthographic' ? 1 : 0);
      updateArray(
        globalUniforms.uViewOffset,
        camera.viewOffset.enabled
          ? Vec2.set(viewOffset, camera.viewOffset.offsetX * 16, camera.viewOffset.offsetY * 16)
          : Vec2.set(viewOffset, 0, 0),
      );
      ValueCell.updateIfChanged(globalUniforms.uModelScale, camera.scale);

      updateArray(globalUniforms.uCameraPosition, Mat4.getTranslation(cameraPosition, invView));
      Vec3.scale(cameraTarget, camera.state.target, camera.scale);
      Vec3.normalize(cameraDir, Vec3.sub(cameraDir, cameraTarget, cameraPosition));
      updateArray(globalUniforms.uCameraDir, cameraDir);

      ValueCell.updateIfChanged(globalUniforms.uFar, camera.far);
      ValueCell.updateIfChanged(globalUniforms.uNear, camera.near);
      ValueCell.updateIfChanged(globalUniforms.uFog, camera.state.fog > 0);
      ValueCell.updateIfChanged(globalUniforms.uFogFar, camera.fogFar);
      ValueCell.updateIfChanged(globalUniforms.uFogNear, camera.fogNear);
      ValueCell.updateIfChanged(globalUniforms.uTransparentBackground, transparentBackground);

      // views drawn together culled as one (stereo's eyes): their joint frustum and camera plane (`RenderCamera.cullView`)
      const { cullView } = camera;
      if (cullView) Frustum3D.copy(frustum, cullView.frustum);
      else Frustum3D.fromProjectionMatrix(frustum, camera.projectionView);
      ValueCell.updateIfChanged(globalUniforms.uCssHeight, camera.viewport.height / ctx.pixelRatio);
      // levels of detail by size on screen: in CSS pixels, of the camera's viewport (not a pass's target)
      LodView.set(
        lodView,
        camera.projection,
        camera.viewport.height / ctx.pixelRatio,
        camera.state.mode === 'orthographic',
      );

      if (cullView) {
        Plane3D.copy(cameraPlane, cullView.cameraPlane);
      } else {
        Plane3D.copy(cameraPlane, frustum[Frustum3D.PlaneIndex.Near]);
        cameraPlane.constant -= Plane3D.distanceToPoint(cameraPlane, cameraPosition);
      }
      Plane3D.copy(cullPlane, cameraPlane);
      Frustum3D.copy(cullFrustum, frustum);
      if (camera.scale !== 1) {
        const s = camera.scale;
        cullPlane.constant /= s;
        for (let i = 0; i < 6; ++i) cullFrustum[i as Frustum3D.PlaneIndex].constant /= s;
      }
      // compared at the first cull after this update: updates of other cameras between (picking) don't cull
      viewChecked = false;
      setViewPlanes();
      updateArray(globalUniforms.uCameraPlane, Plane3D.toArray(cameraPlane, globalUniforms.uCameraPlane.ref.value, 0));

      ValueCell.updateIfChanged(globalUniforms.uMarkerAverage, scene.markerAverage);

      const hasHeadRotation = !Mat4.isZero(camera.headRotation);
      if (hasHeadRotation) {
        ValueCell.updateIfChanged(globalUniforms.uHasHeadRotation, true);
        updateArray(globalUniforms.uInvHeadRotation, Mat4.invert(invHeadRotation, camera.headRotation));
        updateArray(globalUniforms.uLightDirection, getTransformedLightDirection(light, invHeadRotation));
      } else {
        ValueCell.updateIfChanged(globalUniforms.uHasHeadRotation, false);
        ValueCell.updateIfChanged(globalUniforms.uInvHeadRotation, Mat4.id);
        updateArray(globalUniforms.uLightDirection, light.direction);
      }
    };

    const splatModelView = Mat4();
    /**
     * The gaussian splats back to front for `camera`, before any pass: a draw's order (`tSplatOrder`), sorted again where
     * the view turned (the model's uniform scale changes no order)
     */
    /** the GPU's sort where compute shaders run (WebGPU), once its kernels compiled (in the background) */
    let splatSorter: SplatOrderGpu | undefined;
    let splatKernelsPreparing = false;
    const getSplatSorter = () => {
      const { compute } = ctx;
      if (splatSorter || !compute || !isSplatSortSupported(compute)) return splatSorter;
      if (compute.hasKernels(SplatSortKernels)) {
        splatSorter = { compute, sorter: createSplatSorter(compute) };
      } else if (!splatKernelsPreparing) {
        splatKernelsPreparing = true;
        compute.prepareKernels(SplatSortKernels).finally(() => (splatKernelsPreparing = false));
      }
      return splatSorter;
    };

    const sortGaussianSplats = (scene: Scene, camera: RenderCamera) => {
      Mat4.mul(splatModelView, camera.view, scene.primitives.view);
      const gpu = getSplatSorter();
      const { renderables } = scene.primitives;
      for (let i = 0, il = renderables.length; i < il; ++i) {
        const r = renderables[i];
        if (r.state.disposed || !r.state.visible || r.values.dGeometryType.ref.value !== 'gaussianSplats') continue;
        if (updateSplatOrder(r.values as unknown as GaussianSplatsValues, splatModelView, gpu)) r.update();
      }
    };

    const updateInternal = (
      group: Scene.Group,
      camera: RenderCamera,
      depthTexture: Texture | null,
      depthAlphaTexture: Texture | null,
      renderMask: Mask,
      markingDepthTest: boolean,
    ) => {
      arrayMapUpsert(sharedTexturesList, 'tDepth', depthTexture || emptyDepthTexture);
      arrayMapUpsert(sharedTexturesList, 'tDepthAlpha', depthAlphaTexture || emptyDepthAlphaTexture);

      modelScale = camera.scale;

      updateArray(globalUniforms.uModel, Mat4.scaleUniformly(model, group.view, modelScale));
      updateArray(globalUniforms.uModelView, Mat4.mul(modelView, camera.view, model));
      updateArray(globalUniforms.uInvModelView, Mat4.invert(invModelView, modelView));
      updateArray(globalUniforms.uModelViewProjection, Mat4.mul(modelViewProjection, modelView, camera.projection));
      updateArray(globalUniforms.uInvModelViewProjection, Mat4.invert(invModelViewProjection, modelViewProjection));

      ValueCell.updateIfChanged(globalUniforms.uRenderMask, renderMask);
      ValueCell.updateIfChanged(globalUniforms.uMarkingDepthTest, markingDepthTest);

      const hasEyeCamera = !Mat4.isZero(camera.viewEye);
      if (hasEyeCamera) {
        ValueCell.updateIfChanged(globalUniforms.uHasEyeCamera, true);
        updateArray(globalUniforms.uModelViewEye, Mat4.mul(modelViewEye, camera.viewEye, model));
        updateArray(globalUniforms.uInvModelViewEye, Mat4.invert(invModelViewEye, modelViewEye));
      } else {
        ValueCell.updateIfChanged(globalUniforms.uHasEyeCamera, false);
        ValueCell.updateIfChanged(globalUniforms.uModelViewEye, Mat4.id);
        ValueCell.updateIfChanged(globalUniforms.uInvModelViewEye, Mat4.id);
      }

      const { x, y, width, height } = viewport;
      ctx.setViewport(x, y, width, height);
      ctx.setScissor(x, y, width, height);

      globals.invalidate();
    };

    /** draws the objects of `group` as `pass` describes: from its table's rows where it has one */
    const renderPass = (pass: PassDraws, group: Scene.Group) => {
      if (group.table) renderPassRows(pass, group.table);
      else renderPassObjects(pass, group);
    };

    /** `renderPass` from a table's rows (of `part`): visibility per view, the checks on bits */
    const renderPassRows = (pass: PassDraws, table: DrawTable, part: PassPart = WholePass) => {
      const { variant, select, check, draws, caps } = pass;
      const visible = getRowVisibility(table);
      const inlineCaps = !!caps?.inline && canRenderSolidInteriorCaps();
      const { order, bits, objects } = table;
      const begin = part.begin;
      const end = Math.min(part.end, table.count);
      for (let k = part.objects ? begin : end; k < end; ++k) {
        const slot = order[k];
        if (!visible[slot]) continue;
        const o = slot * DrawRow.Stride;
        const facts = bits[o + DrawRow.Facts];
        const state = bits[o + DrawRow.State];
        if (!isPassDrawn(pass, state) || (select && !select(facts, state)) || !check(facts, state)) continue;
        const r = objects[slot]!;
        setRowBoundingSphere(table, slot);
        const objectDraws = draws(facts);
        for (let j = 0, jl = objectDraws.length; j < jl; ++j) {
          drawObject(r, variant, objectDraws[j].flag, objectDraws[j].passState);
        }
        if (inlineCaps && hasSolidInteriorCap(facts) && caps!.check(facts, state)) {
          drawSolidInteriorCaps(r, variant, caps!.mode, caps!.passState);
        }
      }
      if (part.caps && caps && !caps.inline && canRenderSolidInteriorCaps()) {
        for (let k = begin; k < end; ++k) {
          const slot = order[k];
          if (!visible[slot]) continue;
          const o = slot * DrawRow.Stride;
          const facts = bits[o + DrawRow.Facts];
          const state = bits[o + DrawRow.State];
          if (!isPassDrawn(pass, state) || !hasSolidInteriorCap(facts) || !caps.check(facts, state)) continue;
          const r = objects[slot]!;
          setRowBoundingSphere(table, slot);
          drawSolidInteriorCaps(r, variant, caps.mode, caps.passState);
        }
      }
    };

    /** `renderPass` for groups without a table (not a scene's) */
    const renderPassObjects = (pass: PassDraws, group: Scene.Group) => {
      const { variant, select, check, draws, caps } = pass;
      const inlineCaps = !!caps?.inline && canRenderSolidInteriorCaps();
      const { renderables } = group;
      for (let i = 0, il = renderables.length; i < il; ++i) {
        const r = renderables[i];
        const facts = getObjectFacts(r);
        const state = getDrawState(r.state);
        if (!isPassDrawn(pass, state) || (select && !select(facts, state)) || !check(facts, state)) continue;
        const objectDraws = draws(facts);
        for (let j = 0, jl = objectDraws.length; j < jl; ++j) {
          renderObject(r, variant, objectDraws[j].flag, objectDraws[j].passState);
        }
        if (inlineCaps && hasSolidInteriorCap(facts) && caps!.check(facts, state) && isVisible(r, variant)) {
          drawSolidInteriorCaps(r, variant, caps!.mode, caps!.passState);
        }
      }
      if (caps && !caps.inline) renderSolidInteriorCaps(renderables, caps.check, variant, caps.mode, caps.passState);
    };

    const drawLists = createRendererDrawLists({
      ctx,
      globals,
      viewPlanes,
      capSegments,
      get modelScale() {
        return modelScale;
      },
      get near() {
        return globalUniforms.uNear.ref.value;
      },
      get cullRecords() {
        return culling.records;
      },
      get gpuCulledCount() {
        return culling.culledCount;
      },
      get occlusionTest() {
        return isOccluded !== null;
      },
      get drawingMarked() {
        return drawingMarked;
      },
      renderPassRows,
      getRowVisibility,
      setRowBoundingSphere,
      cullObject,
      isRasterDrawn: culling.isRasterDrawn,
      isCapsDrawn: () => solidInteriorCapSupported && canRenderSolidInteriorCaps(),
      hasSolidInteriorCap,
      getSolidInteriorCaps,
      isSameView,
      setDrawingFromList: (value) => {
        drawingFromList = value;
      },
    });
    const collect = (pass: Pass, group: Scene.Group, target: TargetDescriptor, options?: { all?: boolean }) => {
      return collectFrom(pass, group.renderables, target, !!options?.all);
    };

    /**
     * stand-ins of the signatures' objects, with the values the pass draws and object states read and those the passes
     * check (alpha, transparency, emission): one per way the object may draw (`ObjectSignature.opaque`, ...), so that the
     * passes' own checks select what to request
     */
    const signatureStandIns = new WeakMap<ProgramItem, Map<string, GraphicsRenderable[]>>();
    const getSignatureStandIns = (signature: ObjectSignature, item: ProgramItem) => {
      const { opaque, transparent, emissive, gpuCulled } = signature;
      const flags = `${opaque}|${transparent}|${emissive}|${gpuCulled}`;
      let cache = signatureStandIns.get(item);
      if (!cache) signatureStandIns.set(item, (cache = new Map()));
      let standIns = cache.get(flags);
      if (!standIns) {
        const { defines } = signature;
        const fixed: { [k: string]: unknown } = {
          dGeometryType: defines.dGeometryType,
          dSolidInterior: !!defines.dSolidInterior,
          dFlipSided: !!defines.dFlipSided,
          dTransparentBackfaces: defines.dTransparentBackfaces,
          dXrayShaded: defines.dXrayShaded,
          dPointStyle: defines.dPointStyle,
          uDoubleSided: signature.doubleSided,
          hasReflection: signature.hasReflection,
          transparencyAverage: 0,
          markerAverage: 0,
          emissiveAverage: 0,
          uEmissive: emissive ? 1 : 0,
        };
        // only the values of the kind: the object state reads a missing one apart from a false one (`uDoubleSided`)
        const schema = getRenderableProgramSource(signature)?.schema;
        const alphas = [...(opaque ? [1] : []), ...(transparent ? [0.5] : [])];
        standIns = (alphas.length ? alphas : [0]).map((alpha) => {
          const values: { [k: string]: ValueCell<unknown> } = {
            alpha: ValueCell.create(alpha),
            uAlpha: ValueCell.create(alpha),
          };
          for (const k of Object.keys(fixed)) if (!schema || schema[k]) values[k] = ValueCell.create(fixed[k]);
          return {
            values,
            renderItem: item,
            state: { disposed: false, colorOnly: false },
            // whether the sphere raster's `auto` may draw it (`isSphereRasterCandidate`)
            gpuCulled,
          } as unknown as GraphicsRenderable;
        });
        cache.set(flags, standIns);
      }
      return standIns;
    };

    const collectSignatures = (
      pass: Pass,
      items: readonly { readonly signature: ObjectSignature; readonly item: ProgramItem }[],
      target: TargetDescriptor,
      options?: { all?: boolean },
    ) => {
      const standIns = items.flatMap(({ signature, item }) => getSignatureStandIns(signature, item));
      return collectFrom(pass, standIns, target, !!options?.all);
    };

    const collectFrom = (
      pass: Pass,
      renderables: readonly GraphicsRenderable[],
      target: TargetDescriptor,
      all: boolean,
    ) => {
      const requests: PipelineRequest[] = [];
      const { variant, select, check, draws, caps } = getPassDraws(pass, transparentBackground);
      const capsDrawn = solidInteriorCapSupported && target.hasStencil;
      const request = (r: GraphicsRenderable, flag: Flag, passState: RenderState) => {
        const state = getObjectState(getObjectFacts(r), variant, flag, passState);
        if (state) requests.push({ item: r.renderItem, variant, state, target });
      };
      // the opaque pass with the raster on: the resolves of its objects (`renderOpaque`), drawn before the impostors
      const resolves = pass.name === 'opaque' && culling.isRasterAvailable();
      for (const r of renderables) {
        const facts = getObjectFacts(r);
        const state = getDrawState(r.state);
        if (r.state.disposed || (select && !select(facts, state)) || (!all && !check(facts, state))) continue;
        for (const { flag, passState } of draws(facts)) request(r, flag, passState);
        if (resolves && culling.isRasterCandidate(r)) {
          requests.push({ item: r.renderItem, variant: 'color', state: OpaqueState, target, define: 'dSphereResolve' });
        }
        if (caps && capsDrawn && hasSolidInteriorCap(facts) && (all || caps.check(facts, state))) {
          request(r, Flag.SolidInteriorMark, getSolidInteriorMarkState(caps.passState));
          const { resolve, fill } = getSolidInteriorFillStates(caps.passState, caps.mode);
          if (resolve) request(r, Flag.SolidInteriorFill, resolve);
          request(r, Flag.SolidInteriorFill, fill);
        }
      }
      return requests;
    };

    /** `pass` replayed or recorded with the draw lists (`drawLists`) where it can be, drawn directly otherwise */
    const drawPass = (pass: PassDraws, group: Scene.Group) => {
      if (p.drawLists && drawLists.render(pass, group)) return;
      renderPass(pass, group);
    };

    const setOccluder = (occluder: Renderer.Occluder | null | undefined) => {
      if (!occluder) return null;
      ValueCell.updateIfChanged(globalUniforms.uOccluderNear, occluder.near);
      ValueCell.updateIfChanged(globalUniforms.uOccluderFar, occluder.far);
      return occluder.depth;
    };

    const renderPick = (
      group: Scene.Group,
      camera: RenderCamera,
      variant: 'pick' | 'depth',
      pickType: PickType,
      occluder?: Renderer.Occluder | null,
    ) => {
      updateInternal(group, camera, setOccluder(occluder), occluder?.depthAlpha ?? null, Mask.All, false);
      ValueCell.updateIfChanged(globalUniforms.uPickType, pickType);
      drawPass(variant === 'pick' ? PassDraws.pick : PassDraws.pickDepth, group);
    };

    const renderDepth = (group: Scene.Group, camera: RenderCamera) => {
      updateInternal(group, camera, null, null, Mask.All, false);
      drawPass(PassDraws.depth, group);
    };

    const renderDepthOpaque = (group: Scene.Group, camera: RenderCamera) => {
      updateInternal(group, camera, null, null, Mask.Opaque, false);
      drawPass(PassDraws.depthOpaque, group);
    };

    const renderDepthOpaqueBack = (group: Scene.Group, camera: RenderCamera) => {
      updateInternal(group, camera, null, null, Mask.Opaque, false);
      ValueCell.updateIfChanged(globalUniforms.uDepthBack, true);
      renderPass(PassDraws.depthOpaqueBack, group);
      ValueCell.updateIfChanged(globalUniforms.uDepthBack, false);
    };

    const renderDepthTransparent = (group: Scene.Group, camera: RenderCamera, depthTexture: Texture) => {
      updateInternal(group, camera, depthTexture, null, Mask.Transparent, false);
      drawPass(PassDraws.depthTransparent, group);
    };

    const renderClipGhost = (group: Scene.Group, camera: RenderCamera, scale = 1) => {
      updateInternal(group, camera, null, null, Mask.All, false);
      ValueCell.updateIfChanged(globalUniforms.uClipGhost, true);
      const pixelRatio = globalUniforms.uPixelRatio.ref.value;
      ValueCell.updateIfChanged(globalUniforms.uPixelRatio, pixelRatio * scale);
      clipGhost = true;
      try {
        // not from draw lists: their recordings follow the visible pass's rows
        renderPass(PassDraws.clipGhost, group);
      } finally {
        clipGhost = false;
        ValueCell.updateIfChanged(globalUniforms.uClipGhost, false);
        ValueCell.updateIfChanged(globalUniforms.uPixelRatio, pixelRatio);
      }
    };

    const renderMarkingDepth = (group: Scene.Group, camera: RenderCamera) => {
      updateInternal(group, camera, null, null, Mask.All, false);
      ValueCell.updateIfChanged(globalUniforms.uMarkingType, MarkingType.Depth);
      drawPass(PassDraws.markingDepth, group);
    };

    const renderMarkingMask = (group: Scene.Group, camera: RenderCamera, depthTexture: Texture | null) => {
      updateInternal(group, camera, null, depthTexture, Mask.All, !!depthTexture);
      ValueCell.updateIfChanged(globalUniforms.uMarkingType, MarkingType.Mask);
      drawingMarked = true;
      drawPass(PassDraws.markingMask, group);
      drawingMarked = false;
    };

    const renderEmissiveOpaque = (
      group: Scene.Group,
      camera: RenderCamera,
      depthTexture: Texture,
      occludeWithOpaqueDepth: boolean,
    ) => {
      // tDepthAlpha = front transparent depth+alpha; emissive shader dims emitters behind it
      updateInternal(group, camera, null, depthTexture, Mask.Opaque, false);
      drawPass(occludeWithOpaqueDepth ? PassDraws.emissiveOpaqueOccluded : PassDraws.emissiveOpaque, group);
    };

    const renderEmissiveTransparent = (group: Scene.Group, camera: RenderCamera, depthTexture: Texture) => {
      // tDepthAlpha = front transparent depth+alpha; emissive shader dims emitters behind it
      updateInternal(group, camera, null, depthTexture, Mask.Transparent, false);
      drawPass(PassDraws.emissiveTransparent, group);
    };

    const renderTracing = (group: Scene.Group, camera: RenderCamera) => {
      updateInternal(group, camera, null, null, Mask.Opaque, false);
      renderPass(PassDraws.tracing, group);
    };

    const renderTracingSecond = (group: Scene.Group, camera: RenderCamera) => {
      updateInternal(group, camera, null, null, Mask.Opaque, false);
      // the pyramid's passes set theirs
      ctx.setViewport(viewport.x, viewport.y, viewport.width, viewport.height);
      ctx.setScissor(viewport.x, viewport.y, viewport.width, viewport.height);
      drawingSecond = true;
      renderPass(PassDraws.tracingSecond, group);
      drawingSecond = false;
      culling.secondDone();
    };

    const renderBlended = (scene: Scene, camera: RenderCamera) => {
      if (scene.hasOpaque) {
        renderOpaque(scene, camera);
      }
      if (scene.opacityAverage < 1) {
        renderBlendedTransparent(scene, camera);
      }
      if (scene.hasGaussianSplats) {
        renderGaussianSplats(scene, camera);
      }
    };

    const renderOverlay = (group: Scene.Group, camera: RenderCamera, occluder: Renderer.Occluder | null) => {
      const depth = setOccluder(occluder);
      const depthAlpha = occluder?.depthAlpha ?? null;
      updateInternal(group, camera, depth, depthAlpha, Mask.Opaque, false);
      drawPass(PassDraws.opaque, group);
      updateInternal(group, camera, depth, depthAlpha, Mask.Transparent, false);
      drawPass(transparentBackground ? PassDraws.blendedTransparentBackground : PassDraws.blendedTransparent, group);
    };

    const renderOpaque = (group: Scene.Group, camera: RenderCamera) => {
      updateInternal(group, camera, null, null, Mask.Opaque, false);
      // the raster's spheres first, the impostors leave them out (`prepareSphereRaster`)
      culling.resolveRaster();
      // the objects' raster values change per view (set by `prepareSphereRaster`, reset by `finish`): their draws are
      // walked and their data flushed before a recording replays (`flushDataChanges`)
      culling.beginOpaque();
      drawPass(PassDraws.opaque, group);
      culling.endOpaque();
      // the second part draws them too
      if (culling.secondPending === 0) culling.finishRaster();
    };

    const renderOpaqueSecond = (group: Scene.Group, camera: RenderCamera) => {
      updateInternal(group, camera, null, null, Mask.Opaque, false);
      // the pyramid's passes set theirs
      ctx.setViewport(viewport.x, viewport.y, viewport.width, viewport.height);
      ctx.setScissor(viewport.x, viewport.y, viewport.width, viewport.height);
      culling.resolveRaster();
      culling.beginOpaque();
      drawingSecond = true;
      renderPass(PassDraws.opaqueSecond, group);
      drawingSecond = false;
      culling.endOpaque();
      culling.finishRaster();
      culling.secondDone();
    };

    const renderBlendedTransparent = (group: Scene.Group, camera: RenderCamera) => {
      updateInternal(group, camera, null, null, Mask.Transparent, false);
      drawPass(transparentBackground ? PassDraws.blendedTransparentBackground : PassDraws.blendedTransparent, group);
    };

    const renderGaussianSplats = (group: Scene.Group, camera: RenderCamera) => {
      updateInternal(group, camera, null, null, Mask.Transparent, false);
      drawPass(transparentBackground ? PassDraws.gaussianSplatsBackground : PassDraws.gaussianSplats, group);
    };

    const renderVolume = (group: Scene.Group, camera: RenderCamera, depthTexture: Texture) => {
      updateInternal(group, camera, depthTexture, null, Mask.Transparent, false);
      drawPass(PassDraws.volume, group);
    };

    const renderWboitTransparent = (group: Scene.Group, camera: RenderCamera, depthTexture: Texture | null) => {
      updateInternal(group, camera, depthTexture, null, Mask.Transparent, false);
      drawPass(PassDraws.wboitTransparent, group);
    };

    const renderDpoitTransparent = (
      group: Scene.Group,
      camera: RenderCamera,
      depthTexture: Texture,
      dpoitTextures: { depth: Texture; frontColor: Texture },
    ) => {
      arrayMapUpsert(sharedTexturesList, 'tDpoitDepth', dpoitTextures.depth);
      arrayMapUpsert(sharedTexturesList, 'tDpoitFrontColor', dpoitTextures.frontColor);

      updateInternal(group, camera, depthTexture, null, Mask.Transparent, false);
      drawPass(PassDraws.dpoitTransparent, group);
    };

    return {
      clear: (toBackgroundColor: boolean, ignoreTransparentBackground?: boolean, forceToTransparency?: boolean) => {
        ctx.setScissor(viewport.x, viewport.y, viewport.width, viewport.height);
        if (forceToTransparency || (transparentBackground && !ignoreTransparentBackground)) {
          ctx.clear({ color: [0, 0, 0, 0], depth: 0 });
        } else if (toBackgroundColor) {
          ctx.clear({ color: [bgColor[0], bgColor[1], bgColor[2], 1], depth: 0 });
        } else {
          ctx.clear({ color: [1, 1, 1, 1], depth: 0 });
        }
      },
      clearDepth: (packed = false) => {
        ctx.setScissor(viewport.x, viewport.y, viewport.width, viewport.height);
        // the far plane: depth 0, packed as zeros
        if (packed) {
          ctx.clear({ color: [0, 0, 0, 0], depth: 0 });
        } else {
          ctx.clear({ depth: 0 });
        }
      },
      update,
      prepareCull: culling.prepare,
      prepareGaussianSplats: (scene: Scene, camera: RenderCamera) => {
        if (scene.hasGaussianSplats) sortGaussianSplats(scene, camera);
      },
      prepareCullView: culling.prepareView,
      prepareSphereRaster: culling.prepareRaster,
      get cullSecondPending() {
        return culling.secondPending;
      },
      prepareCullSecond: culling.prepareSecond,
      renderOpaqueSecond,
      get gpuCulledCount() {
        return culling.culledCount;
      },

      renderPick,
      renderDepth,
      renderDepthOpaque,
      renderDepthOpaqueBack,
      renderClipGhost,
      renderDepthTransparent,
      renderMarkingDepth,
      renderMarkingMask,
      renderEmissiveOpaque,
      renderEmissiveTransparent,
      renderTracing,
      renderTracingSecond,
      renderBlended,
      renderOverlay,
      renderOpaque,
      renderBlendedTransparent,
      renderGaussianSplats,
      renderVolume,
      renderWboitTransparent,
      renderDpoitTransparent,

      collect,
      collectSignatures,

      setTime: (time: number) => {
        ValueCell.updateIfChanged(globalUniforms.uTime, time);
      },
      setProps: (props: Partial<RendererProps>) => {
        if (props.backgroundColor !== undefined && props.backgroundColor !== p.backgroundColor) {
          p.backgroundColor = props.backgroundColor;
          Color.toVec3Normalized(bgColor, p.backgroundColor);
          ValueCell.update(globalUniforms.uFogColor, Vec3.copy(globalUniforms.uFogColor.ref.value, bgColor));
        }

        if (props.pickingAlphaThreshold !== undefined && props.pickingAlphaThreshold !== p.pickingAlphaThreshold) {
          p.pickingAlphaThreshold = props.pickingAlphaThreshold;
          ValueCell.update(globalUniforms.uPickingAlphaThreshold, p.pickingAlphaThreshold);
        }

        if (props.drawLists !== undefined && props.drawLists !== p.drawLists) {
          p.drawLists = props.drawLists;
          drawLists.reset();
        }
        if (props.sphereRaster !== undefined && props.sphereRaster !== p.sphereRaster) {
          p.sphereRaster = props.sphereRaster;
          if (p.sphereRaster === 'off') culling.rasterOff();
        }
        if (props.sphereRasterRadius !== undefined) p.sphereRasterRadius = props.sphereRasterRadius;
        if (props.sphereRasterObjects !== undefined) p.sphereRasterObjects = props.sphereRasterObjects;

        if (props.colorMarker !== undefined && props.colorMarker !== p.colorMarker) {
          p.colorMarker = props.colorMarker;
        }
        if (props.highlightColor !== undefined && props.highlightColor !== p.highlightColor) {
          p.highlightColor = props.highlightColor;
          ValueCell.update(
            globalUniforms.uHighlightColor,
            Color.toVec3Normalized(globalUniforms.uHighlightColor.ref.value, p.highlightColor),
          );
        }
        if (props.selectColor !== undefined && props.selectColor !== p.selectColor) {
          p.selectColor = props.selectColor;
          ValueCell.update(
            globalUniforms.uSelectColor,
            Color.toVec3Normalized(globalUniforms.uSelectColor.ref.value, p.selectColor),
          );
        }
        if (props.dimColor !== undefined && props.dimColor !== p.dimColor) {
          p.dimColor = props.dimColor;
          ValueCell.update(
            globalUniforms.uDimColor,
            Color.toVec3Normalized(globalUniforms.uDimColor.ref.value, p.dimColor),
          );
        }
        if (props.highlightStrength !== undefined && props.highlightStrength !== p.highlightStrength) {
          p.highlightStrength = props.highlightStrength;
          ValueCell.update(globalUniforms.uHighlightStrength, p.highlightStrength);
        }
        if (props.selectStrength !== undefined && props.selectStrength !== p.selectStrength) {
          p.selectStrength = props.selectStrength;
          ValueCell.update(globalUniforms.uSelectStrength, p.selectStrength);
        }
        if (props.dimStrength !== undefined && props.dimStrength !== p.dimStrength) {
          p.dimStrength = props.dimStrength;
          ValueCell.update(globalUniforms.uDimStrength, p.dimStrength);
        }
        if (props.markerPriority !== undefined && props.markerPriority !== p.markerPriority) {
          p.markerPriority = props.markerPriority;
          ValueCell.update(globalUniforms.uMarkerPriority, p.markerPriority);
        }

        if (props.xrayEdgeFalloff !== undefined && props.xrayEdgeFalloff !== p.xrayEdgeFalloff) {
          p.xrayEdgeFalloff = props.xrayEdgeFalloff;
          ValueCell.update(globalUniforms.uXrayEdgeFalloff, p.xrayEdgeFalloff);
        }

        if (props.celSteps !== undefined && props.celSteps !== p.celSteps) {
          p.celSteps = props.celSteps;
          ValueCell.update(globalUniforms.uCelSteps, p.celSteps);
        }

        if (props.exposure !== undefined && props.exposure !== p.exposure) {
          p.exposure = props.exposure;
          ValueCell.update(globalUniforms.uExposure, p.exposure);
        }

        if (props.light !== undefined && !deepEqual(props.light, p.light)) {
          p.light = props.light;
          Object.assign(light, getLight(props.light, light));
          ValueCell.update(globalUniforms.uLightDirection, light.direction);
          ValueCell.update(globalUniforms.uLightColor, light.color);
        }
        if (props.ambientColor !== undefined && props.ambientColor !== p.ambientColor) {
          p.ambientColor = props.ambientColor;
          Vec3.scale(ambientColor, Color.toArrayNormalized(p.ambientColor, ambientColor, 0), p.ambientIntensity);
          ValueCell.update(globalUniforms.uAmbientColor, ambientColor);
        }
        if (props.ambientIntensity !== undefined && props.ambientIntensity !== p.ambientIntensity) {
          p.ambientIntensity = props.ambientIntensity;
          Vec3.scale(ambientColor, Color.toArrayNormalized(p.ambientColor, ambientColor, 0), p.ambientIntensity);
          ValueCell.update(globalUniforms.uAmbientColor, ambientColor);
        }

        if (props.enableAnimation !== undefined && props.enableAnimation !== p.enableAnimation) {
          p.enableAnimation = props.enableAnimation;
          ValueCell.update(globalUniforms.uEnableAnimation, p.enableAnimation);
        }
      },
      setViewport: (x: number, y: number, width: number, height: number) => {
        ctx.setViewport(x, y, width, height);
        ctx.setScissor(x, y, width, height);
        if (x !== viewport.x || y !== viewport.y || width !== viewport.width || height !== viewport.height) {
          Viewport.set(viewport, x, y, width, height);
          ValueCell.update(globalUniforms.uViewport, Vec4.set(globalUniforms.uViewport.ref.value, x, y, width, height));
        }
      },
      setTransparentBackground: (value: boolean) => {
        transparentBackground = value;
      },
      setDrawingBufferSize: (width: number, height: number) => {
        if (width !== drawingBufferSize[0] || height !== drawingBufferSize[1]) {
          ValueCell.update(globalUniforms.uDrawingBufferSize, Vec2.set(drawingBufferSize, width, height));
        }
      },
      setPixelRatio: (value: number) => {
        ValueCell.update(globalUniforms.uPixelRatio, value);
      },
      setOcclusionTest: (f: ((s: Sphere3D) => boolean) | null) => {
        isOccluded = f;
      },
      get occlusionTest() {
        return isOccluded;
      },
      setOcclusionPyramid: culling.setOcclusionPyramid,

      props: p,
      get stats(): RendererStats {
        return {
          programCount: ctx.stats.resourceCounts.program,
          shaderCount: ctx.stats.resourceCounts.shader,

          attributeCount: ctx.stats.resourceCounts.attribute,
          elementsCount: ctx.stats.resourceCounts.elements,
          framebufferCount: ctx.stats.resourceCounts.framebuffer,
          renderbufferCount: ctx.stats.resourceCounts.renderbuffer,
          textureCount: ctx.stats.resourceCounts.texture,
          vertexArrayCount: ctx.stats.resourceCounts.vertexArray,

          drawCount: stats.drawCount,
          instanceCount: stats.instanceCount,
          instancedDrawCount: stats.instancedDrawCount,
        };
      },
      get light(): Light {
        return light;
      },
      get ambientColor(): Vec3 {
        return globalUniforms.uAmbientColor.ref.value;
      },
      dispose: () => {
        culling.destroy();
        splatSorter?.sorter.destroy();
        globals.destroy();
        if (emptyDepthTexture !== emptyDepthAlphaTexture) emptyDepthTexture.destroy();
        emptyDepthAlphaTexture.destroy();
        // TODO
      },
    };
  }
}

export { Renderer };
