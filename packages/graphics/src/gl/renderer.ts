/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import { Viewport } from '@molstar/graphics/canvas3d/camera/util';
import type { ICamera } from '@molstar/graphics/canvas3d/camera';
import type { Scene } from './scene.js';
import type { GraphicsContext } from './context.js';
import { Mat4, Vec3, Vec4, Vec2 } from '@molstar/core/math/linear-algebra';
import { type GraphicsRenderable, type Frame, createFrame } from './renderable.js';
import { Color } from '@molstar/core/util/color';
import { ValueCell, deepEqual } from '@molstar/core/util';
import { GlobalUniformSchema, type GlobalUniformValues, MaxLightCount } from './renderable/schema.js';
import type { GraphicsRenderVariant } from './render-item.js';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { degToRad } from '@molstar/core/math/misc';
import type { Texture, Textures } from './texture.js';
import { arrayMapUpsert } from '@molstar/core/util/array';
import { clamp } from '@molstar/core/math/interpolate';
import { Frustum3D } from '@molstar/core/math/geometry/primitives/frustum3d';
import { Plane3D } from '@molstar/core/math/geometry/primitives/plane3d';
import { Sphere3D } from '@molstar/core/math/geometry';
import { Clip } from '@molstar/core/util/clip';
import {
  BlendState,
  type CullMode,
  type DepthState,
  type FrontFace,
  RenderState,
  type StencilState,
} from './render-state.js';

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
  update: (camera: ICamera, scene: Scene, frame: Frame) => void;
  setTime: (time: number) => void;

  renderPick: (group: Scene.Group, camera: ICamera, variant: 'pick' | 'depth', pickType: PickType) => void;
  renderDepth: (group: Scene.Group, camera: ICamera) => void;
  renderDepthOpaque: (group: Scene.Group, camera: ICamera) => void;
  renderDepthOpaqueBack: (group: Scene.Group, camera: ICamera) => void;
  renderDepthTransparent: (group: Scene.Group, camera: ICamera, depthTexture: Texture) => void;
  renderMarkingDepth: (group: Scene.Group, camera: ICamera) => void;
  renderMarkingMask: (group: Scene.Group, camera: ICamera, depthTexture: Texture | null) => void;
  renderEmissiveOpaque: (
    group: Scene.Group,
    camera: ICamera,
    depthTexture: Texture,
    occludeWithOpaqueDepth: boolean,
  ) => void;
  renderEmissiveTransparent: (group: Scene.Group, camera: ICamera, depthTexture: Texture) => void;
  renderTracing: (group: Scene.Group, camera: ICamera) => void;
  renderBlended: (group: Scene, camera: ICamera) => void;
  renderOpaque: (group: Scene.Group, camera: ICamera) => void;
  renderBlendedTransparent: (group: Scene.Group, camera: ICamera) => void;
  renderVolume: (group: Scene.Group, camera: ICamera, depthTexture: Texture) => void;
  renderWboitTransparent: (group: Scene.Group, camera: ICamera, depthTexture: Texture) => void;
  renderDpoitTransparent: (
    group: Scene.Group,
    camera: ICamera,
    depthTexture: Texture,
    dpoitTextures: { depth: Texture; frontColor: Texture; backColor: Texture },
  ) => void;

  setProps: (props: Partial<RendererProps>) => void;
  setViewport: (x: number, y: number, width: number, height: number) => void;
  setTransparentBackground: (value: boolean) => void;
  setDrawingBufferSize: (width: number, height: number) => void;
  setPixelRatio: (value: number) => void;
  setOcclusionTest: (f: ((s: Sphere3D) => boolean) | null) => void;

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
};
export type RendererProps = PD.Values<typeof RendererParams>;

export type Light = {
  count: number;
  direction: number[];
  color: number[];
};

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
  const enum Flag {
    None = 0,
    BlendedFront = 1,
    BlendedBack = 2,
    DepthBack = 3,
    SolidInteriorMark = 4,
    SolidInteriorFill = 5,
  }

  type SolidInteriorMode = 'opaque' | 'back' | 'blended' | 'oit';

  const enum Mask {
    All = 0,
    Opaque = 1,
    Transparent = 2,
  }

  const DepthLess: DepthState = { compare: 'less', write: true };
  const DepthLessReadOnly: DepthState = { compare: 'less', write: false };

  const OpaqueState: RenderState = { depth: DepthLess };
  const OpaqueBackState: RenderState = { depth: { compare: 'greater', write: true } };
  /** test against caller-attached opaque depth, read-only + less-equal so an emitter at its own depth still passes */
  const EmissiveOccludedState: RenderState = { depth: { compare: 'less-equal', write: false } };
  /** max blending so overlapping faces don't accumulate */
  const EmissiveTransparentState: RenderState = { blend: BlendState.Max, depth: DepthLessReadOnly };
  const BlendedTransparentState: RenderState = { blend: BlendState.Alpha, depth: DepthLessReadOnly };
  const BlendedTransparentBackgroundState: RenderState = {
    blend: BlendState.Premultiplied,
    depth: DepthLessReadOnly,
  };
  /** depth test done manually in shader against `depthTexture` */
  const VolumeState: RenderState = RenderState.Premultiplied;
  const WboitState: RenderState = {
    blend: BlendState(
      BlendState.Component('add', 'one', 'one'),
      BlendState.Component('add', 'zero', 'one-minus-src-alpha'),
    ),
  };
  const DpoitState: RenderState = { blend: BlendState.Max };

  const SolidInteriorMarkStencil: StencilState = {
    front: { compare: 'always', failOp: 'keep', depthFailOp: 'keep', passOp: 'increment-wrap' },
    back: { compare: 'always', failOp: 'keep', depthFailOp: 'keep', passOp: 'decrement-wrap' },
    readMask: 0xff,
    writeMask: 0xff,
    reference: 0,
  };
  const SolidInteriorFillStencilKeep: StencilState = {
    front: { compare: 'not-equal', failOp: 'keep', depthFailOp: 'keep', passOp: 'keep' },
    back: { compare: 'not-equal', failOp: 'keep', depthFailOp: 'keep', passOp: 'keep' },
    readMask: 0xff,
    writeMask: 0xff,
    reference: 0,
  };
  const SolidInteriorFillStencilZero: StencilState = {
    front: { compare: 'not-equal', failOp: 'keep', depthFailOp: 'keep', passOp: 'zero' },
    back: { compare: 'not-equal', failOp: 'keep', depthFailOp: 'keep', passOp: 'zero' },
    readMask: 0xff,
    writeMask: 0xff,
    reference: 0,
  };
  const NoColorWrite = [false, false, false, false] as const;

  const CullModes: CullMode[] = ['none', 'front', 'back'];
  /** per-object cull and front face variants of the pass states */
  const FaceStates = new WeakMap<RenderState, RenderState[]>();
  function getFaceState(base: RenderState, cull: CullMode, frontFace: FrontFace) {
    let states = FaceStates.get(base);
    if (!states) {
      states = [];
      FaceStates.set(base, states);
    }
    const i = CullModes.indexOf(cull) * 2 + (frontFace === 'cw' ? 1 : 0);
    return (states[i] ??= { ...base, cull, frontFace });
  }

  /** 1x1 depth texture cleared to the far plane */
  function createFarDepthTexture(ctx: GraphicsContext) {
    const texture = ctx.resources.texture('image-depth', 'depth', 'ushort', 'nearest');
    texture.define(1, 1);
    const framebuffer = ctx.resources.framebuffer();
    texture.attachFramebuffer(framebuffer, 'depth');
    ctx.beginRenderPass(framebuffer, { clear: { depth: 1 } });
    ctx.endRenderPass();
    texture.detachFramebuffer(framebuffer, 'depth');
    framebuffer.destroy();
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

    // placeholders at the far plane, packed values of 1 decode to (almost) 1
    const emptyDepthAlphaTexture = ctx.resources.texture('image-uint8', 'rgba', 'ubyte', 'nearest');
    emptyDepthAlphaTexture.define(1, 1);
    emptyDepthAlphaTexture.load({ array: new Uint8Array([255, 255, 255, 255]), width: 1, height: 1 });
    const emptyDepthTexture = ctx.capabilities.depthTexture ? createFarDepthTexture(ctx) : emptyDepthAlphaTexture;
    const sharedTexturesList: Textures = [
      ['tDepth', emptyDepthTexture],
      ['tDepthAlpha', emptyDepthAlphaTexture],
    ];

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

    let modelScale = 1;
    const boundingSphere = Sphere3D();
    let currentFrame: Frame = createFrame();

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
      uIsAsymmetricProjection: ValueCell.create(false),

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
      uCelSteps: ValueCell.create(p.celSteps),
      uExposure: ValueCell.create(p.exposure),

      uTime: ValueCell.create(0),
      uEnableAnimation: ValueCell.create(p.enableAnimation),
    };
    const globals = ctx.createGlobalBindings(GlobalUniformSchema, Object.entries(globalUniforms), sharedTexturesList);

    /** also sets `boundingSphere` to the scaled bounding sphere of `r` */
    const isVisible = (r: GraphicsRenderable, variant: GraphicsRenderVariant) => {
      if (r.state.disposed || !r.state.visible || (!r.state.pickable && variant === 'pick')) {
        return false;
      }

      if (!r.values.drawCount.ref.value) {
        return false;
      }

      Sphere3D.scaleNX(boundingSphere, r.values.boundingSphere.ref.value, modelScale);

      if (!Frustum3D.intersectsSphere3D(frustum, boundingSphere)) {
        return false;
      }

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

    /** assumes `isVisible` returned true, which also set `boundingSphere` */
    const drawObject = (r: GraphicsRenderable, variant: GraphicsRenderVariant, flag: Flag, passState: RenderState) => {
      const unscaled = modelScale === 1;
      if (unscaled) {
        const hasInstanceGrid = r.values.instanceGrid.ref.value.cellSize > 0;
        const hasMultipleInstances = r.values.uInstanceCount.ref.value > 1;
        if (hasInstanceGrid && (hasMultipleInstances || r.values.lodLevels)) {
          r.cull(cameraPlane, frustum, isOccluded, ctx.stats, currentFrame);
        } else {
          r.uncull();
        }
      } else {
        if (r.values.lodLevels) {
          const { center, radius } = boundingSphere;
          const d = Plane3D.distanceToPoint(cameraPlane, center);
          r.cullSimple(d, radius, modelScale);
        } else {
          r.uncull();
        }
      }

      const flipSided = !!r.values.dFlipSided?.ref.value;
      let cull: CullMode;
      let frontFace: FrontFace = flipSided ? 'cw' : 'ccw';
      if (r.values.dGeometryType.ref.value === 'directVolume') {
        if (variant !== 'color') {
          return; // only color supported
        }

        // culling done in fragment shader
        cull = 'none';
        frontFace = 'ccw';
      } else if (flag === Flag.BlendedFront) {
        cull = flipSided ? 'front' : 'back';
      } else if (
        flag === Flag.DepthBack ||
        flag === Flag.SolidInteriorMark ||
        (flag === Flag.SolidInteriorFill && r.values.hasReflection.ref.value)
      ) {
        // reflected instances flip the winding, so their exit faces are not back-facing
        cull = 'none';
      } else if (flag === Flag.BlendedBack || flag === Flag.SolidInteriorFill) {
        cull = flipSided ? 'back' : 'front';
      } else if (r.values.uDoubleSided && !r.values.uDoubleSided.ref.value && !r.values.hasReflection.ref.value) {
        cull = flipSided ? 'front' : 'back';
      } else {
        cull = 'none';
      }
      ctx.setRenderState(getFaceState(passState, cull, frontFace));

      r.render(variant, globals);
    };

    const renderObject = (
      r: GraphicsRenderable,
      variant: GraphicsRenderVariant,
      flag: Flag,
      passState: RenderState,
    ) => {
      if (isVisible(r, variant)) drawObject(r, variant, flag, passState);
    };

    const solidInteriorCapSupported = ctx.capabilities.fragDepth;
    const canRenderSolidInteriorCaps = () => {
      return solidInteriorCapSupported && ctx.hasRenderPassStencil();
    };
    const hasSolidInteriorCap = (r: GraphicsRenderable) => {
      const geomType = r.values.dGeometryType.ref.value;
      return (
        solidInteriorCapSupported &&
        (geomType === 'mesh' || geomType === 'textureMesh') &&
        !!r.values.dSolidInterior?.ref.value
      );
    };

    const setSolidInteriorPass = (pass: number, clipIndex: number) => {
      ValueCell.updateIfChanged(globalUniforms.uSolidInteriorPass, pass);
      ValueCell.updateIfChanged(globalUniforms.uSolidInteriorClip, clipIndex);
      globals.invalidate();
    };

    const renderSolidInteriorMark = (
      r: GraphicsRenderable,
      variant: GraphicsRenderVariant,
      clipIndex: number,
      passState: RenderState,
    ) => {
      ctx.clear({ stencil: 0 });
      setSolidInteriorPass(2, clipIndex);
      const markState: RenderState = {
        ...passState,
        colorWrite: NoColorWrite,
        depth: undefined,
        stencil: SolidInteriorMarkStencil,
      };
      drawObject(r, variant, Flag.SolidInteriorMark, markState);
    };

    const renderSolidInteriorFill = (
      r: GraphicsRenderable,
      variant: GraphicsRenderVariant,
      mode: SolidInteriorMode,
      clipIndex: number,
      passState: RenderState,
    ) => {
      const writeDepth = mode === 'opaque' || mode === 'back';
      setSolidInteriorPass(1, clipIndex);
      let depth: DepthState;
      if (mode === 'oit') {
        // the OIT depth attachment is otherwise unused: resolve the nearest back face, then only draw that one
        ctx.clear({ depth: 1 });
        const resolveState: RenderState = {
          ...passState,
          colorWrite: NoColorWrite,
          depth: { compare: 'less-equal', write: true },
          stencil: SolidInteriorFillStencilKeep,
        };
        drawObject(r, variant, Flag.SolidInteriorFill, resolveState);
        depth = { compare: 'equal', write: false };
      } else {
        depth = { compare: mode === 'back' ? 'greater-equal' : 'less-equal', write: writeDepth };
      }
      const fillState: RenderState = {
        ...passState,
        depth,
        stencil: writeDepth ? SolidInteriorFillStencilKeep : SolidInteriorFillStencilZero,
      };
      drawObject(r, variant, Flag.SolidInteriorFill, fillState);
    };

    const solidInteriorEye = Vec3();

    const renderSolidInteriorCap = (
      r: GraphicsRenderable,
      variant: GraphicsRenderVariant,
      mode: SolidInteriorMode,
      passState: RenderState,
    ) => {
      if (!isVisible(r, variant)) return;

      const back = mode === 'back';
      const near = globalUniforms.uNear.ref.value * 1.0001;
      if (
        !back &&
        Math.abs(Plane3D.distanceToPoint(cameraPlane, boundingSphere.center) - near) <= boundingSphere.radius
      ) {
        renderSolidInteriorMark(r, variant, -1, passState);
        renderSolidInteriorFill(r, variant, mode, -1, passState);
      }

      const clipInfo = r.getClipInfo();
      if (clipInfo && clipInfo.capIndices.length > 0) {
        const { objects, planes, capIndices } = clipInfo;
        Vec3.scale(solidInteriorEye, cameraPosition, 1 / modelScale);
        for (let j = 0, jl = capIndices.length; j < jl; ++j) {
          const i = capIndices[j];
          if (objects.type[i] === Clip.Type.plane) {
            const d = Plane3D.distanceToPoint(planes[i], solidInteriorEye) * modelScale;
            if (back ? d >= -1e-4 : d <= 1e-4) continue;
          }
          renderSolidInteriorMark(r, variant, i, passState);
          renderSolidInteriorFill(r, variant, mode, i, passState);
        }
      }
      setSolidInteriorPass(0, -1);
    };

    const renderSolidInteriorCaps = (
      renderables: ReadonlyArray<GraphicsRenderable>,
      check: (r: GraphicsRenderable) => boolean,
      variant: GraphicsRenderVariant,
      mode: SolidInteriorMode,
      passState: RenderState,
    ) => {
      if (!solidInteriorCapSupported) return;
      let hasCaps = false;
      for (let i = 0, il = renderables.length; i < il; ++i) {
        const r = renderables[i];
        if (!hasSolidInteriorCap(r) || !check(r)) continue;
        if (!hasCaps) {
          if (!canRenderSolidInteriorCaps()) return;
          hasCaps = true;
        }
        renderSolidInteriorCap(r, variant, mode, passState);
      }
    };

    const update = (camera: ICamera, scene: Scene, frame: Frame) => {
      currentFrame = frame;

      ValueCell.update(globalUniforms.uView, camera.view);
      ValueCell.update(globalUniforms.uInvView, Mat4.invert(invView, camera.view));
      ValueCell.update(globalUniforms.uProjection, camera.projection);
      ValueCell.update(globalUniforms.uInvProjection, Mat4.invert(invProjection, camera.projection));

      ValueCell.updateIfChanged(globalUniforms.uIsOrtho, camera.state.mode === 'orthographic' ? 1 : 0);
      ValueCell.update(
        globalUniforms.uViewOffset,
        camera.viewOffset.enabled
          ? Vec2.set(viewOffset, camera.viewOffset.offsetX * 16, camera.viewOffset.offsetY * 16)
          : Vec2.set(viewOffset, 0, 0),
      );
      ValueCell.updateIfChanged(globalUniforms.uModelScale, camera.scale);

      ValueCell.update(globalUniforms.uCameraPosition, Mat4.getTranslation(cameraPosition, invView));
      Vec3.scale(cameraTarget, camera.state.target, camera.scale);
      Vec3.normalize(cameraDir, Vec3.sub(cameraDir, cameraTarget, cameraPosition));
      ValueCell.update(globalUniforms.uCameraDir, cameraDir);

      ValueCell.updateIfChanged(globalUniforms.uFar, camera.far);
      ValueCell.updateIfChanged(globalUniforms.uNear, camera.near);
      ValueCell.updateIfChanged(globalUniforms.uFog, camera.state.fog > 0);
      ValueCell.updateIfChanged(globalUniforms.uFogFar, camera.fogFar);
      ValueCell.updateIfChanged(globalUniforms.uFogNear, camera.fogNear);
      ValueCell.updateIfChanged(globalUniforms.uTransparentBackground, transparentBackground);

      Frustum3D.fromProjectionMatrix(frustum, camera.projectionView);

      Plane3D.copy(cameraPlane, frustum[Frustum3D.PlaneIndex.Near]);
      cameraPlane.constant -= Plane3D.distanceToPoint(cameraPlane, cameraPosition);
      ValueCell.update(
        globalUniforms.uCameraPlane,
        Plane3D.toArray(cameraPlane, globalUniforms.uCameraPlane.ref.value, 0),
      );

      ValueCell.updateIfChanged(globalUniforms.uMarkerAverage, scene.markerAverage);

      const hasHeadRotation = !Mat4.isZero(camera.headRotation);
      if (hasHeadRotation) {
        ValueCell.updateIfChanged(globalUniforms.uHasHeadRotation, true);
        ValueCell.update(globalUniforms.uInvHeadRotation, Mat4.invert(invHeadRotation, camera.headRotation));
        ValueCell.update(globalUniforms.uLightDirection, getTransformedLightDirection(light, invHeadRotation));
      } else {
        ValueCell.updateIfChanged(globalUniforms.uHasHeadRotation, false);
        ValueCell.updateIfChanged(globalUniforms.uInvHeadRotation, Mat4.id);
        ValueCell.update(globalUniforms.uLightDirection, light.direction);
      }

      ValueCell.update(globalUniforms.uIsAsymmetricProjection, camera.isAsymmetricProjection);
    };

    const updateInternal = (
      group: Scene.Group,
      camera: ICamera,
      depthTexture: Texture | null,
      depthAlphaTexture: Texture | null,
      renderMask: Mask,
      markingDepthTest: boolean,
    ) => {
      arrayMapUpsert(sharedTexturesList, 'tDepth', depthTexture || emptyDepthTexture);
      arrayMapUpsert(sharedTexturesList, 'tDepthAlpha', depthAlphaTexture || emptyDepthAlphaTexture);

      modelScale = camera.scale;

      ValueCell.update(globalUniforms.uModel, Mat4.scaleUniformly(model, group.view, modelScale));
      ValueCell.update(globalUniforms.uModelView, Mat4.mul(modelView, camera.view, model));
      ValueCell.update(globalUniforms.uInvModelView, Mat4.invert(invModelView, modelView));
      ValueCell.update(
        globalUniforms.uModelViewProjection,
        Mat4.mul(modelViewProjection, modelView, camera.projection),
      );
      ValueCell.update(
        globalUniforms.uInvModelViewProjection,
        Mat4.invert(invModelViewProjection, modelViewProjection),
      );

      ValueCell.updateIfChanged(globalUniforms.uRenderMask, renderMask);
      ValueCell.updateIfChanged(globalUniforms.uMarkingDepthTest, markingDepthTest);

      const hasEyeCamera = !Mat4.isZero(camera.viewEye);
      if (hasEyeCamera) {
        ValueCell.updateIfChanged(globalUniforms.uHasEyeCamera, true);
        ValueCell.update(globalUniforms.uModelViewEye, Mat4.mul(modelViewEye, camera.viewEye, model));
        ValueCell.update(globalUniforms.uInvModelViewEye, Mat4.invert(invModelViewEye, modelViewEye));
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

    const checkOpaque = function (r: GraphicsRenderable) {
      // uAlpha is updated in `r.render` so we need to recompute it here
      const alpha = clamp(r.values.alpha.ref.value * r.state.alphaFactor, 0, 1);
      const xrayShaded = r.values.dXrayShaded?.ref.value === 'on' || r.values.dXrayShaded?.ref.value === 'inverted';
      return (
        (alpha === 1 &&
          r.values.transparencyAverage.ref.value !== 1 &&
          r.values.dGeometryType.ref.value !== 'directVolume' &&
          r.values.dPointStyle?.ref.value !== 'fuzzy' &&
          !xrayShaded) ||
        r.values.dTransparentBackfaces?.ref.value === 'opaque'
      );
    };

    const checkTransparent = function (r: GraphicsRenderable) {
      // uAlpha is updated in `r.render` so we need to recompute it here
      const alpha = clamp(r.values.alpha.ref.value * r.state.alphaFactor, 0, 1);
      const xrayShaded = r.values.dXrayShaded?.ref.value === 'on' || r.values.dXrayShaded?.ref.value === 'inverted';
      return (
        alpha !== 0 &&
        (alpha < 1 ||
          r.values.transparencyAverage.ref.value > 0 ||
          r.values.dGeometryType.ref.value === 'directVolume' ||
          r.values.dPointStyle?.ref.value === 'fuzzy' ||
          r.values.dGeometryType.ref.value === 'text' ||
          r.values.dGeometryType.ref.value === 'image' ||
          xrayShaded)
      );
    };

    const checkPickable = function (r: GraphicsRenderable) {
      return !r.state.colorOnly;
    };

    const checkMarkingDepth = function (r: GraphicsRenderable) {
      const alpha = clamp(r.values.alpha.ref.value * r.state.alphaFactor, 0, 1);
      return alpha !== 0 && r.values.transparencyAverage.ref.value !== 1 && r.values.markerAverage.ref.value !== 1;
    };

    const checkMarkingMask = function (r: GraphicsRenderable) {
      return r.values.markerAverage.ref.value > 0;
    };

    const checkEmissive = function (r: GraphicsRenderable) {
      return r.values.emissiveAverage.ref.value + r.values.uEmissive.ref.value > 0;
    };

    /** caps are drawn from back faces, which transparent passes discard when `transparentBackfaces` is 'off' */
    const checkTransparentCap = function (r: GraphicsRenderable) {
      return checkTransparent(r) && r.values.dTransparentBackfaces?.ref.value !== 'off';
    };

    const renderPick = (group: Scene.Group, camera: ICamera, variant: 'pick' | 'depth', pickType: PickType) => {
      updateInternal(group, camera, null, null, Mask.All, false);
      ValueCell.updateIfChanged(globalUniforms.uPickType, pickType);

      const { renderables } = group;
      for (let i = 0, il = renderables.length; i < il; ++i) {
        const r = renderables[i];
        if (checkPickable(r)) {
          renderObject(r, variant, Flag.None, OpaqueState);
        }
      }
      renderSolidInteriorCaps(renderables, checkPickable, variant, 'opaque', OpaqueState);
    };

    const renderDepth = (group: Scene.Group, camera: ICamera) => {
      updateInternal(group, camera, null, null, Mask.All, false);

      const { renderables } = group;
      for (let i = 0, il = renderables.length; i < il; ++i) {
        renderObject(renderables[i], 'depth', Flag.None, OpaqueState);
      }
      renderSolidInteriorCaps(renderables, () => true, 'depth', 'opaque', OpaqueState);
    };

    const renderDepthOpaque = (group: Scene.Group, camera: ICamera) => {
      updateInternal(group, camera, null, null, Mask.Opaque, false);

      const { renderables } = group;
      for (let i = 0, il = renderables.length; i < il; ++i) {
        const r = renderables[i];
        if (checkOpaque(r)) {
          renderObject(r, 'depth', Flag.None, OpaqueState);
        }
      }
      renderSolidInteriorCaps(renderables, checkOpaque, 'depth', 'opaque', OpaqueState);
    };

    const renderDepthOpaqueBack = (group: Scene.Group, camera: ICamera) => {
      updateInternal(group, camera, null, null, Mask.Opaque, false);
      ValueCell.updateIfChanged(globalUniforms.uDepthBack, true);

      const { renderables } = group;
      for (let i = 0, il = renderables.length; i < il; ++i) {
        const r = renderables[i];
        if (checkOpaque(r)) {
          renderObject(r, 'depth', Flag.DepthBack, OpaqueBackState);
        }
      }
      renderSolidInteriorCaps(renderables, checkOpaque, 'depth', 'back', OpaqueBackState);
      ValueCell.updateIfChanged(globalUniforms.uDepthBack, false);
    };

    const renderDepthTransparent = (group: Scene.Group, camera: ICamera, depthTexture: Texture) => {
      updateInternal(group, camera, depthTexture, null, Mask.Transparent, false);

      const { renderables } = group;
      for (let i = 0, il = renderables.length; i < il; ++i) {
        const r = renderables[i];
        if (checkTransparent(r)) {
          renderObject(r, 'depth', Flag.None, OpaqueState);
        }
      }
      renderSolidInteriorCaps(renderables, checkTransparentCap, 'depth', 'opaque', OpaqueState);
    };

    const renderMarkingDepth = (group: Scene.Group, camera: ICamera) => {
      updateInternal(group, camera, null, null, Mask.All, false);
      ValueCell.updateIfChanged(globalUniforms.uMarkingType, MarkingType.Depth);

      const { renderables } = group;
      for (let i = 0, il = renderables.length; i < il; ++i) {
        const r = renderables[i];
        if (checkMarkingDepth(r)) {
          renderObject(r, 'marking', Flag.None, OpaqueState);
        }
      }
      renderSolidInteriorCaps(renderables, checkMarkingDepth, 'marking', 'opaque', OpaqueState);
    };

    const renderMarkingMask = (group: Scene.Group, camera: ICamera, depthTexture: Texture | null) => {
      updateInternal(group, camera, null, depthTexture, Mask.All, !!depthTexture);
      ValueCell.updateIfChanged(globalUniforms.uMarkingType, MarkingType.Mask);

      const { renderables } = group;
      for (let i = 0, il = renderables.length; i < il; ++i) {
        const r = renderables[i];
        if (checkMarkingMask(r)) {
          renderObject(r, 'marking', Flag.None, OpaqueState);
        }
      }
      renderSolidInteriorCaps(renderables, checkMarkingMask, 'marking', 'opaque', OpaqueState);
    };

    const renderEmissiveOpaque = (
      group: Scene.Group,
      camera: ICamera,
      depthTexture: Texture,
      occludeWithOpaqueDepth: boolean,
    ) => {
      // without occluding opaque depth (packed depth): build own occluding depth from all opaque geometry
      const passState = occludeWithOpaqueDepth ? EmissiveOccludedState : OpaqueState;

      // tDepthAlpha = front transparent depth+alpha; emissive shader dims emitters behind it
      updateInternal(group, camera, null, depthTexture, Mask.Opaque, false);

      const { renderables } = group;
      for (let i = 0, il = renderables.length; i < il; ++i) {
        const r = renderables[i];
        if (checkOpaque(r) && (!occludeWithOpaqueDepth || checkEmissive(r))) {
          renderObject(r, 'emissive', Flag.None, passState);
        }
      }
      renderSolidInteriorCaps(
        renderables,
        (r) => checkOpaque(r) && (!occludeWithOpaqueDepth || checkEmissive(r)),
        'emissive',
        occludeWithOpaqueDepth ? 'blended' : 'opaque',
        passState,
      );
    };

    const renderEmissiveTransparent = (group: Scene.Group, camera: ICamera, depthTexture: Texture) => {
      // tDepthAlpha = front transparent depth+alpha; emissive shader dims emitters behind it
      updateInternal(group, camera, null, depthTexture, Mask.Transparent, false);

      const { renderables } = group;
      for (let i = 0, il = renderables.length; i < il; ++i) {
        const r = renderables[i];
        if (checkTransparent(r) && r.values.dGeometryType.ref.value !== 'directVolume' && checkEmissive(r)) {
          renderObject(r, 'emissive', Flag.None, EmissiveTransparentState);
        }
      }
      renderSolidInteriorCaps(
        renderables,
        (r) => checkTransparentCap(r) && checkEmissive(r),
        'emissive',
        'blended',
        EmissiveTransparentState,
      );
    };

    const renderTracing = (group: Scene.Group, camera: ICamera) => {
      updateInternal(group, camera, null, null, Mask.Opaque, false);

      const { renderables } = group;
      for (let i = 0, il = renderables.length; i < il; ++i) {
        const r = renderables[i];
        if (checkOpaque(r)) {
          renderObject(r, 'tracing', Flag.None, OpaqueState);
        }
      }
      renderSolidInteriorCaps(renderables, checkOpaque, 'tracing', 'opaque', OpaqueState);
    };

    const renderBlended = (scene: Scene, camera: ICamera) => {
      if (scene.hasOpaque) {
        renderOpaque(scene, camera);
      }
      if (scene.opacityAverage < 1) {
        renderBlendedTransparent(scene, camera);
      }
    };

    const renderOpaque = (group: Scene.Group, camera: ICamera) => {
      updateInternal(group, camera, null, null, Mask.Opaque, false);

      const { renderables } = group;
      for (let i = 0, il = renderables.length; i < il; ++i) {
        const r = renderables[i];
        if (checkOpaque(r)) {
          renderObject(r, 'color', Flag.None, OpaqueState);
        }
      }
      renderSolidInteriorCaps(renderables, checkOpaque, 'color', 'opaque', OpaqueState);
    };

    const renderBlendedTransparent = (group: Scene.Group, camera: ICamera) => {
      const passState = transparentBackground ? BlendedTransparentBackgroundState : BlendedTransparentState;
      updateInternal(group, camera, null, null, Mask.Transparent, false);

      const renderCaps = canRenderSolidInteriorCaps();
      const { renderables } = group;
      for (let i = 0, il = renderables.length; i < il; ++i) {
        const r = renderables[i];
        if (checkTransparent(r)) {
          if (r.values.uDoubleSided?.ref.value) {
            // render frontfaces and backfaces separately to avoid artefacts
            if (r.values.dTransparentBackfaces?.ref.value !== 'opaque') {
              renderObject(r, 'color', Flag.BlendedBack, passState);
            }
            renderObject(r, 'color', Flag.BlendedFront, passState);
          } else {
            renderObject(r, 'color', Flag.None, passState);
          }
          if (renderCaps && hasSolidInteriorCap(r) && checkTransparentCap(r)) {
            renderSolidInteriorCap(r, 'color', 'blended', passState);
          }
        }
      }
    };

    const renderVolume = (group: Scene.Group, camera: ICamera, depthTexture: Texture) => {
      updateInternal(group, camera, depthTexture, null, Mask.Transparent, false);

      const { renderables } = group;
      for (let i = 0, il = renderables.length; i < il; ++i) {
        const r = renderables[i];
        if (r.values.dGeometryType.ref.value === 'directVolume') {
          renderObject(r, 'color', Flag.None, VolumeState);
        }
      }
    };

    const renderWboitTransparent = (group: Scene.Group, camera: ICamera, depthTexture: Texture | null) => {
      updateInternal(group, camera, depthTexture, null, Mask.Transparent, false);

      const { renderables } = group;
      for (let i = 0, il = renderables.length; i < il; ++i) {
        const r = renderables[i];
        if (checkTransparent(r)) {
          renderObject(r, 'color', Flag.None, WboitState);
        }
      }
      renderSolidInteriorCaps(renderables, checkTransparentCap, 'color', 'oit', WboitState);
    };

    const renderDpoitTransparent = (
      group: Scene.Group,
      camera: ICamera,
      depthTexture: Texture,
      dpoitTextures: { depth: Texture; frontColor: Texture; backColor: Texture },
    ) => {
      arrayMapUpsert(sharedTexturesList, 'tDpoitDepth', dpoitTextures.depth);
      arrayMapUpsert(sharedTexturesList, 'tDpoitFrontColor', dpoitTextures.frontColor);
      arrayMapUpsert(sharedTexturesList, 'tDpoitBackColor', dpoitTextures.backColor);

      updateInternal(group, camera, depthTexture, null, Mask.Transparent, false);

      const { renderables } = group;

      for (let i = 0, il = renderables.length; i < il; ++i) {
        const r = renderables[i];
        if (checkTransparent(r)) {
          renderObject(r, 'color', Flag.None, DpoitState);
        }
      }
      renderSolidInteriorCaps(renderables, checkTransparentCap, 'color', 'oit', DpoitState);
    };

    return {
      clear: (toBackgroundColor: boolean, ignoreTransparentBackground?: boolean, forceToTransparency?: boolean) => {
        ctx.setScissor(viewport.x, viewport.y, viewport.width, viewport.height);
        if (forceToTransparency || (transparentBackground && !ignoreTransparentBackground)) {
          ctx.clear({ color: [0, 0, 0, 0], depth: 1 });
        } else if (toBackgroundColor) {
          ctx.clear({ color: [bgColor[0], bgColor[1], bgColor[2], 1], depth: 1 });
        } else {
          ctx.clear({ color: [1, 1, 1, 1], depth: 1 });
        }
      },
      clearDepth: (packed = false) => {
        ctx.setScissor(viewport.x, viewport.y, viewport.width, viewport.height);
        if (packed) {
          ctx.clear({ color: [1, 1, 1, 1], depth: 1 });
        } else {
          ctx.clear({ depth: 1 });
        }
      },
      update,

      renderPick,
      renderDepth,
      renderDepthOpaque,
      renderDepthOpaqueBack,
      renderDepthTransparent,
      renderMarkingDepth,
      renderMarkingMask,
      renderEmissiveOpaque,
      renderEmissiveTransparent,
      renderTracing,
      renderBlended,
      renderOpaque,
      renderBlendedTransparent,
      renderVolume,
      renderWboitTransparent,
      renderDpoitTransparent,

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
        globals.destroy();
        if (emptyDepthTexture !== emptyDepthAlphaTexture) emptyDepthTexture.destroy();
        emptyDepthAlphaTexture.destroy();
        // TODO
      },
    };
  }
}

export { Renderer };
