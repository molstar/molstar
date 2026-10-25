/**
 * Copyright (c) 2024-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import type { PassBuilder } from './pass-builder.js';
import { QuadSchema, QuadValues } from '@molstar/graphics/gl/compute/util';
import {
  DefineSpec,
  MaxLightCount,
  TextureSpec,
  UniformSpec,
  type Values,
} from '@molstar/graphics/gl/renderable/schema';
import type { Texture } from '@molstar/graphics/gl/texture';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import { BlendState, RenderState } from '@molstar/graphics/gl/render-state';
import { ValueCell } from '@molstar/core/util';
import { isTimingMode } from '@molstar/core/util/debug';
import type { Renderer } from '@molstar/graphics/gl/renderer';
import type { Camera, ICamera } from '../camera.js';
import type { Scene } from '@molstar/graphics/gl/scene';
import type { RenderTarget } from '@molstar/graphics/gl/render-target';
import { ShaderCode } from '@molstar/graphics/gl/shader-code';
import { type ComputeRenderable, createComputeRenderable } from '@molstar/graphics/gl/renderable';
import { Vec2 } from '@molstar/core/math/linear-algebra/3d/vec2';
import { Mat4 } from '@molstar/core/math/linear-algebra/3d/mat4';
import { Vec4 } from '@molstar/core/math/linear-algebra/3d/vec4';
import { Vec3 } from '@molstar/core/math/linear-algebra/3d/vec3';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Color } from '@molstar/core/util/color/color';
import type { Framebuffer } from '@molstar/graphics/gl/framebuffer';
import { now } from '@molstar/core/util/now';
import { clamp } from '@molstar/core/math/interpolate';
import type { DrawPass } from './draw.js';
import type { RenderContext } from '../util.js';

export const TracingParams = {
  rendersPerFrame: PD.Interval(
    [1, 16],
    { min: 1, max: 64, step: 1 },
    {
      description:
        'Number of rays per pixel each frame. May be adjusted to reach targetFps but will stay within given interval.',
    },
  ),
  targetFps: PD.Numeric(
    30,
    { min: 0, max: 120, step: 0.1 },
    { description: 'Target FPS per frame. If observed FPS is lower or higher, some parameters may get adjusted.' },
  ),
  steps: PD.Numeric(32, { min: 1, max: 1024, step: 1 }),
  refineSteps: PD.Numeric(
    4,
    { min: 0, max: 8, step: 1 },
    { description: 'Number of refine steps per ray hit. May be lower to reach targetFps.' },
  ),
  rayDistance: PD.Numeric(
    256,
    { min: 1, max: 8192, step: 1 },
    { description: 'Maximum distance a ray can travel (in world units).' },
  ),
  thicknessMode: PD.Select('auto', PD.arrayToOptions(['auto', 'fixed'] as const)),
  minThickness: PD.Numeric(0.5, { min: 0.1, max: 16, step: 0.1 }, { hideIf: (p) => p.thicknessMode === 'fixed' }),
  thicknessFactor: PD.Numeric(1, { min: 0.1, max: 2, step: 0.05 }, { hideIf: (p) => p.thicknessMode === 'fixed' }),
  thickness: PD.Numeric(4, { min: 0.1, max: 512, step: 0.1 }, { hideIf: (p) => p.thicknessMode === 'auto' }),
  bounces: PD.Numeric(4, { min: 1, max: 32, step: 1 }, { description: 'Number of bounces for each ray.' }),
  shadowEnable: PD.Boolean(false),
  shadowSoftness: PD.Numeric(0.1, { min: 0.01, max: 1.0, step: 0.01 }),
  shadowThickness: PD.Numeric(
    0.5,
    { min: 0.0, max: 32, step: 0.1 },
    { description: 'Thickness of the shadow casting geometry. Set to 0.0 for automatic estimation.' },
  ),
};
export type TracingProps = PD.Values<typeof TracingParams>;

export class TracingPass {
  private readonly framebuffer: Framebuffer;

  readonly colorTextureOpaque: Texture;
  readonly normalTextureOpaque: Texture;
  readonly shadedTextureOpaque: Texture;

  private readonly thicknessTarget: RenderTarget;
  readonly accumulateTarget: RenderTarget;
  readonly composeTarget: RenderTarget;

  private readonly traceRenderable: TraceRenderable;

  constructor(
    private readonly graphics: GraphicsContext,
    private readonly drawPass: DrawPass,
  ) {
    const {
      capabilities: { mixedAttachmentFormats },
    } = graphics;
    const { targets } = drawPass;

    const { depthTextureOpaque } = drawPass;
    const width = depthTextureOpaque.getWidth();
    const height = depthTextureOpaque.getHeight();

    if (mixedAttachmentFormats) {
      this.shadedTextureOpaque = targets.texture('image-uint8', 'rgba', 'ubyte', 'nearest');
      this.shadedTextureOpaque.define(width, height);

      // octahedral-encoded normals with the emissive strength, `packNormalToRGBA`
      this.normalTextureOpaque = targets.texture('image-uint8', 'rgba', 'ubyte', 'nearest');
      this.normalTextureOpaque.define(width, height);

      this.colorTextureOpaque = targets.texture('image-uint8', 'rgba', 'ubyte', 'nearest');
      this.colorTextureOpaque.define(width, height);
    } else {
      // webgl1 requires consistent bit plane counts

      this.shadedTextureOpaque = targets.texture('image-float32', 'rgba', 'float', 'nearest');
      this.shadedTextureOpaque.define(width, height);

      this.normalTextureOpaque = targets.texture('image-float32', 'rgba', 'float', 'nearest');
      this.normalTextureOpaque.define(width, height);

      this.colorTextureOpaque = targets.texture('image-float32', 'rgba', 'float', 'nearest');
      this.colorTextureOpaque.define(width, height);
    }

    this.framebuffer = graphics.resources.framebuffer();
    this.framebuffer.setDrawBuffers(3);

    this.shadedTextureOpaque.attachFramebuffer(this.framebuffer, 'color0');
    this.normalTextureOpaque.attachFramebuffer(this.framebuffer, 'color1');
    this.colorTextureOpaque.attachFramebuffer(this.framebuffer, 'color2');

    // the farthest back faces, drawn with the transparency passes' depth-stencil buffer (depth cleared by the pass,
    // stencil per solid interior cap, not read after): the trace reads the packed depth only
    this.thicknessTarget = targets.createRenderTarget(width, height, 'none', 'uint8', 'nearest');
    drawPass.depthTargetTransparent.depthRenderbuffer!.attachFramebuffer(this.thicknessTarget.framebuffer);
    this.accumulateTarget = targets.createRenderTarget(width, height, 'none', 'float32');
    this.composeTarget = targets.createRenderTarget(width, height, 'none', 'uint8', 'linear');

    this.traceRenderable = getTraceRenderable(
      graphics,
      this.colorTextureOpaque,
      this.normalTextureOpaque,
      this.shadedTextureOpaque,
      this.thicknessTarget.texture,
      this.drawPass.depthTextureOpaque,
    );
  }

  getByteCount() {
    return (
      this.thicknessTarget.getByteCount() +
      this.accumulateTarget.getByteCount() +
      this.composeTarget.getByteCount() +
      this.colorTextureOpaque.getByteCount() +
      this.normalTextureOpaque.getByteCount() +
      this.shadedTextureOpaque.getByteCount()
    );
  }

  private renderInput(b: PassBuilder, renderer: Renderer, camera: ICamera, scene: Scene, props: TracingProps) {
    const timing = isTimingMode && b.drawing;
    if (timing) this.graphics.timer.group('TracePass.renderInput');

    this.drawPass.depthTextureOpaque.attachFramebuffer(this.framebuffer, 'depth');
    b.begin(this.framebuffer, { label: 'tracing input' });
    if (b.drawing) renderer.clear(true);
    b.renderer({ name: 'tracing' }, scene.primitives, () => renderer.renderTracing(scene.primitives, camera));

    //

    if (props.thicknessMode === 'auto') {
      b.begin(this.thicknessTarget, {
        label: 'thickness',
        // the farthest back faces (`renderDepthOpaqueBack`), depth reversed: cleared to the near plane (packed as ones)
        clear: { color: [1, 1, 1, 1], depth: 1 },
      });
      b.renderer({ name: 'depthOpaqueBack' }, scene.primitives, () =>
        renderer.renderDepthOpaqueBack(scene.primitives, camera),
      );
    }
    if (timing) this.graphics.timer.groupEnd('TracePass.renderInput');
  }

  setSize(width: number, height: number) {
    const w = this.composeTarget.getWidth();
    const h = this.composeTarget.getHeight();

    if (width !== w || height !== h) {
      this.thicknessTarget.setSize(width, height);
      this.accumulateTarget.setSize(width, height);
      this.composeTarget.setSize(width, height);

      this.colorTextureOpaque.define(width, height);
      this.normalTextureOpaque.define(width, height);
      this.shadedTextureOpaque.define(width, height);

      ValueCell.update(
        this.traceRenderable.values.uTexSize,
        Vec2.set(this.traceRenderable.values.uTexSize.ref.value, width, height),
      );
    }
  }

  private clearAdjustedProps = true;

  restart(clearAdjustedProps = false) {
    this.graphics.beginRenderPass(this.accumulateTarget, { label: 'trace clear', clear: { color: [0, 0, 0, 0] } });
    this.graphics.beginRenderPass(this.composeTarget, { label: 'trace clear', clear: { color: [0, 0, 0, 0] } });

    if (clearAdjustedProps) {
      this.prevTime = 0;
      this.currTime = 0;
      this.clearAdjustedProps = true;
    }
  }

  private prevTime = 0;
  private currTime = 0;
  private rendersPerFrame = 1;
  private refineSteps = 1;
  private steps = 16;

  private increaseAdjustedProps(props: TracingProps) {
    this.steps += 1;
    if (this.steps > props.steps) {
      this.refineSteps += 1;
    }
    if (this.refineSteps > props.refineSteps) {
      this.rendersPerFrame += 1;
    }
  }

  private decreaseAdjustedProps(props: TracingProps) {
    const minRefineSteps = Math.min(1, props.refineSteps);
    this.rendersPerFrame -= 1;
    if (this.rendersPerFrame < 1) {
      this.refineSteps -= 1;
    }
    if (this.refineSteps < minRefineSteps) {
      this.steps -= 1;
    }
  }

  private getAdjustedProps(props: TracingProps, iteration: number) {
    this.currTime = now();
    const minRefineSteps = Math.min(1, props.refineSteps);
    const minSteps = Math.round(props.steps / 2);

    if (this.clearAdjustedProps) {
      this.rendersPerFrame = props.rendersPerFrame[0];
      this.refineSteps = minRefineSteps;
      this.steps = minSteps;
      this.clearAdjustedProps = false;
    }

    if (iteration > 0) {
      const targetTimeMs = 1000 / props.targetFps;
      const deltaTime = this.currTime - this.prevTime;
      let f = Math.round(deltaTime / targetTimeMs);
      if (f >= 2) {
        while (f > 0) {
          this.decreaseAdjustedProps(props);
          f -= 1;
        }
      } else if (deltaTime < targetTimeMs) {
        this.increaseAdjustedProps(props);
      } else if (deltaTime > targetTimeMs + 0.5) {
        this.decreaseAdjustedProps(props);
      }
    }

    this.prevTime = this.currTime;
    this.rendersPerFrame = clamp(this.rendersPerFrame, props.rendersPerFrame[0], props.rendersPerFrame[1]);
    this.refineSteps = clamp(this.refineSteps, minRefineSteps, props.refineSteps);
    this.steps = clamp(this.steps, minSteps, props.steps);

    return {
      rendersPerFrame: iteration === 0 ? Math.ceil(this.rendersPerFrame / 2) : this.rendersPerFrame,
      refineSteps: iteration === 0 ? minRefineSteps : this.refineSteps,
      steps: iteration === 0 ? minSteps : this.steps,
    };
  }

  /**
   * Traces a frame (`b`, see `PassBuilder`); preparing keeps the quality adapted to the frame times so far, which only
   * changes uniforms (defines on WebGL1)
   */
  render(
    b: PassBuilder,
    ctx: RenderContext<Camera>,
    transparentBackground: boolean,
    props: TracingProps,
    iteration: number,
    forceRenderInput: boolean,
  ) {
    const { values } = this.traceRenderable;
    const { rendersPerFrame, refineSteps, steps } = b.drawing
      ? this.getAdjustedProps(props, iteration)
      : {
          rendersPerFrame: values.uRendersPerFrame.ref.value,
          refineSteps: values.uRefineSteps.ref.value,
          steps: values.uSteps.ref.value,
        };

    if (isTimingMode && b.drawing) {
      this.graphics.timer.group('TracePass.render', {
        note: `${rendersPerFrame} rendersPerFrame, ${refineSteps} refineSteps, ${steps} steps`,
      });
    }

    const { renderer, camera, scene, frame } = ctx;
    const { x, y, width, height } = camera.viewport;

    if (iteration === 0 || forceRenderInput) {
      // render color & depth
      renderer.setTransparentBackground(transparentBackground);
      if (b.drawing) {
        renderer.setDrawingBufferSize(this.composeTarget.getWidth(), this.composeTarget.getHeight());
        renderer.setPixelRatio(this.graphics.pixelRatio);
        renderer.setViewport(x, y, width, height);
        renderer.update(camera, scene, frame);
      }
      this.renderInput(b, renderer, camera, scene, props);
    }

    if (b.drawing) {
      this.graphics.setRenderState(RenderState.Default);
      this.graphics.setViewport(x, y, width, height);
      this.graphics.setScissor(x, y, width, height);
    }

    const invProjection = Mat4.identity();
    Mat4.invert(invProjection, camera.projection);
    const orthographic = camera.state.mode === 'orthographic' ? 1 : 0;
    const [w, h] = this.traceRenderable.values.uTexSize.ref.value;
    const v = camera.viewport;

    const ambientColor = Vec3();
    Vec3.scale(
      ambientColor,
      Color.toArrayNormalized(renderer.props.ambientColor, ambientColor, 0),
      renderer.props.ambientIntensity,
    );

    // trace, blended into the running average: the new frame weighs 1 / (iteration + 1), so the first replaces
    b.begin(this.accumulateTarget, { label: 'trace' });
    let needsUpdateTrace = false;
    ValueCell.update(this.traceRenderable.values.uFrameNo, iteration);
    // the adapted quality (see `getAdjustedProps`) is loop bound uniforms, so adapting it compiles nothing (but on
    // WebGL1, where they are defines: `update` takes them, see `needsUpdateTrace`)
    if (this.traceRenderable.values.uRendersPerFrame.ref.value !== rendersPerFrame) {
      ValueCell.update(this.traceRenderable.values.uRendersPerFrame, rendersPerFrame);
      needsUpdateTrace = true;
    }
    ValueCell.update(this.traceRenderable.values.uProjection, camera.projection);
    ValueCell.update(this.traceRenderable.values.uInvProjection, invProjection);
    Vec4.set(
      this.traceRenderable.values.uBounds.ref.value,
      v.x / w,
      v.y / h,
      (v.x + v.width) / w,
      (v.y + v.height) / h,
    );
    ValueCell.update(this.traceRenderable.values.uBounds, this.traceRenderable.values.uBounds.ref.value);
    ValueCell.updateIfChanged(this.traceRenderable.values.uNear, camera.near);
    ValueCell.updateIfChanged(this.traceRenderable.values.uFar, camera.far);
    ValueCell.updateIfChanged(this.traceRenderable.values.uFogFar, camera.fogFar);
    ValueCell.updateIfChanged(this.traceRenderable.values.uFogNear, camera.fogNear);
    ValueCell.update(
      this.traceRenderable.values.uFogColor,
      Color.toVec3Normalized(this.traceRenderable.values.uFogColor.ref.value, renderer.props.backgroundColor),
    );
    if (this.traceRenderable.values.dOrthographic.ref.value !== orthographic) {
      ValueCell.update(this.traceRenderable.values.dOrthographic, orthographic);
      needsUpdateTrace = true;
    }
    ValueCell.update(this.traceRenderable.values.uLightDirection, renderer.light.direction);
    ValueCell.update(this.traceRenderable.values.uLightColor, renderer.light.color);
    if (this.traceRenderable.values.dLightCount.ref.value !== renderer.light.count) {
      ValueCell.update(this.traceRenderable.values.dLightCount, renderer.light.count);
      needsUpdateTrace = true;
    }
    ValueCell.update(this.traceRenderable.values.uAmbientColor, ambientColor);
    ValueCell.updateIfChanged(this.traceRenderable.values.uExposure, renderer.props.exposure);
    if (this.traceRenderable.values.dBounces.ref.value !== props.bounces) {
      ValueCell.update(this.traceRenderable.values.dBounces, props.bounces);
      needsUpdateTrace = true;
    }
    if (this.traceRenderable.values.uSteps.ref.value !== steps) {
      ValueCell.update(this.traceRenderable.values.uSteps, steps);
      needsUpdateTrace = true;
    }
    if (this.traceRenderable.values.uRefineSteps.ref.value !== refineSteps) {
      ValueCell.update(this.traceRenderable.values.uRefineSteps, refineSteps);
      needsUpdateTrace = true;
    }
    ValueCell.updateIfChanged(this.traceRenderable.values.uRayDistance, props.rayDistance);
    if (this.traceRenderable.values.dThicknessMode.ref.value !== props.thicknessMode) {
      ValueCell.update(this.traceRenderable.values.dThicknessMode, props.thicknessMode);
      needsUpdateTrace = true;
    }
    ValueCell.updateIfChanged(this.traceRenderable.values.uMinThickness, props.minThickness);
    ValueCell.updateIfChanged(this.traceRenderable.values.uThicknessFactor, props.thicknessFactor);
    ValueCell.updateIfChanged(this.traceRenderable.values.uThickness, props.thickness);
    if (this.traceRenderable.values.dShadowEnable.ref.value !== props.shadowEnable) {
      ValueCell.update(this.traceRenderable.values.dShadowEnable, props.shadowEnable);
      needsUpdateTrace = true;
    }
    ValueCell.updateIfChanged(this.traceRenderable.values.uShadowSoftness, props.shadowSoftness);
    ValueCell.updateIfChanged(this.traceRenderable.values.uShadowThickness, props.shadowThickness);
    if (needsUpdateTrace) this.traceRenderable.update();
    if (b.drawing) {
      const weight = 1 / (iteration + 1);
      this.graphics.setBlendConstant(weight, weight, weight, weight);
    }
    b.quad(this.traceRenderable, TraceRenderState);
    if (isTimingMode && b.drawing) this.graphics.timer.groupEnd('TracePass.render');
  }
}

/**
 * The color blended with the weight of the frame (the blend constant), the alpha replaced: the shader writes the weight,
 * which the denoiser reads
 */
const TraceRenderState: RenderState = {
  blend: BlendState(
    BlendState.Component('add', 'constant', 'one-minus-constant'),
    BlendState.Component('add', 'one', 'zero'),
  ),
};

//

export const TraceSchema = {
  ...QuadSchema,
  tColor: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
  tNormal: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
  tShaded: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
  tThickness: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
  tDepth: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
  uTexSize: UniformSpec('v2'),

  dOrthographic: DefineSpec('number'),
  uNear: UniformSpec('f'),
  uFar: UniformSpec('f'),
  uFogNear: UniformSpec('f'),
  uFogFar: UniformSpec('f'),
  uFogColor: UniformSpec('v3'),

  uProjection: { ...UniformSpec('m4'), depthRange: 'projection' as const },
  uInvProjection: { ...UniformSpec('m4'), depthRange: 'inverse' as const },
  uBounds: { ...UniformSpec('v4'), windowRect: 'bounds' as const },

  uLightDirection: UniformSpec('v3[]', undefined, undefined, MaxLightCount, 'dLightCount'),
  uLightColor: UniformSpec('v3[]', undefined, undefined, MaxLightCount, 'dLightCount'),
  dLightCount: DefineSpec('number'),
  uAmbientColor: UniformSpec('v3'),
  uExposure: UniformSpec('f'),

  uFrameNo: UniformSpec('i'),
  /** loop bounds: uniforms where loops take them, defines on WebGL1 */
  uRendersPerFrame: { ...UniformSpec('i'), loopBound: true },

  dBounces: DefineSpec('number'),
  uSteps: { ...UniformSpec('i'), loopBound: true },
  uRefineSteps: { ...UniformSpec('i'), loopBound: true },
  uRayDistance: UniformSpec('f'),

  dThicknessMode: DefineSpec('string', ['auto', 'fixed']),
  uMinThickness: UniformSpec('f'),
  uThicknessFactor: UniformSpec('f'),
  uThickness: UniformSpec('f'),

  dShadowEnable: DefineSpec('boolean'),
  uShadowSoftness: UniformSpec('f'),
  uShadowThickness: UniformSpec('f'),
};
export const TraceShaderCode = ShaderCode('trace');
type TraceRenderable = ComputeRenderable<Values<typeof TraceSchema>>;

function getTraceRenderable(
  ctx: GraphicsContext,
  colorTexture: Texture,
  normalTexture: Texture,
  shadedTexture: Texture,
  thicknessTexture: Texture,
  depthTexture: Texture,
): TraceRenderable {
  const values: Values<typeof TraceSchema> = {
    ...QuadValues,
    tColor: ValueCell.create(colorTexture),
    tNormal: ValueCell.create(normalTexture),
    tShaded: ValueCell.create(shadedTexture),
    tThickness: ValueCell.create(thicknessTexture),
    tDepth: ValueCell.create(depthTexture),
    uTexSize: ValueCell.create(Vec2.create(colorTexture.getWidth(), colorTexture.getHeight())),

    dOrthographic: ValueCell.create(0),
    uNear: ValueCell.create(1),
    uFar: ValueCell.create(10000),
    uFogNear: ValueCell.create(10000),
    uFogFar: ValueCell.create(10000),
    uFogColor: ValueCell.create(Vec3.create(1, 1, 1)),

    uProjection: ValueCell.create(Mat4.identity()),
    uInvProjection: ValueCell.create(Mat4.identity()),
    uBounds: ValueCell.create(Vec4()),

    uLightDirection: ValueCell.create([]),
    uLightColor: ValueCell.create([]),
    dLightCount: ValueCell.create(0),
    uAmbientColor: ValueCell.create(Vec3()),
    uExposure: ValueCell.create(1),

    uFrameNo: ValueCell.create(0),
    uRendersPerFrame: ValueCell.create(1),

    dBounces: ValueCell.create(4),
    uSteps: ValueCell.create(32),
    uRefineSteps: ValueCell.create(4),
    uRayDistance: ValueCell.create(256),

    dThicknessMode: ValueCell.create('auto'),
    uMinThickness: ValueCell.create(0.5),
    uThicknessFactor: ValueCell.create(1),
    uThickness: ValueCell.create(4),

    dShadowEnable: ValueCell.create(false),
    uShadowSoftness: ValueCell.create(0.1),
    uShadowThickness: ValueCell.create(0.1),
  };

  const schema = { ...TraceSchema };
  const renderItem = ctx.createComputeRenderItem('triangles', TraceShaderCode, schema, values);

  return createComputeRenderable(renderItem, values);
}
