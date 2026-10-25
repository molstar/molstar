/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { type PassBuilder } from './pass-builder.js';
import { QuadSchema, QuadValues } from '@molstar/graphics/gl/compute/util';
import { TextureSpec, UniformSpec, type Values } from '@molstar/graphics/gl/renderable/schema';
import type { Texture } from '@molstar/graphics/gl/texture';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import { RenderState } from '@molstar/graphics/gl/render-state';
import { ValueCell } from '@molstar/core/util';
import { Vec2 } from '@molstar/core/math/linear-algebra';
import { ShaderCode } from '@molstar/graphics/gl/shader-code';
import { createComputeRenderable, type ComputeRenderable } from '@molstar/graphics/gl/renderable';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { RenderTarget } from '@molstar/graphics/gl/render-target';
import { Camera } from '../camera.js';
import type { PostprocessingProps } from './postprocessing.js';
import type { DrawPass } from './draw.js';
import type { RenderContext } from '../util.js';
import type { MarkingProps } from './marking.js';
import { isTimingMode } from '@molstar/core/util/debug';
import { clearJitter, getJitterOffsets, getTemporalSamplesPerFrame, setJitter } from './jitter.js';
import { StereoCamera } from '../camera/stereo.js';
import type { Viewport } from '../camera/util.js';

export const ComposeSchema = {
  ...QuadSchema,
  tColor: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
  uTexSize: UniformSpec('v2'),
  uWeight: UniformSpec('f'),
};
export const ComposeShaderCode = ShaderCode('multi-sample-compose');
type ComposeRenderable = ComputeRenderable<Values<typeof ComposeSchema>>;

function getComposeRenderable(ctx: GraphicsContext, colorTexture: Texture): ComposeRenderable {
  const values: Values<typeof ComposeSchema> = {
    ...QuadValues,
    tColor: ValueCell.create(colorTexture),
    uTexSize: ValueCell.create(Vec2.create(colorTexture.getWidth(), colorTexture.getHeight())),
    uWeight: ValueCell.create(1.0),
  };

  const schema = { ...ComposeSchema };
  const renderItem = ctx.createComputeRenderItem('triangles', ComposeShaderCode, schema, values);

  return createComputeRenderable(renderItem, values);
}

export const MultiSampleParams = {
  mode: PD.Select('temporal', [
    ['off', 'Off'],
    ['on', 'On'],
    ['temporal', 'Temporal'],
  ]),
  sampleLevel: PD.Numeric(2, { min: 0, max: 5, step: 1 }, { description: 'Take level^2 samples.' }),
  reduceFlicker: PD.Boolean(true, {
    description: 'Reduce flicker in "temporal" mode when marking is disabled (i.e. drawn as part of the scene).',
  }),
  reuseOcclusion: PD.Boolean(true, { description: 'Reuse occlusion data. It is faster but has some artefacts.' }),
};
export type MultiSampleProps = PD.Values<typeof MultiSampleParams>;

type Props = {
  multiSample: MultiSampleProps;
  postprocessing: PostprocessingProps;
  marking: MarkingProps;
  transparentBackground: boolean;
  dpoitIterations: number;
  volumeResolutionScale: number;
};

/** A stereo view can't be marked from a single camera, there marking stays part of the sample rendering. */
function isSeparateMarking(camera: Camera | StereoCamera): camera is Camera {
  return !StereoCamera.is(camera);
}

export class MultiSamplePass {
  static isEnabled(props: MultiSampleProps) {
    return props.mode !== 'off';
  }

  /**
   * accumulated samples, marking is only blended into it when not rendering to the drawing buffer; rendering to the
   * drawing buffer, the compose target holds them where it is a float target (`floatCompose`)
   */
  readonly colorTarget: RenderTarget;

  private composeTarget: RenderTarget;
  /** the compose target holds negative values unclamped, see `planTemporalMultiSample` */
  private readonly floatCompose: boolean;
  private holdTarget: RenderTarget;
  private compose: ComposeRenderable;
  private readonly renderBuilder: PassBuilder;

  constructor(
    private graphics: GraphicsContext,
    private drawPass: DrawPass,
  ) {
    const { floatBlend, halfFloatRenderTarget } = graphics.capabilities;
    const width = drawPass.colorTarget.getWidth();
    const height = drawPass.colorTarget.getHeight();
    this.colorTarget = drawPass.targets.createRenderTarget(width, height, 'none');
    const type = halfFloatRenderTarget ? 'fp16' : floatBlend ? 'float32' : 'uint8';
    this.composeTarget = drawPass.targets.createRenderTarget(width, height, 'none', type);
    this.floatCompose = type !== 'uint8';
    this.holdTarget = drawPass.targets.createRenderTarget(width, height, 'none');
    this.compose = getComposeRenderable(graphics, drawPass.colorTarget.texture);
    this.renderBuilder = drawPass.targets.createRenderBuilder();
  }

  getByteCount() {
    return this.colorTarget.getByteCount() + this.composeTarget.getByteCount() + this.holdTarget.getByteCount();
  }

  syncSize() {
    const width = this.drawPass.colorTarget.getWidth();
    const height = this.drawPass.colorTarget.getHeight();

    const [w, h] = this.compose.values.uTexSize.ref.value;
    if (width !== w || height !== h) {
      this.colorTarget.setSize(width, height);
      this.composeTarget.setSize(width, height);
      this.holdTarget.setSize(width, height);
      ValueCell.update(this.compose.values.uTexSize, Vec2.set(this.compose.values.uTexSize.ref.value, width, height));
    }
  }

  render(sampleIndex: number, ctx: RenderContext, props: Props, toDrawingBuffer: boolean, forceOn: boolean) {
    return this.plan(this.renderBuilder, sampleIndex, ctx, props, toDrawingBuffer, forceOn);
  }

  /**
   * What `render` draws, see `PassBuilder`; preparing covers all samples (and both the first and the following frames
   * of temporal multi-sampling) and returns `sampleIndex` unchanged
   */
  plan(
    b: PassBuilder,
    sampleIndex: number,
    ctx: RenderContext,
    props: Props,
    toDrawingBuffer: boolean,
    forceOn: boolean,
  ) {
    if (props.multiSample.mode === 'temporal' && !forceOn) {
      if (!b.drawing) {
        this.planTemporalMultiSample(b, -1, ctx, props, toDrawingBuffer);
        this.planTemporalMultiSample(b, 0, ctx, props, toDrawingBuffer);
        return sampleIndex;
      }
      return this.planTemporalMultiSample(b, sampleIndex, ctx, props, toDrawingBuffer);
    } else {
      this.planMultiSample(b, ctx, props, toDrawingBuffer);
      return -2;
    }
  }

  private planMultiSample(b: PassBuilder, ctx: RenderContext, props: Props, toDrawingBuffer: boolean) {
    const { camera } = ctx;
    const { compose, composeTarget, drawPass, graphics } = this;
    const timing = isTimingMode && b.drawing;
    if (timing) graphics.timer.group('MultiSamplePass.renderMultiSample');

    // based on the Multisample Anti-Aliasing Render Pass
    // contributed to three.js by bhouston / http://clara.io/
    //
    // This manual approach to MSAA re-renders the scene once for
    // each sample with camera jitter and accumulates the results.
    const offsetList = getJitterOffsets(props.multiSample.sampleLevel);

    const { x, y, width, height } = camera.viewport;
    const baseSampleWeight = 1.0 / offsetList.length;
    const roundingRange = 1 / 32;

    ValueCell.update(compose.values.tColor, drawPass.getColorTarget(ctx.scene, props.postprocessing).texture);
    compose.update();

    // render the scene multiple times, each slightly jitter offset
    // from the last and accumulate the results (preparing: the first sample draws as all of them)
    const sampleCount = b.drawing ? offsetList.length : 1;
    for (let i = 0; i < sampleCount; ++i) {
      const offset = offsetList[i];
      if (b.drawing) setJitter(camera, offset);

      // the theory is that equal weights for each sample lead to an accumulation of rounding
      // errors. The following equation varies the sampleWeight per sample so that it is uniformly
      // distributed across a range of values whose rounding errors cancel each other out.
      const uniformCenteredDistribution = -0.5 + (i + 0.5) / offsetList.length;
      const sampleWeight = baseSampleWeight + roundingRange * uniformCenteredDistribution;
      ValueCell.update(compose.values.uWeight, sampleWeight);

      // render scene
      if (i === 0 || !props.multiSample.reuseOcclusion) {
        drawPass.postprocessing.setOcclusionOffset(0, 0);
      } else {
        drawPass.postprocessing.setOcclusionOffset(offset[0] / width, -offset[1] / height);
      }
      drawPass.plan(b, ctx, props, false, isSeparateMarking(camera));

      // compose rendered scene with compose target
      b.begin(composeTarget, { label: 'multi-sample compose' });
      if (b.drawing) {
        graphics.setRenderState(RenderState.Additive);
        graphics.setViewport(x, y, width, height);
        graphics.setScissor(x, y, width, height);
        if (i === 0) graphics.clear({ color: [0, 0, 0, 0] });
      }
      b.quad(compose, RenderState.Additive);
    }

    drawPass.postprocessing.setOcclusionOffset(0, 0);

    // the drawing buffer gets a copy of the compose target (`MarkingPass.present`); in place, marking needs a copy of
    // its own
    if (!toDrawingBuffer) this.copyToColorTarget(b, camera.viewport);

    if (b.drawing) clearJitter(camera);

    this.drawPass.marking.present(b, ctx, props, {
      base: toDrawingBuffer ? composeTarget : this.colorTarget,
      toDrawingBuffer,
      offsets: offsetList,
      restart: true,
      samples: offsetList.length,
      shading: drawPass.getMarkingShading(props.postprocessing),
    });
    if (timing) graphics.timer.groupEnd('MultiSamplePass.renderMultiSample');
  }

  private planTemporalMultiSample(
    b: PassBuilder,
    sampleIndex: number,
    ctx: RenderContext,
    props: Props,
    toDrawingBuffer: boolean,
  ) {
    const { camera } = ctx;
    const { compose, composeTarget, holdTarget, drawPass, graphics, floatCompose } = this;
    const timing = isTimingMode && b.drawing;
    if (timing) graphics.timer.group('MultiSamplePass.renderTemporalMultiSample');

    // based on the Multisample Anti-Aliasing Render Pass
    // contributed to three.js by bhouston / http://clara.io/
    //
    // This manual approach to MSAA re-renders the scene once for
    // each sample with camera jitter and accumulates the results.
    const offsetList = getJitterOffsets(props.multiSample.sampleLevel);

    if (sampleIndex === -2 || sampleIndex >= offsetList.length) return -2;

    const { x, y, width, height } = camera.viewport;
    const sampleWeight = 1.0 / offsetList.length;
    // preparing: one sample draws as all of them
    const samplesPerFrame = b.drawing ? getTemporalSamplesPerFrame(props.multiSample.sampleLevel) : 1;
    // the unjittered first frame is also the one rendered while the camera moves
    const firstFrame = sampleIndex === -1;

    if (firstFrame) {
      drawPass.plan(b, ctx, props, false, isSeparateMarking(camera));
      ValueCell.update(compose.values.uWeight, 1.0);
      ValueCell.update(compose.values.tColor, drawPass.getColorTarget(ctx.scene, props.postprocessing).texture);
      compose.update();

      b.begin(holdTarget, { label: 'multi-sample hold' });
      if (b.drawing) {
        graphics.setRenderState(RenderState.Default);
        graphics.setViewport(x, y, width, height);
        graphics.setScissor(x, y, width, height);
      }
      b.quad(compose, RenderState.Default);
      if (floatCompose) {
        b.begin(composeTarget, { label: 'multi-sample compose' });
        if (b.drawing) {
          graphics.setViewport(x, y, width, height);
          graphics.setScissor(x, y, width, height);
        }
        b.quad(compose, RenderState.Default);
      }
      sampleIndex += 1;
    } else {
      // render the scene multiple times, each slightly jitter offset
      // from the last and accumulate the results.
      for (let i = 0; i < samplesPerFrame; ++i) {
        const offset = offsetList[sampleIndex];
        if (b.drawing) setJitter(camera, offset);

        // render scene
        if (sampleIndex === 0 || !props.multiSample.reuseOcclusion) {
          drawPass.postprocessing.setOcclusionOffset(0, 0);
        } else {
          drawPass.postprocessing.setOcclusionOffset(offset[0] / width, -offset[1] / height);
        }
        drawPass.plan(b, ctx, props, false, isSeparateMarking(camera));

        // compose rendered scene with compose target
        b.begin(composeTarget, { label: 'multi-sample compose' });
        if (b.drawing) {
          graphics.setRenderState(RenderState.Additive);
          graphics.setViewport(x, y, width, height);
          graphics.setScissor(x, y, width, height);
          if (sampleIndex === 0 && !floatCompose) graphics.clear({ color: [0, 0, 0, 0] });
        }
        ValueCell.update(compose.values.tColor, drawPass.getColorTarget(ctx.scene, props.postprocessing).texture);
        ValueCell.update(compose.values.uWeight, sampleWeight);
        compose.update();
        b.quad(compose, RenderState.Additive);
        if (floatCompose) {
          // the first sample gives this one its weight (unclamped)
          ValueCell.update(compose.values.tColor, holdTarget.texture);
          ValueCell.update(compose.values.uWeight, -sampleWeight);
          compose.update();
          b.quad(compose, RenderState.Additive);
        }

        sampleIndex += 1;
        if (sampleIndex >= offsetList.length) break;
      }
    }

    drawPass.postprocessing.setOcclusionOffset(0, 0);

    // the compose target holds `w * (s_1 + ... + s_k) + (1 - k * w) * s_0` after sample `k` (samples `s`, weight `w`):
    // the first sample, the hold, until every following one has its weight; the drawing buffer gets a copy of it
    // (`MarkingPass.present`), in place, marking needs a copy of its own. An 8 bit compose target clamps: it holds the
    // sum only, the output adds the first sample.
    if (!floatCompose) this.planTemporalOutput(b, camera.viewport, sampleIndex * sampleWeight);
    else if (!toDrawingBuffer) this.copyToColorTarget(b, camera.viewport);

    if (b.drawing) clearJitter(camera);

    this.drawPass.marking.present(b, ctx, props, {
      base: floatCompose && toDrawingBuffer ? composeTarget : this.colorTarget,
      toDrawingBuffer,
      offsets: offsetList,
      restart: firstFrame,
      samples: firstFrame ? 1 : samplesPerFrame,
      shading: drawPass.getMarkingShading(props.postprocessing),
    });
    if (timing) graphics.timer.groupEnd('MultiSamplePass.renderTemporalMultiSample');

    return sampleIndex >= offsetList.length ? -2 : sampleIndex;
  }

  /** The compose target copied into the color target */
  private copyToColorTarget(b: PassBuilder, viewport: Viewport) {
    const { compose, composeTarget, graphics } = this;
    b.begin(this.colorTarget, { label: 'multi-sample output' });
    if (b.drawing) {
      graphics.setViewport(viewport.x, viewport.y, viewport.width, viewport.height);
      graphics.setScissor(viewport.x, viewport.y, viewport.width, viewport.height);
    }
    ValueCell.update(compose.values.uWeight, 1.0);
    ValueCell.update(compose.values.tColor, composeTarget.texture);
    compose.update();
    b.quad(compose, RenderState.Default);
  }

  /** The color target as the sum of the samples so far plus the first sample with the weight of the missing ones */
  private planTemporalOutput(b: PassBuilder, viewport: Viewport, accumulationWeight: number) {
    const { compose, composeTarget, holdTarget, graphics } = this;
    b.begin(this.colorTarget, { label: 'multi-sample output' });
    if (b.drawing) {
      graphics.setViewport(viewport.x, viewport.y, viewport.width, viewport.height);
      graphics.setScissor(viewport.x, viewport.y, viewport.width, viewport.height);
    }
    if (accumulationWeight > 0) {
      ValueCell.update(compose.values.uWeight, 1.0);
      ValueCell.update(compose.values.tColor, composeTarget.texture);
      compose.update();
      b.quad(compose, RenderState.Default);
    }
    if (accumulationWeight < 1.0) {
      ValueCell.update(compose.values.uWeight, 1.0 - accumulationWeight);
      ValueCell.update(compose.values.tColor, holdTarget.texture);
      compose.update();
      b.quad(compose, accumulationWeight === 0 ? RenderState.Default : RenderState.Additive);
    }
  }
}

export class MultiSampleHelper {
  private sampleIndex = -2;

  update(changed: boolean, props: MultiSampleProps) {
    if (changed) this.sampleIndex = -1;
    return props.mode === 'temporal' ? this.sampleIndex !== -2 : false;
  }

  /** Return `true` while more samples are needed */
  render(ctx: RenderContext, props: Props, toDrawingBuffer: boolean, forceOn?: boolean) {
    this.sampleIndex = this.multiSamplePass.render(this.sampleIndex, ctx, props, toDrawingBuffer, !!forceOn);
    return this.sampleIndex < 0;
  }

  /** What `render` draws, see `PassBuilder` */
  plan(b: PassBuilder, ctx: RenderContext, props: Props, toDrawingBuffer: boolean, forceOn?: boolean) {
    this.multiSamplePass.plan(b, this.sampleIndex, ctx, props, toDrawingBuffer, !!forceOn);
  }

  constructor(private multiSamplePass: MultiSamplePass) {}
}
