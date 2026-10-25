/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { GraphicsContext } from '@molstar/graphics/gl/context';
import type { ShaderCode } from '@molstar/graphics/gl/shader-code';
import { SkyboxShaderCode } from './background.js';
import { BlurShaderCode, CompositeShaderCode, LuminosityShaderCode } from './bloom.js';
import { CasShaderCode } from './cas.js';
import { DofShaderCode } from './dof.js';
import { FxaaShaderCode } from './fxaa.js';
import { DimShaderCode, type MarkingProps, OverlayShaderCode } from './marking.js';
import { ComposeShaderCode, type MultiSampleProps } from './multi-sample.js';
import { OutlinesShaderCode } from './outline.js';
import type { PostprocessingProps } from './postprocessing.js';
import { ShadowsShaderCode } from './shadow.js';
import { BlendShaderCode, EdgesShaderCode, WeightsShaderCode } from './smaa.js';

export type SupportedProps = {
  marking: MarkingProps;
  multiSample: MultiSampleProps;
  postprocessing: PostprocessingProps;
};

const Off = { name: 'off', params: {} } as const;

function supports(ctx: GraphicsContext, ...shaderCodes: ShaderCode[]) {
  return shaderCodes.every((c) => ctx.isProgramSupported(c));
}

/**
 * Turns off the features of `props` whose programs `ctx` doesn't implement, replacing the changed sub-objects.
 * Returns the names of the turned off features.
 */
export function limitToSupported(ctx: GraphicsContext, props: SupportedProps): string[] {
  const off: string[] = [];

  if (props.marking.enabled && !supports(ctx, OverlayShaderCode, DimShaderCode)) {
    props.marking = { ...props.marking, enabled: false };
    off.push('marking');
  }
  if (props.multiSample.mode !== 'off' && !supports(ctx, ComposeShaderCode)) {
    props.multiSample = { ...props.multiSample, mode: 'off' };
    off.push('multi-sampling');
  }

  const pp = props.postprocessing;
  const ppOff: Partial<PostprocessingProps> = {};
  if (pp.shadow.name !== 'off' && !supports(ctx, ShadowsShaderCode)) {
    ppOff.shadow = Off;
    off.push('shadows');
  }
  if (pp.outline.name !== 'off' && !supports(ctx, OutlinesShaderCode)) {
    ppOff.outline = Off;
    off.push('outlines');
  }
  if (pp.dof.name !== 'off' && !supports(ctx, DofShaderCode)) {
    ppOff.dof = Off;
    off.push('depth of field');
  }
  if (pp.bloom.name !== 'off' && !supports(ctx, LuminosityShaderCode, BlurShaderCode, CompositeShaderCode)) {
    ppOff.bloom = Off;
    off.push('bloom');
  }
  if (
    (pp.antialiasing.name === 'fxaa' && !supports(ctx, FxaaShaderCode)) ||
    (pp.antialiasing.name === 'smaa' && !supports(ctx, EdgesShaderCode, WeightsShaderCode, BlendShaderCode))
  ) {
    off.push(pp.antialiasing.name);
    ppOff.antialiasing = Off;
  }
  if (pp.sharpening.name !== 'off' && !supports(ctx, CasShaderCode)) {
    ppOff.sharpening = Off;
    off.push('sharpening');
  }
  if (pp.background.variant.name !== 'off' && !supports(ctx, SkyboxShaderCode)) {
    ppOff.background = { ...pp.background, variant: Off };
    off.push('background');
  }
  if (Object.keys(ppOff).length > 0) props.postprocessing = { ...pp, ...ppOff };

  return off;
}

/** `limitToSupported` with a console note on the features turned off for the first time */
export function createSupportedPropsLimiter(ctx: GraphicsContext) {
  const noted = new Set<string>();
  return (props: SupportedProps) => {
    const off = limitToSupported(ctx, props).filter((f) => !noted.has(f));
    if (off.length === 0) return;
    for (const f of off) noted.add(f);
    console.info(`Not supported by the graphics backend, turned off: ${off.join(', ')}`);
  };
}
