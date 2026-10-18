/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { Scene } from '@molstar/graphics/gl/scene';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import type { GraphicsRenderVariant } from '@molstar/graphics/gl/render-item';
import { BloomPass } from '../passes/bloom.js';
import { IlluminationPass, type IlluminationProps } from '../passes/illumination.js';
import { MarkingPass, type MarkingProps } from '../passes/marking.js';
import { PostprocessingPass, type PostprocessingProps } from '../passes/postprocessing.js';

export type ShaderManagerProps = {
  marking: MarkingProps;
  postprocessing: PostprocessingProps;
  illumination: IlluminationProps;
};

export class ShaderManager {
  static ensureRequired(graphics: GraphicsContext, scene: Scene, p: ShaderManagerProps) {
    const sm = new ShaderManager(graphics, scene);
    sm.updateRequired(p);
    sm.finalizeRequired(true);
  }

  private readonly required: GraphicsRenderVariant[] = [];

  constructor(
    private readonly graphics: GraphicsContext,
    private readonly scene: Scene,
  ) {}

  /** The variants frames with `p` draw the scene with (marking only while something is marked) */
  static getRequired(graphics: GraphicsContext, scene: Scene, p: ShaderManagerProps): GraphicsRenderVariant[] {
    const required: GraphicsRenderVariant[] = ['color'];
    if (IlluminationPass.isEnabled(graphics, p.illumination)) {
      required.push('tracing');
    }
    if (MarkingPass.isEnabled(p.marking) && scene.markerAverage > 0) {
      required.push('marking');
    }
    if (BloomPass.isEnabled(p.postprocessing) && scene.emissiveAverage > 0) {
      required.push('emissive');
    }
    if (
      (PostprocessingPass.isTransparentDepthRequired(scene, p.postprocessing) && scene.opacityAverage < 1) ||
      !graphics.capabilities.multipleRenderTargets ||
      !graphics.capabilities.depthTexture ||
      IlluminationPass.isEnabled(graphics, p.illumination)
    ) {
      required.push('depth');
    }
    return required;
  }

  updateRequired(p: ShaderManagerProps) {
    this.required.length = 0;
    this.required.push(...ShaderManager.getRequired(this.graphics, this.scene, p));
    this.graphics.resources.linkPrograms(this.required);
  }

  finalizeRequired(isSynchronous?: boolean) {
    return this.finalize(this.required, isSynchronous);
  }

  finalize(variants?: GraphicsRenderVariant[], isSynchronous?: boolean) {
    return this.graphics.resources.finalizePrograms(variants, isSynchronous);
  }
}
