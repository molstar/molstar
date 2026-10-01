/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { Scene } from '@molstar/graphics/gl/scene';
import type { WebGLContext } from '@molstar/graphics/gl/webgl/context';
import type { GraphicsRenderVariant } from '@molstar/graphics/gl/webgl/render-item';
import { BloomPass } from '../passes/bloom.js';
import { IlluminationPass, type IlluminationProps } from '../passes/illumination.js';
import { MarkingPass, type MarkingProps } from '../passes/marking.js';
import { PostprocessingPass, type PostprocessingProps } from '../passes/postprocessing.js';

export type ShaderManagerProps = {
    marking: MarkingProps
    postprocessing: PostprocessingProps
    illumination: IlluminationProps
}

export class ShaderManager {
    static ensureRequired(webgl: WebGLContext, scene: Scene, p: ShaderManagerProps) {
        const sm = new ShaderManager(webgl, scene);
        sm.updateRequired(p);
        sm.finalizeRequired(true);
    }

    private readonly required: GraphicsRenderVariant[] = [];

    constructor(private readonly webgl: WebGLContext, private readonly scene: Scene) { }

    updateRequired(p: ShaderManagerProps) {
        this.required.length = 0;
        this.required.push('color');
        if (IlluminationPass.isEnabled(this.webgl, p.illumination)) {
            this.required.push('tracing');
        }
        if (MarkingPass.isEnabled(p.marking) && this.scene.markerAverage > 0) {
            this.required.push('marking');
        }
        if (BloomPass.isEnabled(p.postprocessing) && this.scene.emissiveAverage > 0) {
            this.required.push('emissive');
        }
        if (PostprocessingPass.isTransparentDepthRequired(this.scene, p.postprocessing) || !this.webgl.extensions.drawBuffers || !this.webgl.extensions.depthTexture || IlluminationPass.isEnabled(this.webgl, p.illumination)) {
            this.required.push('depth');
        }
        this.webgl.resources.linkPrograms(this.required);
    }

    finalizeRequired(isSynchronous?: boolean) {
        return this.finalize(this.required, isSynchronous);
    }

    finalize(variants?: GraphicsRenderVariant[], isSynchronous?: boolean) {
        return this.webgl.resources.finalizePrograms(variants, isSynchronous);
    }
}
