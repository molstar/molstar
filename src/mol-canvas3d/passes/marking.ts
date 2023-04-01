/**
 * Copyright (c) 2021-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { QuadSchema, QuadValues } from '../../mol-gl/compute/util';
import { ComputeRenderable, createComputeRenderable } from '../../mol-gl/renderable';
import { DefineSpec, TextureSpec, UniformSpec, Values } from '../../mol-gl/renderable/schema';
import { ShaderCode } from '../../mol-gl/shader-code';
import { WebGLContext } from '../../mol-gl/webgl/context';
import { createComputeRenderItem } from '../../mol-gl/webgl/render-item';
import { Texture } from '../../mol-gl/webgl/texture';
import { Vec2, Vec3 } from '../../mol-math/linear-algebra';
import { ValueCell } from '../../mol-util';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { quad_vert } from '../../mol-gl/shader/quad.vert';
import { overlay_frag } from '../../mol-gl/shader/marking/overlay.frag';
import { Viewport } from '../camera/util';
import { RenderTarget } from '../../mol-gl/webgl/render-target';
import { Color } from '../../mol-util/color';
import { edge_frag } from '../../mol-gl/shader/marking/edge.frag';
import { isTimingMode } from '../../mol-util/debug';

export const MarkingParams = {
    enabled: PD.Boolean(true),
    highlightEdgeColor: PD.Color(Color.darken(Color.fromNormalizedRgb(1.0, 0.4, 0.6), 1.0)),
    selectEdgeColor: PD.Color(Color.darken(Color.fromNormalizedRgb(0.2, 1.0, 0.1), 1.0)),
    edgeScale: PD.Numeric(1, { min: 1, max: 3, step: 1 }, { description: 'Thickness of the edge.' }),
    highlightEdgeStrength: PD.Numeric(1.0, { min: 0, max: 1, step: 0.1 }),
    selectEdgeStrength: PD.Numeric(1.0, { min: 0, max: 1, step: 0.1 }),
    ghostEdgeStrength: PD.Numeric(0.3, { min: 0, max: 1, step: 0.1 }, { description: 'Opacity of the hidden edges that are covered by other geometry. When set to 1, one less geometry render pass is done.' }),
    innerEdgeFactor: PD.Numeric(1.5, { min: 0, max: 3, step: 0.1 }, { description: 'Factor to multiply the inner edge color with - for added contrast.' }),
};
export type MarkingProps = PD.Values<typeof MarkingParams>

export class MarkingPass {
    static isEnabled(props: MarkingProps) {
        return props.enabled;
    }

    readonly depthTarget: RenderTarget;
    readonly maskTarget: RenderTarget;
    private readonly edgesTarget: RenderTarget;

    private readonly edge: EdgeRenderable;
    private readonly overlay: OverlayRenderable;

    constructor(private webgl: WebGLContext, width: number, height: number) {
        this.depthTarget = webgl.createRenderTarget(width, height);
        this.maskTarget = webgl.createRenderTarget(width, height);
        this.edgesTarget = webgl.createRenderTarget(width, height);

        this.edge = getEdgeRenderable(webgl, this.maskTarget.texture);
        this.overlay = getOverlayRenderable(webgl, this.edgesTarget.texture);
    }

    private setEdgeState(viewport: Viewport) {
        const { gl, state } = this.webgl;

        state.enable(gl.SCISSOR_TEST);
        state.enable(gl.BLEND);
        state.blendFunc(gl.ONE, gl.ONE);
        state.blendEquation(gl.FUNC_ADD);
        state.disable(gl.DEPTH_TEST);
        state.depthMask(false);

        const { x, y, width, height } = viewport;
        state.viewport(x, y, width, height);
        state.scissor(x, y, width, height);

        state.clearColor(0, 0, 0, 0);
        gl.clear(gl.COLOR_BUFFER_BIT);
    }

    private setOverlayState(viewport: Viewport) {
        const { gl, state } = this.webgl;

        state.enable(gl.SCISSOR_TEST);
        state.enable(gl.BLEND);
        state.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA);
        state.blendEquation(gl.FUNC_ADD);
        state.disable(gl.DEPTH_TEST);
        state.depthMask(false);

        const { x, y, width, height } = viewport;
        state.viewport(x, y, width, height);
        state.scissor(x, y, width, height);
    }

    setSize(width: number, height: number) {
        const w = this.depthTarget.getWidth();
        const h = this.depthTarget.getHeight();

        if (width !== w || height !== h) {
            this.depthTarget.setSize(width, height);
            this.maskTarget.setSize(width, height);
            this.edgesTarget.setSize(width, height);

            ValueCell.update(this.edge.values.uTexSizeInv, Vec2.set(this.edge.values.uTexSizeInv.ref.value, 1 / width, 1 / height));
            ValueCell.update(this.overlay.values.uTexSizeInv, Vec2.set(this.overlay.values.uTexSizeInv.ref.value, 1 / width, 1 / height));
        }
    }

    update(props: MarkingProps) {
        const { highlightEdgeColor, selectEdgeColor, edgeScale, innerEdgeFactor, ghostEdgeStrength, highlightEdgeStrength, selectEdgeStrength } = props;

        const { values: edgeValues } = this.edge;
        const _edgeScale = Math.max(1, Math.round(edgeScale * this.webgl.pixelRatio));
        if (edgeValues.dEdgeScale.ref.value !== _edgeScale) {
            ValueCell.update(edgeValues.dEdgeScale, _edgeScale);
            this.edge.update();
        }

        const { values: overlayValues } = this.overlay;
        ValueCell.update(overlayValues.uHighlightEdgeColor, Color.toVec3Normalized(overlayValues.uHighlightEdgeColor.ref.value, highlightEdgeColor));
        ValueCell.update(overlayValues.uSelectEdgeColor, Color.toVec3Normalized(overlayValues.uSelectEdgeColor.ref.value, selectEdgeColor));
        ValueCell.updateIfChanged(overlayValues.uInnerEdgeFactor, innerEdgeFactor);
        ValueCell.updateIfChanged(overlayValues.uGhostEdgeStrength, ghostEdgeStrength);
        ValueCell.updateIfChanged(overlayValues.uHighlightEdgeStrength, highlightEdgeStrength);
        ValueCell.updateIfChanged(overlayValues.uSelectEdgeStrength, selectEdgeStrength);
    }

    render(viewport: Viewport, target: RenderTarget | undefined) {
        if (isTimingMode) this.webgl.timer.mark('MarkingPass.render');
        this.edgesTarget.bind();
        this.setEdgeState(viewport);
        this.edge.render();

        if (target) {
            target.bind();
        } else {
            this.webgl.unbindFramebuffer();
        }
        this.setOverlayState(viewport);
        this.overlay.render();
        if (isTimingMode) this.webgl.timer.markEnd('MarkingPass.render');
    }
}

//

const EdgeSchema = {
    ...QuadSchema,
    tMaskTexture: TextureSpec('texture', 'rgba', 'ubyte', 'linear'),
    uTexSizeInv: UniformSpec('v2'),
    dEdgeScale: DefineSpec('number'),
};
const EdgeShaderCode = ShaderCode('edge', quad_vert, edge_frag);
type EdgeRenderable = ComputeRenderable<Values<typeof EdgeSchema>>

function getEdgeRenderable(ctx: WebGLContext, maskTexture: Texture): EdgeRenderable {
    const width = maskTexture.getWidth();
    const height = maskTexture.getHeight();

    const values: Values<typeof EdgeSchema> = {
        ...QuadValues,
        tMaskTexture: ValueCell.create(maskTexture),
        uTexSizeInv: ValueCell.create(Vec2.create(1 / width, 1 / height)),
        dEdgeScale: ValueCell.create(1),
    };

    const schema = { ...EdgeSchema };
    const renderItem = createComputeRenderItem(ctx, 'triangles', EdgeShaderCode, schema, values);

    return createComputeRenderable(renderItem, values);
}

//

const OverlaySchema = {
    ...QuadSchema,
    tEdgeTexture: TextureSpec('texture', 'rgba', 'ubyte', 'linear'),
    uTexSizeInv: UniformSpec('v2'),
    uHighlightEdgeColor: UniformSpec('v3'),
    uSelectEdgeColor: UniformSpec('v3'),
    uHighlightEdgeStrength: UniformSpec('f'),
    uSelectEdgeStrength: UniformSpec('f'),
    uGhostEdgeStrength: UniformSpec('f'),
    uInnerEdgeFactor: UniformSpec('f'),
};
const OverlayShaderCode = ShaderCode('overlay', quad_vert, overlay_frag);
type OverlayRenderable = ComputeRenderable<Values<typeof OverlaySchema>>

function getOverlayRenderable(ctx: WebGLContext, edgeTexture: Texture): OverlayRenderable {
    const width = edgeTexture.getWidth();
    const height = edgeTexture.getHeight();

    const values: Values<typeof OverlaySchema> = {
        ...QuadValues,
        tEdgeTexture: ValueCell.create(edgeTexture),
        uTexSizeInv: ValueCell.create(Vec2.create(1 / width, 1 / height)),
        uHighlightEdgeColor: ValueCell.create(Vec3()),
        uSelectEdgeColor: ValueCell.create(Vec3()),
        uHighlightEdgeStrength: ValueCell.create(1),
        uSelectEdgeStrength: ValueCell.create(1),
        uGhostEdgeStrength: ValueCell.create(0),
        uInnerEdgeFactor: ValueCell.create(0),
    };

    const schema = { ...OverlaySchema };
    const renderItem = createComputeRenderItem(ctx, 'triangles', OverlayShaderCode, schema, values);

    return createComputeRenderable(renderItem, values);
}