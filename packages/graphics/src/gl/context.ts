/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { Subject } from 'rxjs';
import type { ValueCell } from '@molstar/core/util';
import type { PendingWork } from './pending-work.js';
import type { now } from '@molstar/core/util/now';
import type { GraphicsCapabilities } from './capabilities.js';
import type { Framebuffer } from './framebuffer.js';
import type { RenderTarget } from './render-target.js';
import type {
    ComputeRenderItem,
    DrawMode,
    GlobalBindings,
    GraphicsRenderItem,
    ProgramVariant,
    Transparency,
} from './render-item.js';
import type { ComputeRenderable } from './renderable.js';
import type { RenderableSchema, RenderableValues } from './renderable/schema.js';
import type { Renderbuffer, RenderbufferAttachment, RenderbufferFormat } from './renderbuffer.js';
import type { ClearValues, RenderState } from './render-state.js';
import type { ShaderCode } from './shader-code.js';
import type { GraphicsStats, ResourceByteCounts } from './stats.js';
import type {
    CubeFaces,
    Texture,
    TextureFilter,
    TextureFormat,
    TextureKind,
    Textures,
    TextureType,
} from './texture.js';
import type { GraphicsTimer } from './timer.js';
import type { UniformBuffer, UniformsList } from './uniform.js';
import type { UniformBlockLayout } from './uniform-layout.js';

export interface GraphicsLimits {
    readonly maxTextureSize: number;
    /** 0 without `texture3d` */
    readonly max3dTextureSize: number;
    readonly maxRenderbufferSize: number;
    /** 0 without `multipleRenderTargets` */
    readonly maxDrawBuffers: number;
    readonly maxTextureImageUnits: number;
}

export interface GraphicsResources {
    readonly limits: GraphicsLimits;

    texture: (kind: TextureKind, format: TextureFormat, type: TextureType, filter: TextureFilter) => Texture;
    cubeTexture: (faces: CubeFaces, mipmaps: boolean, onload?: () => void) => Texture;
    framebuffer: () => Framebuffer;
    renderbuffer: (
        format: RenderbufferFormat,
        attachment: RenderbufferAttachment,
        width: number,
        height: number,
    ) => Renderbuffer;
    uniformBuffer: (layout: UniformBlockLayout) => UniformBuffer;
    /** writes a new snapshot on every change instead of updating in place */
    uniformRingBuffer: (layout: UniformBlockLayout) => UniformBuffer;
    /** cached per `materialId`, shared by render items with the same material */
    materialUniformBuffer: (materialId: number, layout: UniformBlockLayout) => UniformBuffer;

    getByteCounts: () => ResourceByteCounts;

    linkPrograms: (variants?: ProgramVariant[]) => void;
    finalizePrograms: (variants?: ProgramVariant[], isSynchronous?: boolean) => boolean;

    /** Asynchronous work frames may depend on; `ImagePass` waits for it, `Canvas3D` redraws when it settles */
    readonly pending: PendingWork;
    /**
     * Loads the image at `source` (a URL or data URL) into the rgba `ubyte` `texture` with the raw pixel values (no
     * color space conversion), then updates `cell`; tracked in `pending`
     */
    loadImage: (texture: Texture, source: string, cell?: ValueCell<Texture>, label?: string) => Promise<void>;
}

export type RenderTargetDepthStencil = 'none' | 'depth' | 'depth-stencil';
export type RenderTargetType = 'uint8' | 'float32' | 'fp16';
export type RenderTargetFormat = 'rgba' | 'r';

/** Target of a render pass, `drawing-buffer` is the canvas or XR layer */
export type RenderPassTarget = RenderTarget | Framebuffer | 'drawing-buffer';

export interface RenderPassOptions {
    /** Clears the attachments at the start of the pass, like `clear` right after `beginRenderPass` (current scissor) */
    readonly clear?: ClearValues;
    /** Name of the pass in timer results, derived from the target if not given */
    readonly label?: string;
}

export interface GraphicsXR {
    readonly session?: XRSession;
    readonly changed: Subject<void>;
    set: (session: XRSession | undefined, options?: { resolutionScale?: number }) => Promise<void>;
    end: () => Promise<void>;
}

/** Backend-neutral graphics context, implemented by `WebGLContext` */
export interface GraphicsContext {
    readonly capabilities: GraphicsCapabilities;
    readonly resources: GraphicsResources;
    readonly stats: GraphicsStats;
    readonly timer: GraphicsTimer;

    readonly pixelRatio: number;
    setPixelScale: (value: number) => void;

    readonly isContextLost: boolean;
    /** Emits when the context is lost (`setContextLost`), e.g. on a WebGL context loss event or a lost WebGPU device */
    readonly contextLost: Subject<now.Timestamp>;
    readonly contextRestored: Subject<now.Timestamp>;
    setContextLost: () => void;
    /** Recreates the resources after the platform restored the context (WebGL), then `extraResets` */
    handleContextRestored: (extraResets?: () => void) => void;
    /**
     * Recovers from a loss where the backend replaces the context itself (WebGPU: a new device, the resources recreated on
     * it), then `extraResets`; resolves whether the context was restored. Undefined where the platform restores (WebGL).
     */
    readonly restore?: (extraResets?: () => void) => Promise<boolean>;

    readonly xr: GraphicsXR;

    /** Cache for frambuffers, managed by consumers */
    readonly namedFramebuffers: { [name: string]: Framebuffer };
    /** Cache for textures, managed by consumers */
    readonly namedTextures: { [name: string]: Texture };
    /** Cache for compute renderables, managed by consumers */
    readonly namedComputeRenderables: { [name: string]: ComputeRenderable<any> };

    /** Whether the backend implements `shaderCode`; passes whose programs it doesn't are unsupported and fall back */
    isProgramSupported: (shaderCode: ShaderCode) => boolean;
    createGraphicsRenderItem: (
        drawMode: DrawMode,
        shaderCode: ShaderCode,
        schema: RenderableSchema,
        values: RenderableValues,
        materialId: number,
        transparency: Transparency,
    ) => GraphicsRenderItem;
    createComputeRenderItem: (
        drawMode: DrawMode,
        shaderCode: ShaderCode,
        schema: RenderableSchema,
        values: RenderableValues,
        materialId?: number,
    ) => ComputeRenderItem;
    /** Uniforms (`schema` and `uniforms`) and textures shared by the render items of a pass */
    createGlobalBindings: (schema: RenderableSchema, uniforms: UniformsList, textures: Textures) => GlobalBindings;

    createRenderTarget: (
        width: number,
        height: number,
        depthStencil?: RenderTargetDepthStencil,
        type?: RenderTargetType,
        filter?: TextureFilter,
        format?: RenderTargetFormat,
    ) => RenderTarget;
    /** Render target of the drawing buffer (canvas or XR layer) */
    createDrawTarget: () => RenderTarget;
    getDrawingBufferSize: () => { width: number; height: number };

    /**
     * Starts a render pass on `target`, ending the open one. Draws and clears need an open pass.
     * Copies, readback and attachment changes end the open pass (WebGPU runs them outside of passes).
     */
    beginRenderPass: (target: RenderPassTarget, options?: RenderPassOptions) => void;
    endRenderPass: () => void;
    /** Whether the target of the open render pass has a stencil buffer */
    hasRenderPassStencil: () => boolean;

    /** Fixed-function state for the following draws, unset fields are reset to their defaults */
    setRenderState: (state: RenderState) => void;
    /** Used by the `constant` blend factors */
    setBlendConstant: (red: number, green: number, blue: number, alpha: number) => void;
    setViewport: (x: number, y: number, width: number, height: number) => void;
    /** Limits draws and clears to the rectangle, until `disableScissor` */
    setScissor: (x: number, y: number, width: number, height: number) => void;
    /** Draws and clears cover the whole target */
    disableScissor: () => void;
    /**
     * Clears the attachments of the target of the open render pass, within the scissor rectangle and ignoring the
     * write masks. On WebGPU, clears of the whole target before the first draw of a pass are load operations, all
     * others are draws.
     */
    clear: (values: ClearValues) => void;
    /** Ends the open render pass and submits pending commands */
    flush: () => void;
    /** Resolves when the GPU finished the submitted commands */
    waitForGpuCommandsComplete: () => Promise<void>;
    /**
     * Reads the color attachment of `target` as rgba8 in a rectangle with a bottom-left origin, rows top to bottom.
     * Ends the open render pass.
     */
    readRenderTarget: (target: RenderTarget, x: number, y: number, width: number, height: number) => Promise<Uint8Array>;
    /**
     * Reads an rgba8 2d texture in a rectangle with a bottom-left origin, rows top to bottom; the read is ordered after
     * the commands recorded so far. Ends the open render pass.
     */
    readTextureAsync: (texture: Texture, x: number, y: number, width: number, height: number) => Promise<Uint8Array>;

    destroy: () => void;
}
