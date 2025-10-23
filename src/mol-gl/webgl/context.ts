/**
 * Copyright (c) 2018-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { GLRenderingContext, isWebGL2 } from './compat';
import { checkFramebufferStatus, createNullFramebuffer, Framebuffer } from './framebuffer';
import { Scheduler } from '../../mol-task';
import { isDebugMode } from '../../mol-util/debug';
import { createExtensions, resetExtensions, WebGLExtensions } from './extensions';
import { WebGLState, createState } from './state';
import { WebGLResources, createResources } from './resources';
import { RenderTarget, createRenderTarget } from './render-target';
import { Subject } from 'rxjs';
import { now } from '../../mol-util/now';
import { createNullTexture, Texture, TextureFilter } from './texture';
import { ComputeRenderable } from '../renderable';
import { createTimer, WebGLTimer } from './timer';

export function getGLContext(canvas: HTMLCanvasElement, attribs?: WebGLContextAttributes & { preferWebGl1?: boolean }): GLRenderingContext | null {
    function get(id: 'webgl' | 'experimental-webgl' | 'webgl2') {
        try {
            return canvas.getContext(id, attribs) as GLRenderingContext | null;
        } catch (e) {
            return null;
        }
    }
    const gl = (attribs?.preferWebGl1 ? null : get('webgl2')) || get('webgl') || get('experimental-webgl');
    if (isDebugMode) console.log(`isWebgl2: ${isWebGL2(gl)}`);
    return gl;
}

export function getErrorDescription(gl: GLRenderingContext, error: number) {
    switch (error) {
        case gl.NO_ERROR: return 'no error';
        case gl.INVALID_ENUM: return 'invalid enum';
        case gl.INVALID_VALUE: return 'invalid value';
        case gl.INVALID_OPERATION: return 'invalid operation';
        case gl.INVALID_FRAMEBUFFER_OPERATION: return 'invalid framebuffer operation';
        case gl.OUT_OF_MEMORY: return 'out of memory';
        case gl.CONTEXT_LOST_WEBGL: return 'context lost';
    }
    return 'unknown error';
}

export function checkError(gl: GLRenderingContext, message?: string) {
    const error = gl.getError();
    if (error !== gl.NO_ERROR) {
        console.log(`WebGL error: '${getErrorDescription(gl, error)}'${message ? ` (${message})` : ''}`);
        // throw new Error(`WebGL error: '${getErrorDescription(gl, error)}'${message ? ` (${message})` : ''}`);
    }
}

export function glEnumToString(gl: GLRenderingContext, value: number) {
    const keys: string[] = [];
    for (const key in gl) {
        if ((gl as any)[key] === value) {
            keys.push(key);
        }
    }
    return keys.length ? keys.join(' | ') : `0x${value.toString(16)}`;
}

function unbindResources(gl: GLRenderingContext) {
    // bind null to all texture units
    const maxTextureImageUnits = gl.getParameter(gl.MAX_TEXTURE_IMAGE_UNITS);
    for (let i = 0; i < maxTextureImageUnits; ++i) {
        gl.activeTexture(gl.TEXTURE0 + i);
        gl.bindTexture(gl.TEXTURE_2D, null);
        gl.bindTexture(gl.TEXTURE_CUBE_MAP, null);
        if (isWebGL2(gl)) {
            gl.bindTexture(gl.TEXTURE_2D_ARRAY, null);
            gl.bindTexture(gl.TEXTURE_3D, null);
        }
    }

    // assign the smallest possible buffer to all attributes
    const buf = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, buf);
    const maxVertexAttribs = gl.getParameter(gl.MAX_VERTEX_ATTRIBS);
    for (let i = 0; i < maxVertexAttribs; ++i) {
        gl.vertexAttribPointer(i, 1, gl.FLOAT, false, 0, 0);
    }

    // bind null to all buffers
    gl.bindBuffer(gl.ARRAY_BUFFER, null);
    gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, null);
    gl.bindRenderbuffer(gl.RENDERBUFFER, null);
    gl.bindFramebuffer(gl.FRAMEBUFFER, null);
    if (isWebGL2(gl)) {
        gl.bindBuffer(gl.UNIFORM_BUFFER, null);
        gl.bindBuffer(gl.PIXEL_PACK_BUFFER, null);
    }
}

const tmpPixel = new Uint8Array(1 * 4);

function checkSync(gl: WebGL2RenderingContext, sync: WebGLSync, resolve: () => void) {
    if (gl.getSyncParameter(sync, gl.SYNC_STATUS) === gl.SIGNALED) {
        gl.deleteSync(sync);
        resolve();
    } else {
        Scheduler.setImmediate(checkSync, gl, sync, resolve);
    }
}

function fence(gl: WebGL2RenderingContext, resolve: () => void) {
    const sync = gl.fenceSync(gl.SYNC_GPU_COMMANDS_COMPLETE, 0);
    if (!sync) {
        console.warn('Could not create a WebGLSync object');
        gl.readPixels(0, 0, 1, 1, gl.RGBA, gl.UNSIGNED_BYTE, tmpPixel);
        resolve();
    } else {
        Scheduler.setImmediate(checkSync, gl, sync, resolve);
    }
}

let SentWebglSyncObjectNotSupportedInWebglMessage = false;
function waitForGpuCommandsComplete(gl: GLRenderingContext): Promise<void> {
    return new Promise(resolve => {
        if (isWebGL2(gl)) {
            fence(gl, resolve);
        } else {
            if (!SentWebglSyncObjectNotSupportedInWebglMessage) {
                console.info('Sync object not supported in WebGL');
                SentWebglSyncObjectNotSupportedInWebglMessage = true;
            }
            waitForGpuCommandsCompleteSync(gl);
            resolve();
        }
    });
}

function waitForGpuCommandsCompleteSync(gl: GLRenderingContext): void {
    gl.bindFramebuffer(gl.FRAMEBUFFER, null);
    gl.readPixels(0, 0, 1, 1, gl.RGBA, gl.UNSIGNED_BYTE, tmpPixel);
}

export function readPixels(gl: GLRenderingContext, x: number, y: number, width: number, height: number, buffer: Uint8Array | Float32Array | Int32Array) {
    if (isDebugMode) checkFramebufferStatus(gl);
    if (buffer instanceof Uint8Array) {
        gl.readPixels(x, y, width, height, gl.RGBA, gl.UNSIGNED_BYTE, buffer);
    } else if (buffer instanceof Float32Array) {
        gl.readPixels(x, y, width, height, gl.RGBA, gl.FLOAT, buffer);
    } else if (buffer instanceof Int32Array && isWebGL2(gl)) {
        gl.readPixels(x, y, width, height, gl.RGBA_INTEGER, gl.INT, buffer);
    } else {
        throw new Error('unsupported readPixels buffer type');
    }
    if (isDebugMode) checkError(gl);
}

function bindDrawingBuffer(gl: GLRenderingContext, xrLayer?: XRWebGLLayer) {
    if (xrLayer) {
        gl.bindFramebuffer(gl.FRAMEBUFFER, xrLayer.framebuffer);
    } else {
        gl.bindFramebuffer(gl.FRAMEBUFFER, null);
    }
}

function getDrawingBufferSize(gl: GLRenderingContext, xrLayer?: XRWebGLLayer, xrInteractionMode?: XRInteractionMode) {
    let width = xrLayer?.framebufferWidth ?? gl.drawingBufferWidth;
    if (xrInteractionMode === 'screen-space') {
        // workaround so XR with a single view behaves simlar to two views
        width *= 2;
    }
    const height = xrLayer?.framebufferHeight ?? gl.drawingBufferHeight;
    return { width, height };
}

function getShaderPrecisionFormat(gl: GLRenderingContext, shader: 'vertex' | 'fragment', precision: 'low' | 'medium' | 'high', type: 'float' | 'int') {
    const glShader = shader === 'vertex' ? gl.VERTEX_SHADER : gl.FRAGMENT_SHADER;
    const glPrecisionType = gl[`${precision.toUpperCase()}_${type.toUpperCase()}` as 'LOW_FLOAT' | 'MEDIUM_FLOAT' | 'HIGH_FLOAT' | 'LOW_INT' | 'MEDIUM_INT' | 'HIGH_INT'];
    return gl.getShaderPrecisionFormat(glShader, glPrecisionType);
}

function getShaderPrecisionFormats(gl: GLRenderingContext, shader: 'vertex' | 'fragment') {
    return {
        lowFloat: getShaderPrecisionFormat(gl, shader, 'low', 'float'),
        mediumFloat: getShaderPrecisionFormat(gl, shader, 'medium', 'float'),
        highFloat: getShaderPrecisionFormat(gl, shader, 'high', 'float'),
        lowInt: getShaderPrecisionFormat(gl, shader, 'low', 'int'),
        mediumInt: getShaderPrecisionFormat(gl, shader, 'medium', 'int'),
        highInt: getShaderPrecisionFormat(gl, shader, 'high', 'int'),
    };
}

type WebGLShaderPrecisionFormats = ReturnType<typeof getShaderPrecisionFormats>

//

function createStats() {
    const stats = {
        resourceCounts: {
            attribute: 0,
            elements: 0,
            pixelPack: 0,
            framebuffer: 0,
            program: 0,
            renderbuffer: 0,
            shader: 0,
            texture: 0,
            cubeTexture: 0,
            vertexArray: 0,
        },

        drawCount: 0,
        instanceCount: 0,
        instancedDrawCount: 0,

        calls: {
            drawInstanced: 0,
            drawInstancedBase: 0,
            multiDrawInstancedBase: 0,
            counts: 0,
        },

        culled: {
            lod: 0,
            frustum: 0,
            occlusion: 0,
        },
    };
    return stats;
}

export type WebGLStats = ReturnType<typeof createStats>

//

function createParameters(gl: GLRenderingContext, extensions: WebGLExtensions) {
    return {
        maxTextureSize: gl.getParameter(gl.MAX_TEXTURE_SIZE) as number,
        max3dTextureSize: isWebGL2(gl) ? gl.getParameter(gl.MAX_3D_TEXTURE_SIZE) as number : 0,
        maxRenderbufferSize: gl.getParameter(gl.MAX_RENDERBUFFER_SIZE) as number,
        maxDrawBuffers: extensions.drawBuffers ? gl.getParameter(extensions.drawBuffers.MAX_DRAW_BUFFERS) as number : 0,
        maxTextureImageUnits: gl.getParameter(gl.MAX_TEXTURE_IMAGE_UNITS) as number,
        maxVertexTextureImageUnits: gl.getParameter(gl.MAX_VERTEX_TEXTURE_IMAGE_UNITS) as number,
    };
}

export type WebGLParameters = ReturnType<typeof createParameters>

//

/** A WebGL context object, including the rendering context, resource caches and counts */
export interface WebGLContext {
    readonly gl: GLRenderingContext
    readonly isWebGL2: boolean
    readonly pixelRatio: number

    readonly extensions: WebGLExtensions
    readonly state: WebGLState
    readonly stats: WebGLStats
    readonly resources: WebGLResources
    readonly timer: WebGLTimer

    readonly maxTextureSize: number
    readonly max3dTextureSize: number
    readonly maxRenderbufferSize: number
    readonly maxDrawBuffers: number
    readonly maxTextureImageUnits: number
    readonly shaderPrecisionFormats: { vertex: WebGLShaderPrecisionFormats, fragment: WebGLShaderPrecisionFormats }

    readonly isContextLost: boolean
    readonly contextRestored: Subject<now.Timestamp>
    setContextLost: () => void
    handleContextRestored: (extraResets?: () => void) => void

    readonly xr: {
        readonly session?: XRSession
        readonly changed: Subject<void>
        set: (session: XRSession | undefined, options?: { resolutionScale?: number }) => Promise<void>
        end: () => Promise<void>
    }

    setPixelScale: (value: number) => void

    /** Cache for compute renderables, managed by consumers */
    readonly namedComputeRenderables: { [name: string]: ComputeRenderable<any> }
    /** Cache for frambuffers, managed by consumers */
    readonly namedFramebuffers: { [name: string]: Framebuffer }
    /** Cache for textures, managed by consumers */
    readonly namedTextures: { [name: string]: Texture }

    createRenderTarget: (width: number, height: number, depth?: boolean, type?: 'uint8' | 'float32' | 'fp16', filter?: TextureFilter, format?: 'rgba' | 'alpha') => RenderTarget
    createDrawTarget: () => RenderTarget
    bindDrawingBuffer: () => void
    getDrawingBufferSize: () => { width: number, height: number }
    readPixels: (x: number, y: number, width: number, height: number, buffer: Uint8Array | Float32Array | Int32Array) => void
    waitForGpuCommandsComplete: () => Promise<void>
    waitForGpuCommandsCompleteSync: () => void
    getFenceSync: () => WebGLSync | null
    checkSyncStatus: (sync: WebGLSync) => boolean
    deleteSync: (sync: WebGLSync) => void
    clear: (red: number, green: number, blue: number, alpha: number) => void
    checkError: (message?: string) => void
    checkFramebufferStatus: (message?: string) => void
    destroy: (options?: Partial<{ doNotForceWebGLContextLoss: boolean }>) => void
}

export function createContext(gl: GLRenderingContext, props: Partial<{ pixelScale: number }> = {}): WebGLContext {
    const extensions = createExtensions(gl);
    const state = createState(gl, extensions);
    const stats = createStats();
    const parameters = createParameters(gl, extensions);
    const resources = createResources(gl, state, stats, extensions, parameters);
    const timer = createTimer(gl, extensions, stats);

    if (parameters.maxVertexTextureImageUnits < 8) {
        throw new Error('Need "MAX_VERTEX_TEXTURE_IMAGE_UNITS" >= 8');
    }

    const shaderPrecisionFormats = {
        vertex: getShaderPrecisionFormats(gl, 'vertex'),
        fragment: getShaderPrecisionFormats(gl, 'fragment'),
    };

    if (isDebugMode) {
        console.log({ parameters, shaderPrecisionFormats });
    }

    // optimize assuming flats first and last data are same or differences don't matter
    // extension is only available when `FIRST_VERTEX_CONVENTION` is more efficient
    const epv = extensions.provokingVertex;
    epv?.provokingVertex(epv.FIRST_VERTEX_CONVENTION);

    let isContextLost = false;
    const contextRestored = new Subject<now.Timestamp>();

    let pixelScale = props.pixelScale || 1;

    const xr = {
        session: undefined as XRSession | undefined,
        layer: undefined as XRWebGLLayer | undefined,
        changed: new Subject<void>(),
        clear: () => {
            xr.layer = undefined;
            xr.session = undefined;
            xr.changed.next();
        }
    };

    const renderTargets = new Set<RenderTarget>();

    return {
        gl,
        isWebGL2: isWebGL2(gl),
        get pixelRatio() {
            const dpr = (typeof window !== 'undefined') ? (window.devicePixelRatio || 1) : 1;
            return dpr * (pixelScale || 1);
        },

        extensions,
        state,
        stats,
        resources,
        timer,

        get maxTextureSize() { return parameters.maxTextureSize; },
        get max3dTextureSize() { return parameters.max3dTextureSize; },
        get maxRenderbufferSize() { return parameters.maxRenderbufferSize; },
        get maxDrawBuffers() { return parameters.maxDrawBuffers; },
        get maxTextureImageUnits() { return parameters.maxTextureImageUnits; },
        get shaderPrecisionFormats() { return shaderPrecisionFormats; },

        namedComputeRenderables: Object.create(null),
        namedFramebuffers: Object.create(null),
        namedTextures: Object.create(null),

        get isContextLost() {
            return isContextLost || gl.isContextLost();
        },
        contextRestored,
        setContextLost: () => {
            isContextLost = true;
            timer.clear();
        },
        handleContextRestored: (extraResets?: () => void) => {
            resetExtensions(gl, extensions);

            state.reset();
            state.currentMaterialId = -1;
            state.currentProgramId = -1;
            state.currentRenderItemId = -1;

            resources.reset();
            renderTargets.forEach(rt => rt.reset());
            extraResets?.();

            isContextLost = false;
            contextRestored.next(now());
        },

        xr: {
            get session() {
                return xr.session;
            },
            changed: xr.changed,
            set: async (session: XRSession | undefined, options?: { resolutionScale?: number }) => {
                if (xr.session === session) return;

                await xr.session?.end();
                if (session === undefined) return;

                try {
                    await gl.makeXRCompatible();
                    xr.session = session;
                    xr.layer = new XRWebGLLayer(xr.session, gl, {
                        antialias: true,
                        alpha: true,
                        depth: true,
                        framebufferScaleFactor: pixelScale * (options?.resolutionScale ?? 1),
                    });
                    await xr.session.updateRenderState({ baseLayer: xr.layer });

                    xr.session.addEventListener('end', xr.clear);
                    xr.changed.next();
                } catch (err) {
                    if (session) {
                        await session.end();
                    } else {
                        xr.layer = undefined;
                        xr.session = undefined;
                    }
                    throw err;
                }
            },
            end: async () => {
                return xr.session?.end();
            },
        },

        setPixelScale: (value: number) => {
            pixelScale = value;
        },

        createRenderTarget: (width: number, height: number, depth?: boolean, type?: 'uint8' | 'float32' | 'fp16', filter?: TextureFilter, format?: 'rgba' | 'alpha') => {
            const renderTarget = createRenderTarget(gl, resources, width, height, depth, type, filter, format);
            renderTargets.add(renderTarget);
            return {
                ...renderTarget,
                destroy: () => {
                    renderTarget.destroy();
                    renderTargets.delete(renderTarget);
                }
            };
        },
        createDrawTarget: () => {
            return {
                id: -1,
                texture: createNullTexture(gl),
                framebuffer: createNullFramebuffer(),
                depthRenderbuffer: null,

                getByteCount: () => 0,

                getWidth: () => getDrawingBufferSize(gl, xr.layer, xr.session?.interactionMode).width,
                getHeight: () => getDrawingBufferSize(gl, xr.layer, xr.session?.interactionMode).height,
                bind: () => {
                    bindDrawingBuffer(gl, xr.layer);
                },
                setSize: () => {},
                reset: () => {},
                destroy: () => {}
            };
        },
        bindDrawingBuffer: () => bindDrawingBuffer(gl, xr.layer),
        getDrawingBufferSize: () => getDrawingBufferSize(gl, xr.layer, xr.session?.interactionMode),
        readPixels: (x: number, y: number, width: number, height: number, buffer: Uint8Array | Float32Array | Int32Array) => {
            readPixels(gl, x, y, width, height, buffer);
        },
        waitForGpuCommandsComplete: () => waitForGpuCommandsComplete(gl),
        waitForGpuCommandsCompleteSync: () => waitForGpuCommandsCompleteSync(gl),
        getFenceSync: () => {
            return isWebGL2(gl) ? gl.fenceSync(gl.SYNC_GPU_COMMANDS_COMPLETE, 0) : null;
        },
        checkSyncStatus: (sync: WebGLSync) => {
            if (!isWebGL2(gl)) return true;

            if (gl.getSyncParameter(sync, gl.SYNC_STATUS) === gl.SIGNALED) {
                gl.deleteSync(sync);
                return true;
            } else {
                return false;
            }
        },
        deleteSync: (sync: WebGLSync) => {
            if (isWebGL2(gl)) gl.deleteSync(sync);
        },
        clear: (red: number, green: number, blue: number, alpha: number) => {
            const drs = getDrawingBufferSize(gl, xr.layer);
            bindDrawingBuffer(gl, xr.layer);
            state.enable(gl.SCISSOR_TEST);
            state.depthMask(true);
            state.colorMask(true, true, true, true);
            state.clearColor(red, green, blue, alpha);
            state.viewport(0, 0, drs.width, drs.height);
            state.scissor(0, 0, drs.width, drs.height);
            gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
        },

        checkError: (message?: string) => {
            checkError(gl, message);
        },
        checkFramebufferStatus: (message?: string) => {
            checkFramebufferStatus(gl, message);
        },

        destroy: (options?: Partial<{ doNotForceWebGLContextLoss: boolean }>) => {
            resources.destroy();
            unbindResources(gl);

            xr.session?.removeEventListener('end', xr.clear);
            xr.session?.end();

            contextRestored.complete();
            xr.changed.complete();

            // to aid GC
            if (!options?.doNotForceWebGLContextLoss) {
                gl.getExtension('WEBGL_lose_context')?.loseContext();
                gl.getExtension('STACKGL_destroy_context')?.destroy();
            }
        }
    };
}