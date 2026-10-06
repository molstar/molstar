/** Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info. */
import { Camera } from '../../mol-canvas3d/camera';
import { GPUTextureUsage } from './compat';
import { WebGPUContext } from './context';

/** Opaque illumination inputs. Transparent surfaces and volumes retain their normal draw path. */
export class WebGPUTracingInput {
    private shaded?: GPUTexture;
    private normal?: GPUTexture;
    private albedo?: GPUTexture;
    private depth?: GPUTexture;
    private backDepth?: GPUTexture;
    private readonly pipelines = new Map<string, GPURenderPipeline>();

    constructor(private readonly context: WebGPUContext, private readonly module: GPUShaderModule, private readonly layout: GPUPipelineLayout, private readonly sphereModule = module) { }

    get textures() {
        if (!this.shaded) throw new Error('Render illumination inputs before accessing their textures.');
        return { shaded: this.shaded, normal: this.normal!, albedo: this.albedo!, depth: this.depth!, backDepth: this.backDepth! };
    }

    pipeline(back: boolean, cull: boolean, flip: boolean, sphere = false) {
        const key = `${back}/${cull}/${flip}/${sphere}`;
        const module = sphere ? this.sphereModule : this.module;
        let pipeline = this.pipelines.get(key);
        if (!pipeline) {
            pipeline = this.context.device.createRenderPipeline({
                label: `molstar-tracing-input-${key}`, layout: this.layout,
                vertex: { module, entryPoint: 'vs' },
                fragment: { module, entryPoint: back ? 'tracingBackDepth' : 'tracing', targets: back ? [] : Array.from({ length: 3 }, () => ({ format: 'rgba16float' as const })) },
                primitive: { topology: 'triangle-list', cullMode: back || !cull ? 'none' : flip ? 'front' : 'back' },
                depthStencil: { format: 'depth32float', depthWriteEnabled: true, depthCompare: back ? 'greater' : 'less-equal' },
            });
            this.pipelines.set(key, pipeline);
        }
        return pipeline;
    }

    private setSize(width: number, height: number) {
        if (this.shaded?.width === width && this.shaded.height === height) return;
        this.destroyTextures();
        const usage = GPUTextureUsage.RENDER_ATTACHMENT | GPUTextureUsage.TEXTURE_BINDING | GPUTextureUsage.COPY_SRC;
        const size = [width, height];
        this.shaded = this.context.device.createTexture({ label: 'molstar-tracing-shaded', size, format: 'rgba16float', usage });
        this.normal = this.context.device.createTexture({ label: 'molstar-tracing-normal', size, format: 'rgba16float', usage });
        this.albedo = this.context.device.createTexture({ label: 'molstar-tracing-albedo-density', size, format: 'rgba16float', usage });
        this.depth = this.context.device.createTexture({ label: 'molstar-tracing-depth', size, format: 'depth32float', usage });
        this.backDepth = this.context.device.createTexture({ label: 'molstar-tracing-back-depth', size, format: 'depth32float', usage });
    }

    render(encoder: GPUCommandEncoder, camera: Camera, width: number, height: number, draw: (pass: GPURenderPassEncoder, back: boolean) => void, renderBackDepth = true) {
        this.setSize(width, height);
        const { shaded, normal, albedo, depth, backDepth } = this.textures;
        const viewport = camera.viewport;
        const setViewport = (pass: GPURenderPassEncoder) => pass.setViewport(viewport.x, height - viewport.y - viewport.height, viewport.width, viewport.height, 0, 1);
        const pass = encoder.beginRenderPass({ label: 'molstar-tracing-gbuffer',
            colorAttachments: [shaded, normal, albedo].map(texture => ({ view: texture.createView(), loadOp: 'clear' as const, storeOp: 'store' as const })),
            depthStencilAttachment: { view: depth.createView(), depthClearValue: 1, depthLoadOp: 'clear', depthStoreOp: 'store' },
        });
        setViewport(pass); draw(pass, false); pass.end();
        if (!renderBackDepth) return;
        const back = encoder.beginRenderPass({ label: 'molstar-tracing-thickness', colorAttachments: [],
            depthStencilAttachment: { view: backDepth.createView(), depthClearValue: 0, depthLoadOp: 'clear', depthStoreOp: 'store' },
        });
        setViewport(back); draw(back, true); back.end();
    }

    private destroyTextures() { for (const texture of [this.shaded, this.normal, this.albedo, this.depth, this.backDepth]) texture?.destroy(); }
    dispose() { this.destroyTextures(); this.pipelines.clear(); }
}
