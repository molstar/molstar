/** Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info. */

import { Camera } from '../../mol-canvas3d/camera';
import { createWebGPUHandleHelper } from '../../mol-canvas3d/helper/handle-helper';
import { createWebGPUPointerHelper } from '../../mol-canvas3d/helper/pointer-helper';
import { PickData } from '../../mol-canvas3d/passes/pick';
import { Ray3D } from '../../mol-math/geometry/primitives/ray3d';
import { Mat4, Vec3 } from '../../mol-math/linear-algebra';
import { degToRad, spiral2d } from '../../mol-math/misc';
import { deepClone } from '../../mol-util/object';
import { GraphicsRenderObject } from '../render-object';
import { RendererProps } from '../renderer';
import { copyWebGPUCameraState } from './camera';
import { WebGPUContext } from './context';
import { WebGPURenderer } from './renderer';

/** Ray-aligned native picking owns its targets and leaves the displayed frame intact. */
export class WebGPURayPick {
    private renderer: Promise<WebGPURenderer> | undefined;
    private queue: Promise<unknown> = Promise.resolve();
    private disposed = false;

    constructor(private readonly context: WebGPUContext, private readonly helpers?: { handle: ReturnType<typeof createWebGPUHandleHelper>, pointer: ReturnType<typeof createWebGPUPointerHelper> }) { }

    pick(ray: Ray3D, source: Camera, objects: readonly GraphicsRenderObject[], props: RendererProps, padding: number, time = 0): Promise<PickData | undefined> {
        if (this.disposed || ![...ray.origin, ...ray.direction, source.scale].every(Number.isFinite) || source.scale <= 0 || Vec3.magnitude(ray.direction) === 0) return Promise.resolve(undefined);
        const radius = Math.max(0, Math.min(10, Math.round(padding))), size = radius * 2 + 1;
        const camera = new Camera(source.getSnapshot(), { x: 0, y: 0, width: size, height: size });
        copyWebGPUCameraState(camera, source);
        const direction = Vec3.normalize(Vec3(), ray.direction), target = Vec3.add(Vec3(), ray.origin, direction);
        const up = Vec3.clone(source.up);
        if (Math.abs(Vec3.dot(direction, up)) > 0.999) Vec3.copy(up, Math.abs(direction[1]) < 0.9 ? Vec3.unitY : Vec3.unitX);
        camera.setState({ mode: 'orthographic', position: Vec3.scale(Vec3(), ray.origin, 1 / source.scale), target: Vec3.scale(Vec3(), target, 1 / source.scale), up }, 0);
        camera.near = source.near; camera.far = source.far;
        camera.fogNear = source.fogNear; camera.fogFar = source.fogFar;
        // Match the narrow ray frustum used by the existing RayHelper.
        const half = Math.tan(degToRad(0.1) / 2) * Vec3.distance(source.position, source.target) * source.scale;
        camera.zoom = size / (2 * half);
        Mat4.ortho(camera.projection, -half, half, half, -half, camera.near, camera.far);
        Mat4.lookAt(camera.view, ray.origin, target, up);
        Mat4.mul(camera.projectionView, camera.projection, camera.view);
        Mat4.tryInvert(camera.inverseProjectionView, camera.projectionView);
        Mat4.copy(camera.viewEye, source.view);
        const rendererProps = deepClone(props);
        const run = this.queue.then(async () => {
            if (this.disposed) return;
            const renderer = await (this.renderer ??= WebGPURenderer.create(this.context));
            if (this.disposed) return;
            renderer.setWorldHelpers(this.helpers?.handle, this.helpers?.pointer);
            renderer.setTime(time);
            renderer.render(objects, camera, rendererProps, false, 1, { width: size, height: size, present: false });
            // Read the small ID/depth target once; inspect padding in spiral order.
            const { array: ids } = await renderer.readPickingPixels();
            const depths = new Float32Array(ids.buffer, ids.byteOffset, ids.length);
            for (const [dx, dy] of spiral2d(radius)) {
                const x = radius + dx, y = radius + dy, index = (y * size + x) * 4;
                if (!ids[index]) continue;
                return { id: { objectId: ids[index] - 1, instanceId: ids[index + 1], groupId: ids[index + 2] },
                    position: Vec3.scale(Vec3(), camera.unproject(Vec3(), Vec3.create(x + 0.5, size - y - 0.5, depths[index + 3])), 1 / camera.scale) };
            }
        });
        this.queue = run.catch(() => {});
        return run;
    }

    dispose() {
        this.disposed = true;
        if (this.renderer) void this.renderer.then(renderer => renderer.dispose(), () => {});
    }
}
