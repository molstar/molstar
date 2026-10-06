/** Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info. */
import { TracingParams } from '../../../mol-canvas3d/passes/tracing';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { WebGPUTracingQuality } from '../tracing-quality';

function props() { return { ...PD.getDefaultValues(TracingParams), steps: 8, refineSteps: 2, rendersPerFrame: [1, 4] as [number, number], targetFps: 30 }; }

describe('native illumination workload controller', () => {
    it('increases workload with spare frame time and reduces it after a stall', () => {
        const controller = new WebGPUTracingQuality(), p = props();
        const first = controller.update(p, 0, 0);
        let current = first;
        for (let i = 1; i <= 30; i++) current = controller.update(p, i, i * 10);
        expect(current).toEqual({ rendersPerFrame: 4, steps: 8, refineSteps: 2 });
        const slow = controller.update(p, 31, 10000);
        expect(slow.rendersPerFrame * slow.steps).toBeLessThan(current.rendersPerFrame * current.steps);
        expect(slow).toEqual({ rendersPerFrame: 1, steps: 4, refineSteps: 1 });
    });
    it('respects user quality bounds over varying frame times', () => {
        const controller = new WebGPUTracingQuality(), p = props();
        p.rendersPerFrame = [4, 8];
        expect(controller.update(p, 0, 0).rendersPerFrame).toBe(2); // Existing cheap first-frame behavior.
        let time = 0;
        for (let i = 1; i < 200; i++) {
            time += [1, 16, 34, 100, 1000][i % 5];
            const quality = controller.update(p, i, time);
            expect(quality.rendersPerFrame).toBeGreaterThanOrEqual(4); expect(quality.rendersPerFrame).toBeLessThanOrEqual(8);
            expect(quality.steps).toBeGreaterThanOrEqual(4); expect(quality.steps).toBeLessThanOrEqual(8);
            expect(quality.refineSteps).toBeGreaterThanOrEqual(1); expect(quality.refineSteps).toBeLessThanOrEqual(2);
        }
    });
    it('resets changed workload settings while preserving quality across scene restarts', () => {
        const controller = new WebGPUTracingQuality(), p = props();
        for (let i = 0; i < 30; i++) controller.update(p, i, i);
        expect(controller.update(p, 0, 31)).toEqual({ rendersPerFrame: 2, steps: 4, refineSteps: 1 });
        p.steps = 16; p.refineSteps = 0; p.rendersPerFrame = [2, 2];
        expect(controller.update(p, 0, 32)).toEqual({ rendersPerFrame: 1, steps: 8, refineSteps: 0 });
    });
    it('ramps toward full quality when no frame-rate target is set', () => {
        const controller = new WebGPUTracingQuality(), p = props(); p.targetFps = 0;
        let quality = controller.update(p, 0, 0);
        for (let i = 1; i < 30; i++) quality = controller.update(p, i, i * 1000);
        expect(quality).toEqual({ rendersPerFrame: 4, steps: 8, refineSteps: 2 });
    });
    it('keeps minimum and maximum legal settings finite after long pauses', () => {
        for (const [steps, refineSteps, rendersPerFrame] of [[1, 0, 1], [1024, 8, 64]]) {
            const controller = new WebGPUTracingQuality(), p = { ...props(), steps, refineSteps, rendersPerFrame: [rendersPerFrame, rendersPerFrame] as [number, number] };
            controller.update(p, 0, 0);
            const quality = controller.update(p, 1, 1e12);
            expect(Object.values(quality).every(Number.isFinite)).toBe(true);
            expect(quality.steps).toBeGreaterThanOrEqual(1); expect(quality.refineSteps).toBeGreaterThanOrEqual(0);
            expect(quality.rendersPerFrame).toBe(rendersPerFrame);
        }
    });
});
