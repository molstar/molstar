/** Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info. */
import { TracingProps } from '../../mol-canvas3d/passes/tracing';

/** Matches the existing illumination controller's quality ramp and initial cheap frame. */
export class WebGPUTracingQuality {
    private key = '';
    private previousTime = 0;
    private rays = 1;
    private steps = 1;
    private refine = 0;

    update(props: TracingProps, iteration: number, time = performance.now()) {
        const key = JSON.stringify([props.rendersPerFrame, props.targetFps, props.steps, props.refineSteps]);
        const minRays = Math.max(1, Math.min(64, Math.round(props.rendersPerFrame[0])));
        const maxRays = Math.max(minRays, Math.min(64, Math.round(props.rendersPerFrame[1])));
        const maxSteps = Math.max(1, Math.min(1024, Math.round(props.steps))), minSteps = Math.max(1, Math.round(maxSteps / 2));
        const maxRefine = Math.max(0, Math.min(8, Math.round(props.refineSteps))), minRefine = Math.min(1, maxRefine);
        if (key !== this.key) {
            this.key = key; this.rays = minRays; this.steps = minSteps; this.refine = minRefine;
        } else if (iteration > 0) {
            const target = props.targetFps > 0 ? 1000 / props.targetFps : Infinity;
            const delta = Math.max(0, time - this.previousTime), missedFrames = Math.round(delta / target);
            const decrease = () => {
                this.rays--;
                if (this.rays < 1) this.refine--;
                if (this.refine < minRefine) this.steps--;
            };
            if (missedFrames >= 2) {
                // Once every setting is below its minimum, additional decreases cannot change the result.
                for (let i = 0; i < Math.min(missedFrames, maxRays + maxRefine + maxSteps + 1); i++) decrease();
            } else if (delta < target) {
                this.steps++;
                if (this.steps > maxSteps) this.refine++;
                if (this.refine > maxRefine) this.rays++;
            } else if (delta > target + 0.5) decrease();
        }
        this.previousTime = time;
        this.rays = Math.min(maxRays, Math.max(minRays, this.rays));
        this.steps = Math.min(maxSteps, Math.max(minSteps, this.steps));
        this.refine = Math.min(maxRefine, Math.max(minRefine, this.refine));
        return {
            rendersPerFrame: iteration === 0 ? Math.ceil(this.rays / 2) : this.rays,
            steps: iteration === 0 ? minSteps : this.steps,
            refineSteps: iteration === 0 ? minRefine : this.refine,
        };
    }
}
