/**
 * Copyright (c) 2022-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { now } from '../../mol-util/now';
import { GLRenderingContext } from './compat';
import { WebGLStats } from './context';
import { WebGLExtensions } from './extensions';

function movingAverage(avg: number, sample: number, count: number) {
    avg -= avg / count;
    avg += sample / count;
    return avg;
}

class MovingAverage {
    private readonly avgs = new Map<string, number>();

    add(label: string, sample: number) {
        let avg = this.avgs.get(label) || sample;
        avg = movingAverage(avg, sample, this.count);
        this.avgs.set(label, avg);
        return avg;
    }

    get(label: string) {
        return this.avgs.get(label);
    }

    stats() {
        return Object.fromEntries(this.avgs.entries());
    }

    clear() {
        this.avgs.clear();
    }

    constructor(private count: number) { }
}

function clearStats(stats: WebGLStats) {
    stats.calls.drawInstanced = 0;
    stats.calls.drawInstancedBase = 0;
    stats.calls.multiDrawInstancedBase = 0;
    stats.calls.counts = 0;

    stats.culled.lod = 0;
    stats.culled.frustum = 0;
    stats.culled.occlusion = 0;
}

export type TimerResult = {
    readonly label: string
    readonly gpuElapsed: number
    readonly gpuAvg: number
    readonly cpuElapsed: number
    readonly cpuAvg: number
    readonly children: TimerResult[]
    readonly calls?: Calls
    readonly note?: string
}

function getQuery(extensions: WebGLExtensions) {
    return extensions.disjointTimerQuery ? extensions.disjointTimerQuery.createQuery() : null;
}

type WebGLTimerOptions = { captureStats?: boolean, note?: string }

export type WebGLTimer = {
    /** Check with GPU for finished timers. */
    resolve: () => TimerResult[]
    mark: (label: string, options?: WebGLTimerOptions) => void
    markEnd: (label: string) => void
    stats: () => { gpu: Record<string, number>, cpu: Record<string, number> }
    formatedStats: () => Record<string, string>
    clear: () => void
    destroy: () => void
}

type Calls = {
    drawInstanced: number,
    counts: number,
}

type Measure = {
    label: string,
    queries: WebGLQuery[],
    children: Measure[],
    root: boolean,
    cpu: { start: number, end: number },
    captureStats: boolean,
    timeElapsed?: number,
    calls?: Calls,
    note?: string,
};

type QueryResult = { timeElapsed?: number, refCount: number };

export function createTimer(gl: GLRenderingContext, extensions: WebGLExtensions, stats: WebGLStats, options?: { avgCount: number }): WebGLTimer {
    const dtq = extensions.disjointTimerQuery;
    const avgCount = options?.avgCount ?? 30;

    const queries = new Map<WebGLQuery, QueryResult>();
    const pending = new Map<string, Measure>();
    const stack: Measure[] = [];
    const gpuAvgs = new MovingAverage(avgCount);
    const cpuAvgs = new MovingAverage(avgCount);

    let measures: Measure[] = [];
    let current: WebGLQuery | null = null;
    let capturingStats = false;

    const clear = () => {
        pending.clear();
        stack.length = 0;
        gpuAvgs.clear();
        cpuAvgs.clear();

        measures = [];
        current = null;
        capturingStats = false;

        if (dtq) {
            queries.forEach((_, query) => {
                dtq.deleteQuery(query);
            });
        }
        queries.clear();
    };

    const add = () => {
        if (!dtq) return;

        const query = getQuery(extensions);
        if (!query) return;

        dtq.beginQuery(dtq.TIME_ELAPSED, query);
        pending.forEach((measure, _) => {
            measure.queries.push(query);
        });
        queries.set(query, { refCount: pending.size });
        current = query;
    };

    return {
        resolve: () => {
            const results: TimerResult[] = [];
            if (!dtq || !measures.length || capturingStats) return results;

            // console.log('resolve');
            queries.forEach((result, query) => {
                if (result.timeElapsed !== undefined) return;

                const available = dtq.getQueryParameter(query, dtq.QUERY_RESULT_AVAILABLE);
                const disjoint = gl.getParameter(dtq.GPU_DISJOINT);

                if (available && !disjoint) {
                    const timeElapsed = dtq.getQueryParameter(query, dtq.QUERY_RESULT) as number;
                    result.timeElapsed = timeElapsed;
                    // console.log('timeElapsed', result.timeElapsed);
                }

                if (available || disjoint) {
                    dtq.deleteQuery(query);
                }
            });

            const unresolved: Measure[] = [];
            for (const measure of measures) {
                if (measure.queries.every(q => queries.get(q)?.timeElapsed !== undefined)) {
                    let timeElapsed = 0;
                    for (const query of measure.queries) {
                        const result = queries.get(query)!;
                        timeElapsed += result.timeElapsed!;
                        result.refCount -= 1;
                    }
                    measure.timeElapsed = timeElapsed;
                    if (measure.root) {
                        const children: TimerResult[] = [];
                        const add = (measures: Measure[], children: TimerResult[]) => {
                            for (const measure of measures) {
                                const timeElapsed = measure.timeElapsed!;
                                const cpuElapsed = measure.cpu.end - measure.cpu.start;
                                const result: TimerResult = {
                                    label: measure.label,
                                    gpuElapsed: timeElapsed,
                                    gpuAvg: gpuAvgs.add(measure.label, timeElapsed),
                                    cpuElapsed,
                                    cpuAvg: cpuAvgs.add(measure.label, cpuElapsed),
                                    children: [],
                                    calls: measure.calls,
                                    note: measure.note,
                                };
                                children.push(result);
                                add(measure.children, result.children);
                            }
                        };
                        add(measure.children, children);
                        const cpuElapsed = measure.cpu.end - measure.cpu.start;
                        results.push({
                            label: measure.label,
                            gpuElapsed: timeElapsed,
                            gpuAvg: gpuAvgs.add(measure.label, timeElapsed),
                            cpuElapsed,
                            cpuAvg: cpuAvgs.add(measure.label, cpuElapsed),
                            children,
                            calls: measure.calls,
                            note: measure.note,
                        });
                    }
                } else {
                    unresolved.push(measure);
                }
            }
            measures = unresolved;

            queries.forEach((result, query) => {
                if (result.refCount === 0) {
                    queries.delete(query);
                }
            });

            return results;
        },
        mark: (label: string, options?: WebGLTimerOptions) => {
            if (!dtq) return;

            if (pending.has(label)) {
                throw new Error(`Timer mark for '${label}' already exists`);
            }

            const captureStats = options?.captureStats ?? false;

            if (current !== null) {
                dtq.endQuery(dtq.TIME_ELAPSED);
            }
            const measure: Measure = {
                label,
                queries: [],
                children: [],
                root: current === null,
                cpu: { start: now(), end: -1 },
                captureStats,
            };
            if (options?.note) measure.note = options.note;
            pending.set(label, measure);

            if (stack.length) {
                stack[stack.length - 1].children.push(measure);
            }
            stack.push(measure);

            if (captureStats) {
                if (capturingStats) {
                    throw new Error('Already capturing stats');
                }
                clearStats(stats);
                capturingStats = true;
            }

            add();
        },
        markEnd: (label: string) => {
            if (!dtq) return;

            const measure = pending.get(label);
            if (!measure) {
                throw new Error(`Timer mark for '${label}' does not exist`);
            }

            if (stack.pop()?.label !== label) {
                throw new Error(`Timer mark for '${label}' has pending nested mark`);
            }

            dtq.endQuery(dtq.TIME_ELAPSED);
            pending.delete(label);

            measure.cpu.end = now();
            if (measure.captureStats) {
                measure.calls = { ...stats.calls };
                capturingStats = false;
            }

            measures.push(measure);

            if (pending.size > 0) {
                add();
            } else {
                current = null;
            }
        },
        stats: () => {
            return {
                gpu: gpuAvgs.stats(),
                cpu: cpuAvgs.stats(),
            };
        },
        formatedStats: () => {
            const stats: Record<string, string> = {};
            const gpu = gpuAvgs.stats();
            const cpu = cpuAvgs.stats();
            for (const l of Object.keys(gpu)) {
                const g = `${(gpu[l] / 1000 / 1000).toFixed(2)}`;
                const c = `${cpu[l].toFixed(2)}`;
                stats[l] = `${g} ms | CPU: ${c} ms`;
            }
            return stats;
        },
        clear,
        destroy: () => {
            clear();
        }
    };
}

function formatTimerResult(result: TimerResult) {
    const gpu = `${(result.gpuElapsed / 1000 / 1000).toFixed(2)}`;
    const gpuAvg = `${(result.gpuAvg / 1000 / 1000).toFixed(2)}`;
    const cpu = `${result.cpuElapsed.toFixed(2)}`;
    const cpuAvg = `${result.cpuAvg.toFixed(2)}`;
    return `${result.label} ${gpu} ms (avg. ${gpuAvg} ms) | CPU: ${cpu} ms (avg. ${cpuAvg} ms)`;
}

export function printTimerResults(results: TimerResult[]) {
    results.map(r => {
        const f = formatTimerResult(r);
        if (r.children.length || r.calls || r.note) {
            console.groupCollapsed(f);
            if (r.calls) console.log(r.calls);
            if (r.note) console.log(r.note);
            printTimerResults(r.children);
            console.groupEnd();
        } else {
            console.log(f);
        }
    });
}
