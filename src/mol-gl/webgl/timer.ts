/**
 * Copyright (c) 2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { GLRenderingContext } from './compat';
import { WebGLExtensions } from './extensions';

export type TimerResult = {
    readonly label: string
    readonly timeElapsed: number
    readonly children: TimerResult[]
}

function getQuery(extensions: WebGLExtensions) {
    return extensions.disjointTimerQuery ? extensions.disjointTimerQuery.createQuery() : null;
}

export type WebGLTimer = {
    /** Check with GPU for finished timers. */
    resolve: () => TimerResult[]
    mark: (label: string) => void
    markEnd: (label: string) => void
    clear: () => void
    destroy: () => void
}

type Measure = { label: string, queries: WebGLQuery[], children: Measure[], root: boolean, timeElapsed?: number };
type QueryResult = { timeElapsed?: number, refCount: number };

export function createTimer(gl: GLRenderingContext, extensions: WebGLExtensions): WebGLTimer {
    const dtq = extensions.disjointTimerQuery;

    const queries = new Map<WebGLQuery, QueryResult>();
    const pending = new Map<string, Measure>();
    const stack: Measure[] = [];

    let measures: Measure[] = [];
    let current: WebGLQuery | null = null;

    const clear = () => {
        if (!dtq) return;

        queries.forEach((_, query) => {
            dtq.deleteQuery(query);
        });
        pending.clear();
        measures = [];
        current = null;
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
            if (!dtq || !measures.length) return results;
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
                                const result: TimerResult = {
                                    label: measure.label,
                                    timeElapsed: measure.timeElapsed!,
                                    children: []
                                };
                                children.push(result);
                                add(measure.children, result.children);
                            }
                        };
                        add(measure.children, children);
                        results.push({ label: measure.label, timeElapsed, children });
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
        mark: (label: string) => {
            if (!dtq) return;

            if (pending.has(label)) {
                throw new Error(`Timer mark for '${label}' already exists`);
            }

            if (current !== null) {
                dtq.endQuery(dtq.TIME_ELAPSED);
            }
            const measure: Measure = { label, queries: [], children: [], root: current === null };
            pending.set(label, measure);

            if (stack.length) {
                stack[stack.length - 1].children.push(measure);
            }
            stack.push(measure);

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
            measures.push(measure);

            if (pending.size > 0) {
                add();
            } else {
                current = null;
            }
        },
        clear,
        destroy: () => {
            clear();
        }
    };
}

function formatTimerResult(result: TimerResult) {
    const timeElapsed = result.timeElapsed / 1000 / 1000;
    return `${result.label} ${timeElapsed.toFixed(2)}ms`;
}

export function printTimerResults(results: TimerResult[]) {
    return results.map(r => {
        const f = formatTimerResult(r);
        if (r.children.length) {
            console.groupCollapsed(f);
            printTimerResults(r.children);
            console.groupEnd();
        } else {
            console.log(f);
        }
    });
}
