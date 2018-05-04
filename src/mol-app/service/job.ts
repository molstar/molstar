/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Adapted from LiteMol
 * Copyright (c) 2016 - now David Sehnal, licensed under Apache 2.0, See LICENSE file for more info.
 */

import { Context } from '../context/context'
import { JobEvents } from '../event/basic';
import { PerformanceMonitor } from 'mol-util/performance-monitor';
import { formatProgress } from 'mol-util';
import { Progress, Task, Run } from 'mol-task';

export class Job<T> {
    private info: Job.Info;
    get id() { return this.info.id; }
    get reportTime() { return this.info.reportTime; }

    run(context: Context) {
        return this.runWithContext(context).result;
    }

    runWithContext(context: Context): Job.Running<T> {
        return new Job.Running(context, this.task, this.info);
    }

    setReportTime(report: boolean) {
        this.info.reportTime = report;
        return this;
    }

    constructor(public name: string, public type: Job.Type, private task: Task<T>) {
        this.info = {
            id: serialJobId++,
            name,
            type,
            reportTime: false
        };
    }
}

let serialJobId = 0;
export namespace Job {
    export let __DEBUG_MODE__ = false;

    export type Type = 'Normal' | 'Background' | 'Silent';

    export interface Info {
        id: number,
        type: Type,
        name: string,
        reportTime: boolean
    }

    export class Running<T> {
        result: Promise<T>;

        private progressUpdated(progress: Progress) {
            JobEvents.StateUpdated.dispatch(this.context, {
                jobId: this.info.id,
                type: this.info.type,
                name: this.info.name,
                message: formatProgress(progress),
                abort: progress.requestAbort
            });
        }

        private resolved() {
            try {
                this.context.performance.end('job' + this.info.id);
                if (this.info.reportTime) {
                    let time = this.context.performance.time('job' + this.info.id);
                    if (this.info.type !== 'Silent') this.context.logger.info(`${this.info.name} finished in ${PerformanceMonitor.format(time)}.`)
                }
            } finally {
                JobEvents.Completed.dispatch(this.context, this.info.id);
            }
        }

        private rejected(err: any) {
            this.context.performance.end('job' + this.info.id);
            this.context.performance.formatTime('job' + this.info.id);

            if (__DEBUG_MODE__) {
                console.error(err);
            }

            try {
                if (this.info.type === 'Silent') {
                    if (err.warn)  this.context.logger.warning(`Warning (${this.info.name}): ${err.message}`);
                    else console.error(`Error (${this.info.name})`, err);
                } else {
                    if (err.warn) {
                        this.context.logger.warning(`Warning (${this.info.name}): ${err.message}`);
                    } else {
                        let e = '' + err;
                        if (e.indexOf('Aborted') >= 0) this.context.logger.info(`${this.info.name}: Aborted.`);
                        else this.context.logger.error(`Error (${this.info.name}): ${err}`);
                    }
                }
            } catch (e) {
                console.error(e);
            } finally {
                JobEvents.Completed.dispatch(this.context, this.info.id);
            }
        }

        private run() {
            JobEvents.Started.dispatch(this.context, this.info);
            this.context.performance.start('job' + this.info.id);

            this.result = Run(this.task, (p: Progress) => this.progressUpdated(p), 250)
            this.result.then(() => this.resolved()).catch(e => this.rejected(e));
        }

        constructor(private context: Context, private task: Task<T>, private info: Info) {
            this.run();
        }
    }

    export interface State {
        jobId: number,
        type: Type,
        name: string,
        message: string,
        abort?: () => void
    }

    export function create<T>(name: string, type: Type, task: Task<T>) {
        return new Job<T>(name, type, task);
    }
}