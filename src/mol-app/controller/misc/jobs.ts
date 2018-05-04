/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Adapted from LiteMol
 * Copyright (c) 2016 - now David Sehnal, licensed under Apache 2.0, See LICENSE file for more info.
 */

import { Map } from 'immutable'
import { filter } from 'rxjs/operators';

import { Controller } from '../controller'
import { JobEvents } from '../../event/basic';
import { Context } from '../../context/context';
import { Job } from '../../service/job';


export interface JobInfo {
    name: string;
    message: string;
    abort?: () => void
}

export interface JobsState {
    jobs: Map<number, JobInfo>
}

export class JobsController extends Controller<JobsState> {
    private updated(state: Job.State) {
        let isWatched = state.type === this.type;
        let jobs = this.latestState.jobs!;

        if (!isWatched) {
            if (jobs.has(state.jobId)) {
                jobs = jobs.delete(state.jobId);
                this.setState({ jobs });
            }
            return;
        }

        jobs = jobs.set(state.jobId, {
            name: state.name,
            message: state.message,
            abort: state.abort
        });
        this.setState({ jobs });
    }

    private started(job: Job.Info) {
        this.setState({
            jobs: this.latestState.jobs!.set(job.id, { name: job.name, message: 'Running...' })
        });
    }

    private completed(taskId: number) {
        if (!this.latestState.jobs!.has(taskId)) return;

        this.setState({
            jobs: this.latestState.jobs!.delete(taskId)
        });
    }

    constructor(context: Context, private type: Job.Type) {
        super(context, {
            jobs: Map<number, JobInfo>()
        });

        JobEvents.StateUpdated.getStream(this.context)
            .subscribe(e => this.updated(e.data));

        JobEvents.Started.getStream(this.context).pipe(
            filter(e => e.data.type === type))
            .subscribe(e => this.started(e.data));

        JobEvents.Completed.getStream(this.context)
            .subscribe(e => this.completed(e.data));
    }
}