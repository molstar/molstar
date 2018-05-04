/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Adapted from LiteMol
 * Copyright (c) 2016 - now David Sehnal, licensed under Apache 2.0, See LICENSE file for more info.
 */

import * as React from 'react'
import { JobInfo, JobsController } from '../../controller/misc/jobs';
import { Button } from '../controls/common';
import { View } from '../view';

class JobState extends React.Component<{ info: JobInfo, isSmall?: boolean }, {}> {
    render() {
        const info = this.props.info;
        return <div className='molstar-task-state'>
            <div>
                { info.abort ? <Button onClick={() => info.abort!.call(null) } style='remove'
                    icon='abort' title='Abort' customClass='molstar-btn-icon'
                /> : void 0 }
                <div>
                    {info.name}: {info.message}
                </div>
            </div>
        </div>;
    }
}

export class Overlay extends View<JobsController, {}, {}> {
    render() {
        const state = this.controller.latestState;

        if (!state.jobs!.count()) return <div className='molstar-empty-control' />

        const jobs: any[] = [];
        state.jobs!.forEach((t, k) => jobs.push(<JobState key={k} info={t!} />));

        return <div className='molstar-overlay'>
            <div className='molstar-overlay-background' />
            <div className='molstar-overlay-content-wrap'>
                <div className='molstar-overlay-content'>
                    <div>
                        {jobs}
                    </div>
                </div>
            </div>
        </div>;
    }
}

export class BackgroundJobs extends View<JobsController, {}, {}> {
    render() {
        const state = this.controller.latestState;

        if (!state.jobs!.count()) return <div className='molstar-empty-control' />

        const jobs: any[] = [];
        state.jobs!.forEach((t, k) => jobs.push(<JobState key={k} info={t!} isSmall={true} />));

        return <div className='molstar-background-jobs'>
            {jobs}
        </div>;
    }
}