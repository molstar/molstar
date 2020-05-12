/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react';
import { PluginUIComponent } from './base';
import { OrderedMap } from 'immutable';
import { TaskManager } from '../mol-plugin/util/task-manager';
import { filter } from 'rxjs/operators';
import { Progress } from '../mol-task';
import { IconButton } from './controls/common';
import { CancelSvg } from './controls/icons';

export class BackgroundTaskProgress extends PluginUIComponent<{ }, { tracked: OrderedMap<number, TaskManager.ProgressEvent> }> {
    componentDidMount() {
        this.subscribe(this.plugin.events.task.progress.pipe(filter(e => e.level !== 'none')), e => {
            this.setState({ tracked: this.state.tracked.set(e.id, e) });
        });
        this.subscribe(this.plugin.events.task.finished, ({ id }) => {
            this.setState({ tracked: this.state.tracked.delete(id) });
        });
    }

    state = { tracked: OrderedMap<number, TaskManager.ProgressEvent>() };

    render() {
        return <div className='msp-background-tasks'>
            {this.state.tracked.valueSeq().map(e => <ProgressEntry key={e!.id} event={e!} />)}
        </div>;
    }
}

class ProgressEntry extends PluginUIComponent<{ event: TaskManager.ProgressEvent }> {
    abort = () => {
        this.plugin.requestTaskAbort(this.props.event.progress, 'User Request');
    }

    render() {
        const root = this.props.event.progress.root;
        const subtaskCount = countSubtasks(this.props.event.progress.root) - 1;
        const pr = root.progress.isIndeterminate
            ? void 0
            : <>[{root.progress.current}/{root.progress.max}]</>;
        const subtasks = subtaskCount > 0
            ? <>[{subtaskCount} subtask(s)]</>
            : void 0;

        return <div className='msp-task-state'>
            <div>
                {root.progress.canAbort && <IconButton svg={CancelSvg} onClick={this.abort} title='Abort' />}
                <div>
                    {root.progress.message} {pr} {subtasks}
                </div>
            </div>
        </div>;
    }
}

function countSubtasks(progress: Progress.Node) {
    if (progress.children.length === 0) return 1;
    let sum = 0;
    for (const c of progress.children) sum += countSubtasks(c);
    return sum;
}