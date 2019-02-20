/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react';
import { PluginUIComponent } from './base';
import { OrderedMap } from 'immutable';
import { TaskManager } from 'mol-plugin/util/task-manager';
import { filter } from 'rxjs/operators';
import { Progress } from 'mol-task';

export class BackgroundTaskProgress extends PluginUIComponent<{ }, { tracked: OrderedMap<number, TaskManager.ProgressEvent> }> {
    componentDidMount() {
        this.subscribe(this.plugin.events.task.progress.pipe(filter(e => e.level !== 'none')), e => {
            this.setState({ tracked: this.state.tracked.set(e.id, e) })
        });
        this.subscribe(this.plugin.events.task.finished, ({ id }) => {
            this.setState({ tracked: this.state.tracked.delete(id) })
        })
    }

    state = { tracked: OrderedMap<number, TaskManager.ProgressEvent>() };

    render() {
        return <div>
            {this.state.tracked.valueSeq().map(e => <ProgressEntry key={e!.id} event={e!} />)}
        </div>;
    }
}

class ProgressEntry extends PluginUIComponent<{ event: TaskManager.ProgressEvent }> {
    render() {
        const root = this.props.event.progress.root;
        const subtaskCount = countSubtasks(this.props.event.progress.root) - 1;
        const pr = root.progress.isIndeterminate
            ? void 0
            : <>[{root.progress.current}/{root.progress.max}]</>;
        const subtasks = subtaskCount > 0
            ? <>[{subtaskCount} subtask(s)]</>
            : void 0
        return <div>
            {root.progress.message} {pr} {subtasks}
        </div>;
    }
}

function countSubtasks(progress: Progress.Node) {
    if (progress.children.length === 0) return 1;
    let sum = 0;
    for (const c of progress.children) sum += countSubtasks(c);
    return sum;
}