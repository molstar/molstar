/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginReactContext, PluginUIComponent } from './base';
import { OrderedMap } from 'immutable';
import { TaskManager } from '../mol-plugin/util/task-manager';
import { Progress } from '../mol-task';
import { IconButton } from './controls/common';
import { CancelSvg } from './controls/icons';
import { useContext, useEffect, useState } from 'react';
import { useBehavior } from './hooks/use-behavior';

export function BackgroundTaskProgress() {
    const plugin = useContext(PluginReactContext);
    const [tracked, setTracked] = useState<OrderedMap<number, TaskManager.ProgressEvent>>(OrderedMap());

    useEffect(() => {
        const started = plugin.events.task.progress.subscribe(e => {
            const hideOverlay = !!plugin.spec.components?.hideTaskOverlay;
            if (e.level === 'background' && (hideOverlay || !e.useOverlay)) {
                setTracked(tracked => tracked.set(e.id, e));
            }
        });

        const finished = plugin.events.task.finished.subscribe(({ id }) => {
            setTracked(tracked => tracked.delete(id));
        });

        return () => {
            started.unsubscribe();
            finished.unsubscribe();
        };
    }, [plugin]);

    return <div className='msp-background-tasks'>
        {tracked.valueSeq().map(e => <ProgressEntry key={e!.id} event={e!} />)}
        <CanvasCommitState />
    </div>;
}

function CanvasCommitState() {
    const plugin = useContext(PluginReactContext);
    const queueSize = useBehavior(plugin.canvas3d?.commitQueueSize);

    if (!queueSize) return null;

    return <div className='msp-task-state'>
        <div>
            <div>
                Commiting renderables... {queueSize} remaining
            </div>
        </div>
    </div>;
}

class ProgressEntry extends PluginUIComponent<{ event: TaskManager.ProgressEvent }> {
    abort = () => {
        this.plugin.managers.task.requestAbort(this.props.event.progress.root.progress.taskId, 'User Request');
    };

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

export function OverlayTaskProgress() {
    const plugin = useContext(PluginReactContext);
    const [tracked, setTracked] = useState<OrderedMap<number, TaskManager.ProgressEvent>>(OrderedMap());

    useEffect(() => {
        const started = plugin.events.task.progress.subscribe(e => {
            if (!!e.useOverlay) {
                setTracked(tracked => tracked.set(e.id, e));
            }
        });

        const finished = plugin.events.task.finished.subscribe(({ id }) => {
            setTracked(tracked => tracked.delete(id));
        });

        return () => {
            started.unsubscribe();
            finished.unsubscribe();
        };
    }, [plugin]);

    if (tracked.size === 0) return null;

    return <div className='msp-overlay-tasks'>
        {tracked.valueSeq().map(e => <ProgressEntry key={e!.id} event={e!} />)}
    </div>;
}
