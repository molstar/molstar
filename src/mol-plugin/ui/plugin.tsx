/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react';
import { PluginContext } from '../context';
import { StateTree } from './state-tree';
import { Viewport, ViewportControls } from './viewport';
import { Controls, _test_UpdateTransform, _test_ApplyAction, TrajectoryControls } from './controls';
import { PluginComponent, PluginReactContext } from './base';
import { merge } from 'rxjs';
import { State } from 'mol-state';
import { CameraSnapshots } from './camera';
import { StateSnapshots } from './state';
import { List } from 'immutable';
import { LogEntry } from 'mol-util/log-entry';
import { formatTime } from 'mol-util';
import { BackgroundTaskProgress } from './task';

export class Plugin extends React.Component<{ plugin: PluginContext }, {}> {
    render() {
        return <PluginReactContext.Provider value={this.props.plugin}>
            <div style={{ position: 'absolute', width: '100%', height: '100%', fontFamily: 'monospace' }}>
                <div style={{ position: 'absolute', width: '350px', height: '100%', overflowY: 'scroll', padding: '10px' }}>
                    <StateTree state={this.props.plugin.state.data} />
                    <h3>Behaviors</h3>
                    <StateTree state={this.props.plugin.state.behavior} />
                </div>
                <div style={{ position: 'absolute', left: '350px', right: '300px', top: '0', bottom: '100px' }}>
                    <Viewport />
                    <div style={{ position: 'absolute', left: '10px', top: '10px', height: '100%', color: 'white' }}>
                        <TrajectoryControls />
                    </div>
                    <ViewportControls />
                    <div style={{ position: 'absolute', left: '10px', bottom: '10px', color: 'white' }}>
                        <BackgroundTaskProgress />
                    </div>
                </div>
                <div style={{ position: 'absolute', width: '300px', right: '0', top: '0', padding: '10px', overflowY: 'scroll' }}>
                    <CurrentObject />
                    <hr />
                    <Controls />
                    <hr />
                    <CameraSnapshots />
                    <hr />
                    <StateSnapshots />
                </div>
                <div style={{ position: 'absolute', right: '300px', left: '350px', bottom: '0', height: '100px', overflow: 'hidden' }}>
                    <Log />
                </div>
            </div>
        </PluginReactContext.Provider>;
    }
}

export class Log extends PluginComponent<{}, { entries: List<LogEntry> }> {
    private wrapper = React.createRef<HTMLDivElement>();

    componentDidMount() {
        this.subscribe(this.plugin.events.log, e => this.setState({ entries: this.state.entries.push(e) }));
    }

    componentDidUpdate() {
        this.scrollToBottom();
    }

    state = { entries: List<LogEntry>() };

    private scrollToBottom() {
        const log = this.wrapper.current;
        if (log) log.scrollTop = log.scrollHeight - log.clientHeight - 1;
    }

    render() {
        return <div ref={this.wrapper} style={{ position: 'absolute', top: '0', right: '0', bottom: '0', left: '0', padding: '10px', overflowY: 'scroll' }}>
            <ul style={{ listStyle: 'none' }}>
                {this.state.entries.map((e, i) => <li key={i} style={{ borderBottom: '1px solid #999', padding: '3px' }}>
                    [{e!.type}] [{formatTime(e!.timestamp)}] {e!.message}
                </li>)}
            </ul>
        </div>;
    }
}

export class CurrentObject extends PluginComponent {
    componentDidMount() {
        let current: State.ObjectEvent | undefined = void 0;
        this.subscribe(merge(this.plugin.behaviors.state.data.currentObject, this.plugin.behaviors.state.behavior.currentObject), o => {
            current = o;
            this.forceUpdate()
        });

        this.subscribe(this.plugin.events.state.data.object.updated, ({ ref, state }) => {
            if (!current || current.ref !== ref && current.state !== state) return;
            this.forceUpdate();
        });
    }

    render() {
        const current = this.plugin.behaviors.state.data.currentObject.value;

        const ref = current.ref;
        // const n = this.props.plugin.state.data.tree.nodes.get(ref)!;
        const obj = this.plugin.state.data.cells.get(ref)!;

        const type = obj && obj.obj ? obj.obj.type : void 0;

        const actions = type
            ? current.state.actions.fromType(type)
            : []
        return <div>
            <hr />
            <h3>Update {obj.obj ? obj.obj.label : ref}</h3>
            <_test_UpdateTransform key={`${ref} update`} state={current.state} nodeRef={ref} />
            <hr />
            <h3>Create</h3>
            {
                actions.map((act, i) => <_test_ApplyAction key={`${act.id} ${ref} ${i}`}
                    state={current.state} action={act} nodeRef={ref} />)
            }
        </div>;
    }
}