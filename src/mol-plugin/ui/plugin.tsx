/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react';
import { PluginContext } from '../context';
import { StateTree } from './state-tree';
import { Viewport, ViewportControls } from './viewport';
import { Controls, TrajectoryControls, LociLabelControl } from './controls';
import { PluginComponent, PluginReactContext } from './base';
import { CameraSnapshots } from './camera';
import { StateSnapshots } from './state';
import { List } from 'immutable';
import { LogEntry } from 'mol-util/log-entry';
import { formatTime } from 'mol-util';
import { BackgroundTaskProgress } from './task';
import { ApplyActionContol } from './state/apply-action';
import { PluginState } from 'mol-plugin/state';
import { UpdateTransformContol } from './state/update-transform';
import { StateObjectCell } from 'mol-state';

export class Plugin extends React.Component<{ plugin: PluginContext }, {}> {

    region(kind: 'left' | 'right' | 'bottom' | 'main', element: JSX.Element) {
        return <div className={`msp-layout-region msp-layout-${kind}`}>
            <div className='msp-layout-static'>
                {element}
            </div>
        </div>
    }

    render() {
        return <PluginReactContext.Provider value={this.props.plugin}>
            <div className='msp-plugin'>
                <div className='msp-plugin-content msp-layout-expanded'>
                    <div className='msp-layout-hide-top'>
                        {this.region('main', <ViewportWrapper />)}
                        {this.region('left', <State />)}
                        {this.region('right', <div className='msp-scrollable-container msp-right-controls'>
                            <CurrentObject />
                            <Controls />
                            <CameraSnapshots />
                            <StateSnapshots />
                        </div>)}
                        {this.region('bottom', <Log />)}
                    </div>
                </div>
            </div>
        </PluginReactContext.Provider>;
    }
}

export class ViewportWrapper extends PluginComponent {
    render() {
        return <>
            <Viewport />
            <div style={{ position: 'absolute', left: '10px', top: '10px', height: '100%', color: 'white' }}>
                <TrajectoryControls />
            </div>
            <ViewportControls />
            <div style={{ position: 'absolute', left: '10px', bottom: '10px' }}>
                <BackgroundTaskProgress />
            </div>
            <div style={{ position: 'absolute', right: '10px', bottom: '10px' }}>
                <LociLabelControl />
            </div>
        </>;
    }
}

export class State extends PluginComponent {
    componentDidMount() {
        this.subscribe(this.plugin.state.behavior.kind, () => this.forceUpdate());
    }

    set(kind: PluginState.Kind) {
        // TODO: do command for this?
        this.plugin.state.setKind(kind);
    }

    render() {
        const kind = this.plugin.state.behavior.kind.value;
        return <div className='msp-scrollable-container'>
            <div className='msp-btn-row-group msp-data-beh'>
                <button className='msp-btn msp-btn-block msp-form-control' onClick={() => this.set('data')} style={{ fontWeight: kind === 'data' ? 'bold' : 'normal'}}>Data</button>
                <button className='msp-btn msp-btn-block msp-form-control' onClick={() => this.set('behavior')} style={{ fontWeight: kind === 'behavior' ? 'bold' : 'normal'}}>Behavior</button>
            </div>
            <StateTree state={kind === 'data' ? this.plugin.state.dataState : this.plugin.state.behaviorState} />
        </div>
    }
}

export class Log extends PluginComponent<{}, { entries: List<LogEntry> }> {
    private wrapper = React.createRef<HTMLDivElement>();

    componentDidMount() {
        // TODO: only show last 100 entries.
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
        return <div ref={this.wrapper} className='msp-log' style={{ position: 'absolute', top: '0', right: '0', bottom: '0', left: '0', overflowY: 'auto' }}>
            <ul className='msp-list-unstyled'>
                {this.state.entries.map((e, i) => <li key={i}>
                    <div className={'msp-log-entry-badge msp-log-entry-' + e!.type} />
                    <div className='msp-log-timestamp'>{formatTime(e!.timestamp)}</div>
                    <div className='msp-log-entry'>{e!.message}</div>
                </li>)}
            </ul>
        </div>;
    }
}

export class CurrentObject extends PluginComponent {
    get current() {
        return this.plugin.state.behavior.currentObject.value;
    }

    componentDidMount() {
        this.subscribe(this.plugin.state.behavior.currentObject, o => {
            this.forceUpdate();
        });

        this.subscribe(this.plugin.events.state.object.updated, ({ ref, state }) => {
            const current = this.current;
            if (current.ref !== ref || current.state !== state) return;
            this.forceUpdate();
        });
    }

    render() {
        const current = this.current;
        const ref = current.ref;
        const cell = current.state.cells.get(ref)!;
        const parent: StateObjectCell | undefined = (cell.sourceRef && current.state.cells.get(cell.sourceRef)!) || void 0;

        const type = cell && cell.obj ? cell.obj.type : void 0;
        const transform = cell.transform;
        const def = transform.transformer.definition;

        const actions = type ? current.state.actions.fromType(type) : [];
        return <>
            <div className='msp-section-header'>
                {cell.obj ? cell.obj.label : (def.display && def.display.name) || def.name}
            </div>
            { (parent && parent.status === 'ok') && <UpdateTransformContol state={current.state} transform={transform} /> }
            {cell.status === 'ok' &&
                actions.map((act, i) => <ApplyActionContol plugin={this.plugin} key={`${act.id}`} state={current.state} action={act} nodeRef={ref} />)
            }
        </>;
    }
}