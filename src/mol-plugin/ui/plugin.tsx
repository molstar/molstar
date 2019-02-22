/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react';
import { PluginContext } from '../context';
import { StateTree } from './state/tree';
import { Viewport, ViewportControls } from './viewport';
import { Controls, TrajectoryControls, LociLabelControl } from './controls';
import { PluginUIComponent, PluginReactContext } from './base';
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
import { AnimationControls } from './state/animation';

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
            <Layout />
        </PluginReactContext.Provider>;
    }
}

class Layout extends PluginUIComponent {
    componentDidMount() {
        this.subscribe(this.plugin.layout.events.updated, () => this.forceUpdate());
    }

    region(kind: 'left' | 'right' | 'bottom' | 'main', element: JSX.Element) {
        return <div className={`msp-layout-region msp-layout-${kind}`}>
            <div className='msp-layout-static'>
                {element}
            </div>
        </div>;
    }

    render() {
        const layout = this.plugin.layout.state;
        return <div className='msp-plugin'>
            <div className={`msp-plugin-content ${layout.isExpanded ? 'msp-layout-expanded' : 'msp-layout-standard msp-layout-standard-outside'}`}>
                <div className={layout.showControls ? 'msp-layout-hide-top' : 'msp-layout-hide-top msp-layout-hide-right msp-layout-hide-bottom msp-layout-hide-left'}>
                    {this.region('main', <ViewportWrapper />)}
                    {layout.showControls && this.region('left', <State />)}
                    {layout.showControls && this.region('right', <div className='msp-scrollable-container msp-right-controls'>
                        <CurrentObject />
                        <Controls />
                        <AnimationControls />
                        <CameraSnapshots />
                        <StateSnapshots />
                    </div>)}
                    {layout.showControls && this.region('bottom', <Log />)}
                </div>
            </div>
        </div>;
    }
}

export class ViewportWrapper extends PluginUIComponent {
    render() {
        return <>
            <Viewport />
            <div style={{ position: 'absolute', left: '10px', top: '10px', color: 'white' }}>
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

export class State extends PluginUIComponent {
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

export class Log extends PluginUIComponent<{}, { entries: List<LogEntry> }> {
    private wrapper = React.createRef<HTMLDivElement>();

    componentDidMount() {
        this.subscribe(this.plugin.events.log, () => this.setState({ entries: this.plugin.log.entries }));
    }

    componentDidUpdate() {
        this.scrollToBottom();
    }

    state = { entries: this.plugin.log.entries };

    private scrollToBottom() {
        const log = this.wrapper.current;
        if (log) log.scrollTop = log.scrollHeight - log.clientHeight - 1;
    }

    render() {
        // TODO: ability to show full log
        // showing more entries dramatically slows animations.
        const maxEntries = 10;
        const xs = this.state.entries, l = xs.size;
        const entries: JSX.Element[] = [];
        for (let i = Math.max(0, l - maxEntries), o = 0; i < l; i++) {
            const e = xs.get(i);
            entries.push(<li key={o++}>
                <div className={'msp-log-entry-badge msp-log-entry-' + e!.type} />
                <div className='msp-log-timestamp'>{formatTime(e!.timestamp)}</div>
                <div className='msp-log-entry'>{e!.message}</div>
            </li>);
        }
        return <div ref={this.wrapper} className='msp-log' style={{ position: 'absolute', top: '0', right: '0', bottom: '0', left: '0', overflowY: 'auto' }}>
            <ul className='msp-list-unstyled'>{entries}</ul>
        </div>;
    }
}

export class CurrentObject extends PluginUIComponent {
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

        const transform = cell.transform;
        const def = transform.transformer.definition;

        const actions =  current.state.actions.fromCell(cell, this.plugin);
        return <>
            <div className='msp-current-header'>
                {cell.obj ? cell.obj.label : (def.display && def.display.name) || def.name}
            </div>
            {(parent && parent.status === 'ok') && <UpdateTransformContol state={current.state} transform={transform} />}
            {cell.status === 'ok' && <>
                <div className='msp-section-header'>Actions</div>
                {actions.map((act, i) => <ApplyActionContol plugin={this.plugin} key={`${act.id}`} state={current.state} action={act} nodeRef={ref} />)}
                </>}
        </>;
    }
}