/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { List } from 'immutable';
import { PluginState } from '../../mol-plugin/state';
import { formatTime } from '../../mol-util';
import { LogEntry } from '../../mol-util/log-entry';
import * as React from 'react';
import { PluginContext } from '../context';
import { PluginReactContext, PluginUIComponent } from './base';
import { LociLabels, TrajectoryViewportControls, StateSnapshotViewportControls, AnimationViewportControls, StructureToolsWrapper } from './controls';
import { StateSnapshots } from './state';
import { StateObjectActions } from './state/actions';
import { StateTree } from './state/tree';
import { BackgroundTaskProgress } from './task';
import { Viewport, ViewportControls } from './viewport';
import { StateTransform } from '../../mol-state';
import { UpdateTransformControl } from './state/update-transform';
import { SequenceView } from './sequence';
import { Toasts } from './toast';
import { ImageControls } from './image';
import { Icon } from './controls/common';

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

export class PluginContextContainer extends React.Component<{ plugin: PluginContext }> {
    render() {
        return <PluginReactContext.Provider value={this.props.plugin}>
            <div className='msp-plugin'>
                {this.props.children}
            </div>
        </PluginReactContext.Provider>;
    }
}

type RegionKind = 'top' | 'left' | 'right' | 'bottom' | 'main'

class Layout extends PluginUIComponent {
    componentDidMount() {
        this.subscribe(this.plugin.layout.events.updated, () => this.forceUpdate());
    }

    region(kind: RegionKind, Element?: React.ComponentClass) {
        return <div className={`msp-layout-region msp-layout-${kind}`}>
            <div className='msp-layout-static'>
                {Element ? <Element /> : null}
            </div>
        </div>;
    }

    get layoutVisibilityClassName() {
        const layout = this.plugin.layout.state;
        const controls = (this.plugin.spec.layout && this.plugin.spec.layout.controls) || { };

        const classList: string[] = []
        if (controls.top === 'none' || !layout.showControls) {
            classList.push('msp-layout-hide-top')
        }
        if (controls.left === 'none' || !layout.showControls) {
            classList.push('msp-layout-hide-left')
        }
        if (controls.right === 'none' || !layout.showControls) {
            classList.push('msp-layout-hide-right')
        }
        if (controls.bottom === 'none' || !layout.showControls) {
            classList.push('msp-layout-hide-bottom')
        }

        return classList.join(' ')
    }

    get layoutClassName() {
        const layout = this.plugin.layout.state;

        const classList: string[] = ['msp-plugin-content']
        if (layout.isExpanded) {
            classList.push('msp-layout-expanded')
        } else {
            classList.push('msp-layout-standard', `msp-layout-standard-${layout.controlsDisplay}`)
        }

        return classList.join(' ')
    }

    render() {
        const layout = this.plugin.layout.state;
        const controls = (this.plugin.spec.layout && this.plugin.spec.layout.controls) || { };
        const viewport = (this.plugin.spec.layout && this.plugin.spec.layout.viewport) || ViewportWrapper;

        return <div className='msp-plugin'>
            <div className={this.layoutClassName}>
                <div className={this.layoutVisibilityClassName}>
                    {this.region('main', viewport)}
                    {layout.showControls && controls.top !== 'none' && this.region('top', controls.top || SequenceView)}
                    {layout.showControls && controls.left !== 'none' && this.region('left', controls.left || State)}
                    {layout.showControls && controls.right !== 'none' && this.region('right', controls.right || ControlsWrapper)}
                    {layout.showControls && controls.bottom !== 'none' && this.region('bottom', controls.bottom || Log)}
                </div>
            </div>
        </div>;
    }
}

export class ControlsWrapper extends PluginUIComponent {
    render() {
        return <div className='msp-scrollable-container msp-right-controls'>
            <CurrentObject />
            {/* <AnimationControlsWrapper /> */}
            {/* <CameraSnapshots /> */}
            <StructureToolsWrapper />
            <ImageControls />
            <StateSnapshots />
        </div>;
    }
}

export class ViewportWrapper extends PluginUIComponent {
    render() {
        return <>
            <Viewport />
            <div className='msp-viewport-top-left-controls'>
                <AnimationViewportControls />
                <TrajectoryViewportControls />
                <StateSnapshotViewportControls />
            </div>
            <ViewportControls />
            <BackgroundTaskProgress />
            <div className='msp-highlight-toast-wrapper'>
                <LociLabels />
                <Toasts />
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
                <button className='msp-btn msp-btn-block msp-form-control' onClick={() => this.set('data')} style={{ fontWeight: kind === 'data' ? 'bold' : 'normal' }}>
                    <span className='msp-icon msp-icon-database' /> Data
                </button>
                <button className='msp-btn msp-btn-block msp-form-control' onClick={() => this.set('behavior')} style={{ fontWeight: kind === 'behavior' ? 'bold' : 'normal' }}>
                    <span className='msp-icon msp-icon-tools' /> Behavior
                </button>
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
        const transform = cell.transform;

        let showActions = true;
        if (ref === StateTransform.RootRef) {
            const children = current.state.tree.children.get(ref);
            showActions = children.size !== 0;
        }

        if (!showActions) return null;

        return <>
            {(cell.status === 'ok' || cell.status === 'error') && <>
                <div className='msp-section-header' style={{ margin: '0 0 -1px 0' }}>
                    <Icon name='flow-cascade' />
                    {`${cell.obj?.label || transform.transformer.definition.display.name}`} <small>{transform.transformer.definition.display.name}</small>
                </div>
                <UpdateTransformControl state={current.state} transform={transform} customHeader='none' />
            </> }
            {cell.status === 'ok' && <StateObjectActions state={current.state} nodeRef={ref} initiallyCollapsed />}
        </>;
    }
}