/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { List } from 'immutable';
import * as React from 'react';
import { PluginContext } from '../mol-plugin/context';
import { formatTime } from '../mol-util';
import { LogEntry } from '../mol-util/log-entry';
import { PluginReactContext, PluginUIComponent } from './base';
import { AnimationViewportControls, DefaultStructureTools, LociLabels, StateSnapshotViewportControls, TrajectoryViewportControls, SelectionViewportControls } from './controls';
import { LeftPanelControls } from './left-panel';
import { SequenceView } from './sequence';
import { BackgroundTaskProgress } from './task';
import { Toasts } from './toast';
import { Viewport, ViewportControls } from './viewport';
import { PluginCommands } from '../mol-plugin/commands';

export class Plugin extends React.Component<{ plugin: PluginContext }, {}> {
    region(kind: 'left' | 'right' | 'bottom' | 'main', element: JSX.Element) {
        return <div className={`msp-layout-region msp-layout-${kind}`}>
            <div className='msp-layout-static'>
                {element}
            </div>
        </div>;
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
        const controls = (this.plugin.spec.layout && this.plugin.spec.layout.controls) || {};

        const classList: string[] = [];
        if (controls.top === 'none' || !layout.showControls || layout.regionState.top === 'hidden') {
            classList.push('msp-layout-hide-top');
        }

        if (controls.left === 'none' || !layout.showControls || layout.regionState.left === 'hidden') {
            classList.push('msp-layout-hide-left');
        } else if (layout.regionState.left === 'collapsed') {
            classList.push('msp-layout-collapse-left');
        }

        if (controls.right === 'none' || !layout.showControls || layout.regionState.right === 'hidden') {
            classList.push('msp-layout-hide-right');
        }

        if (controls.bottom === 'none' || !layout.showControls || layout.regionState.bottom === 'hidden') {
            classList.push('msp-layout-hide-bottom');
        }

        return classList.join(' ');
    }

    get layoutClassName() {
        const layout = this.plugin.layout.state;

        const classList: string[] = ['msp-plugin-content'];
        if (layout.isExpanded) {
            classList.push('msp-layout-expanded');
        } else {
            classList.push('msp-layout-standard', `msp-layout-standard-${layout.controlsDisplay}`);
        }

        return classList.join(' ');
    }

    onDrop = (ev: React.DragEvent<HTMLDivElement>) => {
        ev.preventDefault();

        let file: File | undefined | null;
        if (ev.dataTransfer.items) {
            // Use DataTransferItemList interface to access the file(s)
            for (let i = 0; i < ev.dataTransfer.items.length; i++) {
                if (ev.dataTransfer.items[i].kind !== 'file') continue;
                file = ev.dataTransfer.items[i].getAsFile();
                break;
            }
        } else {
            file = ev.dataTransfer.files[0];
        }
        if (!file) return;

        const fn = file?.name.toLowerCase() || '';
        if (fn.endsWith('molx') || fn.endsWith('molj')) {
            PluginCommands.State.Snapshots.OpenFile(this.plugin, { file });
        }
    }

    onDragOver = (ev: React.DragEvent<HTMLDivElement>) => {
        ev.preventDefault();
    }

    render() {
        const layout = this.plugin.layout.state;
        const controls = this.plugin.spec.layout?.controls || {};
        const viewport = this.plugin.spec.components?.viewport?.view || DefaultViewport;

        return <div className='msp-plugin' onDrop={this.onDrop} onDragOver={this.onDragOver}>
            <div className={this.layoutClassName}>
                <div className={this.layoutVisibilityClassName}>
                    {this.region('main', viewport)}
                    {layout.showControls && controls.top !== 'none' && this.region('top', controls.top || SequenceView)}
                    {layout.showControls && controls.left !== 'none' && this.region('left', controls.left || LeftPanelControls)}
                    {layout.showControls && controls.right !== 'none' && this.region('right', controls.right || ControlsWrapper)}
                    {layout.showControls && controls.bottom !== 'none' && this.region('bottom', controls.bottom || Log)}
                </div>
            </div>
        </div>;
    }
}

export class ControlsWrapper extends PluginUIComponent {
    render() {
        const StructureTools = this.plugin.spec.components?.structureTools || DefaultStructureTools;
        return <div className='msp-scrollable-container'>
            {/* <CurrentObject /> */}
            <StructureTools />
        </div>;
    }
}

export class DefaultViewport extends PluginUIComponent {
    render() {
        const VPControls = this.plugin.spec.components?.viewport?.controls || ViewportControls;

        return <>
            <Viewport />
            <div className='msp-viewport-top-left-controls'>
                <AnimationViewportControls />
                <TrajectoryViewportControls />
                <StateSnapshotViewportControls />
            </div>
            <SelectionViewportControls />
            <VPControls />
            <BackgroundTaskProgress />
            <div className='msp-highlight-toast-wrapper'>
                <LociLabels />
                <Toasts />
            </div>
        </>;
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