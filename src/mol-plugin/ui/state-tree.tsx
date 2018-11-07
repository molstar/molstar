/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react';
import { PluginContext } from '../context';
import { PluginStateObject } from 'mol-plugin/state/base';
import { StateObject } from 'mol-state'
import { PluginCommands } from 'mol-plugin/command';

export class StateTree extends React.Component<{ plugin: PluginContext }, { }> {
    componentWillMount() {
        this.props.plugin.events.state.data.updated.subscribe(() => this.forceUpdate());
    }
    render() {
        // const n = this.props.plugin.state.data.tree.nodes.get(this.props.plugin.state.data.tree.rootRef)!;
        const n = this.props.plugin.state.data.tree.rootRef;
        return <div>
            <StateTreeNode plugin={this.props.plugin} nodeRef={n} key={n} />
            { /* n.children.map(c => <StateTreeNode plugin={this.props.plugin} nodeRef={c!} key={c} />) */}
        </div>;
    }
}

export class StateTreeNode extends React.Component<{ plugin: PluginContext, nodeRef: string }, { }> {
    render() {
        const n = this.props.plugin.state.data.tree.nodes.get(this.props.nodeRef)!;
        const obj = this.props.plugin.state.data.objects.get(this.props.nodeRef)!;
        if (!obj.obj) {
            return <div style={{ borderLeft: '1px solid black', paddingLeft: '5px' }}>
                {StateObject.StateType[obj.state]} {obj.errorText}
            </div>;
        }
        const props = obj.obj!.props as PluginStateObject.Props;
        return <div style={{ borderLeft: '1px solid black', paddingLeft: '5px' }}>
            <a href='#' onClick={e => {
                e.preventDefault();
                PluginCommands.Data.SetCurrentObject.dispatch(this.props.plugin, { ref: this.props.nodeRef });
            }}>{props.label}</a>
            {props.description ? <small>{props.description}</small> : void 0}
            {n.children.size === 0
                ? void 0
                : <div style={{ marginLeft: '10px' }}>{n.children.map(c => <StateTreeNode plugin={this.props.plugin} nodeRef={c!} key={c} />)}</div>
            }
        </div>;
    }
}