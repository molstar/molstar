/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react';
import { PluginContext } from '../context';
import { PluginStateObject } from 'mol-plugin/state/base';
import { StateObject } from 'mol-state'

export class Tree extends React.Component<{ plugin: PluginContext }, { }> {

    componentWillMount() {
        this.props.plugin.events.stateUpdated.subscribe(() => this.forceUpdate());
    }
    render() {
        const n = this.props.plugin.state.data.tree.nodes.get(this.props.plugin.state.data.tree.rootRef)!;
        return <div>
            {n.children.map(c => <TreeNode plugin={this.props.plugin} nodeRef={c!} key={c} />)}
        </div>;
    }
}

export class TreeNode extends React.Component<{ plugin: PluginContext, nodeRef: string }, { }> {
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
            {props.label} {props.description ? <small>{props.description}</small> : void 0}
            {n.children.size === 0
                ? void 0
                : <div style={{ marginLeft: '10px' }}>{n.children.map(c => <TreeNode plugin={this.props.plugin} nodeRef={c!} key={c} />)}</div>
            }
        </div>;
    }
}