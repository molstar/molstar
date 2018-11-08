/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react';
import { PluginContext } from '../context';
import { PluginStateObject } from 'mol-plugin/state/base';
import { State } from 'mol-state'
import { PluginCommands } from 'mol-plugin/command';

export class StateTree extends React.Component<{ plugin: PluginContext, state: State }, { }> {
    componentDidMount() {
        // TODO: move to constructor?
        this.props.state.context.events.updated.subscribe(() => this.forceUpdate());
    }
    render() {
        // const n = this.props.plugin.state.data.tree.nodes.get(this.props.plugin.state.data.tree.rootRef)!;
        const n = this.props.state.tree.rootRef;
        return <div>
            <StateTreeNode plugin={this.props.plugin} state={this.props.state} nodeRef={n} key={n} />
            { /* n.children.map(c => <StateTreeNode plugin={this.props.plugin} nodeRef={c!} key={c} />) */}
        </div>;
    }
}

export class StateTreeNode extends React.Component<{ plugin: PluginContext, nodeRef: string, state: State }, { }> {
    render() {
        const n = this.props.state.tree.nodes.get(this.props.nodeRef)!;
        const cell = this.props.state.cells.get(this.props.nodeRef)!;

        const remove = <>[<a href='#' onClick={e => {
            e.preventDefault();
            PluginCommands.Data.RemoveObject.dispatch(this.props.plugin, { ref: this.props.nodeRef });
        }}>X</a>]</>

        let label: any;
        if (cell.status !== 'ok' || !cell.obj) {
            const name = (n.value.transformer.definition.display && n.value.transformer.definition.display.name) || n.value.transformer.definition.name;
            label = <><b>{cell.status}</b> <a href='#' onClick={e => {
                e.preventDefault();
                PluginCommands.Data.SetCurrentObject.dispatch(this.props.plugin, { ref: this.props.nodeRef });
            }}>{name}</a>: <i>{cell.errorText}</i></>;
        } else {
            const obj = cell.obj as PluginStateObject.Any;
            const props = obj.props;
            const type = obj.type;
            label = <>[<span title={type.description}>{ type.shortName }</span>] <a href='#' onClick={e => {
                e.preventDefault();
                PluginCommands.Data.SetCurrentObject.dispatch(this.props.plugin, { ref: this.props.nodeRef });
            }}>{props.label}</a> {props.description ? <small>{props.description}</small> : void 0}</>;
        }

        return <div>
            {remove}{label}
            {n.children.size === 0
                ? void 0
                : <div style={{ marginLeft: '7px', paddingLeft: '3px', borderLeft: '1px solid #999' }}>{n.children.map(c => <StateTreeNode plugin={this.props.plugin} state={this.props.state} nodeRef={c!} key={c} />)}</div>
            }
        </div>;
    }
}