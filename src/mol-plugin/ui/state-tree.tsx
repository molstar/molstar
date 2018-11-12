/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react';
import { PluginStateObject } from 'mol-plugin/state/objects';
import { State } from 'mol-state'
import { PluginCommands } from 'mol-plugin/command';
import { PluginComponent } from './base';
import { merge } from 'rxjs';

export class StateTree extends PluginComponent<{ state: State }, { }> {
    componentDidMount() {
        this.subscribe(this.props.state.events.changed, () => this.forceUpdate());
    }

    render() {
        // const n = this.props.plugin.state.data.tree.nodes.get(this.props.plugin.state.data.tree.rootRef)!;
        const n = this.props.state.tree.root.ref;
        return <div>
            <StateTreeNode state={this.props.state} nodeRef={n} key={n} />
            { /* n.children.map(c => <StateTreeNode plugin={this.props.plugin} nodeRef={c!} key={c} />) */}
        </div>;
    }
}

export class StateTreeNode extends PluginComponent<{ nodeRef: string, state: State }, { }> {
    componentDidMount() {
        this.subscribe(merge(this.context.events.state.data.object.cellState, this.context.events.state.behavior.object.cellState), o => {
            if (o.ref === this.props.nodeRef && o.state === this.props.state) this.forceUpdate();
        });
    }

    render() {
        const n = this.props.state.tree.nodes.get(this.props.nodeRef)!;
        const cell = this.props.state.cells.get(this.props.nodeRef)!;

        const remove = <>[<a href='#' onClick={e => {
            e.preventDefault();
            PluginCommands.State.RemoveObject.dispatch(this.context, { state: this.props.state, ref: this.props.nodeRef });
        }}>X</a>]</>

        let label: any;
        if (cell.status !== 'ok' || !cell.obj) {
            const name = (n.transformer.definition.display && n.transformer.definition.display.name) || n.transformer.definition.name;
            label = <><b>{cell.status}</b> <a href='#' onClick={e => {
                e.preventDefault();
                PluginCommands.State.SetCurrentObject.dispatch(this.context, { state: this.props.state, ref: this.props.nodeRef });
            }}>{name}</a>: <i>{cell.errorText}</i></>;
        } else {
            const obj = cell.obj as PluginStateObject.Any;
            label = <><a href='#' onClick={e => {
                e.preventDefault();
                PluginCommands.State.SetCurrentObject.dispatch(this.context, { state: this.props.state, ref: this.props.nodeRef });
            }}>{obj.label}</a> {obj.description ? <small>{obj.description}</small> : void 0}</>;
        }

        const expander = <>
            [<a href='#' onClick={e => {
                e.preventDefault();
                PluginCommands.State.ToggleExpanded.dispatch(this.context, { state: this.props.state, ref: this.props.nodeRef });
            }}>{cell.transform.cellState.isCollapsed ? '+' : '-'}</a>]
        </>;

        const children = this.props.state.tree.children.get(this.props.nodeRef);
        return <div>
            {remove}{children.size === 0 ? void 0 : expander} {label}
            {cell.transform.cellState.isCollapsed || children.size === 0
                ? void 0
                : <div style={{ marginLeft: '7px', paddingLeft: '3px', borderLeft: '1px solid #999' }}>{children.map(c => <StateTreeNode state={this.props.state} nodeRef={c!} key={c} />)}</div>
            }
        </div>;
    }
}