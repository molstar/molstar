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

export class StateTree extends PluginComponent<{ state: State }, { }> {
    componentDidMount() {
        // this.subscribe(this.props.state.events.changed, () => {
        //     this.forceUpdate()
        // });
    }

    render() {
        // const n = this.props.plugin.state.data.tree.nodes.get(this.props.plugin.state.data.tree.rootRef)!;
        const n = this.props.state.tree.root.ref;
        return <div>
            <StateTreeNode state={this.props.state} nodeRef={n} />
            {/* n.children.map(c => <StateTreeNode plugin={this.props.plugin} nodeRef={c!} key={c} />) */}
        </div>;
    }
}

class StateTreeNode extends PluginComponent<{ nodeRef: string, state: State }, { }> {
    is(e: State.ObjectEvent) {
        return e.ref === this.props.nodeRef && e.state === this.props.state;
    }

    get cellState() {
        return this.props.state.tree.cellStates.get(this.props.nodeRef);
    }

    componentDidMount() {
        let isCollapsed = this.cellState.isCollapsed;
        this.subscribe(this.plugin.events.state.cell.stateUpdated, e => {
            if (this.is(e) && isCollapsed !== e.cellState.isCollapsed) {
                isCollapsed = e.cellState.isCollapsed;
                this.forceUpdate();
            }
        });

        this.subscribe(this.plugin.events.state.cell.created, e => {
            if (this.props.state === e.state && this.props.nodeRef === e.cell.transform.parent) {
                this.forceUpdate();
            }
        });

        this.subscribe(this.plugin.events.state.cell.removed, e => {
            if (this.props.state === e.state && this.props.nodeRef === e.parent) {
                this.forceUpdate();
            }
        });
    }

    toggleExpanded = (e: React.MouseEvent<HTMLElement>) => {
        e.preventDefault();
        PluginCommands.State.ToggleExpanded.dispatch(this.plugin, { state: this.props.state, ref: this.props.nodeRef });
    }

    render() {
        const cellState = this.cellState;

        const expander = <>
            [<a href='#' onClick={this.toggleExpanded}>{cellState.isCollapsed ? '+' : '-'}</a>]
        </>;

        const children = this.props.state.tree.children.get(this.props.nodeRef);
        return <div>
            {children.size === 0 ? void 0 : expander} <StateTreeNodeLabel nodeRef={this.props.nodeRef} state={this.props.state} />
            {children.size === 0
                ? void 0
                : <div style={{ marginLeft: '7px', paddingLeft: '3px', borderLeft: '1px solid #999', display: cellState.isCollapsed ? 'none' : 'block' }}>
                    {children.map(c => <StateTreeNode state={this.props.state} nodeRef={c!} key={c} />)}
                </div>
            }
        </div>;
    }
}

class StateTreeNodeLabel extends PluginComponent<{ nodeRef: string, state: State }> {
    is(e: State.ObjectEvent) {
        return e.ref === this.props.nodeRef && e.state === this.props.state;
    }

    componentDidMount() {
        this.subscribe(this.plugin.events.state.cell.stateUpdated, e => {
            if (this.is(e)) this.forceUpdate();
        });

        let isCurrent = this.is(this.props.state.behaviors.currentObject.value);

        this.subscribe(this.plugin.state.behavior.currentObject, e => {
            if (this.is(e)) {
                if (!isCurrent) {
                    isCurrent = true;
                    this.forceUpdate();
                }
            } else if (isCurrent) {
                isCurrent = false;
                // have to check the node wasn't remove
                if (e.state.tree.transforms.has(this.props.nodeRef)) this.forceUpdate();
            }
        });
    }

    setCurrent = (e: React.MouseEvent<HTMLElement>) => {
        e.preventDefault();
        PluginCommands.State.SetCurrentObject.dispatch(this.plugin, { state: this.props.state, ref: this.props.nodeRef });
    }

    remove = (e: React.MouseEvent<HTMLElement>) => {
        e.preventDefault();
        PluginCommands.State.RemoveObject.dispatch(this.plugin, { state: this.props.state, ref: this.props.nodeRef });
    }

    toggleVisible = (e: React.MouseEvent<HTMLElement>) => {
        e.preventDefault();
        PluginCommands.State.ToggleVisibility.dispatch(this.plugin, { state: this.props.state, ref: this.props.nodeRef });
    }

    render() {
        const n = this.props.state.tree.transforms.get(this.props.nodeRef)!;
        const cell = this.props.state.cells.get(this.props.nodeRef)!;

        const isCurrent = this.is(this.props.state.behaviors.currentObject.value);

        const remove = <>[<a href='#' onClick={this.remove}>X</a>]</>

        let label: any;
        if (cell.status !== 'ok' || !cell.obj) {
            const name = (n.transformer.definition.display && n.transformer.definition.display.name) || n.transformer.definition.name;
            label = <><b>{cell.status}</b> <a href='#' onClick={this.setCurrent}>{name}</a>: <i>{cell.errorText}</i></>;
        } else {
            const obj = cell.obj as PluginStateObject.Any;
            label = <><a href='#' onClick={this.setCurrent}>{obj.label}</a> {obj.description ? <small>{obj.description}</small> : void 0}</>;
        }

        const cellState = this.props.state.tree.cellStates.get(this.props.nodeRef);
        const visibility = <>[<a href='#' onClick={this.toggleVisible}>{cellState.isHidden ? 'H' : 'V'}</a>]</>;

        return <>
            {remove}{visibility} {isCurrent ? <b>{label}</b> : label}
        </>;
    }
}