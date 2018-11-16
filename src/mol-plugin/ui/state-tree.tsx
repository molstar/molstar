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
        return this.props.state.cellStates.get(this.props.nodeRef);
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

    render() {
        const cellState = this.cellState;

        const children = this.props.state.tree.children.get(this.props.nodeRef);
        return <div>
            <StateTreeNodeLabel nodeRef={this.props.nodeRef} state={this.props.state} />
            {children.size === 0
                ? void 0
                : <div className='msp-tree-children' style={{ display: cellState.isCollapsed ? 'none' : 'block' }}>
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
                // have to check the node wasn't removed
                if (e.state.transforms.has(this.props.nodeRef)) this.forceUpdate();
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
        e.currentTarget.blur();
    }

    toggleExpanded = (e: React.MouseEvent<HTMLElement>) => {
        e.preventDefault();
        PluginCommands.State.ToggleExpanded.dispatch(this.plugin, { state: this.props.state, ref: this.props.nodeRef });
        e.currentTarget.blur();
    }

    render() {
        const n = this.props.state.transforms.get(this.props.nodeRef)!;
        const cell = this.props.state.cells.get(this.props.nodeRef)!;

        const isCurrent = this.is(this.props.state.behaviors.currentObject.value);


        let label: any;
        if (cell.status !== 'ok' || !cell.obj) {
            const name = (n.transformer.definition.display && n.transformer.definition.display.name) || n.transformer.definition.name;
            const title = `${cell.errorText}`
            label = <><b>{cell.status}</b> <a title={title} href='#' onClick={this.setCurrent}>{name}</a>: <i>{cell.errorText}</i></>;
        } else {
            const obj = cell.obj as PluginStateObject.Any;
            const title = `${obj.label} ${obj.description ? obj.description : ''}`
            label = <><a title={title} href='#' onClick={this.setCurrent}>{obj.label}</a> {obj.description ? <small>{obj.description}</small> : void 0}</>;
        }

        const children = this.props.state.tree.children.get(this.props.nodeRef);
        const cellState = this.props.state.cellStates.get(this.props.nodeRef);

        const remove = <button onClick={this.remove} className='msp-btn msp-btn-link msp-tree-remove-button'>
            <span className='msp-icon msp-icon-remove' />
        </button>;

        const visibility = <button onClick={this.toggleVisible} className={`msp-btn msp-btn-link msp-tree-visibility${cellState.isHidden ? ' msp-tree-visibility-hidden' : ''}`}>
            <span className='msp-icon msp-icon-visual-visibility' />
        </button>;

        return <div className={`msp-tree-row${isCurrent ? ' msp-tree-row-current' : ''}`}>
            {isCurrent ? <b>{label}</b> : label}
            {children.size > 0 &&  <button onClick={this.toggleExpanded} className='msp-btn msp-btn-link msp-tree-toggle-exp-button'>
                <span className={`msp-icon msp-icon-${cellState.isCollapsed ? 'expand' : 'collapse'}`} />
            </button>}
            {remove}{visibility}
        </div>
    }
}