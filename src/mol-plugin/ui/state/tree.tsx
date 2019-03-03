/**
 * Copyright (c) 2018 - 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react';
import { PluginStateObject } from 'mol-plugin/state/objects';
import { State, StateObject, StateTransform } from 'mol-state'
import { PluginCommands } from 'mol-plugin/command';
import { PluginUIComponent } from '../base';
import { StateObjectActions } from './actions';

export class StateTree extends PluginUIComponent<{ state: State }, { showActions: boolean }> {
    state = { showActions: true };

    componentDidMount() {
        this.subscribe(this.plugin.events.state.cell.created, e => {
            if (e.cell.transform.parent === StateTransform.RootRef) this.forceUpdate();
        });

        this.subscribe(this.plugin.events.state.cell.removed, e => {
            if (e.parent === StateTransform.RootRef) this.forceUpdate();
        });
    }

    static getDerivedStateFromProps(props: { state: State }, state: { showActions: boolean }) {
        const n = props.state.tree.root.ref;
        const children = props.state.tree.children.get(n);
        const showActions = children.size === 0;
        if (state.showActions === showActions) return null;
        return { showActions };
    }

    render() {
        const ref = this.props.state.tree.root.ref;
        if (this.state.showActions) {
            return <StateObjectActions state={this.props.state} nodeRef={ref} hideHeader={true} />
        }
        return <StateTreeNode state={this.props.state} nodeRef={ref} depth={0} />;
    }
}

class StateTreeNode extends PluginUIComponent<{ nodeRef: string, state: State, depth: number }, { state: State, isCollapsed: boolean }> {
    is(e: State.ObjectEvent) {
        return e.ref === this.props.nodeRef && e.state === this.props.state;
    }

    get cellState() {
        return this.props.state.cellStates.get(this.props.nodeRef);
    }

    componentDidMount() {
        this.subscribe(this.plugin.events.state.cell.stateUpdated, e => {
            if (this.is(e) && e.state.transforms.has(this.props.nodeRef)) {
                this.setState({ isCollapsed: e.cellState.isCollapsed });
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

    state = {
        isCollapsed: this.props.state.cellStates.get(this.props.nodeRef).isCollapsed,
        state: this.props.state
    }

    static getDerivedStateFromProps(props: { nodeRef: string, state: State }, state: { state: State, isCollapsed: boolean }) {
        if (props.state === state.state) return null;
        return {
            isCollapsed: props.state.cellStates.get(props.nodeRef).isCollapsed,
            state: props.state
        };
    }

    render() {
        const cell = this.props.state.cells.get(this.props.nodeRef);
        if (!cell || cell.obj === StateObject.Null) return null;

        const cellState = this.cellState;
        const showLabel = cell.status !== 'ok' || !cell.transform.props || !cell.transform.props.isGhost;
        const children = this.props.state.tree.children.get(this.props.nodeRef);
        const newDepth = showLabel ? this.props.depth + 1 : this.props.depth;

        if (!showLabel) {
            if (children.size === 0) return null;
            return <div style={{ display: cellState.isCollapsed ? 'none' : 'block' }}>
                {children.map(c => <StateTreeNode state={this.props.state} nodeRef={c!} key={c} depth={newDepth} />)}
            </div>;
        }

        return <>
            <StateTreeNodeLabel nodeRef={this.props.nodeRef} state={this.props.state} depth={this.props.depth} />
            {children.size === 0
                ? void 0
                : <div style={{ display: cellState.isCollapsed ? 'none' : 'block' }}>
                    {children.map(c => <StateTreeNode state={this.props.state} nodeRef={c!} key={c} depth={newDepth} />)}
                </div>
            }
        </>;
    }
}

class StateTreeNodeLabel extends PluginUIComponent<
    { nodeRef: string, state: State, depth: number },
    { state: State, isCurrent: boolean, isCollapsed: boolean /*, updaterCollapsed: boolean */ }> {

    is(e: State.ObjectEvent) {
        return e.ref === this.props.nodeRef && e.state === this.props.state;
    }

    componentDidMount() {
        this.subscribe(this.plugin.events.state.cell.stateUpdated, e => {
            if (this.is(e)) this.forceUpdate();
        });

        this.subscribe(this.plugin.state.behavior.currentObject, e => {
            if (!this.is(e)) {
                if (this.state.isCurrent && e.state.transforms.has(this.props.nodeRef)) {
                    this.setState({ isCurrent: this.props.state.current === this.props.nodeRef });
                }
                return;
            }

            if (e.state.transforms.has(this.props.nodeRef)) {
                this.setState({
                    isCurrent: this.props.state.current === this.props.nodeRef,
                    isCollapsed: this.props.state.cellStates.get(this.props.nodeRef).isCollapsed
                });
            }
        });
    }

    state = {
        isCurrent: this.props.state.current === this.props.nodeRef,
        isCollapsed: this.props.state.cellStates.get(this.props.nodeRef).isCollapsed,
        state: this.props.state,
        // updaterCollapsed: true
    }

    static getDerivedStateFromProps(props: { nodeRef: string, state: State }, state: { state: State, isCurrent: boolean, isCollapsed: boolean }) {
        if (props.state === state.state) return null;
        return {
            isCurrent: props.state.current === props.nodeRef,
            isCollapsed: props.state.cellStates.get(props.nodeRef).isCollapsed,
            state: props.state,
            updaterCollapsed: true
        };
    }

    setCurrent = (e: React.MouseEvent<HTMLElement>) => {
        e.preventDefault();
        e.currentTarget.blur();
        PluginCommands.State.SetCurrentObject.dispatch(this.plugin, { state: this.props.state, ref: this.props.nodeRef });
    }

    remove = (e: React.MouseEvent<HTMLElement>) => {
        e.preventDefault();
        PluginCommands.State.RemoveObject.dispatch(this.plugin, { state: this.props.state, ref: this.props.nodeRef, removeParentGhosts: true });
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

    highlight = (e: React.MouseEvent<HTMLElement>) => {
        e.preventDefault();
        PluginCommands.State.Highlight.dispatch(this.plugin, { state: this.props.state, ref: this.props.nodeRef });
        e.currentTarget.blur();
    }

    clearHighlight = (e: React.MouseEvent<HTMLElement>) => {
        e.preventDefault();
        PluginCommands.State.ClearHighlight.dispatch(this.plugin, { state: this.props.state, ref: this.props.nodeRef });
        e.currentTarget.blur();
    }

    // private toggleUpdaterObs = new Subject();
    // toggleUpdater = (e: React.MouseEvent<HTMLAnchorElement>) => {
    //     e.preventDefault();
    //     e.currentTarget.blur();
    //     this.toggleUpdaterObs.next();
    // }

    render() {
        const n = this.props.state.transforms.get(this.props.nodeRef)!;
        const cell = this.props.state.cells.get(this.props.nodeRef);
        if (!cell) return null;

        const isCurrent = this.is(this.props.state.behaviors.currentObject.value);


        let label: any;
        if (cell.status === 'pending' || cell.status === 'processing') {
            const name = n.transformer.definition.display.name;
            label = <><b>[{cell.status}]</b> <span title={name}>{name}</span></>;
        } else if (cell.status !== 'ok' || !cell.obj) {
            const name = n.transformer.definition.display.name;
            const title = `${cell.errorText}`
            label = <><b>[{cell.status}]</b> <a title={title} href='#' onClick={this.setCurrent}>{name}</a>: <i>{cell.errorText}</i></>;
        } else {
            const obj = cell.obj as PluginStateObject.Any;
            const title = `${obj.label} ${obj.description ? obj.description : ''}`
            if (this.state.isCurrent) {
                label = <><a title={title} href='#'><b>{obj.label}</b> {obj.description ? <small>{obj.description}</small> : void 0}</a></>;
            } else {
                label = <><a title={title} href='#' onClick={this.setCurrent}>{obj.label} {obj.description ? <small>{obj.description}</small> : void 0}</a></>;
            }
        }

        const children = this.props.state.tree.children.get(this.props.nodeRef);
        const cellState = this.props.state.cellStates.get(this.props.nodeRef);

        const visibility = <button onClick={this.toggleVisible} className={`msp-btn msp-btn-link msp-tree-visibility${cellState.isHidden ? ' msp-tree-visibility-hidden' : ''}`}>
            <span className='msp-icon msp-icon-visual-visibility' />
        </button>;

        const style: React.HTMLAttributes<HTMLDivElement>['style'] = {
            marginLeft: /* this.state.isCurrent ? void 0 :*/ `${this.props.depth * 10}px`,
            // paddingLeft: !this.state.isCurrent ? void 0 : `${this.props.depth * 10}px`,
            borderLeft: /* isCurrent || */ this.props.depth === 0 ? 'none' : void 0
        }

        const row = <div className={`msp-tree-row${isCurrent ? ' msp-tree-row-current' : ''}`} onMouseEnter={this.highlight} onMouseLeave={this.clearHighlight} style={style}>
            {label}
            {children.size > 0 &&  <button onClick={this.toggleExpanded} className='msp-btn msp-btn-link msp-tree-toggle-exp-button'>
                <span className={`msp-icon msp-icon-${cellState.isCollapsed ? 'expand' : 'collapse'}`} />
            </button>}
            {!cell.transform.props.isLocked && <button onClick={this.remove} className='msp-btn msp-btn-link msp-tree-remove-button'>
                <span className='msp-icon msp-icon-remove' />
            </button>}{visibility}
        </div>;

        // if (this.state.isCurrent) {
        //     return <>
        //         {row}
        //         <StateTreeNodeTransform {...this.props} toggleCollapsed={this.toggleUpdaterObs} />
        //     </>
        // }

        return row;
    }
}

// class StateTreeNodeTransform extends PluginUIComponent<{ nodeRef: string, state: State, depth: number, toggleCollapsed?: Observable<any> }> {
//     componentDidMount() {
//         // this.subscribe(this.plugin.events.state.object.updated, ({ ref, state }) => {
//         //     if (this.props.nodeRef !== ref || this.props.state !== state) return;
//         //     this.forceUpdate();
//         // });
//     }

//     render() {
//         const ref = this.props.nodeRef;
//         const cell = this.props.state.cells.get(ref)!;
//         const parent: StateObjectCell | undefined = (cell.sourceRef && this.props.state.cells.get(cell.sourceRef)!) || void 0;

//         if (!parent || parent.status !== 'ok') return null;

//         const transform = cell.transform;
//         return <UpdateTransformContol state={this.props.state} transform={transform} initiallyCollapsed={true} toggleCollapsed={this.props.toggleCollapsed} />;
//     }
// }