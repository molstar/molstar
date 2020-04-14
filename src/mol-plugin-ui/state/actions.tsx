/**
 * Copyright (c) 2018 - 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react';
import { State } from '../../mol-state';
import { PluginUIComponent } from '../base';
import { Icon } from '../controls/icons';
import { ApplyActionControl } from './apply-action';

export class StateObjectActions extends PluginUIComponent<{ state: State, nodeRef: string, hideHeader?: boolean, initiallyCollapsed?: boolean, alwaysExpandFirst?: boolean }> {
    get current() {
        return this.props.state.behaviors.currentObject.value;
    }

    componentDidMount() {
        // TODO: handle tree change: some state actions might become invalid
        // this.subscribe(this.props.state.events.changed, o => {
        //     this.setState(createStateObjectActionSelectState(this.props));
        // });

        this.subscribe(this.plugin.events.state.object.updated, ({ ref, state }) => {
            const current = this.current;
            if (current.ref !== ref || current.state !== state) return;
            this.forceUpdate();
        });

        this.subscribe(this.plugin.state.data.actions.events.added, () => this.forceUpdate());
        this.subscribe(this.plugin.state.data.actions.events.removed, () => this.forceUpdate());
    }

    render() {
        const { state, nodeRef: ref } = this.props;
        const cell = state.cells.get(ref)!;
        const actions = state.actions.fromCell(cell, this.plugin);
        if (actions.length === 0) return null;

        const def = cell.transform.transformer.definition;
        const display = cell.obj ? cell.obj.label : (def.display && def.display.name) || def.name;

        return <div className='msp-state-actions'>
            {!this.props.hideHeader && <div className='msp-section-header'><Icon name='code' /> {`Actions (${display})`}</div> }
            {actions.map((act, i) => <ApplyActionControl
                key={`${act.id}`} state={state} action={act} nodeRef={ref}
                initiallyCollapsed={i === 0 ? !this.props.alwaysExpandFirst && this.props.initiallyCollapsed : this.props.initiallyCollapsed} />)}
        </div>;
    }
}

// interface StateObjectActionSelectProps {
//     plugin: PluginContext,
//     state: State,
//     nodeRef: string
// }

// interface StateObjectActionSelectState {
//     state: State,
//     nodeRef: string,
//     version: string,
//     actions: readonly StateAction[],
//     currentActionIndex: number
// }

// function createStateObjectActionSelectState(props: StateObjectActionSelectProps): StateObjectActionSelectState {
//     const cell = props.state.cells.get(props.nodeRef)!;
//     const actions = [...props.state.actions.fromCell(cell, props.plugin)];
//     actions.sort((a, b) => a.definition.display.name < b.definition.display.name ? -1 : a.definition.display.name === b.definition.display.name ? 0 : 1);
//     return {
//         state: props.state,
//         nodeRef: props.nodeRef,
//         version: cell.transform.version,
//         actions,
//         currentActionIndex: -1
//     }
// }

// export class StateObjectActionSelect extends PluginUIComponent<StateObjectActionSelectProps, StateObjectActionSelectState> {
//     state = createStateObjectActionSelectState(this.props);

//     get current() {
//         return this.plugin.state.behavior.currentObject.value;
//     }

//     static getDerivedStateFromProps(props: StateObjectActionSelectProps, state: StateObjectActionSelectState) {
//         if (state.state !== props.state || state.nodeRef !== props.nodeRef) return createStateObjectActionSelectState(props);
//         const cell = props.state.cells.get(props.nodeRef)!;
//         if (cell.transform.version !== state.version) return createStateObjectActionSelectState(props);
//         return null;
//     }

//     componentDidMount() {
//         // TODO: handle tree change: some state actions might become invalid
//         // this.subscribe(this.props.state.events.changed, o => {
//         //     this.setState(createStateObjectActionSelectState(this.props));
//         // });

//         this.subscribe(this.plugin.events.state.object.updated, ({ ref, state }) => {
//             const current = this.current;
//             if (current.ref !== ref || current.state !== state) return;
//             this.setState(createStateObjectActionSelectState(this.props));
//         });
//     }

//     onChange =  (e: React.ChangeEvent<HTMLSelectElement>) => {
//         this.setState({ currentActionIndex: parseInt(e.target.value, 10) });
//     }

//     render() {
//         const actions = this.state.actions;
//         if (actions.length === 0) return null;

//         const current = this.state.currentActionIndex >= 0 && actions[this.state.currentActionIndex];
//         const title = current ? current.definition.display.description : 'Select Action';

//         return <>
//             <div className='msp-contol-row msp-action-select'>
//                 <select className='msp-form-control' title={title} value={this.state.currentActionIndex} onChange={this.onChange} style={{ fontWeight: 'bold' }}>
//                     <option key={-1} value={-1} style={{ color: '#999' }}>[ Select Action ]</option>
//                     {actions.map((a, i) => <option key={i} value={i}>{a.definition.display.name}</option>)}
//                 </select>
//                 <Icon name='flow-tree' />
//             </div>
//             {current && <ApplyActionControl key={current.id} plugin={this.plugin} state={this.props.state} action={current} nodeRef={this.props.nodeRef} hideHeader />}
//         </>;
//     }
// }