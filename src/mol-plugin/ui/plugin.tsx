/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react';
import { PluginContext } from '../context';
import { StateTree } from './state-tree';
import { Viewport } from './viewport';
import { Controls, _test_UpdateTransform, _test_ApplyAction, _test_TrajectoryControls } from './controls';
import { PluginComponent, PluginReactContext } from './base';

export class Plugin extends React.Component<{ plugin: PluginContext }, {}> {
    render() {
        return <PluginReactContext.Provider value={this.props.plugin}>
            <div style={{ position: 'absolute', width: '100%', height: '100%', fontFamily: 'monospace' }}>
                <div style={{ position: 'absolute', width: '350px', height: '100%', overflowY: 'scroll', padding: '10px' }}>
                    <StateTree state={this.props.plugin.state.data} />
                    <h3>Behaviors</h3>
                    <StateTree state={this.props.plugin.state.behavior} />
                </div>
                <div style={{ position: 'absolute', left: '350px', right: '300px', height: '100%' }}>
                    <Viewport />
                    <div style={{ position: 'absolute', left: '10px', top: '10px', height: '100%', color: 'white' }}>
                        <_test_TrajectoryControls />
                    </div>
                </div>
                <div style={{ position: 'absolute', width: '300px', right: '0', height: '100%', padding: '10px' }}>
                    <_test_CurrentObject />
                    <hr />
                    <Controls />
                </div>
            </div>
        </PluginReactContext.Provider>;
    }
}

export class _test_CurrentObject extends PluginComponent {
    init() {
        this.subscribe(this.context.behaviors.state.data.currentObject, () => this.forceUpdate());
    }

    render() {
        const current = this.context.behaviors.state.data.currentObject.value;

        const ref = current.ref;
        // const n = this.props.plugin.state.data.tree.nodes.get(ref)!;
        const obj = this.context.state.data.cells.get(ref)!;

        const type = obj && obj.obj ? obj.obj.type : void 0;

        const actions = type
            ? current.state.actions.fromType(type)
            : []
        return <div>
            <hr />
            <h3>Update {ref}</h3>
            <_test_UpdateTransform key={`${ref} update`} state={current.state} nodeRef={ref} />
            <hr />
            <h3>Create</h3>
            {
                actions.map((act, i) => <_test_ApplyAction key={`${act.id} ${ref} ${i}`}
                    state={current.state} action={act} nodeRef={ref} />)
            }
        </div>;
    }
}