/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react';
import { PluginContext } from '../context';
import { StateTree } from './state-tree';
import { Viewport } from './viewport';
import { Controls, _test_CreateTransform, _test_UpdateTransform } from './controls';
import { Transformer } from 'mol-state';

// TODO: base object with subscribe helpers, separate behavior list etc

export class Plugin extends React.Component<{ plugin: PluginContext }, { }> {
    render() {
        return <div style={{ position: 'absolute', width: '100%', height: '100%', fontFamily: 'monospace' }}>
            <div style={{ position: 'absolute', width: '350px', height: '100%', overflowY: 'scroll' }}>
                <h3>Data</h3>
                <StateTree plugin={this.props.plugin} state={this.props.plugin.state.data} />
                <hr />
                <_test_CurrentObject plugin={this.props.plugin} />
                <h3>Behaviors</h3>
                <StateTree plugin={this.props.plugin} state={this.props.plugin.state.behavior} />
            </div>
            <div style={{ position: 'absolute', left: '350px', right: '250px', height: '100%' }}>
                <Viewport plugin={this.props.plugin} />
            </div>
            <div style={{ position: 'absolute', width: '250px', right: '0', height: '100%' }}>
                <Controls plugin={this.props.plugin} />
            </div>
        </div>;
    }
}

export class _test_CurrentObject extends React.Component<{ plugin: PluginContext }, { }> {
    componentDidMount() {
        // TODO: move to constructor?
        this.props.plugin.behaviors.state.data.currentObject.subscribe(() => this.forceUpdate());
    }
    render() {
        const current = this.props.plugin.behaviors.state.data.currentObject.value;
        const ref = current.ref;
        // const n = this.props.plugin.state.data.tree.nodes.get(ref)!;
        const obj = this.props.plugin.state.data.cells.get(ref)!;

        const type = obj && obj.obj ? obj.obj.type : void 0;

        const transforms = type
            ? Transformer.fromType(type)
            : []
        return <div>
            Current Ref: {ref}
            <hr />
            <h3>Update</h3>
            <_test_UpdateTransform key={`${ref} update`} plugin={this.props.plugin} state={current.state} nodeRef={ref} />
            <hr />
            <h3>Create</h3>
            {
                transforms.map((t, i) => <_test_CreateTransform key={`${t.id} ${ref} ${i}`} plugin={this.props.plugin} state={current.state} transformer={t} nodeRef={ref} />)
            }
        </div>;
    }
}