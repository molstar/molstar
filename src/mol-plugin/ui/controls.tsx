/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react';
import { PluginCommands } from 'mol-plugin/command';
import { UpdateTrajectory } from 'mol-plugin/state/actions/basic';
import { PluginComponent } from './base';

export class Controls extends PluginComponent<{ }, { }> {
    render() {
        return <>

        </>;
    }
}

export class TrajectoryControls extends PluginComponent {
    render() {
        return <div>
            <b>Trajectory: </b>
            <button onClick={() => PluginCommands.State.ApplyAction.dispatch(this.plugin, {
                state: this.plugin.state.dataState,
                action: UpdateTrajectory.create({ action: 'advance', by: -1 })
            })}>&lt;&lt;</button>
            <button onClick={() => PluginCommands.State.ApplyAction.dispatch(this.plugin, {
                state: this.plugin.state.dataState,
                action: UpdateTrajectory.create({ action: 'reset' })
            })}>Reset</button>
            <button onClick={() => PluginCommands.State.ApplyAction.dispatch(this.plugin, {
                state: this.plugin.state.dataState,
                action: UpdateTrajectory.create({ action: 'advance', by: +1 })
            })}>&gt;&gt;</button><br />
        </div>
    }
}