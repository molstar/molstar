/**
 * Copyright (c) 2018 - 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react';
import { PluginUIComponent } from '../base';
import { ApplyActionControl } from './apply-action';
import { State } from '../../../mol-state';
import { Icon } from '../controls/common';

export class StateObjectActions extends PluginUIComponent<{ state: State, nodeRef: string, hideHeader?: boolean, initiallyColapsed?: boolean }> {
    get current() {
        return this.plugin.state.behavior.currentObject.value;
    }

    componentDidMount() {
        this.subscribe(this.plugin.state.behavior.currentObject, o => {
            this.forceUpdate();
        });

        this.subscribe(this.plugin.events.state.object.updated, ({ ref, state }) => {
            const current = this.current;
            if (current.ref !== ref || current.state !== state) return;
            this.forceUpdate();
        });
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
            {actions.map((act, i) => <ApplyActionControl plugin={this.plugin} key={`${act.id}`} state={state} action={act} nodeRef={ref} initiallyCollapsed={this.props.initiallyColapsed} />)}
        </div>;
    }
}