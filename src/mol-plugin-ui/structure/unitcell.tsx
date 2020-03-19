/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react';
import { PurePluginUIComponent } from '../base';
import { StateTransformer, StateTransform, StateObjectCell } from '../../mol-state';
import { IconButton } from '../controls/common';
import { UpdateTransformControl } from '../state/update-transform';
import { PluginStateObject } from '../../mol-plugin-state/objects';
import { PluginCommands } from '../../mol-plugin/commands';

export type UnitcellCell = StateObjectCell<PluginStateObject.Shape.Representation3D, StateTransform<StateTransformer<PluginStateObject.Molecule.Model, PluginStateObject.Shape.Representation3D, any>>>

export class UnitcellEntry extends PurePluginUIComponent<{ cell: UnitcellCell }, { showOptions: boolean }> {
    state = { showOptions: false }

    componentDidMount() {
        this.subscribe(this.plugin.events.state.cell.stateUpdated, e => {
            this.forceUpdate();
        });
    }

    toggleVisibility = (e: React.MouseEvent<HTMLElement>) => {
        e.preventDefault();
        PluginCommands.State.ToggleVisibility(this.plugin, { state: this.props.cell.parent, ref: this.props.cell.transform.ref });
        e.currentTarget.blur();
    }

    highlight = (e: React.MouseEvent<HTMLElement>) => {
        e.preventDefault();
        if (!this.props.cell.state.isHidden) {
            PluginCommands.Interactivity.Object.Highlight(this.plugin, { state: this.props.cell.parent, ref: this.props.cell.transform.ref });
        }
    }

    clearHighlight = (e: React.MouseEvent<HTMLElement>) => {
        e.preventDefault();
        PluginCommands.Interactivity.ClearHighlights(this.plugin);
    }

    focus = (e: React.MouseEvent<HTMLElement>) => {
        e.preventDefault();
        if (!this.props.cell.state.isHidden) {
            const loci = this.props.cell.obj?.data.repr.getLoci()
            if (loci) this.plugin.managers.camera.focusLoci(loci);
        }
    }

    toggleOptions = () => this.setState({ showOptions: !this.state.showOptions })

    render() {
        const { cell } = this.props;
        const { obj } = cell;
        if (!obj) return null;

        return <>
            <div className='msp-btn-row-group' style={{ marginTop: '6px' }}>
                <button className='msp-form-control msp-control-button-label' title={`Unitcell. Click to focus.`} onClick={this.focus} onMouseEnter={this.highlight} onMouseLeave={this.clearHighlight} style={{ textAlign: 'left' }}>
                    Unitcell
                </button>
                <IconButton customClass='msp-form-control' onClick={this.toggleVisibility} icon='visual-visibility' toggleState={!cell.state.isHidden} title={`${cell.state.isHidden ? 'Show' : 'Hide'}`} small style={{ flex: '0 0 32px' }} />
                <IconButton customClass='msp-form-control' onClick={this.toggleOptions} icon='dot-3' title='Options' toggleState={this.state.showOptions} style={{ flex: '0 0 32px', padding: '0px' }} />
            </div>
            {this.state.showOptions && <>
                <div className='msp-control-offset'>
                    <UpdateTransformControl state={cell.parent} transform={cell.transform} customHeader='none' autoHideApply />
                </div>
            </>}
        </>;
    }
}