/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react';
import { ModelUnitcellRef } from '../../mol-plugin-state/manager/structure/hierarchy-state';
import { PluginCommands } from '../../mol-plugin/commands';
import { State } from '../../mol-state';
import { PurePluginUIComponent } from '../base';
import { IconButton } from '../controls/common';
import { UpdateTransformControl } from '../state/update-transform';

export class UnitcellEntry extends PurePluginUIComponent<{ refs: ModelUnitcellRef[] }, { showOptions: boolean }> {
    state = { showOptions: false }

    componentDidMount() {
        this.subscribe(this.plugin.events.state.cell.stateUpdated, e => {
            if (State.ObjectEvent.isCell(e, this.pivot?.cell)) this.forceUpdate();
        });
    }

    get pivot() { return this.props.refs[0]; }

    toggleVisibility = (e: React.MouseEvent<HTMLElement>) => {
        e.preventDefault();
        this.plugin.managers.structure.hierarchy.toggleVisibility(this.props.refs);
        e.currentTarget.blur();
    }

    highlight = (e: React.MouseEvent<HTMLElement>) => {
        e.preventDefault();
        PluginCommands.Interactivity.Object.Highlight(this.plugin, {
            state: this.pivot.cell.parent,
            ref: this.props.refs.map(c => c.cell.transform.ref)
        });
    }

    clearHighlight = (e: React.MouseEvent<HTMLElement>) => {
        e.preventDefault();
        PluginCommands.Interactivity.ClearHighlights(this.plugin);
    }

    focus = (e: React.MouseEvent<HTMLElement>) => {
        e.preventDefault();

        const loci = [];
        for (const uc of this.props.refs) {
            if (uc.cell.state.isHidden) continue;

            const l = uc.cell.obj?.data.repr.getLoci()
            if (l) loci.push(l);
        }
        this.plugin.managers.camera.focusLoci(loci);
    }

    toggleOptions = () => this.setState({ showOptions: !this.state.showOptions })

    render() {
        const { refs } = this.props;
        if (refs.length === 0) return null;

        const pivot = refs[0];

        let label, description;
        if (refs.length === 1) {
            const { obj } = pivot.cell;
            if (!obj) return null;
            label = obj?.label;
            description = obj?.description;
        } else {
            label = 'Unitcells';
        }

        return <>
            <div className='msp-btn-row-group' style={{ marginTop: '6px' }}>
                <button className='msp-form-control msp-control-button-label' title={`${label}. Click to focus.`} onClick={this.focus} onMouseEnter={this.highlight} onMouseLeave={this.clearHighlight} style={{ textAlign: 'left' }}>
                    {label} <small>{description}</small>
                </button>
                <IconButton customClass='msp-form-control' onClick={this.toggleVisibility} icon='visual-visibility' toggleState={!pivot.cell.state.isHidden} title={`${pivot.cell.state.isHidden ? 'Show' : 'Hide'}`} small style={{ flex: '0 0 32px' }} />
                {refs.length === 1 && <IconButton customClass='msp-form-control' onClick={this.toggleOptions} icon='dot-3' title='Options' toggleState={this.state.showOptions} style={{ flex: '0 0 32px', padding: '0px' }} />}
            </div>
            {(refs.length === 1 && this.state.showOptions) && <>
                <div className='msp-control-offset'>
                    <UpdateTransformControl state={pivot.cell.parent} transform={pivot.cell.transform} customHeader='none' autoHideApply />
                </div>
            </>}
        </>;
    }
}