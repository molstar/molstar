/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react';
import { CollapsableControls, CollapsableState } from '../base';
import { StateTransformer, StateTransform } from '../../mol-state';
import { ModelRef } from '../../mol-plugin-state/manager/structure/hierarchy-state';
import { StateTransforms } from '../../mol-plugin-state/transforms';
import { StructureBuilderTags } from '../../mol-plugin-state/builder/structure';
import { IconButton } from '../controls/common';

interface ObjectControlState extends CollapsableState {
    isBusy: boolean,
    showOptions: boolean,
}

export class ObjectControls extends CollapsableControls<{}, ObjectControlState> {
    protected defaultState(): ObjectControlState {
        return {
            header: 'Objects',
            isCollapsed: false,
            isBusy: false,
            showOptions: false
        };
    }

    componentDidMount() {
        this.subscribe(this.plugin.managers.structure.hierarchy.behaviors.selection, () => this.forceUpdate());
        this.subscribe(this.plugin.behaviors.state.isBusy, v => {
            this.setState({ isBusy: v })
        });
    }

    getUnitcell(model: ModelRef) {
        return model.genericRepresentations?.filter(r => {
            return r.cell.transform.transformer === StateTransforms.Representation.ModelUnitcell3D
        })[0]
    }

    toggleVisible = (e: React.MouseEvent<HTMLElement>) => {
        e.preventDefault();
        e.currentTarget.blur();

        for (const model of this.plugin.managers.structure.hierarchy.current.models) {
            const unitcell = this.getUnitcell(model)
            if (unitcell) {
                this.plugin.state.data.updateCellState(unitcell.cell.transform.ref, { isHidden: !unitcell.cell.state.isHidden });
            }
        }
    }

    isVisible() {
        for (const model of this.plugin.managers.structure.hierarchy.current.models) {
            const unitcell = this.getUnitcell(model)
            if (unitcell && !unitcell.cell.state.isHidden) return true
        }
        return false
    }

    async createUnitcell(model: ModelRef, params?: StateTransformer.Params<StateTransforms['Representation']['ModelUnitcell3D']>, initialState?: Partial<StateTransform.State>) {
        const state = this.plugin.state.data;
        const unitcell = state.build().to(model.cell)
            .apply(StateTransforms.Representation.ModelUnitcell3D, params, { tags: StructureBuilderTags.ModelGenericRepresentation, state: initialState });

        await this.plugin.updateDataState(unitcell, { revertOnError: true });
        return unitcell.selector;
    }

    ensureUnitcell = async () => {
        for (const model of this.plugin.managers.structure.hierarchy.current.models) {
            if (!this.getUnitcell(model)) await this.createUnitcell(model)
        }
    }

    highlight = (e: React.MouseEvent<HTMLElement>) => {
        e.preventDefault();
        // TODO
        // PluginCommands.Interactivity.Object.Highlight(this.plugin, { state: this.props.group[0].cell.parent, ref: this.props.group.map(c => c.cell.transform.ref) });
    }

    clearHighlight = (e: React.MouseEvent<HTMLElement>) => {
        e.preventDefault();
        // TODO
        // PluginCommands.Interactivity.ClearHighlights(this.plugin);
    }

    focus = () => {
        // TODO
        // const sphere = this.pivot.cell.obj?.data.boundary.sphere;
        // if (sphere) this.plugin.managers.camera.focusSphere(sphere);
    }

    toggleOptions = () => {
        // TODO
        this.setState({ showOptions: !this.state.showOptions })
    }

    renderControls() {
        const label = 'Unitcell'
        const isVisible = this.isVisible()
        return <>
            <div className='msp-control-row'>
                <button className='msp-control-button-label' title={`${label}. Click to focus.`} onClick={this.focus} onMouseEnter={this.highlight} onMouseLeave={this.clearHighlight} style={{ textAlign: 'left' }}>
                    {label}
                </button>
                <div className='msp-select-row'>
                    <IconButton onClick={this.toggleVisible} icon='visual-visibility' toggleState={isVisible} title={`${isVisible ? 'Hide' : 'Show'} component`} disabled={this.state.isBusy} style={{ flex: '0 0 40px' }} />
                    <IconButton onClick={this.toggleOptions} icon='cog' title='Options' toggleState={this.state.showOptions} disabled={this.state.isBusy} style={{ flex: '0 0 40px' }} />
                </div>
            </div>
            <div className='msp-control-group-header msp-flex-row' style={{ marginTop: '1px' }}>
                <button className='msp-btn msp-form-control msp-flex-item msp-no-overflow' onClick={this.ensureUnitcell}>
                    Unitcell
                </button>
            </div>
        </>;
    }
}