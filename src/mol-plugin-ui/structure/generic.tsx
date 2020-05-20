/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react';
import { StructureHierarchyRef } from '../../mol-plugin-state/manager/structure/hierarchy-state';
import { PluginCommands } from '../../mol-plugin/commands';
import { State } from '../../mol-state';
import { PurePluginUIComponent } from '../base';
import { IconButton } from '../controls/common';
import { UpdateTransformControl } from '../state/update-transform';
import { VisibilityOffOutlinedSvg, VisibilityOutlinedSvg, MoreHorizSvg } from '../controls/icons';

export class GenericEntryListControls extends PurePluginUIComponent {
    get current() {
        return this.plugin.managers.structure.hierarchy.behaviors.selection;
    }

    componentDidMount() {
        this.subscribe(this.current, () => this.forceUpdate());
    }

    get unitcell() {
        const { selection } = this.plugin.managers.structure.hierarchy;
        if (selection.structures.length === 0) return null;

        const refs = [];
        for (const s of selection.structures) {
            const model = s.model;
            if (model?.unitcell && model.unitcell?.cell.obj) refs.push(model.unitcell);
        }
        if (refs.length === 0) return null;

        return <GenericEntry refs={refs} labelMultiple='Unit Cells' />;
    }

    get customControls(): JSX.Element[] | null {
        const controls: JSX.Element[] = [];
        this.plugin.genericRepresentationControls.forEach((provider, key) => {
            const [refs, labelMultiple] = provider(this.plugin.managers.structure.hierarchy.selection);
            if (refs.length > 0) {
                controls.push(<div key={key}>
                    <GenericEntry refs={refs} labelMultiple={labelMultiple} />
                </div>);
            }
        });
        return controls.length > 0 ? controls : null;
    }

    render() {
        return <>
            <div style={{ marginTop: '6px' }}>
                {this.unitcell}
                {this.customControls}
            </div>
        </>;
    }
}

export class GenericEntry<T extends StructureHierarchyRef> extends PurePluginUIComponent<{ refs: T[], labelMultiple?: string }, { showOptions: boolean }> {
    state = { showOptions: false }

    componentDidMount() {
        this.subscribe(this.plugin.state.events.cell.stateUpdated, e => {
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
        if (!this.pivot.cell.parent) return;
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

        let allHidden = true;
        for (const uc of this.props.refs) {
            if (!uc.cell.state.isHidden) {
                allHidden = false;
                break;
            }
        }

        if (allHidden) {
            this.plugin.managers.structure.hierarchy.toggleVisibility(this.props.refs, 'show');
        }

        const loci = [];
        for (const uc of this.props.refs) {
            if (uc.cell.state.isHidden) {
                continue;
            }

            const l = uc.cell.obj?.data.repr.getLoci();
            if (l) loci.push(l);
        }
        this.plugin.managers.camera.focusLoci(loci);
    }

    toggleOptions = () => this.setState({ showOptions: !this.state.showOptions })

    render() {
        const { refs, labelMultiple } = this.props;
        if (refs.length === 0) return null;

        const pivot = refs[0];

        let label, description;
        if (refs.length === 1) {
            const { obj } = pivot.cell;
            if (!obj) return null;
            label = obj?.label;
            description = obj?.description;
        } else {
            label = `${refs.length} ${labelMultiple || 'Objects'}`;
        }

        return <>
            <div className='msp-flex-row'>
                <button className='msp-form-control msp-control-button-label' title={`${label}. Click to focus.`} onClick={this.focus} onMouseEnter={this.highlight} onMouseLeave={this.clearHighlight} style={{ textAlign: 'left' }}>
                    {label} <small>{description}</small>
                </button>
                <IconButton svg={pivot.cell.state.isHidden ? VisibilityOffOutlinedSvg : VisibilityOutlinedSvg} toggleState={false} className='msp-form-control' onClick={this.toggleVisibility} title={`${pivot.cell.state.isHidden ? 'Show' : 'Hide'}`} small flex />
                {refs.length === 1 && <IconButton svg={MoreHorizSvg} className='msp-form-control' onClick={this.toggleOptions} title='Options' toggleState={this.state.showOptions} flex />}
            </div>
            {(refs.length === 1 && this.state.showOptions && pivot.cell.parent) && <>
                <div className='msp-control-offset'>
                    <UpdateTransformControl state={pivot.cell.parent} transform={pivot.cell.transform} customHeader='none' autoHideApply />
                </div>
            </>}
        </>;
    }
}