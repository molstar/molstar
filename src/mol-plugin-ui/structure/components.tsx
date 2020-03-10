/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react';
import { CollapsableControls, CollapsableState, PurePluginUIComponent } from '../base';
import { StructureHierarchyManager } from '../../mol-plugin-state/manager/structure';
import { StructureComponentRef, StructureRepresentationRef } from '../../mol-plugin-state/manager/structure/hierarchy';
import { Icon } from '../controls/icons';
import { State, StateAction } from '../../mol-state';
import { PluginCommands } from '../../mol-plugin/commands';
import { ExpandGroup, IconButton } from '../controls/common';
import { UpdateTransformControl } from '../state/update-transform';
import { ActionMenu } from '../controls/action-menu';
import { ApplyActionControl } from '../state/apply-action';
import { StateTransforms } from '../../mol-plugin-state/transforms';

interface StructureComponentControlState extends CollapsableState {
    isDisabled: boolean
}

export class StructureComponentControls extends CollapsableControls<{}, StructureComponentControlState> {
    protected defaultState(): StructureComponentControlState {
        return { header: 'Structure Components', isCollapsed: false, isDisabled: false };
    }

    get currentModels() {
        return this.plugin.managers.structure.hierarchy.behaviors.currentModels;
    }

    componentDidMount() {
        this.subscribe(this.currentModels, () => this.forceUpdate());
    }

    renderControls() {
        const components = StructureHierarchyManager.getCommonComponentPivots(this.currentModels.value)
        return <>
            {components.map(c => <StructureComponentEntry key={c.cell.transform.ref} component={c} />)}
        </>;
    }
}

const createRepr = StateAction.fromTransformer(StateTransforms.Representation.StructureRepresentation3D);
class StructureComponentEntry extends PurePluginUIComponent<{ component: StructureComponentRef }, { showActions: boolean, showAddRepr: boolean }> {
    state = { showActions: false, showAddRepr: false }

    is(e: State.ObjectEvent) {
        return e.ref === this.ref && e.state === this.props.component.cell.parent;
    }

    get ref() {
        return this.props.component.cell.transform.ref;
    }

    componentDidMount() {
        this.subscribe(this.plugin.events.state.cell.stateUpdated, e => {
            if (this.is(e)) this.forceUpdate();
        });
    }

    toggleVisible = (e: React.MouseEvent<HTMLElement>) => {
        e.preventDefault();
        PluginCommands.State.ToggleVisibility(this.plugin, { state: this.props.component.cell.parent, ref: this.ref });
        e.currentTarget.blur();
    }

    remove(ref: string) {
        return () => {
            this.setState({ showActions: false });
            PluginCommands.State.RemoveObject(this.plugin, { state: this.props.component.cell.parent, ref, removeParentGhosts: true });
        }
    }

    
    get actions(): ActionMenu.Items {
        const ret = [
            ActionMenu.Item(`${this.state.showAddRepr ? 'Hide ' : ''}Add Representation`, 'plus', this.toggleAddRepr),
            ActionMenu.Item('Remove', 'remove', this.remove(this.ref))
        ];
        for (const repr of this.props.component.representations) {
            ret.push(ActionMenu.Item(`Remove ${repr.cell.obj?.label}`, 'remove', this.remove(repr.cell.transform.ref)))
        }
        return ret;
    }
    
    selectAction: ActionMenu.OnSelect = item => {
        if (!item) return;
        (item?.value as any)();
    }
    
    toggleAddRepr = () => this.setState({ showActions: false, showAddRepr: !this.state.showAddRepr });
    toggleActions = () => this.setState({ showActions: !this.state.showActions });

    render() {
        const component = this.props.component;
        const cell = component.cell;
        const label = cell.obj?.label;
        return <>
            <div className='msp-control-row'>
                <span title={label}>{label}</span>
                <div className='msp-select-row'>
                    <button onClick={this.toggleVisible}><Icon name='visual-visibility' style={{ fontSize: '80%' }} /> {cell.state.isHidden ? 'Show' : 'Hide'}</button>
                    <IconButton onClick={this.toggleActions} icon='menu' style={{ width: '64px' }} toggleState={this.state.showActions} title='Actions' />
                </div>
            </div>
            {this.state.showActions && <ActionMenu items={this.actions} onSelect={this.selectAction} />}
            <div className='msp-control-offset'>
                {this.state.showAddRepr && 
                    <ApplyActionControl plugin={this.plugin} state={cell.parent} action={createRepr} nodeRef={this.ref} hideHeader noMargin onApply={this.toggleAddRepr} applyLabel='Add' />}
                {component.representations.map(r => <StructureRepresentationEntry key={r.cell.transform.ref} representation={r} />)}
            </div>
        </>;
    }
}

class StructureRepresentationEntry extends PurePluginUIComponent<{ representation: StructureRepresentationRef }> {
    render() {
        const repr = this.props.representation.cell;
        return <ExpandGroup header={repr.obj?.label || ''} noOffset>
            <UpdateTransformControl state={repr.parent} transform={repr.transform} customHeader='none' noMargin />
        </ExpandGroup>;
    }
}