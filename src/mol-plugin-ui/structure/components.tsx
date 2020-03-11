/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react';
import { CollapsableControls, CollapsableState, PurePluginUIComponent } from '../base';
import { StructureHierarchyManager } from '../../mol-plugin-state/manager/structure/hierarchy';
import { StructureComponentRef, StructureRepresentationRef } from '../../mol-plugin-state/manager/structure/hierarchy-state';
import { State, StateAction } from '../../mol-state';
import { PluginCommands } from '../../mol-plugin/commands';
import { ExpandGroup, IconButton, ControlGroup, ToggleButton } from '../controls/common';
import { UpdateTransformControl } from '../state/update-transform';
import { ActionMenu } from '../controls/action-menu';
import { ApplyActionControl } from '../state/apply-action';
import { StateTransforms } from '../../mol-plugin-state/transforms';
import { Icon } from '../controls/icons';

interface StructureComponentControlState extends CollapsableState {
    isDisabled: boolean
}

const MeasurementFocusOptions = {
    minRadius: 6,
    extraRadius: 6,
    durationMs: 250,
}

export class StructureComponentControls extends CollapsableControls<{}, StructureComponentControlState> {
    protected defaultState(): StructureComponentControlState {
        return { header: 'Components', isCollapsed: false, isDisabled: false };
    }

    renderControls() {
        return <>
            <ComponentEditorControls />
            <ComponentListControls />
        </>;
    }
}

interface ComponentEditorControlsState {
    action?: 'preset' | 'modify' | 'options',
    isDisabled: boolean
}

class ComponentEditorControls extends PurePluginUIComponent<{}, ComponentEditorControlsState> {
    state: ComponentEditorControlsState = {
        isDisabled: false
    };

    private toggleAction(action: ComponentEditorControlsState['action']) {
        return () => this.setState({ action: this.state.action === action ? void 0 : action });
    }

    togglePreset = this.toggleAction('preset');
    toggleModify = this.toggleAction('modify');
    toggleOptions = this.toggleAction('options');
    hideAction = () => this.setState({ action: void 0 });

    modifyControls = {
        add: () => { },
    }

    modifySelect() {
        return <div className='msp-control-row msp-select-row'>
            <button><Icon name='plus' /> Add</button>
            <button><Icon name='flow-branch' /> Merge</button>
            <button><Icon name='minus' /> Sub</button>
            <button><Icon name='brush' /> Color</button>
        </div>;
    }

    get presetControls() {
        return <ActionMenu items={this.presetActions} onSelect={this.applyPreset} />
    }

    get presetActions() {
        const actions = [
            ActionMenu.Item('Clear', null),
        ];
        // TODO: filter by applicable
        for (const p of this.plugin.builders.structure.representation.providerList) {
            actions.push(ActionMenu.Item(p.display.name, p));
        }
        return actions;
    }

    applyPreset: ActionMenu.OnSelect = item => {
        this.hideAction();

        if (!item) return;
        const mng = this.plugin.managers.structure;

        // TODO: proper list
        const structures = mng.hierarchy.state.currentModels[0].structures;
        if (item.value === null) mng.component.clear(structures);
        else mng.component.applyPreset(this.plugin.managers.structure.hierarchy.state.currentModels[0].structures, item.value as any);
    }

    render() {
        return <>
            <div className='msp-control-row msp-select-row'>
                <ToggleButton icon='bookmarks' label='Preset' toggle={this.togglePreset} isSelected={this.state.action === 'preset'} disabled={this.state.isDisabled} />
                <ToggleButton icon='flow-cascade' label='Modify' toggle={this.toggleModify} isSelected={this.state.action === 'modify'} disabled={this.state.isDisabled} />
                <ToggleButton icon='cog' label='Options' toggle={this.toggleOptions} isSelected={this.state.action === 'options'} disabled={this.state.isDisabled} />
            </div>
            {this.state.action === 'preset' && this.presetControls}
        </>;
    }
}

class ComponentListControls extends PurePluginUIComponent {
    get currentModels() {
        return this.plugin.managers.structure.hierarchy.behaviors.currentModels;
    }

    componentDidMount() {
        this.subscribe(this.currentModels, () => this.forceUpdate());
    }

    render() {
        const components = StructureHierarchyManager.getCommonComponentPivots(this.currentModels.value)
        return components.map(c => <StructureComponentEntry key={c.cell.transform.ref} component={c} />)
    }
}

type StructureComponentEntryActions = 'add-repr' | 'remove' | 'none'

const createRepr = StateAction.fromTransformer(StateTransforms.Representation.StructureRepresentation3D);
class StructureComponentEntry extends PurePluginUIComponent<{ component: StructureComponentRef }, { action: StructureComponentEntryActions }> {
    state = { action: 'none' as StructureComponentEntryActions }

    get ref() {
        return this.props.component.cell.transform.ref;
    }

    componentDidMount() {
        this.subscribe(this.plugin.events.state.cell.stateUpdated, e => {
            if (State.ObjectEvent.isCell(e, this.props.component.cell)) this.forceUpdate();
        });
    }

    toggleVisible = (e: React.MouseEvent<HTMLElement>) => {
        e.preventDefault();
        PluginCommands.State.ToggleVisibility(this.plugin, { state: this.props.component.cell.parent, ref: this.ref });
        e.currentTarget.blur();
    }

    remove(ref: string) {
        return () => {
            this.setState({ action: 'none' });
            PluginCommands.State.RemoveObject(this.plugin, { state: this.props.component.cell.parent, ref, removeParentGhosts: true });
        }
    }

    
    get removeActions(): ActionMenu.Items {
        const ret = [
            ActionMenu.Item('Remove Selection', 'remove', this.remove(this.ref))
        ];
        for (const repr of this.props.component.representations) {
            ret.push(ActionMenu.Item(`Remove ${repr.cell.obj?.label}`, 'remove', this.remove(repr.cell.transform.ref)))
        }
        return ret;
    }
    
    selectRemoveAction: ActionMenu.OnSelect = item => {
        if (!item) return;
        (item?.value as any)();
    }
    
    toggleAddRepr = () => this.setState({ action: this.state.action === 'none' ? 'add-repr' : 'none' });
    toggleRemoveActions = () => this.setState({ action: this.state.action === 'none' ? 'remove' : 'none' });

    highlight = (e: React.MouseEvent<HTMLElement>) => {
        e.preventDefault();
        PluginCommands.State.Highlight(this.plugin, { state: this.props.component.cell.parent, ref: this.ref });
    }

    clearHighlight = (e: React.MouseEvent<HTMLElement>) => {
        e.preventDefault();
        PluginCommands.State.ClearHighlight(this.plugin, { state: this.props.component.cell.parent, ref: this.ref });
    }

    focus = () => {
        const sphere = this.props.component.cell.obj?.data.boundary.sphere;
        if (sphere) {
            const { extraRadius, minRadius, durationMs } = MeasurementFocusOptions;
            const radius = Math.max(sphere.radius + extraRadius, minRadius);
            PluginCommands.Camera.Focus(this.plugin, { center: sphere.center, radius, durationMs });
        }
    }

    render() {
        const component = this.props.component;
        const cell = component.cell;
        const label = cell.obj?.label;
        return <>
            <div className='msp-control-row'>
                <button className='msp-control-button-label' title={`${label}. Click to focus.`} onClick={this.focus} onMouseEnter={this.highlight} onMouseLeave={this.clearHighlight} style={{ textAlign: 'left' }}>
                    {label}
                </button>
                <div className='msp-select-row'>
                    <IconButton onClick={this.toggleVisible} icon='visual-visibility' toggleState={!cell.state.isHidden} title={`${cell.state.isHidden ? 'Show' : 'Hide'} component`} small />
                    <IconButton onClick={this.toggleRemoveActions} icon='remove' title='Remove' small toggleState={this.state.action === 'remove'} />
                    <IconButton onClick={this.toggleAddRepr} icon='plus' title='Add Representation' toggleState={this.state.action === 'add-repr'} />
                </div>
            </div>
            {this.state.action === 'remove' && <ActionMenu items={this.removeActions} onSelect={this.selectRemoveAction} />}
            <div className='msp-control-offset'>
                {this.state.action === 'add-repr' && 
                <ControlGroup header='Add Representation' initialExpanded={true} hideExpander={true} hideOffset={true} onHeaderClick={this.toggleAddRepr} topRightIcon='off'>
                    <ApplyActionControl plugin={this.plugin} state={cell.parent} action={createRepr} nodeRef={this.ref} hideHeader noMargin onApply={this.toggleAddRepr} applyLabel='Add' />
                </ControlGroup>}
                {component.representations.map(r => <StructureRepresentationEntry key={r.cell.transform.ref} representation={r} />)}
            </div>
        </>;
    }
}

class StructureRepresentationEntry extends PurePluginUIComponent<{ representation: StructureRepresentationRef }> {
    render() {
        const repr = this.props.representation.cell;
        return <ExpandGroup header={`${repr.obj?.label || ''} Representation`} noOffset>
            <UpdateTransformControl state={repr.parent} transform={repr.transform} customHeader='none' noMargin />
        </ExpandGroup>;
    }
}