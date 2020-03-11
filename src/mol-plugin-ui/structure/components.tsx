/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react';
import { CollapsableControls, CollapsableState, PurePluginUIComponent } from '../base';
import { StructureComponentRef, StructureRepresentationRef } from '../../mol-plugin-state/manager/structure/hierarchy-state';
import { State, StateAction } from '../../mol-state';
import { PluginCommands } from '../../mol-plugin/commands';
import { ExpandGroup, IconButton, ControlGroup, ToggleButton } from '../controls/common';
import { UpdateTransformControl } from '../state/update-transform';
import { ActionMenu } from '../controls/action-menu';
import { ApplyActionControl } from '../state/apply-action';
import { StateTransforms } from '../../mol-plugin-state/transforms';
import { Icon } from '../controls/icons';
import { StructureComponentManager } from '../../mol-plugin-state/manager/structure/component';
import { ParamDefinition } from '../../mol-util/param-definition';
import { ParameterControls } from '../controls/parameters';

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
        return { header: 'Representation', isCollapsed: false, isDisabled: false };
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

    get current() {
        return this.plugin.managers.structure.hierarchy.behaviors.current;
    }

    componentDidMount() {
        this.subscribe(this.current, () => this.setState({ action: void 0 }));
        this.subscribe(this.plugin.behaviors.state.isBusy, v => {
            this.setState({ isDisabled: v, action: void 0 })
        });
    }

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

        const structures = mng.hierarchy.state.currentStructures;
        if (item.value === null) mng.component.clear(structures);
        else mng.component.applyPreset(structures, item.value as any);
    }

    modifyComponentControls = <div className='msp-control-offset'><ModifyComponentControls onApply={this.hideAction} /></div>
    optionsControls = <div className='msp-control-offset'><ComponentOptionsControls /></div>

    render() {
        return <>
            <div className='msp-control-row msp-select-row'>
                <ToggleButton icon='bookmarks' label='Preset' toggle={this.togglePreset} isSelected={this.state.action === 'preset'} disabled={this.state.isDisabled} />
                <ToggleButton icon='flow-cascade' label='Modify' toggle={this.toggleModify} isSelected={this.state.action === 'modify'} disabled={this.state.isDisabled} />
                <ToggleButton icon='cog' label='Options' toggle={this.toggleOptions} isSelected={this.state.action === 'options'} disabled={this.state.isDisabled} />
            </div>
            {this.state.action === 'preset' && this.presetControls}
            {this.state.action === 'modify' && this.modifyComponentControls}
            {this.state.action === 'options' && this.optionsControls}
        </>;
    }
}

interface ModifyComponentControlsState {
    action?: StructureComponentManager.ActionType,
    actionParams?: ParamDefinition.Params,
    actionParamValues?: StructureComponentManager.ModifyAction
}

class ModifyComponentControls extends PurePluginUIComponent<{ onApply: () => void }, ModifyComponentControlsState> {
    state: ModifyComponentControlsState = { };

    private toggleAction(action: StructureComponentManager.ActionType) {
        return () => {
            if (this.state.action === action) {
                this.setState({ action: void 0, actionParams: void 0, actionParamValues: void 0 });
            } else {
                const actionParams = StructureComponentManager.getActionParams(this.plugin, action) as any;
                const actionParamValues = ParamDefinition.getDefaultValues(actionParams) as StructureComponentManager.ModifyAction;
                this.setState({ action, actionParams, actionParamValues });
            }
        }
    }

    toggleAdd = this.toggleAction('add');
    toggleMerge = this.toggleAction('merge');
    toggleSubtract = this.toggleAction('subtract');
    toggleColor = this.toggleAction('color');

    hideAction = () => this.setState({ action: void 0 });

    apply = () => {
        this.plugin.managers.structure.component.modify(this.state.actionParamValues!);
        this.props.onApply();
    }

    paramsChanged = (actionParamValues: any) => this.setState({ actionParamValues })
    get paramControls() {
        if (!this.state.action) return null;
        return <>
            <ParameterControls params={this.state.actionParams!} values={this.state.actionParamValues!} onChangeObject={this.paramsChanged} />
            <button className={`msp-btn msp-btn-block msp-btn-commit msp-btn-commit-on`} onClick={this.apply} style={{ marginTop: '1px' }}>
                <Icon name='ok' /> Apply
            </button>
        </>
    }

    render() {
        return <>
            <div className='msp-control-row msp-select-row'>
                <ToggleButton icon='plus' label='Add' toggle={this.toggleAdd} isSelected={this.state.action === 'add'} />
                <ToggleButton icon='flow-branch' label='Merge' toggle={this.toggleMerge} isSelected={this.state.action === 'merge'} />
                <ToggleButton icon='minus' label='Sub' toggle={this.toggleSubtract} isSelected={this.state.action === 'subtract'} />
                <ToggleButton icon='brush' label='Color' toggle={this.toggleColor} isSelected={this.state.action === 'color'} />
            </div>
            {this.paramControls}
        </>;
    }
}

class ComponentOptionsControls extends PurePluginUIComponent {
    componentDidMount() {
        this.subscribe(this.plugin.managers.structure.component.events.optionsUpdated, () => this.forceUpdate());
    }

    update = (options: StructureComponentManager.Options) => this.plugin.managers.structure.component.setOptions(options)

    render() {
        return <ParameterControls params={StructureComponentManager.OptionsParams} values={this.plugin.managers.structure.component.state.options} onChangeObject={this.update} />;
    }
}

class ComponentListControls extends PurePluginUIComponent {
    get current() {
        return this.plugin.managers.structure.hierarchy.behaviors.current;
    }

    componentDidMount() {
        this.subscribe(this.current, () => this.forceUpdate());
    }

    render() {
        const componentGroups = this.plugin.managers.structure.hierarchy.componentGroups;
        return <div style={{ marginTop: '8px' }}>
            {componentGroups.map(g => <StructureComponentGroup key={g[0].cell.transform.ref} group={g} />)}
        </div>;
    }
}

type StructureComponentEntryActions = 'add-repr' | 'remove' | 'none'

const createRepr = StateAction.fromTransformer(StateTransforms.Representation.StructureRepresentation3D);
class StructureComponentGroup extends PurePluginUIComponent<{ group: StructureComponentRef[] }, { action: StructureComponentEntryActions }> {
    state = { action: 'none' as StructureComponentEntryActions }

    get pivot() {
        return this.props.group[0];
    }

    componentDidMount() {
        this.subscribe(this.plugin.events.state.cell.stateUpdated, e => {
            if (State.ObjectEvent.isCell(e, this.pivot.cell)) this.forceUpdate();
        });
    }

    toggleVisible = (e: React.MouseEvent<HTMLElement>) => {
        e.preventDefault();
        // TODO: check visibility beforehand to set correct value if user made individual change
        for (const c of this.props.group) {
            PluginCommands.State.ToggleVisibility(this.plugin, { state: c.cell.parent, ref: c.cell.transform.ref });
        }
        e.currentTarget.blur();
    }

    get removeActions(): ActionMenu.Items {
        const ret = [
            ActionMenu.Item('Remove Selection', 'remove', () => this.plugin.managers.structure.hierarchy.remove(this.props.group))
        ];
        for (const repr of this.pivot.representations) {
            ret.push(ActionMenu.Item(`Remove ${repr.cell.obj?.label}`, 'remove', () => this.plugin.managers.structure.component.removeRepresentations(this.props.group, repr)))
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
        for (const c of this.props.group) {
            PluginCommands.State.Highlight(this.plugin, { state: c.cell.parent, ref: c.cell.transform.ref });
        }
    }

    clearHighlight = (e: React.MouseEvent<HTMLElement>) => {
        e.preventDefault();
        for (const c of this.props.group) {
            PluginCommands.State.ClearHighlight(this.plugin, { state: c.cell.parent, ref: c.cell.transform.ref });
        }
    }

    focus = () => {
        const sphere = this.pivot.cell.obj?.data.boundary.sphere;
        if (sphere) {
            const { extraRadius, minRadius, durationMs } = MeasurementFocusOptions;
            const radius = Math.max(sphere.radius + extraRadius, minRadius);
            PluginCommands.Camera.Focus(this.plugin, { center: sphere.center, radius, durationMs });
        }
    }

    render() {
        const component = this.pivot;
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
                    <ApplyActionControl plugin={this.plugin} state={cell.parent} action={createRepr} nodeRef={component.cell.transform.ref} hideHeader noMargin onApply={this.toggleAddRepr} applyLabel='Add' />
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