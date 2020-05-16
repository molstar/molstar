/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react';
import { StructureSelectionQueries, StructureSelectionQuery, getNonStandardResidueQueries, getElementQueries, getPolymerAndBranchedEntityQueries } from '../../mol-plugin-state/helpers/structure-selection-query';
import { InteractivityManager } from '../../mol-plugin-state/manager/interactivity';
import { StructureComponentManager } from '../../mol-plugin-state/manager/structure/component';
import { StructureRef, StructureComponentRef } from '../../mol-plugin-state/manager/structure/hierarchy-state';
import { StructureSelectionModifier } from '../../mol-plugin-state/manager/structure/selection';
import { memoizeLatest } from '../../mol-util/memoize';
import { ParamDefinition } from '../../mol-util/param-definition';
import { stripTags } from '../../mol-util/string';
import { PluginUIComponent, PurePluginUIComponent } from '../base';
import { ActionMenu } from '../controls/action-menu';
import { Button, ControlGroup, IconButton, ToggleButton } from '../controls/common';
import { ParameterControls, ParamOnChange, PureSelectControl } from '../controls/parameters';
import { UnionSvg, SubtractSvg, IntersectSvg, SetSvg, CubeOutlineSvg, Icon, SelectionModeSvg, RemoveSvg, RestoreSvg, HelpOutlineSvg, CancelOutlinedSvg, BrushSvg, CloseSvg } from '../controls/icons';
import { AddComponentControls } from './components';
import Structure from '../../mol-model/structure/structure/structure';
import { ViewportHelpContent, HelpGroup, HelpText } from '../viewport/help';


export class ToggleSelectionModeButton extends PurePluginUIComponent<{ inline?: boolean }> {
    componentDidMount() {
        this.subscribe(this.plugin.events.canvas3d.settingsUpdated, () => this.forceUpdate());
        this.subscribe(this.plugin.layout.events.updated, () => this.forceUpdate());
        this.subscribe(this.plugin.behaviors.interaction.selectionMode, () => this.forceUpdate());
    }

    _toggleSelMode = () => {
        this.plugin.selectionMode = !this.plugin.selectionMode;
    }

    render() {
        const style = this.props.inline
            ? { background: 'transparent', width: 'auto', height: 'auto', lineHeight: 'unset' }
            : { background: 'transparent' };
        return <IconButton svg={SelectionModeSvg} onClick={this._toggleSelMode} title={'Toggle Selection Mode'} style={style} toggleState={this.plugin.selectionMode} />;
    }
}

const StructureSelectionParams = {
    granularity: InteractivityManager.Params.granularity,
};

interface StructureSelectionActionsControlsState {
    isEmpty: boolean,
    isBusy: boolean,
    canUndo: boolean,

    action?: StructureSelectionModifier | 'theme' | 'add-component' | 'help'
}

const ActionHeader = new Map<StructureSelectionModifier, string>([
    ['add', 'Add/Union Selection'],
    ['remove', 'Remove/Subtract Selection'],
    ['intersect', 'Intersect Selection'],
    ['set', 'Set Selection']
] as const);

export class StructureSelectionActionsControls extends PluginUIComponent<{}, StructureSelectionActionsControlsState> {
    state = {
        action: void 0 as StructureSelectionActionsControlsState['action'],

        isEmpty: true,
        isBusy: false,
        canUndo: false,
    }

    componentDidMount() {
        this.subscribe(this.plugin.managers.structure.hierarchy.behaviors.selection, c => {
            const isEmpty = c.hierarchy.structures.length === 0;
            if (this.state.isEmpty !== isEmpty) {
                this.setState({ isEmpty });
            }
            // trigger elementQueries and nonStandardResidueQueries recalculation
            this.queriesVersion = -1;
            this.forceUpdate();
        });

        this.subscribe(this.plugin.behaviors.state.isBusy, v => {
            this.setState({ isBusy: v, action: void 0 });
        });

        this.subscribe(this.plugin.managers.interactivity.events.propsUpdated, () => {
            this.forceUpdate();
        });

        this.subscribe(this.plugin.state.data.events.historyUpdated, ({ state }) => {
            this.setState({ canUndo: state.canUndo });
        });
    }

    get isDisabled() {
        return this.state.isBusy || this.state.isEmpty;
    }

    set = (modifier: StructureSelectionModifier, selectionQuery: StructureSelectionQuery) => {
        this.plugin.managers.structure.selection.fromSelectionQuery(modifier, selectionQuery, false);
    }

    selectQuery: ActionMenu.OnSelect = (item, e) => {
        if (!item || !this.state.action) {
            this.setState({ action: void 0 });
            return;
        }
        const q = this.state.action! as StructureSelectionModifier;
        if (e?.shiftKey) {
            this.set(q, item.value as StructureSelectionQuery);
        } else {
            this.setState({ action: void 0 }, () => {
                this.set(q, item.value as StructureSelectionQuery);
            });
        }
    }

    get structures () {
        const structures: Structure[] = [];
        for (const s of this.plugin.managers.structure.hierarchy.selection.structures) {
            const structure = s.cell.obj?.data;
            if (structure) structures.push(structure);
        }
        return structures;
    }

    private queriesItems: ActionMenu.Items[] = []
    private queriesVersion = -1
    get queries () {
        const { registry } = this.plugin.query.structure;
        if (registry.version !== this.queriesVersion) {
            const structures = this.structures;
            const queries = [
                ...registry.list,
                ...getPolymerAndBranchedEntityQueries(structures),
                ...getNonStandardResidueQueries(structures),
                ...getElementQueries(structures)
            ].sort((a, b) => b.priority - a.priority);
            this.queriesItems = ActionMenu.createItems(queries, {
                filter: q => q !== StructureSelectionQueries.current && !q.isHidden,
                label: q => q.label,
                category: q => q.category,
                description: q => q.description
            });
            this.queriesVersion = registry.version;
        }
        return this.queriesItems;
    }

    private showAction(q: StructureSelectionActionsControlsState['action']) {
        return () => this.setState({ action: this.state.action === q ? void 0 : q });
    }

    toggleAdd = this.showAction('add')
    toggleRemove = this.showAction('remove')
    toggleIntersect = this.showAction('intersect')
    toggleSet = this.showAction('set')
    toggleTheme = this.showAction('theme')
    toggleAddComponent = this.showAction('add-component')
    toggleHelp = this.showAction('help')

    setGranuality: ParamOnChange = ({ value }) => {
        this.plugin.managers.interactivity.setProps({ granularity: value });
    }

    turnOff = () => this.plugin.selectionMode = false;

    undo = () => {
        const task = this.plugin.state.data.undo();
        if (task) this.plugin.runTask(task);
    }

    subtract = () => {
        const sel = this.plugin.managers.structure.hierarchy.getStructuresWithSelection();
        const components: StructureComponentRef[] = [];
        for (const s of sel) components.push(...s.components);
        if (components.length === 0) return;
        this.plugin.managers.structure.component.modifyByCurrentSelection(components, 'subtract');
    }

    render() {
        const granularity = this.plugin.managers.interactivity.props.granularity;
        const undoTitle = this.state.canUndo
            ? `Undo ${this.plugin.state.data.latestUndoLabel}`
            : 'Some mistakes of the past can be undone.';

        return <>
            <div className='msp-flex-row' style={{ background: 'none' }}>
                <PureSelectControl title={`Picking Level for selecting and highlighting`} param={StructureSelectionParams.granularity} name='granularity' value={granularity} onChange={this.setGranuality} isDisabled={this.isDisabled} />
                <ToggleButton icon={UnionSvg} title={`${ActionHeader.get('add')}. Hold shift key to keep menu open.`} toggle={this.toggleAdd} isSelected={this.state.action === 'add'} disabled={this.isDisabled} />
                <ToggleButton icon={SubtractSvg} title={`${ActionHeader.get('remove')}. Hold shift key to keep menu open.`} toggle={this.toggleRemove} isSelected={this.state.action === 'remove'} disabled={this.isDisabled} />
                <ToggleButton icon={IntersectSvg} title={`${ActionHeader.get('intersect')}. Hold shift key to keep menu open.`} toggle={this.toggleIntersect} isSelected={this.state.action === 'intersect'} disabled={this.isDisabled} />
                <ToggleButton icon={SetSvg} title={`${ActionHeader.get('set')}. Hold shift key to keep menu open.`} toggle={this.toggleSet} isSelected={this.state.action === 'set'} disabled={this.isDisabled} />

                <ToggleButton icon={BrushSvg} title='Apply Theme to Selection' toggle={this.toggleTheme} isSelected={this.state.action === 'theme'} disabled={this.isDisabled} style={{ marginLeft: '10px' }}  />
                <ToggleButton icon={CubeOutlineSvg} title='Create Component of Selection with Representation' toggle={this.toggleAddComponent} isSelected={this.state.action === 'add-component'} disabled={this.isDisabled} />
                <IconButton svg={RemoveSvg} title='Remove/subtract Selection from all Components' onClick={this.subtract} disabled={this.isDisabled} />
                <IconButton svg={RestoreSvg} onClick={this.undo} disabled={!this.state.canUndo || this.isDisabled} title={undoTitle} />

                <ToggleButton icon={HelpOutlineSvg} title='Show/hide help' toggle={this.toggleHelp} style={{ marginLeft: '10px' }} isSelected={this.state.action === 'help'} />
                <IconButton svg={CancelOutlinedSvg} title='Turn selection mode off' onClick={this.turnOff} />
            </div>
            {(this.state.action && this.state.action !== 'theme' && this.state.action !== 'add-component' && this.state.action !== 'help') && <div className='msp-selection-viewport-controls-actions'>
                <ActionMenu header={ActionHeader.get(this.state.action as StructureSelectionModifier)} title='Click to close.' items={this.queries} onSelect={this.selectQuery} noOffset />
            </div>}
            {this.state.action === 'theme' && <div className='msp-selection-viewport-controls-actions'>
                <ControlGroup header='Theme' title='Click to close.' initialExpanded={true} hideExpander={true} hideOffset={true} onHeaderClick={this.toggleTheme} topRightIcon={CloseSvg}>
                    <ApplyThemeControls onApply={this.toggleTheme} />
                </ControlGroup>
            </div>}
            {this.state.action === 'add-component' && <div className='msp-selection-viewport-controls-actions'>
                <ControlGroup header='Add Component' title='Click to close.' initialExpanded={true} hideExpander={true} hideOffset={true} onHeaderClick={this.toggleAddComponent} topRightIcon={CloseSvg}>
                    <AddComponentControls onApply={this.toggleAddComponent} forSelection />
                </ControlGroup>
            </div>}
            {this.state.action === 'help' && <div className='msp-selection-viewport-controls-actions'>
                <ControlGroup header='Help' title='Click to close.' initialExpanded={true} hideExpander={true} hideOffset={true} onHeaderClick={this.toggleHelp} topRightIcon={CloseSvg} maxHeight='300px'>
                    <HelpGroup header='Selection Operations'>
                        <HelpText>Use <Icon svg={UnionSvg} inline /> <Icon svg={SubtractSvg} inline /> <Icon svg={IntersectSvg} inline /> <Icon svg={SetSvg} inline /> to modify the selection.</HelpText>
                    </HelpGroup>
                    <HelpGroup header='Representation Operations'>
                        <HelpText>Use <Icon svg={BrushSvg} inline /> <Icon svg={CubeOutlineSvg} inline /> <Icon svg={RemoveSvg} inline /> <Icon svg={RestoreSvg} inline /> to color, create components, remove from components, or undo actions.</HelpText>
                    </HelpGroup>
                    <ViewportHelpContent selectOnly={true} />
                </ControlGroup>
            </div>}
        </>;
    }
}

export class StructureSelectionStatsControls extends PluginUIComponent<{ hideOnEmpty?: boolean }, { isEmpty: boolean, isBusy: boolean }> {
    state = {
        isEmpty: true,
        isBusy: false
    }

    componentDidMount() {
        this.subscribe(this.plugin.managers.structure.selection.events.changed, () => {
            this.forceUpdate();
        });

        this.subscribe(this.plugin.managers.structure.hierarchy.behaviors.selection, c => {
            const isEmpty = c.structures.length === 0;
            if (this.state.isEmpty !== isEmpty) {
                this.setState({ isEmpty });
            }
        });

        this.subscribe(this.plugin.behaviors.state.isBusy, v => {
            this.setState({ isBusy: v });
        });
    }

    get isDisabled() {
        return this.state.isBusy || this.state.isEmpty;
    }

    get stats() {
        const stats = this.plugin.managers.structure.selection.stats;
        if (stats.structureCount === 0 || stats.elementCount === 0) {
            return 'Nothing Selected';
        } else {
            return `${stripTags(stats.label)} Selected`;
        }
    }

    clear = () => this.plugin.managers.interactivity.lociSelects.deselectAll();

    focus = () => {
        if (this.plugin.managers.structure.selection.stats.elementCount === 0) return;
        const { sphere } = this.plugin.managers.structure.selection.getBoundary();
        this.plugin.managers.camera.focusSphere(sphere);
    }

    highlight = (e: React.MouseEvent<HTMLElement>) => {
        this.plugin.managers.interactivity.lociHighlights.clearHighlights();
        this.plugin.managers.structure.selection.entries.forEach(e => {
            this.plugin.managers.interactivity.lociHighlights.highlight({ loci: e.selection }, false);
        });
    }

    clearHighlight = () => {
        this.plugin.managers.interactivity.lociHighlights.clearHighlights();
    }

    render() {
        const stats = this.plugin.managers.structure.selection.stats;
        const empty = stats.structureCount === 0 || stats.elementCount === 0;

        if (empty && this.props.hideOnEmpty) return null;

        return <>
            <div className='msp-flex-row'>
                <Button noOverflow onClick={this.focus} title='Click to Focus Selection' disabled={empty} onMouseEnter={this.highlight} onMouseLeave={this.clearHighlight}
                    style={{ textAlignLast: !empty ? 'left' : void 0 }}>
                    {this.stats}
                </Button>
                {!empty && <IconButton svg={CancelOutlinedSvg} onClick={this.clear} title='Clear' className='msp-form-control' flex />}
            </div>
        </>;
    }
}

interface ApplyThemeControlsState {
    values: StructureComponentManager.ThemeParams
}

interface ApplyThemeControlsProps {
    onApply?: () => void
}

class ApplyThemeControls extends PurePluginUIComponent<ApplyThemeControlsProps, ApplyThemeControlsState> {
    _params = memoizeLatest((pivot: StructureRef | undefined) => StructureComponentManager.getThemeParams(this.plugin, pivot));
    get params() { return this._params(this.plugin.managers.structure.component.pivotStructure); }

    state = { values: ParamDefinition.getDefaultValues(this.params) };

    apply = () => {
        this.plugin.managers.structure.component.applyTheme(this.state.values, this.plugin.managers.structure.hierarchy.current.structures);
        this.props.onApply?.();
    }

    paramsChanged = (values: any) => this.setState({ values })

    render() {
        return <>
            <ParameterControls params={this.params} values={this.state.values} onChangeValues={this.paramsChanged} />
            <Button icon={BrushSvg} className='msp-btn-commit msp-btn-commit-on' onClick={this.apply} style={{ marginTop: '1px' }}>
                Apply Theme
            </Button>
        </>;
    }
}