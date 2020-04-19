/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import Close from '@material-ui/icons/Close';
import Clear from '@material-ui/icons/Clear';
import Brush from '@material-ui/icons/Brush';
import * as React from 'react';
import { StructureSelectionQueries, StructureSelectionQuery } from '../../mol-plugin-state/helpers/structure-selection-query';
import { InteractivityManager } from '../../mol-plugin-state/manager/interactivity';
import { StructureComponentManager } from '../../mol-plugin-state/manager/structure/component';
import { StructureRef } from '../../mol-plugin-state/manager/structure/hierarchy-state';
import { StructureSelectionModifier } from '../../mol-plugin-state/manager/structure/selection';
import { memoizeLatest } from '../../mol-util/memoize';
import { ParamDefinition } from '../../mol-util/param-definition';
import { stripTags } from '../../mol-util/string';
import { PluginUIComponent, PurePluginUIComponent } from '../base';
import { ActionMenu } from '../controls/action-menu';
import { Button, ControlGroup, IconButton, ToggleButton } from '../controls/common';
import { ParameterControls, ParamOnChange, PureSelectControl } from '../controls/parameters';
import { Union, Subtract, Intersect, SetSvg as SetSvg, CubeSvg } from '../controls/icons';
import { AddComponentControls } from './components';

const StructureSelectionParams = {
    granularity: InteractivityManager.Params.granularity,
};

interface StructureSelectionActionsControlsState {
    isEmpty: boolean,
    isBusy: boolean,

    action?: StructureSelectionModifier | 'color' | 'add-repr'
}

const ActionHeader = new Map<StructureSelectionModifier, string>([
    ['add', 'Add/Union'],
    ['remove', 'Remove/Subtract'],
    ['intersect', 'Intersect'],
    ['set', 'Set']
] as const);

export class StructureSelectionActionsControls extends PluginUIComponent<{}, StructureSelectionActionsControlsState> {
    state = {
        action: void 0 as StructureSelectionActionsControlsState['action'],

        isEmpty: true,
        isBusy: false,
    }

    componentDidMount() {
        this.subscribe(this.plugin.managers.structure.hierarchy.behaviors.selection, c => {
            const isEmpty = c.hierarchy.structures.length === 0;
            if (this.state.isEmpty !== isEmpty) {
                this.setState({ isEmpty });
            }
        });

        this.subscribe(this.plugin.behaviors.state.isBusy, v => {
            this.setState({ isBusy: v, action: void 0 });
        });

        this.subscribe(this.plugin.managers.interactivity.events.propsUpdated, () => {
            this.forceUpdate();
        });
    }

    get isDisabled() {
        return this.state.isBusy || this.state.isEmpty;
    }

    set = (modifier: StructureSelectionModifier, selectionQuery: StructureSelectionQuery) => {
        this.plugin.managers.structure.selection.fromSelectionQuery(modifier, selectionQuery, false);
    }

    selectQuery: ActionMenu.OnSelect = item => {
        if (!item || !this.state.action) {
            this.setState({ action: void 0 });
            return;
        }
        const q = this.state.action! as StructureSelectionModifier;
        this.setState({ action: void 0 }, () => {
            this.set(q, item.value as StructureSelectionQuery);
        });
    }

    private queriesItems: ActionMenu.Items[] = []
    private queriesVersion = -1
    get queries () {
        const { registry } = this.plugin.query.structure;
        if (registry.version !== this.queriesVersion) {
            this.queriesItems = ActionMenu.createItems(registry.list, {
                filter: q => q !== StructureSelectionQueries.current,
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
    toggleColor = this.showAction('color')
    toggleAddRepr = this.showAction('add-repr')

    setGranuality: ParamOnChange = ({ value }) => {
        this.plugin.managers.interactivity.setProps({ granularity: value });
    }

    turnOff = () => this.plugin.selectionMode = false;

    render() {
        const granularity = this.plugin.managers.interactivity.props.granularity;
        return <>
            <div className='msp-flex-row'>
                <ToggleButton icon={Union} title={ActionHeader.get('add')} toggle={this.toggleAdd} isSelected={this.state.action === 'add'} disabled={this.isDisabled} />
                <ToggleButton icon={Subtract} title={ActionHeader.get('remove')} toggle={this.toggleRemove} isSelected={this.state.action === 'remove'} disabled={this.isDisabled} />
                <ToggleButton icon={Intersect} title={ActionHeader.get('intersect')} toggle={this.toggleIntersect} isSelected={this.state.action === 'intersect'} disabled={this.isDisabled} />
                <ToggleButton icon={SetSvg} title={ActionHeader.get('set')} toggle={this.toggleSet} isSelected={this.state.action === 'set'} disabled={this.isDisabled} />
                <ToggleButton icon={Brush} title='Color' toggle={this.toggleColor} isSelected={this.state.action === 'color'} disabled={this.isDisabled} />
                <ToggleButton icon={CubeSvg} title='Create Representation' toggle={this.toggleAddRepr} isSelected={this.state.action === 'add-repr'} disabled={this.isDisabled} />
                <PureSelectControl title={`Picking Level`} param={StructureSelectionParams.granularity} name='granularity' value={granularity} onChange={this.setGranuality} isDisabled={this.isDisabled} />
                <IconButton svg={Close} title='Turn selection mode off' onClick={this.turnOff} />
            </div>
            {(this.state.action && this.state.action !== 'color' && this.state.action !== 'add-repr') && <div className='msp-selection-viewport-controls-actions'>
                <ActionMenu header={ActionHeader.get(this.state.action as StructureSelectionModifier)} items={this.queries} onSelect={this.selectQuery} noOffset />
            </div>}
            {this.state.action === 'color' && <div className='msp-selection-viewport-controls-actions'>
                <ControlGroup header='Color' initialExpanded={true} hideExpander={true} hideOffset={true} onHeaderClick={this.toggleColor} topRightIcon={Close}>
                    <ApplyColorControls />
                </ControlGroup>
            </div>}
            {this.state.action === 'add-repr' && <div className='msp-selection-viewport-controls-actions'>
                <ControlGroup header='Add Representation' initialExpanded={true} hideExpander={true} hideOffset={true} onHeaderClick={this.toggleAddRepr} topRightIcon={Close}>
                    <AddComponentControls onApply={this.toggleAddRepr} forSelection />
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
                {!empty && <IconButton svg={Clear} onClick={this.clear} title='Clear' className='msp-form-control' flex />}
            </div>
        </>;
    }
}

interface ApplyColorControlsState {
    values: StructureComponentManager.ColorParams
}

interface ApplyColorControlsProps {
    onApply?: () => void
}

class ApplyColorControls extends PurePluginUIComponent<ApplyColorControlsProps, ApplyColorControlsState> {
    _params = memoizeLatest((pivot: StructureRef | undefined) => StructureComponentManager.getColorParams(this.plugin, pivot));
    get params() { return this._params(this.plugin.managers.structure.component.pivotStructure); }

    state = { values: ParamDefinition.getDefaultValues(this.params) };

    apply = () => {
        this.plugin.managers.structure.component.applyColor(this.state.values);
        this.props.onApply?.();
    }

    paramsChanged = (values: any) => this.setState({ values })

    render() {
        return <>
            <ParameterControls params={this.params} values={this.state.values} onChangeValues={this.paramsChanged} />
            <Button icon={Brush} className='msp-btn-commit msp-btn-commit-on' onClick={this.apply} style={{ marginTop: '1px' }}>
                Apply Coloring
            </Button>
        </>;
    }
}