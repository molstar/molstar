/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react';
import { StructureSelectionQueries, StructureSelectionQuery } from '../../mol-plugin-state/helpers/structure-selection-query';
import { InteractivityManager } from '../../mol-plugin-state/manager/interactivity';
import { StructureComponentManager } from '../../mol-plugin-state/manager/structure/component';
import { StructureRef } from '../../mol-plugin-state/manager/structure/hierarchy-state';
import { StructureSelectionModifier } from '../../mol-plugin-state/manager/structure/selection';
import { memoizeLatest } from '../../mol-util/memoize';
import { ParamDefinition } from '../../mol-util/param-definition';
import { stripTags } from '../../mol-util/string';
import { CollapsableControls, CollapsableState, PurePluginUIComponent } from '../base';
import { ActionMenu } from '../controls/action-menu';
import { ControlGroup, ToggleButton, IconButton } from '../controls/common';
import { Icon } from '../controls/icons';
import { ParameterControls } from '../controls/parameters';
import { StructureMeasurementsControls } from './measurements';

const StructureSelectionParams = {
    granularity: InteractivityManager.Params.granularity,
}

interface StructureSelectionControlsState extends CollapsableState {
    isEmpty: boolean,
    isBusy: boolean,

    action?: StructureSelectionModifier | 'color'
}

const ActionHeader = new Map<StructureSelectionModifier, string>([
    ['add', 'Add/Union'],
    ['remove', 'Remove/Subtract'],
    ['intersect', 'Intersect'],
    ['set', 'Set']
] as const);

export class StructureSelectionControls<P, S extends StructureSelectionControlsState> extends CollapsableControls<P, S> {
    componentDidMount() {
        this.subscribe(this.plugin.managers.structure.selection.events.changed, () => {
            this.forceUpdate()
        });

        this.subscribe(this.plugin.managers.interactivity.events.propsUpdated, () => {
            this.forceUpdate()
        });

        this.subscribe(this.plugin.managers.structure.hierarchy.behaviors.selection, c => {
            const isEmpty = c.structures.length === 0;
            if (this.state.isEmpty !== isEmpty) {
                this.setState({ isEmpty });
            }
        });

        this.subscribe(this.plugin.behaviors.state.isBusy, v => {
            this.setState({ isBusy: v, action: void 0 })
        })
    }

    get isDisabled() {
        return this.state.isBusy || this.state.isEmpty
    }

    get stats() {
        const stats = this.plugin.managers.structure.selection.stats
        if (stats.structureCount === 0 || stats.elementCount === 0) {
            return 'Nothing Selected'
        } else {
            return `${stripTags(stats.label)} Selected`
        }
    }

    clear = () => this.plugin.managers.interactivity.lociSelects.deselectAll();

    focus = () => {
        if (this.plugin.managers.structure.selection.stats.elementCount === 0) return;
        const principalAxes = this.plugin.managers.structure.selection.getPrincipalAxes();
        const { sphere } = this.plugin.managers.structure.selection.getBoundary();
        this.plugin.managers.camera.focusSphere(sphere, { principalAxes });
    }

    setProps = (props: any) => {
        this.plugin.managers.interactivity.setProps(props);
    }

    get values () {
        return {
            granularity: this.plugin.managers.interactivity.props.granularity,
        }
    }

    set = (modifier: StructureSelectionModifier, selectionQuery: StructureSelectionQuery) => {
        this.plugin.managers.structure.selection.fromSelectionQuery(modifier, selectionQuery, false)
    }

    selectQuery: ActionMenu.OnSelect = item => {
        if (!item || !this.state.action) {
            this.setState({ action: void 0 });
            return;
        }
        const q = this.state.action! as StructureSelectionModifier;
        this.setState({ action: void 0 }, () => {
            this.set(q, item.value as StructureSelectionQuery);
        })
    }

    private queriesItems: ActionMenu.Items[] = []
    private queriesVersion = -1
    get queries () {
        const { registry } = this.plugin.query.structure
        if (registry.version !== this.queriesVersion) {
            this.queriesItems = ActionMenu.createItems(registry.list, {
                filter: q => q !== StructureSelectionQueries.current,
                label: q => q.label,
                category: q => q.category
            });
            this.queriesVersion = registry.version
        }
        return this.queriesItems
    }

    private showAction(q: StructureSelectionControlsState['action']) {
        return () => this.setState({ action: this.state.action === q ? void 0 : q });
    }

    toggleAdd = this.showAction('add')
    toggleRemove = this.showAction('remove')
    toggleIntersect = this.showAction('intersect')
    toggleSet = this.showAction('set')
    toggleColor = this.showAction('color')

    get controls() {
        return <>
            <div className='msp-control-row msp-select-row'>
                <ToggleButton icon='union' title={ActionHeader.get('add')} toggle={this.toggleAdd} isSelected={this.state.action === 'add'} disabled={this.isDisabled} />
                <ToggleButton icon='subtract' title={ActionHeader.get('remove')} toggle={this.toggleRemove} isSelected={this.state.action === 'remove'} disabled={this.isDisabled} />
                <ToggleButton icon='intersect' title={ActionHeader.get('intersect')} toggle={this.toggleIntersect} isSelected={this.state.action === 'intersect'} disabled={this.isDisabled} />
                <ToggleButton icon='set' title={ActionHeader.get('set')} toggle={this.toggleSet} isSelected={this.state.action === 'set'} disabled={this.isDisabled} />
                <ToggleButton icon='brush' title='Color' toggle={this.toggleColor} isSelected={this.state.action === 'color'} disabled={this.isDisabled} />
            </div>
            {(this.state.action && this.state.action !== 'color') && <ActionMenu header={ActionHeader.get(this.state.action as StructureSelectionModifier)} items={this.queries} onSelect={this.selectQuery} />}
            {this.state.action === 'color' && <ControlGroup header='Color' initialExpanded={true} hideExpander={true} hideOffset={false} onHeaderClick={this.toggleColor} topRightIcon='off'>
                <ApplyColorControls />
            </ControlGroup>}
        </>
    }

    defaultState() {
        return {
            isCollapsed: false,
            header: 'Selection',

            action: void 0,

            isEmpty: true,
            isBusy: false
        } as S
    }

    renderControls() {
        const stats = this.plugin.managers.structure.selection.stats
        const empty = stats.structureCount === 0 || stats.elementCount === 0;

        return <>
            <ParameterControls params={StructureSelectionParams} values={this.values} onChangeValues={this.setProps} />
            {this.controls}
            <div className='msp-control-row msp-select-row' style={{ margin: '6px 0' }}>
                <button className='msp-btn msp-btn-block msp-no-overflow' onClick={this.focus} title='Click to Focus Selection' disabled={empty}
                    style={{ textAlignLast: !empty ? 'left' : void 0 }}>
                    {this.stats}
                </button>
                {!empty && <IconButton onClick={this.clear} icon='cancel' title='Clear' customClass='msp-form-control' flex />}
            </div>
            <StructureMeasurementsControls />
        </>
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
            <button className={`msp-btn msp-btn-block msp-btn-commit msp-btn-commit-on`} onClick={this.apply} style={{ marginTop: '1px' }}>
                <Icon name='brush' /> Apply Coloring
            </button>
        </>;
    }
}