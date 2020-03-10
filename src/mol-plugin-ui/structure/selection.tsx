/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react';
import { CollapsableControls, CollapsableState } from '../base';
import { StructureSelectionQuery, StructureSelectionQueryList } from '../../mol-plugin-state/helpers/structure-selection-query';
import { PluginCommands } from '../../mol-plugin/commands';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { Interactivity } from '../../mol-plugin/util/interactivity';
import { ParameterControls } from '../controls/parameters';
import { stripTags } from '../../mol-util/string';
import { StructureElement } from '../../mol-model/structure';
import { ActionMenu } from '../controls/action-menu';
import { ToggleButton, ExpandGroup } from '../controls/common';
import { Icon } from '../controls/icons';
import { StructureSelectionModifier } from '../../mol-plugin-state/manager/structure/selection';

export const DefaultQueries = ActionMenu.createItems(StructureSelectionQueryList, {
    label: q => q.label,
    category: q => q.category
});

const StructureSelectionParams = {
    granularity: Interactivity.Params.granularity,
}

interface StructureSelectionControlsState extends CollapsableState {
    minRadius: number,
    extraRadius: number,
    durationMs: number,

    isDisabled: boolean,

    queryAction?: StructureSelectionModifier
}

export class StructureSelectionControls<P, S extends StructureSelectionControlsState> extends CollapsableControls<P, S> {
    componentDidMount() {
        this.subscribe(this.plugin.managers.structure.selection.events.changed, () => {
            this.forceUpdate()
        });

        this.subscribe(this.plugin.events.interactivity.propsUpdated, () => {
            this.forceUpdate()
        });

        this.subscribe(this.plugin.behaviors.state.isBusy, v => {
            this.setState({ isDisabled: v, queryAction: void 0 })
        })
    }

    get stats() {
        const stats = this.plugin.managers.structure.selection.stats
        if (stats.structureCount === 0 || stats.elementCount === 0) {
            return 'Nothing Selected'
        } else {
            return `Selected ${stripTags(stats.label)}`
        }
    }

    focus = () => {
        const { extraRadius, minRadius, durationMs } = this.state
        if (this.plugin.managers.structure.selection.stats.elementCount === 0) return
        const principalAxes = this.plugin.managers.structure.selection.getPrincipalAxes();
        const { origin, dirA, dirC } = principalAxes.boxAxes
        const { sphere } = this.plugin.managers.structure.selection.getBoundary()
        const radius = Math.max(sphere.radius + extraRadius, minRadius);
        this.plugin.canvas3d?.camera.focus(origin, radius, this.plugin.canvas3d.boundingSphere.radius, durationMs, dirA, dirC);
    }

    focusLoci(loci: StructureElement.Loci) {
        return () => {
            const { extraRadius, minRadius, durationMs } = this.state
            if (this.plugin.managers.structure.selection.stats.elementCount === 0) return
            const { sphere } = StructureElement.Loci.getBoundary(loci)
            const radius = Math.max(sphere.radius + extraRadius, minRadius);
            this.plugin.canvas3d?.camera.focus(sphere.center, radius, this.plugin.canvas3d.boundingSphere.radius, durationMs);
        }
    }

    setProps = (p: { param: PD.Base<any>, name: string, value: any }) => {
        if (p.name === 'granularity') {
            PluginCommands.Interactivity.SetProps(this.plugin, { props: { granularity: p.value } });
        }
    }

    get values () {
        return {
            granularity: this.plugin.interactivity.props.granularity,
        }
    }

    set = (modifier: StructureSelectionModifier, selectionQuery: StructureSelectionQuery) => {
        this.plugin.managers.structure.selection.fromSelectionQuery(modifier, selectionQuery, false)
    }

    selectQuery: ActionMenu.OnSelect = item => {
        if (!item || !this.state.queryAction) {
            this.setState({ queryAction: void 0 });
            return;
        }
        const q = this.state.queryAction!;
        this.setState({ queryAction: void 0 }, () => {
            this.set(q, item.value as StructureSelectionQuery);
        })
    }

    queries = DefaultQueries

    private showQueries(q: StructureSelectionModifier) {
        return () => this.setState({ queryAction: this.state.queryAction === q ? void 0 : q });
    }

    toggleAdd = this.showQueries('add')
    toggleRemove = this.showQueries('remove')
    toggleOnly = this.showQueries('set')

    get controls() {
        return <div>
            <div className='msp-control-row msp-select-row'>
                <ToggleButton icon='plus' label='Select' toggle={this.toggleAdd} isSelected={this.state.queryAction === 'add'} disabled={this.state.isDisabled} />
                <ToggleButton icon='minus' label='Deselect' toggle={this.toggleRemove} isSelected={this.state.queryAction === 'remove'} disabled={this.state.isDisabled} />
                <ToggleButton icon='flash' label='Only' toggle={this.toggleOnly} isSelected={this.state.queryAction === 'set'} disabled={this.state.isDisabled} />
            </div>
            {this.state.queryAction && <ActionMenu items={this.queries} onSelect={this.selectQuery} />}
        </div>
    }

    defaultState() {
        return {
            isCollapsed: false,
            header: 'Selection',

            minRadius: 8,
            extraRadius: 4,
            durationMs: 250,

            queryAction: void 0,

            isDisabled: false
        } as S
    }

    renderControls() {
        const history: JSX.Element[] = [];

        const mng = this.plugin.managers.structure.selection;

        // TODO: fix the styles, move them to CSS

        for (let i = 0, _i = Math.min(4, mng.history.length); i < _i; i++) {
            const e = mng.history[i];
            history.push(<li key={e!.label}>
                <button className='msp-btn msp-btn-block msp-form-control' style={{ borderRight: '6px solid transparent', overflow: 'hidden' }}
                    title='Click to focus.' onClick={this.focusLoci(e.loci)}>
                    <span dangerouslySetInnerHTML={{ __html: e.label.split('|').reverse().join(' | ') }} />
                </button>
                {/* <div>
                    <IconButton icon='remove' title='Remove' onClick={() => {}} />
                </div> */}
            </li>)
        }

        return <div>
            <div className='msp-control-row msp-row-text'>
                <button className='msp-btn msp-btn-block' onClick={this.focus}>
                    <Icon name='focus-on-visual' style={{ position: 'absolute', left: '5px' }} />
                    {this.stats}
                </button>
            </div>
            <ParameterControls params={StructureSelectionParams} values={this.values} onChange={this.setProps} isDisabled={this.state.isDisabled} />
            {this.controls}
            {history.length > 0 && <ExpandGroup header='Selection History'>
                <ul style={{ listStyle: 'none', marginTop: '1px', marginBottom: '0' }} className='msp-state-list'>
                    {history}
                </ul>
            </ExpandGroup>}
        </div>
    }
}