/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react';
import { StructureElement } from '../../mol-model/structure';
import { StructureSelectionQueries, StructureSelectionQuery, StructureSelectionQueryList } from '../../mol-plugin-state/helpers/structure-selection-query';
import { InteractivityManager } from '../../mol-plugin-state/manager/interactivity';
import { StructureComponentManager } from '../../mol-plugin-state/manager/structure/component';
import { StructureRef } from '../../mol-plugin-state/manager/structure/hierarchy-state';
import { StructureSelectionModifier } from '../../mol-plugin-state/manager/structure/selection';
import { memoize1 } from '../../mol-util/memoize';
import { ParamDefinition } from '../../mol-util/param-definition';
import { stripTags } from '../../mol-util/string';
import { CollapsableControls, CollapsableState, PurePluginUIComponent } from '../base';
import { ActionMenu } from '../controls/action-menu';
import { ExpandGroup, ToggleButton } from '../controls/common';
import { Icon } from '../controls/icons';
import { ParameterControls } from '../controls/parameters';

export const DefaultQueries = ActionMenu.createItems(StructureSelectionQueryList, {
    filter: q => q !== StructureSelectionQueries.current,
    label: q => q.label,
    category: q => q.category
});

const StructureSelectionParams = {
    granularity: InteractivityManager.Params.granularity,
}

interface StructureSelectionControlsState extends CollapsableState {
    minRadius: number,
    extraRadius: number,
    durationMs: number,

    isEmpty: boolean,
    isBusy: boolean,

    action?: StructureSelectionModifier | 'color'
}

export class StructureSelectionControls<P, S extends StructureSelectionControlsState> extends CollapsableControls<P, S> {
    componentDidMount() {
        this.subscribe(this.plugin.managers.structure.selection.events.changed, () => {
            this.forceUpdate()
        });

        this.subscribe(this.plugin.managers.interactivity.events.propsUpdated, () => {
            this.forceUpdate()
        });

        this.subscribe(this.plugin.managers.structure.hierarchy.behaviors.current, c => {
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

    queries = DefaultQueries

    private showAction(q: StructureSelectionControlsState['action']) {
        return () => this.setState({ action: this.state.action === q ? void 0 : q });
    }

    toggleAdd = this.showAction('add')
    toggleRemove = this.showAction('remove')
    toggleIntersect = this.showAction('intersect')
    toggleSet = this.showAction('set')
    toggleColor = this.showAction('color')

    // TODO better icons
    get controls() {
        return <div>
            <div className='msp-control-row msp-select-row'>
                <ToggleButton icon='plus' title='Add' toggle={this.toggleAdd} isSelected={this.state.action === 'add'} disabled={this.isDisabled} />
                <ToggleButton icon='minus' title='Remove' toggle={this.toggleRemove} isSelected={this.state.action === 'remove'} disabled={this.isDisabled} />
                <ToggleButton icon='star' title='Intersect' toggle={this.toggleIntersect} isSelected={this.state.action === 'intersect'} disabled={this.isDisabled} />
                <ToggleButton icon='flash' title='Set' toggle={this.toggleSet} isSelected={this.state.action === 'set'} disabled={this.isDisabled} />
                <ToggleButton icon='brush' title='Color' toggle={this.toggleColor} isSelected={this.state.action === 'color'} disabled={this.isDisabled} />
            </div>
            {(this.state.action && this.state.action !== 'color') && <ActionMenu items={this.queries} onSelect={this.selectQuery} />}
            {this.state.action === 'color' && <div className='msp-control-offset'><ApplyColorControls /></div>}
        </div>
    }

    defaultState() {
        return {
            isCollapsed: false,
            header: 'Selection',

            minRadius: 8,
            extraRadius: 4,
            durationMs: 250,

            action: void 0,

            isEmpty: true,
            isBusy: false
        } as S
    }

    renderControls() {
        const history: JSX.Element[] = [];

        const mng = this.plugin.managers.structure.selection;

        // TODO: fix the styles, move them to CSS

        for (let i = 0, _i = Math.min(4, mng.history.length); i < _i; i++) {
            const e = mng.history[i];
            history.push(<li key={e!.label}>
                <button className='msp-btn msp-btn-block msp-form-control' style={{ overflow: 'hidden' }}
                    title='Click to focus.' onClick={this.focusLoci(e.loci)}>
                    <span dangerouslySetInnerHTML={{ __html: e.label.split('|').reverse().join(' | ') }} />
                </button>
                {/* <div>
                    <IconButton icon='remove' title='Remove' onClick={() => {}} />
                </div> */}
            </li>)
        }

        return <>
            <ParameterControls params={StructureSelectionParams} values={this.values} onChangeObject={this.setProps} />
            {this.controls}
            <div className='msp-control-row msp-row-text' style={{ marginTop: '6px' }}>
                <button className='msp-btn msp-btn-block' onClick={this.focus}>
                    <Icon name='focus-on-visual' style={{ position: 'absolute', left: '5px' }} />
                    {this.stats}
                </button>
            </div>
            {history.length > 0 && <ExpandGroup header='Selection History'>
                <ul style={{ listStyle: 'none', marginTop: '1px', marginBottom: '0' }} className='msp-state-list'>
                    {history}
                </ul>
            </ExpandGroup>}
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
    _params = memoize1((pivot: StructureRef | undefined) => StructureComponentManager.getColorParams(this.plugin, pivot));
    get params() { return this._params(this.plugin.managers.structure.component.pivotStructure); }

    state = { values: ParamDefinition.getDefaultValues(this.params) };

    apply = () => {
        this.plugin.managers.structure.component.applyColor(this.state.values);
        this.props.onApply?.();
    }

    paramsChanged = (values: any) => this.setState({ values })

    render() {
        return <>
            <ParameterControls params={this.params} values={this.state.values} onChangeObject={this.paramsChanged} />
            <button className={`msp-btn msp-btn-block msp-btn-commit msp-btn-commit-on`} onClick={this.apply} style={{ marginTop: '1px' }}>
                <Icon name='brush' /> Apply Coloring
            </button>
        </>;
    }
}