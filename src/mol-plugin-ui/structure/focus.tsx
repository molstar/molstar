/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react';
import { CollapsableState, CollapsableControls } from '../base';
import { ToggleButton } from '../controls/common';
import { StructureHierarchyManager } from '../../mol-plugin-state/manager/structure/hierarchy';
import { ActionMenu } from '../controls/action-menu';
import { stringToWords } from '../../mol-util/string';
import { StructureElement, StructureProperties } from '../../mol-model/structure';
import { OrderedSet, SortedArray } from '../../mol-data/int';
import { UnitIndex } from '../../mol-model/structure/structure/element/element';
import { FocusEntry } from '../../mol-plugin-state/manager/structure/focus';
import { Icon } from '../controls/icons';

type FocusAction = 'presets' | 'history'

interface StructureComponentControlState extends CollapsableState {
    isBusy: boolean
    action?: FocusAction
}

export class StructureFocusControls extends CollapsableControls<{}, StructureComponentControlState> {
    protected defaultState(): StructureComponentControlState {
        return {
            header: 'Focus',
            isCollapsed: false,
            isBusy: false,
        };
    }

    componentDidMount() {
        this.subscribe(this.plugin.managers.structure.hierarchy.behaviors.selection, c => {
            this.setState({
                description: StructureHierarchyManager.getSelectedStructuresDescription(this.plugin)
            });
            // if setState is called on non-pure component, forceUpdate is reduntant
            // this.forceUpdate();
        });
    }

    get presetsItems() {
        const items: FocusEntry[] = []
        const l = StructureElement.Location.create()
        const { structures } = this.plugin.managers.structure.hierarchy.selection;
        for (const s of structures) {
            const d = s.cell.obj?.data
            if (d) {
                l.structure = d
                for (const ug of d.unitSymmetryGroups) {
                    l.unit = ug.units[0]
                    l.element = ug.elements[0]
                    const et = StructureProperties.entity.type(l)
                    if (et === 'non-polymer') {
                        const idx = SortedArray.indexOf(ug.elements, l.element) as UnitIndex
                        const loci = StructureElement.Loci(d, [{ unit: l.unit, indices: OrderedSet.ofSingleton(idx) }])
                        items.push({
                            label: StructureProperties.entity.pdbx_description(l).join(', '),
                            loci: StructureElement.Loci.extendToWholeResidues(loci)
                        })
                    }
                }
            }
        }

        return items
    }

    get historyItems() {
        return this.plugin.managers.structure.focus.history
    }

    get actionItems() {
        let items: FocusEntry[]
        switch (this.state.action) {
            case 'presets': items = this.presetsItems; break
            case 'history': items = this.historyItems; break
            default: items = []
        }
        return ActionMenu.createItems(items, {
            label: f => f.label,
            category: f => f.category
        })
    }

    selectAction: ActionMenu.OnSelect = item => {
        if (!item || !this.state.action) {
            this.setState({ action: void 0 });
            return;
        }
        this.setState({ action: void 0 }, async () => {
            const f = item.value as FocusEntry
            this.plugin.managers.structure.focus.set(f)
            this.plugin.managers.camera.focusLoci(f.loci, { durationMs: 0 })
        })
    }

    private showAction(a: FocusAction) {
        return () => this.setState({ action: this.state.action === a ? void 0 : a });
    }

    togglePresets = this.showAction('presets')
    toggleHistory = this.showAction('history')

    focus = () => {
        const { current } = this.plugin.managers.structure.focus
        if (current) this.plugin.managers.camera.focusLoci(current.loci);
    }

    highlightCurrent = () => {
        const { current } = this.plugin.managers.structure.focus
        if (current) this.plugin.managers.interactivity.lociHighlights.highlightOnly({ loci: current.loci }, false);
    }

    clearHighlights = () => {
        this.plugin.managers.interactivity.lociHighlights.clearHighlights()
    }

    renderControls() {
        const { current } = this.plugin.managers.structure.focus
        const label = current?.label || 'Nothing Focused'

        return <>
            <div className='msp-control-row msp-select-row'>
                <ToggleButton icon='bookmarks' title='Preset' label='Preset' toggle={this.togglePresets} isSelected={this.state.action === 'presets'} disabled={this.state.isBusy} />
                <ToggleButton icon='clock' title='History' label='History' toggle={this.toggleHistory} isSelected={this.state.action === 'history'} disabled={this.state.isBusy} />
            </div>
            {this.state.action && <ActionMenu header={stringToWords(this.state.action)} items={this.actionItems} onSelect={this.selectAction} />}
            <div className='msp-control-row msp-row-text' style={{ marginTop: '6px' }}>
                <button className='msp-btn msp-btn-block msp-no-overflow' onClick={this.focus} title='Click to Center Focused' disabled={!current} onMouseEnter={this.highlightCurrent} onMouseLeave={this.clearHighlights}>
                    <Icon name='focus-on-visual' style={{ position: 'absolute', left: '5px' }} />
                    {label}
                </button>
            </div>
        </>;
    }
}