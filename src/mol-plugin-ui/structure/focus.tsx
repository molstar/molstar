/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react';
import { PluginUIComponent } from '../base';
import { ToggleButton } from '../controls/common';
import { ActionMenu } from '../controls/action-menu';
import { StructureElement, StructureProperties, Structure } from '../../mol-model/structure';
import { OrderedSet, SortedArray } from '../../mol-data/int';
import { UnitIndex } from '../../mol-model/structure/structure/element/element';
import { FocusEntry } from '../../mol-plugin-state/manager/structure/focus';
import { lociLabel } from '../../mol-theme/label';

interface StructureFocusControlsState {
    isBusy: boolean
    showAction: boolean
}

function getFocusEntries(structure: Structure) {
    const entityEntries = new Map<string, FocusEntry[]>()
    const l = StructureElement.Location.create(structure)

    for (const ug of structure.unitSymmetryGroups) {
        l.unit = ug.units[0]
        l.element = ug.elements[0]
        const et = StructureProperties.entity.type(l)
        if (et === 'non-polymer') {
            for (const u of ug.units) {
                l.unit = u
                const idx = SortedArray.indexOf(u.elements, l.element) as UnitIndex
                const loci = StructureElement.Loci.extendToWholeResidues(
                    StructureElement.Loci(structure, [
                        { unit: l.unit, indices: OrderedSet.ofSingleton(idx) }
                    ])
                )
                let label = lociLabel(loci, { reverse: true, hidePrefix: true, htmlStyling: false })
                if (ug.units.length > 1) {
                    label += ` | ${u.conformation.operator.name}`
                }
                const name = StructureProperties.entity.pdbx_description(l).join(', ')
                const item: FocusEntry = { label, category: name, loci }

                if (entityEntries.has(name)) entityEntries.get(name)!.push(item)
                else entityEntries.set(name, [item])
            }
        } else if (et === 'branched') {
            // TODO split into residues
        }
    }

    const entries: FocusEntry[] = []
    entityEntries.forEach((e, name) => {
        if (e.length === 1) {
            entries.push({ label: name, loci: e[0].loci })
        } else {
            entries.push(...e)
        }
    })

    return entries
}

export class StructureFocusControls extends PluginUIComponent<{}, StructureFocusControlsState> {
    state = { isBusy: false, showAction: false }

    componentDidMount() {
        this.subscribe(this.plugin.managers.structure.focus.events.changed, c => {
            this.forceUpdate();
        });

        this.subscribe(this.plugin.managers.structure.focus.events.historyUpdated, c => {
            this.forceUpdate();
        });

        this.subscribe(this.plugin.behaviors.state.isBusy, v => {
            this.setState({ isBusy: v, showAction: false })
        })
    }

    get isDisabled() {
        return this.state.isBusy || this.actionItems.length === 0
    }

    get actionItems() {
        const historyItems: ActionMenu.Items[] = []
        const { history } = this.plugin.managers.structure.focus
        if (history.length > 0) {
            historyItems.push([
                ActionMenu.Header('History'),
                ...ActionMenu.createItems(history, {
                    label: f => f.label,
                    category: f => f.category
                })
            ])
        }

        const presetItems: ActionMenu.Items[] = []
        const { structures } = this.plugin.managers.structure.hierarchy.selection;
        for (const s of structures) {
            const d = s.cell.obj?.data
            if (d) {
                const entries = getFocusEntries(d)
                if (entries.length > 0) {
                    presetItems.push([
                        ActionMenu.Header(d.label),
                        ...ActionMenu.createItems(entries, {
                            label: f => f.label,
                            category: f => f.category
                        })
                    ])
                }
            }
        }
        if (presetItems.length === 1) {
            const item = presetItems[0] as ActionMenu.Items[]
            const header = item[0] as ActionMenu.Header
            header.initiallyExpanded = true
        }

        const items: ActionMenu.Items[] = []
        if (historyItems.length > 0) items.push(...historyItems)
        if (presetItems.length > 0) items.push(...presetItems)

        return items
    }

    selectAction: ActionMenu.OnSelect = item => {
        if (!item || !this.state.showAction) {
            this.setState({ showAction: false });
            return;
        }
        this.setState({ showAction: false }, () => {
            const f = item.value as FocusEntry
            this.plugin.managers.structure.focus.set(f)
            this.plugin.managers.camera.focusLoci(f.loci, { durationMs: 0 })
        })
    }

    toggleAction = () => this.setState({ showAction: !this.state.showAction })

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

    render() {
        const { current } = this.plugin.managers.structure.focus
        const label = current?.label || 'Nothing Focused'

        return <>
            <div className='msp-control-row msp-select-row'>
                <button className='msp-btn msp-btn-block msp-no-overflow' onClick={this.focus} title='Click to Center Focused' onMouseEnter={this.highlightCurrent} onMouseLeave={this.clearHighlights} disabled={this.isDisabled || !current}>
                    {label}
                </button>
                <ToggleButton icon='target' title='Focus Target' toggle={this.toggleAction} isSelected={this.state.showAction} disabled={this.isDisabled} style={{ flex: '0 0 40px' }} />
            </div>
            {this.state.showAction && <ActionMenu items={this.actionItems} onSelect={this.selectAction} />}
        </>;
    }
}