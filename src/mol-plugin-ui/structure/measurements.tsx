/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react';
import { CollapsableControls, CollapsableState } from '../base';
import { lociLabel, dihedralLabel, angleLabel, distanceLabel } from '../../mol-theme/label';
import { Loci } from '../../mol-model/loci';
import { FiniteArray } from '../../mol-util/type-helpers';
import { StateObjectCell, StateTransform, StateTransformer } from '../../mol-state';
import { PluginStateObject } from '../../mol-plugin-state/objects';
import { ShapeRepresentation } from '../../mol-repr/shape/representation';
import { IconButton } from '../controls/common';
import { PluginCommands } from '../../mol-plugin/commands';

// TODO details, options (e.g. change text for labels)
// TODO better updates on state changes

type MeasurementTransform = StateObjectCell<PluginStateObject.Shape.Representation3D, StateTransform<StateTransformer<PluginStateObject.Molecule.Structure.Selections, PluginStateObject.Shape.Representation3D, any>>>


interface StructureMeasurementsControlsState extends CollapsableState {
    minRadius: number,
    extraRadius: number,
    durationMs: number,
    unitLabel: string,

    isDisabled: boolean,
}

export class StructureMeasurementsControls<P, S extends StructureMeasurementsControlsState> extends CollapsableControls<P, S> {
    componentDidMount() {
        this.subscribe(this.plugin.events.state.object.updated, ({ }) => {
            // TODO
            this.forceUpdate()
        })

        this.subscribe(this.plugin.events.state.object.created, ({ }) => {
            // TODO
            this.forceUpdate()
        })

        this.subscribe(this.plugin.events.state.object.removed, ({ }) => {
            // TODO
            this.forceUpdate()
        })

        this.subscribe(this.plugin.state.dataState.events.isUpdating, v => {
            this.setState({ isDisabled: v })
        })
    }

    focus(selections: PluginStateObject.Molecule.Structure.Selections) {
        return () => {
            const sphere = Loci.getBundleBoundingSphere(toLociBundle(selections.data))
            if (sphere) {
                const { extraRadius, minRadius, durationMs } = this.state
                const radius = Math.max(sphere.radius + extraRadius, minRadius);
                this.plugin.canvas3d?.camera.focus(sphere.center, radius, this.plugin.canvas3d.boundingSphere.radius, durationMs);
            }
        }
    }

    defaultState() {
        return {
            isCollapsed: false,
            header: 'Measurements',

            minRadius: 8,
            extraRadius: 4,
            durationMs: 250,
            unitLabel: '\u212B',

            isDisabled: false
        } as S
    }

    getLabel(selections: PluginStateObject.Molecule.Structure.Selections) {
        switch (selections.data.length) {
            case 1: return lociLabel(selections.data[0].loci, { condensed: true })
            case 2: return distanceLabel(toLociBundle(selections.data), { condensed: true })
            case 3: return angleLabel(toLociBundle(selections.data), { condensed: true })
            case 4: return dihedralLabel(toLociBundle(selections.data), { condensed: true })
        }
        return ''
    }

    highlight(repr: ShapeRepresentation<any, any, any>, selections: PluginStateObject.Molecule.Structure.Selections) {
        return () => {
            this.plugin.interactivity.lociHighlights.clearHighlights();
            for (const d of selections.data) {
                this.plugin.interactivity.lociHighlights.highlight({ loci: d.loci }, false);
            }
            this.plugin.interactivity.lociHighlights.highlight({ loci: repr.getLoci() }, false);
        }
    }

    clearHighlight() {
        return () => {
            this.plugin.interactivity.lociHighlights.clearHighlights();
        }
    }

    delete(cell: MeasurementTransform) {
        return () => {
            PluginCommands.State.RemoveObject(this.plugin, { state: cell.parent, ref: cell.transform.parent, removeParentGhosts: true });
        }
    }

    toggleVisibility(cell: MeasurementTransform) {
        return (e: React.MouseEvent<HTMLElement>) => {
            e.preventDefault();
            PluginCommands.State.ToggleVisibility(this.plugin, { state: cell.parent, ref: cell.transform.parent });
            e.currentTarget.blur();
        }
    }

    getRow(cell: MeasurementTransform) {
        const { obj } = cell
        if (!obj) return null

        const selections = obj.data.source as PluginStateObject.Molecule.Structure.Selections

        return <div className='msp-btn-row-group' key={obj.id} onMouseEnter={this.highlight(obj.data.repr, selections)} onMouseLeave={this.clearHighlight()}>
            <button className='msp-btn msp-btn-block msp-form-control' title='Click to focus. Hover to highlight.' onClick={this.focus(selections)}>
                <span dangerouslySetInnerHTML={{ __html: this.getLabel(selections) }} />
            </button>
            <IconButton isSmall={true} customClass='msp-form-control' onClick={this.delete(cell)} icon='remove' style={{ width: '52px' }} title='Delete' />
            <IconButton isSmall={true} customClass='msp-form-control' onClick={this.toggleVisibility(cell)} icon='eye' style={{ width: '52px' }} title={cell.state.isHidden ? 'Show' : 'Hide'} toggleState={cell.state.isHidden} />
        </div>
    }

    getData() {

    }

    renderControls() {
        const labels: JSX.Element[] = [];
        const distances: JSX.Element[] = [];
        const angles: JSX.Element[] = [];
        const dihedrals: JSX.Element[] = [];

        const measurements = this.plugin.helpers.measurement.getMeasurements();

        for (const d of measurements.labels) {
            const row = this.getRow(d)
            if (row) labels.push(row)
        }

        for (const d of measurements.distances) {
            const row = this.getRow(d)
            if (row) distances.push(row)
        }

        for (const d of measurements.angles) {
            const row = this.getRow(d)
            if (row) angles.push(row)
        }

        for (const d of measurements.dihedrals) {
            const row = this.getRow(d)
            if (row) dihedrals.push(row)
        }

        return <div>
            {labels.length > 0 && <div>
                <div className='msp-control-group-header' style={{ marginTop: '1px' }}><span>Labels</span></div>
                <div className='msp-control-offset'>{labels}</div>
            </div>}
            {distances.length > 0 && <div>
                <div className='msp-control-group-header' style={{ marginTop: '1px' }}><span>Distances</span></div>
                <div className='msp-control-offset'>{distances}</div>
            </div>}
            {angles.length > 0 && <div>
                <div className='msp-control-group-header' style={{ marginTop: '1px' }}><span>Angles</span></div>
                <div className='msp-control-offset'>{angles}</div>
            </div>}
            {dihedrals.length > 0 && <div>
                <div className='msp-control-group-header' style={{ marginTop: '1px' }}><span>Dihedrals</span></div>
                <div className='msp-control-offset'>{dihedrals}</div>
            </div>}
        </div>
    }
}

function toLociBundle(data: FiniteArray<{ loci: Loci }, any>): { loci: FiniteArray<Loci, any> } {
    return { loci: (data.map(d => d.loci) as unknown as FiniteArray<Loci, any>) }
}