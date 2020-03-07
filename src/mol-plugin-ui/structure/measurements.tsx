/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react';
import { CollapsableControls, CollapsableState } from '../base';
import { lociLabel } from '../../mol-theme/label';
import { StructureElement } from '../../mol-model/structure';

// TODO hide/show, delete, details, options (e.g. change text for labels)
// TODO better labels: shorter, include measure
// TODO better updates on state changes

interface StructureMeasurementsControlsState extends CollapsableState {
    minRadius: number,
    extraRadius: number,
    durationMs: number,

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

    focusLoci(loci: StructureElement.Loci) {
        return () => {
            const { extraRadius, minRadius, durationMs } = this.state
            if (this.plugin.helpers.structureSelectionManager.stats.elementCount === 0) return
            const { sphere } = StructureElement.Loci.getBoundary(loci)
            const radius = Math.max(sphere.radius + extraRadius, minRadius);
            this.plugin.canvas3d?.camera.focus(sphere.center, radius, this.plugin.canvas3d.boundingSphere.radius, durationMs);
        }
    }

    defaultState() {
        return {
            isCollapsed: false,
            header: 'Measurements',

            minRadius: 8,
            extraRadius: 4,
            durationMs: 250,

            isDisabled: false
        } as S
    }

    renderControls() {
        const labels: JSX.Element[] = [];
        const distances: JSX.Element[] = [];
        const angles: JSX.Element[] = [];
        const dihedrals: JSX.Element[] = [];

        const measurements = this.plugin.helpers.measurement.getMeasurements();

        for (const d of measurements.labels) {
            const source = d.obj?.data.source
            if (source) {
                const lA = simpleLabel(source.data[0].loci)
                labels.push(<li key={source.id}>
                    <button className='msp-btn msp-btn-block msp-form-control' style={{ borderRight: '6px solid transparent', overflow: 'hidden' }}
                        title='Click to focus.' onClick={this.focusLoci(source.data[0].loci)}>
                        <span dangerouslySetInnerHTML={{ __html: `${lA}` }} />
                    </button>
                </li>)
            }
        }

        for (const d of measurements.distances) {
            const source = d.obj?.data.source
            if (source) {
                const lA = simpleLabel(source.data[0].loci)
                const lB = simpleLabel(source.data[1].loci)
                distances.push(<li key={source.id}>
                    <button className='msp-btn msp-btn-block msp-form-control' style={{ borderRight: '6px solid transparent', overflow: 'hidden' }}
                        title='Click to focus.' onClick={this.focusLoci(source.data[0].loci)}>
                        <span dangerouslySetInnerHTML={{ __html: `${lA} \u2014 ${lB}` }} />
                    </button>
                </li>)
            }
        }

        for (const d of measurements.angles) {
            const source = d.obj?.data.source
            if (source) {
                const lA = simpleLabel(source.data[0].loci)
                const lB = simpleLabel(source.data[1].loci)
                const lC = simpleLabel(source.data[2].loci)
                angles.push(<li key={source.id}>
                    <button className='msp-btn msp-btn-block msp-form-control' style={{ borderRight: '6px solid transparent', overflow: 'hidden' }}
                        title='Click to focus.' onClick={this.focusLoci(source.data[0].loci)}>
                        <span dangerouslySetInnerHTML={{ __html: `${lA} \u2014 ${lB} \u2014 ${lC}` }} />
                    </button>
                </li>)
            }
        }

        for (const d of measurements.dihedrals) {
            const source = d.obj?.data.source
            if (source) {
                const lA = simpleLabel(source.data[0].loci)
                const lB = simpleLabel(source.data[1].loci)
                const lC = simpleLabel(source.data[2].loci)
                const lD = simpleLabel(source.data[3].loci)
                dihedrals.push(<li key={source.id}>
                    <button className='msp-btn msp-btn-block msp-form-control' style={{ borderRight: '6px solid transparent', overflow: 'hidden' }}
                        title='Click to focus.' onClick={this.focusLoci(source.data[0].loci)}>
                        <span dangerouslySetInnerHTML={{ __html: `${lA} \u2014 ${lB} \u2014 ${lC} \u2014 ${lD}` }} />
                    </button>
                </li>)
            }
        }

        return <div>
            {labels.length > 0 && <div>
                <div className='msp-control-group-header' style={{ marginTop: '1px' }}><span>Labels</span></div>
                <ul style={{ listStyle: 'none', marginTop: '1px', marginBottom: '0' }} className='msp-state-list'>
                    {labels}
                </ul>
            </div>}
            {distances.length > 0 && <div>
                <div className='msp-control-group-header' style={{ marginTop: '1px' }}><span>Distances</span></div>
                <ul style={{ listStyle: 'none', marginTop: '1px', marginBottom: '0' }} className='msp-state-list'>
                    {distances}
                </ul>
            </div>}
            {angles.length > 0 && <div>
                <div className='msp-control-group-header' style={{ marginTop: '1px' }}><span>Angles</span></div>
                <ul style={{ listStyle: 'none', marginTop: '1px', marginBottom: '0' }} className='msp-state-list'>
                    {angles}
                </ul>
            </div>}
            {dihedrals.length > 0 && <div>
                <div className='msp-control-group-header' style={{ marginTop: '1px' }}><span>Dihedrals</span></div>
                <ul style={{ listStyle: 'none', marginTop: '1px', marginBottom: '0' }} className='msp-state-list'>
                    {dihedrals}
                </ul>
            </div>}
        </div>
    }
}

function simpleLabel(loci: StructureElement.Loci) {
    return lociLabel(loci, { htmlStyling: false })
        .split('|')
        .reverse()[0]
        .replace(/\[.*\]/g, '')
        .trim()
}