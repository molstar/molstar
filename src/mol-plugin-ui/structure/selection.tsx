/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react';
import { CollapsableControls, CollapsableState } from '../base';
import { StructureSelectionQueries, SelectionModifier } from '../../mol-plugin/util/structure-selection-helper';
import { ButtonSelect, Options } from '../controls/common';
import { PluginCommands } from '../../mol-plugin/command';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { Interactivity } from '../../mol-plugin/util/interactivity';
import { ParameterControls } from '../controls/parameters';
import { stripTags } from '../../mol-util/string';
import { StructureElement } from '../../mol-model/structure';

const SSQ = StructureSelectionQueries
const DefaultQueries: (keyof typeof SSQ)[] = [
    'all', 'polymer', 'trace', 'backbone', 'protein', 'nucleic', 'water', 'branched', 'ligand', 'nonStandardPolymer',
    'surroundings', 'complement', 'bonded'
]

const StructureSelectionParams = {
    granularity: Interactivity.Params.granularity,
}

interface StructureSelectionControlsState extends CollapsableState {
    minRadius: number,
    extraRadius: number,
    durationMs: number,

    isDisabled: boolean
}

export class StructureSelectionControls<P, S extends StructureSelectionControlsState> extends CollapsableControls<P, S> {
    componentDidMount() {
        this.subscribe(this.plugin.events.interactivity.selectionUpdated, () => {
            this.forceUpdate()
        });

        this.subscribe(this.plugin.events.interactivity.propsUpdated, () => {
            this.forceUpdate()
        });

        this.subscribe(this.plugin.state.dataState.events.isUpdating, v => this.setState({ isDisabled: v }))
    }

    get stats() {
        const stats = this.plugin.helpers.structureSelectionManager.stats
        if (stats.structureCount === 0 || stats.elementCount === 0) {
            return 'Selected nothing'
        } else {
            return `Selected ${stripTags(stats.label)}`
        }
    }

    focus = () => {
        const { extraRadius, minRadius, durationMs } = this.state
        if (this.plugin.helpers.structureSelectionManager.stats.elementCount === 0) return
        const principalAxes = this.plugin.helpers.structureSelectionManager.getPrincipalAxes();
        const { origin, dirA, dirC } = principalAxes.boxAxes
        const { sphere } = this.plugin.helpers.structureSelectionManager.getBoundary()
        const radius = Math.max(sphere.radius + extraRadius, minRadius);
        this.plugin.canvas3d?.camera.focus(origin, radius, durationMs, dirA, dirC);
    }

    focusLoci(loci: StructureElement.Loci) {
        return () => {
            const { extraRadius, minRadius, durationMs } = this.state
            if (this.plugin.helpers.structureSelectionManager.stats.elementCount === 0) return
            const { sphere } = StructureElement.Loci.getBoundary(loci)
            const radius = Math.max(sphere.radius + extraRadius, minRadius);
            this.plugin.canvas3d?.camera.focus(sphere.center, radius, durationMs);
        }
    }

    measureDistance = () => {
        const loci = this.plugin.helpers.structureSelectionManager.latestLoci;
        this.plugin.helpers.measurement.addDistance(loci[0].loci, loci[1].loci);
    }

    measureAngle = () => {
        const loci = this.plugin.helpers.structureSelectionManager.latestLoci;
        this.plugin.helpers.measurement.addAngle(loci[0].loci, loci[1].loci, loci[2].loci);
    }

    measureDihedral = () => {
        const loci = this.plugin.helpers.structureSelectionManager.latestLoci;
        this.plugin.helpers.measurement.addDihedral(loci[0].loci, loci[1].loci, loci[2].loci, loci[3].loci);
    }

    addLabel = () => {
        const loci = this.plugin.helpers.structureSelectionManager.latestLoci;
        this.plugin.helpers.measurement.addLabel(loci[0].loci);
    }

    addOrientation = () => {
        const loci = this.plugin.helpers.structureSelectionManager.latestLoci;
        this.plugin.helpers.measurement.addOrientation(loci[0].loci);
    }

    setProps = (p: { param: PD.Base<any>, name: string, value: any }) => {
        if (p.name === 'granularity') {
            PluginCommands.Interactivity.SetProps.dispatch(this.plugin, { props: { granularity: p.value } });
        }
    }

    get values () {
        return {
            granularity: this.plugin.interactivity.props.granularity,
        }
    }

    set = (modifier: SelectionModifier, value: string) => {
        const query = SSQ[value as keyof typeof SSQ]
        this.plugin.helpers.structureSelection.set(modifier, query.query, false)
    }

    add = (value: string) => this.set('add', value)
    remove = (value: string) => this.set('remove', value)
    only = (value: string) => this.set('only', value)

    queries = Options(Object.keys(StructureSelectionQueries)
            .map(name => [name, SSQ[name as keyof typeof SSQ].label] as [string, string])
            .filter(pair => DefaultQueries.includes(pair[0] as keyof typeof SSQ)));

    controls = <div className='msp-control-row'>
        <div className='msp-select-row'>
            <ButtonSelect label='Select' onChange={this.add} disabled={this.state.isDisabled}>
                <optgroup label='Select'>
                    {this.queries}
                </optgroup>
            </ButtonSelect>
            <ButtonSelect label='Deselect' onChange={this.remove} disabled={this.state.isDisabled}>
                <optgroup label='Deselect'>
                    {this.queries}
                </optgroup>
            </ButtonSelect>
            <ButtonSelect label='Only' onChange={this.only} disabled={this.state.isDisabled}>
                <optgroup label='Only'>
                    {this.queries}
                </optgroup>
            </ButtonSelect>
        </div>
    </div>

    defaultState() {
        return {
            isCollapsed: false,
            header: 'Selection',

            minRadius: 8,
            extraRadius: 4,
            durationMs: 250,

            isDisabled: false
        } as S
    }

    renderControls() {
        const latest: JSX.Element[] = [];

        const mng = this.plugin.helpers.structureSelectionManager;

        // TODO: fix the styles, move them to CSS

        for (let i = 0, _i = Math.min(4, mng.latestLoci.length); i < _i; i++) {
            const e = mng.latestLoci[i];
            latest.push(<li key={e!.label}>
                <button className='msp-btn msp-btn-block msp-form-control' style={{ borderRight: '6px solid transparent', overflow: 'hidden' }}
                    title='Click to focus.' onClick={this.focusLoci(e.loci)}>
                    <span dangerouslySetInnerHTML={{ __html: e.label }} />
                </button>
                {/* <div>
                    <IconButton icon='remove' title='Remove' onClick={() => {}} />
                </div> */}
            </li>)
        }

        return <div>
            <div className='msp-control-row msp-row-text'>
                <button className='msp-btn msp-btn-block' onClick={this.focus}>
                    <span className={`msp-icon msp-icon-focus-on-visual`} style={{ position: 'absolute', left: '5px' }} />
                    {this.stats}
                </button>
            </div>
            <ParameterControls params={StructureSelectionParams} values={this.values} onChange={this.setProps} isDisabled={this.state.isDisabled} />
            {this.controls}
            { latest.length > 0 &&
            <>
                <div className='msp-control-group-header' style={{ marginTop: '1px' }}><span>Latest Selections &amp; Measurement</span></div>
                <ul style={{ listStyle: 'none', marginTop: '1px', marginBottom: '0' }} className='msp-state-list'>
                    {latest}
                </ul>
                {latest.length >= 1 && <div className='msp-control-row msp-row-text'>
                    <button className='msp-btn msp-btn-block' onClick={this.addLabel} title='Add label for latest selection'>
                        Add Label
                    </button>
                </div>}
                {latest.length >= 1 && <div className='msp-control-row msp-row-text'>
                    <button className='msp-btn msp-btn-block' onClick={this.addOrientation} title='Add orientation box/axes for latest selection'>
                        Add Orientation
                    </button>
                </div>}
                {latest.length >= 2 && <div className='msp-control-row msp-row-text'>
                    <button className='msp-btn msp-btn-block' onClick={this.measureDistance} title='Measure distance between latest 2 selections'>
                        Measure Distance
                    </button>
                </div>}
                {latest.length >= 3 && <div className='msp-control-row msp-row-text'>
                    <button className='msp-btn msp-btn-block' onClick={this.measureAngle} title='Measure angle between latest 3 selections'>
                        Measure Angle
                    </button>
                </div>}
                {latest.length >= 4 && <div className='msp-control-row msp-row-text'>
                    <button className='msp-btn msp-btn-block' onClick={this.measureDihedral} title='Measure dihedral between latest 4 selections'>
                        Measure Dihedral
                    </button>
                </div>}
            </>}
        </div>
    }
}