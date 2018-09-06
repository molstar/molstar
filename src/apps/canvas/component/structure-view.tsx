/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react'
import { StructureView } from '../structure-view';
import { StructureRepresentation } from 'mol-geo/representation/structure';
import { StructureRepresentationComponent } from './structure-representation';

// export function FileInput (props: {
//     accept: string
//     onChange: (v: FileList | null) => void,
// }) {
//     return <input
//         accept={props.accept || '*.*'}
//         type='file'
//         onChange={e => props.onChange.call(null, e.target.files)}
//     />
// }

export interface StructureViewComponentProps {
    structureView: StructureView
}

export interface StructureViewComponentState {
    label: string
    modelId: number
    modelIds: { id: number, label: string }[]
    assemblyId: string
    assemblyIds: { id: string, label: string }[]
    symmetryFeatureId: number
    symmetryFeatureIds: { id: number, label: string }[]

    active: { [k: string]: boolean }
    structureRepresentations: { [k: string]: StructureRepresentation<any> }
}

export class StructureViewComponent extends React.Component<StructureViewComponentProps, StructureViewComponentState> {
    state = {
        label: this.props.structureView.label,
        modelId: this.props.structureView.modelId,
        modelIds: this.props.structureView.getModelIds(),
        assemblyId: this.props.structureView.assemblyId,
        assemblyIds: this.props.structureView.getAssemblyIds(),
        symmetryFeatureId: this.props.structureView.symmetryFeatureId,
        symmetryFeatureIds: this.props.structureView.getSymmetryFeatureIds(),

        active: this.props.structureView.active,
        structureRepresentations: this.props.structureView.structureRepresentations
    }

    componentWillMount() {
        const sv = this.props.structureView

        this.setState({
            ...this.state,
            label: sv.label,
            modelId: sv.modelId,
            modelIds: sv.getModelIds(),
            assemblyId: sv.assemblyId,
            assemblyIds: sv.getAssemblyIds(),
            symmetryFeatureId: sv.symmetryFeatureId,
            symmetryFeatureIds: sv.getSymmetryFeatureIds(),

            active: sv.active,
            structureRepresentations: sv.structureRepresentations
        })
    }

    componentDidMount() {
        const sv = this.props.structureView

        this.props.structureView.updated.subscribe(() => this.setState({
            symmetryFeatureIds: sv.getSymmetryFeatureIds(),
            structureRepresentations: sv.structureRepresentations
        }))
    }

    async update(state: Partial<StructureViewComponentState>) {
        const sv = this.props.structureView

        console.log(state)

        if (state.modelId !== undefined) await sv.setModel(state.modelId)
        if (state.assemblyId !== undefined) await sv.setAssembly(state.assemblyId)
        if (state.symmetryFeatureId !== undefined) await sv.setSymmetryFeature(state.symmetryFeatureId)

        const newState = {
            ...this.state,
            label: sv.label,
            modelId: sv.modelId,
            modelIds: sv.getModelIds(),
            assemblyId: sv.assemblyId,
            assemblyIds: sv.getAssemblyIds(),
            symmetryFeatureId: sv.symmetryFeatureId,
            symmetryFeatureIds: sv.getSymmetryFeatureIds(),

            active: sv.active,
            structureRepresentations: sv.structureRepresentations
        }
        this.setState(newState)
    }

    render() {
        const { label, modelIds, assemblyIds, symmetryFeatureIds, active, structureRepresentations } = this.state

        const modelIdOptions = modelIds.map(m => {
            return <option key={m.id} value={m.id}>{m.label}</option>
        })
        const assemblyIdOptions = assemblyIds.map(a => {
            return <option key={a.id} value={a.id}>{a.label}</option>
        })
        const symmetryFeatureIdOptions = symmetryFeatureIds.map(f => {
            return <option key={f.id} value={f.id}>{f.label}</option>
        })

        return <div>
            <div>
                <h2>{label}</h2>
            </div>
            <div>
                <div>
                    <span>Model</span>
                    <select
                        value={this.state.modelId}
                        onChange={(e) => {
                            this.update({ modelId: parseInt(e.target.value) })
                        }}
                    >
                        {modelIdOptions}
                    </select>
                    <input type='range'
                        value={this.state.modelId}
                        min={Math.min(...modelIds.map(m => m.id))}
                        max={Math.max(...modelIds.map(m => m.id))}
                        step='1'
                        onInput={(e) => {
                            this.update({ modelId: parseInt(e.currentTarget.value) })
                        }}
                    >
                    </input>
                </div>
                <div>
                    <span>Assembly</span>
                    <select
                        value={this.state.assemblyId}
                        onChange={(e) => {
                            this.update({ assemblyId: e.target.value })
                        }}
                    >
                        {assemblyIdOptions}
                    </select>
                </div>
                <div>
                    <span>Symmetry Feature</span>
                    <select
                        value={this.state.symmetryFeatureId}
                        onChange={(e) => {
                            this.update({ symmetryFeatureId: parseInt(e.target.value) })
                        }}
                    >
                        {symmetryFeatureIdOptions}
                    </select>
                </div>
                <div>
                    <h4>Active</h4>
                    { Object.keys(active).map((k, i) => {
                        return <div key={i}>
                            <input
                                type='checkbox'
                                checked={active[k]}
                                onChange={(e) => {
                                    const sv = this.props.structureView
                                    if (k === 'symmetryAxes') {
                                        sv.setSymmetryAxes(e.target.checked)
                                    } else if (Object.keys(sv.structureRepresentations).includes(k)) {
                                        sv.setStructureRepresentation(k, e.target.checked)
                                    }
                                }}
                            /> {k}
                        </div>
                    } ) }
                </div>
                <div>
                    <h3>Structure Representations</h3>
                    { Object.keys(structureRepresentations).map((k, i) => {
                        if (active[k]) {
                            return <div key={i}>
                                <StructureRepresentationComponent
                                    representation={structureRepresentations[k]}
                                    viewer={this.props.structureView.viewer}
                                />
                            </div>
                        } else {
                            return ''
                        }
                    } ) }
                </div>
            </div>
        </div>;
    }
}