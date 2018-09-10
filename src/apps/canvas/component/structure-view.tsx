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
    structureView: StructureView

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
    state = this.stateFromStructureView(this.props.structureView)

    private stateFromStructureView(sv: StructureView) {
        return {
            structureView: sv,

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
    }

    componentWillMount() {
        this.setState(this.stateFromStructureView(this.props.structureView))
    }

    componentDidMount() {
        const sv = this.props.structureView

        this.props.structureView.updated.subscribe(() => this.setState({
            symmetryFeatureIds: sv.getSymmetryFeatureIds(),
            structureRepresentations: sv.structureRepresentations
        }))
    }

    componentWillReceiveProps(nextProps: StructureViewComponentProps) {
        if (nextProps.structureView !== this.props.structureView) {
            this.setState(this.stateFromStructureView(nextProps.structureView))

            nextProps.structureView.updated.subscribe(() => this.setState({
                symmetryFeatureIds: nextProps.structureView.getSymmetryFeatureIds(),
                structureRepresentations: nextProps.structureView.structureRepresentations
            }))
        }
    }

    async update(state: Partial<StructureViewComponentState>) {
        const sv = this.state.structureView

        if (state.modelId !== undefined) await sv.setModel(state.modelId)
        if (state.assemblyId !== undefined) await sv.setAssembly(state.assemblyId)
        if (state.symmetryFeatureId !== undefined) await sv.setSymmetryFeature(state.symmetryFeatureId)

        this.setState(this.stateFromStructureView(sv))
    }

    render() {
        const { structureView, label, modelIds, assemblyIds, symmetryFeatureIds, active, structureRepresentations } = this.state

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
                    <span>Model </span>
                    <select
                        style={{width: '100px'}}
                        value={this.state.modelId}
                        onChange={(e) => {
                            this.update({ modelId: parseInt(e.target.value) })
                        }}
                    >
                        {modelIdOptions}
                    </select>
                    <span> </span>
                    <input type='range'
                        defaultValue={this.state.modelId.toString()}
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
                    <span>Assembly </span>
                    <select
                        style={{width: '150px'}}
                        value={this.state.assemblyId}
                        onChange={(e) => {
                            this.update({ assemblyId: e.target.value })
                        }}
                    >
                        {assemblyIdOptions}
                    </select>
                </div>
                <div>
                    <span>Symmetry Feature </span>
                    <select
                        style={{width: '150px'}}
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
                                    const sv = structureView
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
                                    viewer={structureView.viewer}
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