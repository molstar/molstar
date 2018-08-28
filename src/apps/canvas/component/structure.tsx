/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react'
import { StructureView } from '../view';

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

export interface StructureComponentProps {
    structureView: StructureView
}

export interface StructureComponentState {
    label: string
    assemblyId: string
    assemblyIds: { id: string, label: string }[]
    symmetryFeatureId: number
    symmetryFeatureIds: { id: number, label: string }[]
}

export class StructureComponent extends React.Component<StructureComponentProps, StructureComponentState> {
    state = {
        label: this.props.structureView.label,
        assemblyId: this.props.structureView.assemblyId,
        assemblyIds: this.props.structureView.getAssemblyIds(),
        symmetryFeatureId: this.props.structureView.symmetryFeatureId,
        symmetryFeatureIds: this.props.structureView.getSymmetryFeatureIds()
    }

    componentWillMount() {
        const sv = this.props.structureView

        this.setState({
            ...this.state,
            label: sv.label,
            assemblyId: sv.assemblyId,
            assemblyIds: sv.getAssemblyIds(),
            symmetryFeatureId: sv.symmetryFeatureId,
            symmetryFeatureIds: sv.getSymmetryFeatureIds()
        })
    }

    async update(state: Partial<StructureComponentState>) {
        const sv = this.props.structureView

        if (state.assemblyId !== undefined) await sv.setAssembly(state.assemblyId)
        if (state.symmetryFeatureId !== undefined) await sv.setSymmetryFeature(state.symmetryFeatureId)

        const newState = {
            ...this.state,
            label: sv.label,
            assemblyId: sv.assemblyId,
            assemblyIds: sv.getAssemblyIds(),
            symmetryFeatureId: sv.symmetryFeatureId,
            symmetryFeatureIds: sv.getSymmetryFeatureIds()
        }
        this.setState(newState)
    }

    render() {
        const { label, assemblyIds, symmetryFeatureIds } = this.state

        const assemblyIdOptions = assemblyIds.map(a => {
            return <option key={a.id} value={a.id}>{a.label}</option>
        })
        const symmetryFeatureIdOptions = symmetryFeatureIds.map(f => {
            return <option key={f.id} value={f.id}>{f.label}</option>
        })

        return <div>
            <div>
                <span>{label}</span>
            </div>
            <div>
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
            </div>
        </div>;
    }
}