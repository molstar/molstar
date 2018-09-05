/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react'
import { StructureRepresentation, StructureProps } from 'mol-geo/representation/structure';
import Viewer from 'mol-view/viewer';

export interface StructureRepresentationComponentProps {
    viewer: Viewer
    representation: StructureRepresentation<StructureProps>
}

export interface StructureRepresentationComponentState {
    label: string
    visible: boolean
}

export class StructureRepresentationComponent extends React.Component<StructureRepresentationComponentProps, StructureRepresentationComponentState> {
    state = {
        label: this.props.representation.label,
        visible: this.props.representation.props.visible,
    }

    componentWillMount() {
        const repr = this.props.representation

        this.setState({
            ...this.state,
            label: repr.label,
            visible: repr.props.visible,
        })
    }

    async update(state: Partial<StructureRepresentationComponentState>) {
        const repr = this.props.representation

        if (state.visible !== undefined) {
            await repr.createOrUpdate({ visible: state.visible }).run()
            // this.props.viewer.add(repr)
            this.props.viewer.requestDraw()
        }

        const newState = {
            ...this.state,
            visible: repr.props.visible,
        }
        this.setState(newState)
    }

    render() {
        const { label, visible } = this.state

        return <div>
            <div>
                <h4>{label}</h4>
            </div>
            <div>
                <div>
                    <span>Visible</span>
                    <button onClick={(e) => this.update({ visible: !visible }) }>
                        {visible ? 'Hide' : 'Show'}
                    </button>
                </div>
            </div>
        </div>;
    }
}