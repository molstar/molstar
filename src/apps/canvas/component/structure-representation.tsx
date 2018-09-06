/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react'
import { StructureRepresentation, StructureProps } from 'mol-geo/representation/structure';
import Viewer from 'mol-view/viewer';
import { VisualQuality } from 'mol-geo/representation/util';

export interface StructureRepresentationComponentProps {
    viewer: Viewer
    representation: StructureRepresentation<StructureProps>
}

export interface StructureRepresentationComponentState {
    label: string
    visible: boolean
    quality: VisualQuality
}

export class StructureRepresentationComponent extends React.Component<StructureRepresentationComponentProps, StructureRepresentationComponentState> {
    state = {
        label: this.props.representation.label,
        visible: this.props.representation.props.visible,
        quality: this.props.representation.props.quality,
    }

    componentWillMount() {
        const repr = this.props.representation

        this.setState({
            ...this.state,
            label: repr.label,
            visible: repr.props.visible,
            quality: repr.props.quality,
        })
    }

    async update(state: Partial<StructureRepresentationComponentState>) {
        const repr = this.props.representation
        const props: Partial<StructureProps> = {}

        if (state.visible !== undefined) props.visible = state.visible
        if (state.quality !== undefined) props.quality = state.quality

        await repr.createOrUpdate(props).run()
        this.props.viewer.add(repr)
        this.props.viewer.requestDraw()

        const newState = {
            ...this.state,
            visible: repr.props.visible,
            quality: repr.props.quality
        }
        this.setState(newState)
    }

    render() {
        const { label, visible, quality } = this.state

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
                <div>
                    <span>Quality</span>
                    <select value={quality} onChange={(e) => this.update({ quality: e.target.value as VisualQuality }) }>
                        <option value='auto'>auto</option>
                        <option value='lowest'>lowest</option>
                        <option value='low'>low</option>
                        <option value='medium'>medium</option>
                        <option value='high'>high</option>
                        <option value='highest'>highest</option>
                    </select>
                </div>
            </div>
        </div>;
    }
}