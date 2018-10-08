/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react'
import { RepresentationComponent } from './representation';
import { Representation } from 'mol-geo/representation';
import { VolumeView } from '../volume-view';
import { VolumeRepresentation } from 'mol-geo/representation/volume';

export interface VolumeViewComponentProps {
    volumeView: VolumeView
}

export interface VolumeViewComponentState {
    volumeView: VolumeView
    label: string
    active: { [k: string]: boolean }
    volumeRepresentations: { [k: string]: VolumeRepresentation<any> }
}

export class VolumeViewComponent extends React.Component<VolumeViewComponentProps, VolumeViewComponentState> {
    state = this.stateFromVolumeView(this.props.volumeView)

    private stateFromVolumeView(vv: VolumeView) {
        return {
            volumeView: vv,
            label: vv.label,
            active: vv.active,
            volumeRepresentations: vv.volumeRepresentations
        }
    }

    componentWillMount() {
        this.setState(this.stateFromVolumeView(this.props.volumeView))
    }

    componentDidMount() {
        const vv = this.props.volumeView

        this.props.volumeView.updated.subscribe(() => this.setState({
            volumeRepresentations: vv.volumeRepresentations
        }))
    }

    componentWillReceiveProps(nextProps: VolumeViewComponentProps) {
        if (nextProps.volumeView !== this.props.volumeView) {
            this.setState(this.stateFromVolumeView(nextProps.volumeView))

            nextProps.volumeView.updated.subscribe(() => this.setState({
                volumeRepresentations: nextProps.volumeView.volumeRepresentations
            }))
        }
    }

    // async update(state: Partial<VolumeViewComponentState>) {
    //     const vv = this.state.volumeView
    //     this.setState(this.stateFromVolumeView(vv))
    // }

    render() {
        const { volumeView, label, active, volumeRepresentations } = this.state

        return <div>
            <div>
                <h2>{label}</h2>
            </div>
            <div>
                <div>
                    <h4>Active</h4>
                    { Object.keys(active).map((k, i) => {
                        return <div key={i}>
                            <input
                                type='checkbox'
                                checked={active[k]}
                                onChange={(e) => {
                                    volumeView.setVolumeRepresentation(k, e.target.checked)
                                }}
                            /> {k}
                        </div>
                    } ) }
                </div>
                <div>
                    <h3>Volume Representations</h3>
                    { Object.keys(volumeRepresentations).map((k, i) => {
                        if (active[k]) {
                            return <div key={i}>
                                <RepresentationComponent
                                    repr={volumeRepresentations[k] as Representation<any>}
                                    viewer={volumeView.viewer}
                                    app={volumeView.app}
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