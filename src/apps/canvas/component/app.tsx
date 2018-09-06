/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react'
import { StructureView } from '../structure-view';
import { App } from '../app';
import { Viewport } from './viewport';
import { StructureViewComponent } from './structure-view';

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

export interface AppProps {
    app: App
}

export interface AppState {
    structureView: StructureView | null
}

export class AppComponent extends React.Component<AppProps, AppState> {
    state = {
        structureView: this.props.app.structureView,
    }

    componentDidMount() {
        this.props.app.pdbIdLoaded.subscribe((structureView) => {
            this.setState({
                structureView: this.props.app.structureView
            })
        })
    }

    render() {
        const { structureView } = this.state

        return <div style={{width: '100%', height: '100%'}}>
            <div style={{float: 'left', width: '70%', height: '100%'}}>
                <Viewport app={this.props.app} />
            </div>

            <div style={{float: 'right', width: '25%', height: '100%'}}>
                <div>
                    <span>Load PDB ID</span>
                    <input
                        type='text'
                        onKeyDown={e => {
                            if (e.keyCode === 13) {
                                const value = e.currentTarget.value.trim()
                                if (value) {
                                    this.props.app.loadPdbId(value)
                                }
                            }
                        }}
                    />
                </div>
                <div>
                    <span>Load CIF file</span>
                    <input
                        accept='*.cif'
                        type='file'
                        onChange={e => {
                            if (e.target.files) this.props.app.loadCifFile(e.target.files[0])
                        }}
                    />
                </div>
                <hr/>
                <div>
                    {structureView ? <StructureViewComponent structureView={structureView} /> : ''}
                </div>
            </div>
        </div>;
    }
}