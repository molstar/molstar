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
import { Examples } from '../examples';
import { VolumeViewComponent } from './volume-view';
import { VolumeView } from '../volume-view';

export interface AppProps {
    app: App
}

export interface AppState {
    structureView: StructureView | null,
    volumeView: VolumeView | null,
    mmcifBinary: boolean,
    volcifBinary: boolean
}

export class AppComponent extends React.Component<AppProps, AppState> {
    state = {
        structureView: this.props.app.structureView,
        volumeView: this.props.app.volumeView,
        mmcifBinary: false,
        volcifBinary: true
    }

    componentDidMount() {
        this.props.app.structureLoaded.subscribe((structureView) => {
            this.setState({ structureView: this.props.app.structureView })
        })
        this.props.app.volumeLoaded.subscribe((volumeView) => {
            this.setState({ volumeView: this.props.app.volumeView })
        })
    }

    render() {
        const { structureView, volumeView } = this.state

        return <div style={{width: '100%', height: '100%'}}>
            <div style={{left: '0px', right: '350px', height: '100%', position: 'absolute'}}>
                <Viewport app={this.props.app} />
            </div>

            <div style={{width: '330px', paddingLeft: '10px', paddingRight: '10px', right: '0px', height: '100%', position: 'absolute', overflow: 'auto'}}>
                <div style={{marginTop: '10px'}}>
                    <span>Load PDB ID or URL</span>&nbsp;&nbsp;
                    <input type='checkbox' checked={this.state.mmcifBinary} onChange={e => this.setState({ mmcifBinary: e.target.checked })} /> Binary<br />
                    <input
                        style={{ width: '100%' }}
                        type='text'
                        onKeyDown={e => {
                            if (e.keyCode === 13) {
                                const value = e.currentTarget.value.trim()
                                if (value) {
                                    this.props.app.loadPdbIdOrMmcifUrl(value, { binary: this.state.mmcifBinary })
                                }
                            }
                        }}
                    />
                </div>
                <div>
                    <span>Load CIF file </span>
                    <input
                        accept='*.cif'
                        type='file'
                        onChange={e => {
                            if (e.target.files) this.props.app.loadMmcifFile(e.target.files[0])
                        }}
                    />
                </div>
                <div>
                    <span>Load CCP4/MRC file </span>
                    <input
                        accept='*.ccp4,*.mrc, *.map'
                        type='file'
                        onChange={e => {
                            if (e.target.files) this.props.app.loadCcp4File(e.target.files[0])
                        }}
                    />
                </div>
                <div style={{marginTop: '10px'}}>
                    <span>Load DensityServer URL</span>&nbsp;&nbsp;
                    <input type='checkbox' checked={this.state.volcifBinary} onChange={e => this.setState({ volcifBinary: e.target.checked })} /> Binary<br />
                    <input
                        style={{ width: '100%' }}
                        type='text'
                        onKeyDown={e => {
                            if (e.keyCode === 13) {
                                const value = e.currentTarget.value.trim()
                                if (value) {
                                    this.props.app.loadVolcifUrl(value, this.state.volcifBinary)
                                }
                            }
                        }}
                    />
                </div>
                <div>
                    <span>Load example </span>
                    <select
                        style={{width: '200px'}}
                        onChange={e => {
                            this.props.app.loadPdbIdOrMmcifUrl(e.target.value)
                        }}
                    >
                        <option value=''></option>
                        {Examples.map(({label, id, description}, i) => {
                            return <option key={i} value={id}>{`${label ? label : id} - ${description}`}</option>
                        })}
                    </select>
                </div>
                <hr/>
                <div style={{marginBottom: '10px'}}>
                    {structureView ? <StructureViewComponent structureView={structureView} /> : ''}
                </div>
                <hr/>
                <div style={{marginBottom: '10px'}}>
                    {volumeView ? <VolumeViewComponent volumeView={volumeView} /> : ''}
                </div>
            </div>
        </div>;
    }
}