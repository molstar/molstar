/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react'
import { StructureView } from '../structure-view';
import { App } from '../app';
import { Viewport } from './viewport';
import { StructureComponent } from './structure';

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
        this.props.app.pdbIdLoaded.subscribe(() => this.setState({
            structureView: this.props.app.structureView
        }))
    }

    render() {
        const { structureView } = this.state

        return <div style={{width: '100%', height: '100%'}}>
            <div style={{float: 'left', width: '70%', height: '100%'}}>
                <Viewport app={this.props.app} />
            </div>

            <div style={{float: 'right', width: '25%', height: '100%'}}>
                {structureView ? <StructureComponent structureView={structureView} /> : ''}
            </div>
        </div>;
    }
}