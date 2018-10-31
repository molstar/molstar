/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react';
import { PluginContext } from '../context';

export class Controls extends React.Component<{ plugin: PluginContext }, { id: string }> {
    state = { id: '1grm' };

    private createState = () => {
        const url = `http://www.ebi.ac.uk/pdbe/static/entry/${this.state.id.toLowerCase()}_updated.cif`;
        // const url = `https://webchem.ncbr.muni.cz/CoordinateServer/${this.state.id.toLowerCase()}/full`
        this.props.plugin._test_createState(url);
    }

    render() {
        return <div>
            <input type='text' defaultValue={this.state.id} onChange={e => this.setState({ id: e.currentTarget.value })} />
            <button onClick={this.createState}>Create State</button><br/>
            <button onClick={() => this.props.plugin._test_centerView()}>Center View</button><br/>
            <button onClick={() => this.props.plugin._test_nextModel()}>Next Model</button><br/>
            <button onClick={() => this.props.plugin._test_playModels()}>Play Models</button><br/>
        </div>;
    }
}