/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react'
import { WithStyles } from 'material-ui/styles';
import TextField from 'material-ui/TextField';
// import FileUpload from '@material-ui/icons/FileUpload';

import State from '../state'
import Observer from './observer';

export default class FileInput extends Observer<{ state: State } & WithStyles, { loading: boolean }> {
    state = { loading: false }

    componentDidMount() {
        this.subscribe(this.props.state.loading, value => {
           this.setState({ loading: value });
        });
    }

    render() {
        const { classes, state } = this.props;

        return <TextField
            label='PDB ID'
            className={classes.textField}
            disabled={this.state.loading}
            margin='normal'
            onChange={(event) => {
                state.pdbId = event.target.value
            }}
            onKeyPress={(event) => {
                if (event.key === 'Enter') state.loadPdbId()
            }}
        />
    }
}