/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react'
import { WithStyles } from 'material-ui/styles';
import TextField from 'material-ui/TextField';
import Button from 'material-ui/Button';

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

        return <div>
            <TextField
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
            <input
                accept='*.cif'
                className={classes.input}
                id='button-file'
                type='file'
                onChange={(event) => {
                    if (event.target.files) {
                        state.loadFile(event.target.files[0])
                    }
                }}
            />
            <label htmlFor='button-file'>
                <Button
                    variant='raised'
                    component='span'
                    className={classes.button}
                >
                    Open CIF
                </Button>
            </label>
        </div>
    }
}