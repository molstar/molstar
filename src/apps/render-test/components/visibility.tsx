/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react'
import { WithStyles } from 'material-ui/styles';
import { FormLabel, FormControl, FormGroup, FormControlLabel } from 'material-ui/Form';
import Checkbox from 'material-ui/Checkbox';

import State from '../state'
import Observer from './observer';

interface VisibilityState {
    loading: boolean
    spacefill: boolean
    point: boolean
}

export default class Visibility extends Observer<{ state: State } & WithStyles, VisibilityState> {
    state = { loading: false, spacefill: true, point: true }

    componentDidMount() {
        this.subscribe(this.props.state.loading, value => {
           this.setState({ loading: value });
        });
        this.subscribe(this.props.state.spacefillVisibility, value => {
            this.setState({ spacefill: value });
        });
        this.subscribe(this.props.state.pointVisibility, value => {
            this.setState({ point: value });
        });
    }

    handleChange = (event: React.ChangeEvent<any>) => {
        switch (event.target.name) {
            case 'point': this.props.state.pointVisibility.next(event.target.checked); break;
            case 'spacefill': this.props.state.spacefillVisibility.next(event.target.checked); break;
        }
    }

    render() {
        const { classes } = this.props

        return <div className={classes.formControl}>
            <FormControl component='fieldset'>
                <FormLabel component='legend'>Visibility</FormLabel>
                <FormGroup>
                    <FormControlLabel
                        control={<Checkbox
                            checked={this.state.point}
                            onChange={this.handleChange}
                            name='point'
                        />}
                        label='Point'
                    />
                    <FormControlLabel
                        control={<Checkbox
                            checked={this.state.spacefill}
                            onChange={this.handleChange}
                            name='spacefill'
                        />}
                        label='Spacefill'
                    />
                </FormGroup>
            </FormControl>
        </div>
    }
}