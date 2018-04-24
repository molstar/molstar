/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react'
import { WithStyles } from 'material-ui/styles';
import { MenuItem } from 'material-ui/Menu';
import { InputLabel } from 'material-ui/Input';
import { FormControl } from 'material-ui/Form';
import Select from 'material-ui/Select';

import State from '../state'
import Observer from './observer';

interface DetailState {
    loading: boolean
    value: number
}

export default class Detail extends Observer<{ state: State } & WithStyles, DetailState> {
    state = { loading: false, value: 2 }

    componentDidMount() {
        this.subscribe(this.props.state.loading, value => {
           this.setState({ loading: value });
        });
        this.subscribe(this.props.state.detail, value => {
            this.setState({ value });
         });
    }

    handleValueChange = (event: React.ChangeEvent<any>) => {
        this.props.state.detail.next(event.target.value)
    }

    render() {
        const { classes } = this.props;

        const items = [0, 1, 2].map((value, idx) => {
            return <MenuItem key={idx} value={value}>{value.toString()}</MenuItem>
        })

        return <FormControl className={classes.formControl}>
            <InputLabel htmlFor='detail-value'>Detail</InputLabel>
            <Select
                className={classes.selectField}
                value={this.state.value}
                onChange={this.handleValueChange}
                inputProps={{
                    name: 'value',
                    id: 'detail-value',
                }}
            >
                {items}
            </Select>
        </FormControl>
    }
}