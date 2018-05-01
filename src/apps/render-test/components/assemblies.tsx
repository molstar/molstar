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
import { Assembly } from 'mol-model/structure/model/properties/symmetry';

interface AssemblyState {
    loading: boolean
    assemblies: ReadonlyArray<Assembly>
    value: string
}

export default class Assemblies extends Observer<{ state: State } & WithStyles, AssemblyState> {
    state: AssemblyState = { loading: false, assemblies: [], value: '' }

    componentDidMount() {
        this.subscribe(this.props.state.loading, value => {
           this.setState({ loading: value });
        });
        this.subscribe(this.props.state.model, value => {
            this.setState({ assemblies: value ? value.symmetry.assemblies : [] });
        });
        this.subscribe(this.props.state.assembly, value => {
            this.setState({ value });
        });
    }

    handleValueChange = (event: React.ChangeEvent<any>) => {
        this.props.state.assembly.next(event.target.value)
    }

    render() {
        const { classes } = this.props;

        const items = this.state.assemblies.map((value, idx) => {
            return <MenuItem key={idx} value={value.id}>{value.details}</MenuItem>
        })

        return <FormControl className={classes.formControl}>
            <InputLabel htmlFor='assembly-value'>Assembly</InputLabel>
            <Select
                className={classes.selectField}
                value={this.state.value}
                onChange={this.handleValueChange}
                inputProps={{
                    name: 'value',
                    id: 'assembly-value',
                }}
            >
                {items}
            </Select>
        </FormControl>
    }
}