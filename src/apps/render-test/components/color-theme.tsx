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

import State, { ColorTheme as _ColorTheme } from '../state'
import Observer from './observer';

interface ColorThemeState {
    loading: boolean
    name: _ColorTheme
}

export default class ColorTheme extends Observer<{ state: State } & WithStyles, ColorThemeState> {
    state = { loading: false, name: 'element-symbol' as _ColorTheme }

    componentDidMount() {
        this.subscribe(this.props.state.loading, value => {
           this.setState({ loading: value });
        });
        this.subscribe(this.props.state.colorTheme, value => {
            this.setState({ name: value });
         });
    }

    handleNameChange = (event: React.ChangeEvent<any>) => {
        this.props.state.colorTheme.next(event.target.value)
    }

    render() {
        const { classes } = this.props;

        const items = Object.keys(_ColorTheme).map((name, idx) => {
            return <MenuItem key={idx} value={name}>{name}</MenuItem>
        })

        return <FormControl className={classes.formControl}>
            <InputLabel htmlFor='color-theme-name'>Color Theme</InputLabel>
            <Select
                className={classes.selectField}
                value={this.state.name}
                onChange={this.handleNameChange}
                inputProps={{
                    name: 'name',
                    id: 'color-theme-name',
                }}
            >
                {items}
            </Select>
        </FormControl>
    }
}