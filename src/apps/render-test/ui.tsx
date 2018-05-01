/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react'
import { withStyles, WithStyles, Theme, StyleRulesCallback } from 'material-ui/styles';
import Typography from 'material-ui/Typography';
import Toolbar from 'material-ui/Toolbar';
import AppBar from 'material-ui/AppBar';
import Drawer from 'material-ui/Drawer';

import State from './state'

import Viewport from './components/viewport'
import FileInput from './components/file-input'
import ColorTheme from './components/color-theme'
import Detail from './components/detail'
import Visibility from './components/visibility'
import Assemblies from './components/assemblies'


const styles: StyleRulesCallback = (theme: Theme) => ({
    root: {
        flexGrow: 1,
        height: 830,
        zIndex: 1,
        overflow: 'hidden',
        position: 'relative',
        display: 'flex',
    },
    appBar: {
        zIndex: theme.zIndex.drawer + 1,
    },
    drawerPaper: {
        position: 'relative',
        width: 240,
    },
    content: {
        flexGrow: 1,
        backgroundColor: theme.palette.background.default,
        padding: theme.spacing.unit * 3,
        minWidth: 0, // So the Typography noWrap works
    },
    toolbar: theme.mixins.toolbar,
    formControl: {
        margin: theme.spacing.unit,
        width: 200,
    },
    textField: {
        marginLeft: theme.spacing.unit,
        marginRight: theme.spacing.unit,
        width: 200,
    },
    button: {
        margin: theme.spacing.unit,
    },
    input: {
        display: 'none',
    },
} as any);

const decorate = withStyles(styles);

interface Props {
    state: State;
};

class UI extends React.Component<{ state: State } & WithStyles, {  }> {
    render() {
        const { classes, state } = this.props;
        return (
            <div className={classes.root}>
                <AppBar position='absolute' className={classes.appBar}>
                    <Toolbar>
                    <Typography variant='title' color='inherit' noWrap>
                        Mol* Render Test
                    </Typography>
                    </Toolbar>
                </AppBar>
                <Drawer variant='permanent' classes={{ paper: classes.drawerPaper, }}>
                    <div className={classes.toolbar} />
                    <FileInput state={state} classes={classes}></FileInput>
                    <form className={classes.root} autoComplete='off'>
                        <div>
                            <Assemblies state={state} classes={classes}></Assemblies>
                            <ColorTheme state={state} classes={classes}></ColorTheme>
                            <Detail state={state} classes={classes}></Detail>
                            <Visibility state={state} classes={classes}></Visibility>
                        </div>
                    </form>
                </Drawer>
                <main className={classes.content}>
                    <div className={classes.toolbar} />
                    <Viewport state={state}></Viewport>
                </main>
            </div>
        );
    }
}

export default decorate<Props>(UI)