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

import Viewport from './components/viewport'
import FileInput from './components/file-input'
import State from './state'

const styles: StyleRulesCallback<any> = (theme: Theme) => ({
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
    textField: {
        marginLeft: theme.spacing.unit,
        marginRight: theme.spacing.unit,
        width: 200,
    },
});

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