/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react';
import { PluginUIComponent } from '@molstar/plugin-ui/base';
import { Button, type ColorAccent } from '@molstar/plugin-ui/controls/common';
import { Icon, ArrowRightSvg, ArrowDropDownSvg } from '@molstar/plugin-ui/controls/icons';

export type CollapsableProps = { initiallyCollapsed?: boolean, header?: string }
export type CollapsableState = {
    isCollapsed: boolean,
    header: string,
    description?: string,
    isHidden?: boolean,
    brand?: { svg?: React.FC, accent: ColorAccent }
}

export abstract class CollapsableControls<P = {}, S = {}, SS = {}> extends PluginUIComponent<P & CollapsableProps, S & CollapsableState, SS> {
    toggleCollapsed() {
        this.setState({ isCollapsed: !this.state.isCollapsed } as (S & CollapsableState));
    };

    componentDidUpdate(prevProps: P & CollapsableProps) {
        if (this.props.initiallyCollapsed !== undefined && prevProps.initiallyCollapsed !== this.props.initiallyCollapsed) {
            this.setState({ isCollapsed: this.props.initiallyCollapsed as any });
        }
    }

    protected abstract defaultState(): (S & CollapsableState)
    protected abstract renderControls(): JSX.Element | null

    render() {
        if (this.state.isHidden) return null;
        const divid = this.state.header.toLowerCase().replace(/\s/g, '');
        const wrapClass = this.state.isCollapsed
            ? 'msp-transform-wrapper msp-transform-wrapper-collapsed'
            : 'msp-transform-wrapper';

        return <div className={wrapClass}>
            <div id={divid} className='msp-transform-header'>
                <Button icon={this.state.brand ? void 0 : this.state.isCollapsed ? ArrowRightSvg : ArrowDropDownSvg} noOverflow onClick={() => this.toggleCollapsed()}
                    className={this.state.brand ? `msp-transform-header-brand msp-transform-header-brand-${this.state.brand.accent}` : void 0} title={`Click to ${this.state.isCollapsed ? 'expand' : 'collapse'}`}>
                    <Icon svg={this.state.brand?.svg} inline />
                    {this.state.header}
                    <small style={{ margin: '0 6px' }}>{this.state.isCollapsed ? '' : this.state.description}</small>
                </Button>
            </div>
            {!this.state.isCollapsed && this.renderControls()}
        </div>;
    }

    constructor(props: P & CollapsableProps, context?: any) {
        super(props, context);

        const state = this.defaultState();
        if (props.initiallyCollapsed !== undefined) state.isCollapsed = props.initiallyCollapsed;
        if (props.header !== undefined) state.header = props.header;
        this.state = state;
    }
}
