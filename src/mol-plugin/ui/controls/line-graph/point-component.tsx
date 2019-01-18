/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Paul Luna <paulluna0215@gmail.com>
 */

import * as React from 'react';

import { Vec2 } from 'mol-math/linear-algebra';

export default class PointComponent extends React.Component<any, {show: boolean}> {
    constructor(props: any){
        super(props);
        this.state = {show: false}
        
        this.handleHover = this.handleHover.bind(this);
        this.handleHoverOff = this.handleHoverOff.bind(this);
    }

    private handleHover() {
        this.setState({show: true});
        const point = Vec2.create(this.props.nX, this.props.nY);
        this.props.onmouseover(point);
    }

    private handleHoverOff(){
        this.setState({show: false});
        this.props.onmouseover(undefined);
    }
    
    public render() {
        return([
            <circle 
                r="10"
                key={`${this.props.id}circle`}
                id={`${this.props.id}`}
                cx={this.props.x} 
                cy={this.props.y} 
                onClick={this.props.onclick}
                onMouseEnter={this.handleHover} 
                onMouseLeave={this.handleHoverOff}
                onMouseDown={this.props.onmousedown}
                fill="black"
            />
        ]);
    }
}