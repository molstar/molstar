
import * as React from 'react';

import { UUID } from '../../../mol-util';
import { Color } from '../../../mol-util/color';
import { ControlPointData } from './line-graph-component';

export class PointComponent extends React.Component<any, {show: boolean, id: UUID, index: number}> {
    constructor(props: any) {
        super(props);
        this.state = { show: false, id: props.idm, index: props.index };

        this.handleHover = this.handleHover.bind(this);
        this.handleHoverOff = this.handleHoverOff.bind(this);
        this.deletePoint = this.deletePoint.bind(this);
    }

    private handleHover() {
        this.setState({ show: true });
        const point: ControlPointData = {
            x: this.props.nX,
            alpha: this.props.nY
        };
        this.props.onmouseover(point);
    }

    private handleHoverOff() {
        this.setState({ show: false });
        this.props.onmouseover(undefined);
    }

    private deletePoint() {
        this.props.delete(this.props.id);
    }

    public render() {
        const rgb = Color.toRgb(this.props.color);
        const fill = `rgb(${rgb})`;
        return ([
            <circle
                r="10"
                key={`${this.props.id}circle`}
                id={`${this.props.id}`}
                cx={this.props.x}
                cy={this.props.y}
                onClick={this.props.onclick}
                onDoubleClick={this.props.delete(this.props.id)}
                onMouseEnter={this.handleHover}
                onMouseLeave={this.handleHoverOff}
                onMouseDown={this.props.onmousedown}
                fill={fill}
            >
            </circle>
        ]);
    }
}