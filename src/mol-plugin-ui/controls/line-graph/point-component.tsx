
import * as React from 'react';

import { UUID } from '../../../mol-util';
import { Color } from '../../../mol-util/color';
import { ControlPointData } from './line-graph-component';
import { forwardRef, useState } from 'react';

// export class PointComponent extends React.Component<any, {show: boolean, id: UUID, index: number}> {
//     constructor(props: any) {
//         super(props);
//         this.state = { show: false, id: props.idm, index: props.index };

//         this.handleHover = this.handleHover.bind(this);
//         this.handleHoverOff = this.handleHoverOff.bind(this);
//     }

//     private handleHover() {
//         this.setState({ show: true });
//         const point: ControlPointData = {
//             x: this.props.nX,
//             alpha: this.props.nY
//         };
//         this.props.onmouseover(point);
//     }

//     private handleHoverOff() {
//         this.setState({ show: false });
//         this.props.onmouseover(undefined);
//     }

//     public render() {
//         const rgb = Color.toRgb(this.props.color);
//         const fill = `rgb(${rgb})`;
//         return ([
//             <circle
//                 r="10"
//                 key={`${this.props.id}circle`}
//                 id={`${this.props.id}`}
//                 cx={this.props.x}
//                 cy={this.props.y}
//                 onClick={this.props.onclick}
//                 // onDoubleClick={this.props.delete(this.props.id)}
//                 onMouseEnter={this.handleHover}
//                 onMouseLeave={this.handleHoverOff}
//                 onMouseDown={this.props.onmousedown}
//                 fill={fill}
//             >
//             </circle>
//         ]);
//     }
// }

interface PointComponentProps {
    index: any
    onmouseover: any
    nX: number
    nY: number
    color: Color
    onmousedown: any
    onclick: any
    x: number
    y: number
    id: UUID
    ref: any
}

export const PointComponent = forwardRef(_PointComponent);

function _PointComponent(props: PointComponentProps, ref: any) {
    const [show, setShow] = useState(false);
    const [id, setID] = useState(props.id);
    const [index, setIndex] = useState(props.index);

    function handleHoverOff() {
        setShow(false);
        // this.setState({ show: false });
        props.onmouseover(undefined);
    }



    function handleHover() {
        setShow(true);
        // this.setState({ show: true });
        const point: ControlPointData = {
            x: props.nX,
            alpha: props.nY
        };
        props.onmouseover(point);
    }


    function render() {
        const rgb = Color.toRgb(props.color);
        const fill = `rgb(${rgb})`;
        return ([
            <circle
                ref={ref}
                r="10"
                key={`${props.id}circle`}
                id={`${props.id}`}
                cx={props.x}
                cy={props.y}
                onClick={props.onclick}
                // onDoubleClick={props.delete(props.id)}
                onMouseEnter={handleHover}
                onMouseLeave={handleHoverOff}
                onMouseDown={props.onmousedown}
                fill={fill}
            >
            </circle>
        ]);
    }
    return render();
}