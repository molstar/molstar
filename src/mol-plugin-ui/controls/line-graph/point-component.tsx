
import * as React from 'react';

import { UUID } from '../../../mol-util';
import { Color } from '../../../mol-util/color';
import { ControlPointData } from './line-graph-component';
import { forwardRef, useState } from 'react';

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