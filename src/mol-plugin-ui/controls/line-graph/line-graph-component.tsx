/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Paul Luna <paulluna0215@gmail.com>
 */
import PointComponent from './point-component';

import * as React from 'react';
import { Vec2 } from '../../../mol-math/linear-algebra';

interface LineGraphComponentState {
    points: Vec2[],
    copyPoint: any,
    canSelectMultiple: boolean,
}

export default class LineGraphComponent extends React.Component<any, LineGraphComponentState> {
    private myRef: any;
    private height: number;
    private width: number;
    private padding: number;
    private updatedX: number;
    private updatedY: number;
    private selected?: number[];
    private ghostPoints: SVGElement[];
    private gElement: SVGElement;
    private namespace: string;

    constructor(props: any) {
        super(props);
        this.myRef = React.createRef();
        this.state = {
            points:[
                Vec2.create(0, 0),
                Vec2.create(1, 0)
            ],
            copyPoint: undefined,
            canSelectMultiple: false,
        };
        this.height = 400;
        this.width = 600;
        this.padding = 70;
        this.selected = undefined;
        this.ghostPoints = [];
        this.namespace = 'http://www.w3.org/2000/svg';

        for (const point of this.props.data){
            this.state.points.push(point);
        }

        this.state.points.sort((a, b) => {
            if(a[0] === b[0]){
                if(a[0] === 0){
                    return a[1] - b[1];
                }
                if(a[1] === 1){
                    return b[1] - a[1];
                }
                return a[1] - b[1];
            }
            return a[0] - b[0];
        });

        this.handleDrag = this.handleDrag.bind(this);
        this.handleMultipleDrag = this.handleMultipleDrag.bind(this);
        this.handleDoubleClick = this.handleDoubleClick.bind(this);
        this.refCallBack = this.refCallBack.bind(this);
        this.handlePointUpdate = this.handlePointUpdate.bind(this);
        this.change = this.change.bind(this);
        this.handleKeyUp = this.handleKeyUp.bind(this);
        this.handleLeave = this.handleLeave.bind(this);
        this.handleEnter = this.handleEnter.bind(this);

    }

    public render() {
        const points = this.renderPoints();
        const lines = this.renderLines();

        return ([
            <div key="LineGraph">
                <svg
                    className="msp-canvas"
                    ref={this.refCallBack}
                    viewBox={`0 0 ${this.width + this.padding} ${this.height + this.padding}`}
                    onMouseMove={this.handleDrag}
                    onMouseUp={this.handlePointUpdate}
                    onMouseLeave={this.handleLeave}
                    onMouseEnter={this.handleEnter}
                    tabIndex={0}
                    onKeyDown={this.handleKeyDown}
                    onKeyUp={this.handleKeyUp}
                    onDoubleClick={this.handleDoubleClick}>

                    <g stroke="black" fill="black">
                        {lines}
                        {points}
                    </g>
                    <g className="ghost-points" stroke="black" fill="black">
                    </g>
                </svg>
            </div>,
            <div key="modal" id="modal-root" />
        ]);
    }

    componentDidMount() {
        this.gElement = document.getElementsByClassName('ghost-points')[0] as SVGElement;
    }

    private change(points: Vec2[]){
        let copyPoints = points.slice();
        copyPoints.shift();
        copyPoints.pop();
        this.props.onChange(copyPoints);
    }

    private handleKeyDown = (event: any) => {
        // TODO: set canSelectMultiple = true
    }

    private handleKeyUp = (event: any) => {
        // TODO: SET canSelectMultiple = fasle
    }

    private handleClick = (id: number) => (event: any) => {
        // TODO: add point to selected array
    }

    private handleMouseDown = (id: number) => (event: any) => {
        if(id === 0 || id === this.state.points.length - 1){
            return;
        }

        if (this.state.canSelectMultiple) {
            return;
        }

        const copyPoint: Vec2 = this.normalizePoint(Vec2.create(this.state.points[id][0], this.state.points[id][1]));
        this.ghostPoints.push(document.createElementNS(this.namespace, 'circle') as SVGElement);
        this.ghostPoints[0].setAttribute('r', '10');
        this.ghostPoints[0].setAttribute('fill', 'orange');
        this.ghostPoints[0].setAttribute('cx', `${copyPoint[0]}`);
        this.ghostPoints[0].setAttribute('cy', `${copyPoint[1]}`);
        this.ghostPoints[0].setAttribute('style', 'display: none');
        this.gElement.appendChild(this.ghostPoints[0]);
        this.updatedX = copyPoint[0];
        this.updatedY = copyPoint[1];
        this.selected = [id];
    }

    private handleDrag(event: any) {
        if(this.selected === undefined){
            return;
        }

        const pt = this.myRef.createSVGPoint();
        let updatedCopyPoint;
        const padding = this.padding / 2;
        pt.x = event.clientX;
        pt.y = event.clientY;
        const svgP = pt.matrixTransform(this.myRef.getScreenCTM().inverse());
        updatedCopyPoint = Vec2.create(svgP.x, svgP.y);

        if ((svgP.x < (padding) || svgP.x > (this.width + (padding))) && (svgP.y > (this.height + (padding)) || svgP.y < (padding))) {
            updatedCopyPoint = Vec2.create(this.updatedX, this.updatedY);
        } else if (svgP.x < padding) {
            updatedCopyPoint = Vec2.create(padding, svgP.y);
        } else if( svgP.x > (this.width + (padding))) {
            updatedCopyPoint = Vec2.create(this.width + padding, svgP.y);
        } else if (svgP.y > (this.height + (padding))) {
            updatedCopyPoint = Vec2.create(svgP.x, this.height + padding);
        } else if (svgP.y < (padding)) {
            updatedCopyPoint = Vec2.create(svgP.x, padding);
        } else {
            updatedCopyPoint = Vec2.create(svgP.x, svgP.y);
        }

        this.updatedX = updatedCopyPoint[0];
        this.updatedY = updatedCopyPoint[1];
        const unNormalizePoint = this.unNormalizePoint(updatedCopyPoint);
        this.ghostPoints[0].setAttribute('style', 'display: visible');
        this.ghostPoints[0].setAttribute('cx', `${updatedCopyPoint[0]}`);
        this.ghostPoints[0].setAttribute('cy', `${updatedCopyPoint[1]}`);


        this.props.onDrag(unNormalizePoint);
    }

    private handleMultipleDrag() {
        // TODO
    }

    private handlePointUpdate(event: any) {
        const selected = this.selected;
        if (this.state.canSelectMultiple) {
            return;
        }

        if(selected === undefined || selected[0] === 0 || selected[0] === this.state.points.length - 1) {
            this.setState({
                copyPoint: undefined,
            });
            return;
        }
        this.selected = undefined;

        const updatedPoint = this.unNormalizePoint(Vec2.create(this.updatedX, this.updatedY));
        const points = this.state.points.filter((_, i) => i !== selected[0]);
        points.push(updatedPoint);;
        points.sort((a, b) => {
            if(a[0] === b[0]){
                if(a[0] === 0){
                    return a[1] - b[1];
                }
                if(a[1] === 1){
                    return b[1] - a[1];
                }
                return a[1] - b[1];
            }
            return a[0] - b[0];
        });
        this.setState({
            points,
        });
        this.change(points);
        this.gElement.innerHTML = '';
        this.ghostPoints = [];
        document.removeEventListener('mousemove', this.handleDrag, true);
        document.removeEventListener('mouseup', this.handlePointUpdate, true);
    }

    private handleDoubleClick(event: any) {
        let newPoint;
        const pt = this.myRef.createSVGPoint();
        pt.x = event.clientX;
        pt.y = event.clientY;
        const svgP = pt.matrixTransform(this.myRef.getScreenCTM().inverse());
        const points = this.state.points;
        const padding = this.padding / 2;

        if( svgP.x < (padding) ||
            svgP.x > (this.width + (padding)) ||
            svgP.y > (this.height + (padding)) ||
            svgP.y < (this.padding / 2)) {
            return;
        }
        newPoint = this.unNormalizePoint(Vec2.create(svgP.x, svgP.y));
        points.push(newPoint);
        points.sort((a, b) => {
            if(a[0] === b[0]){
                if(a[0] === 0){
                    return a[1] - b[1];
                }
                if(a[1] === 1){
                    return b[1] - a[1];
                }
                return a[1] - b[1];
            }
            return a[0] - b[0];
        });
        this.setState({points});
        this.change(points);
    }

    private deletePoint = (i: number) => (event: any) => {
        if(i === 0 || i === this.state.points.length - 1){ return; }
        const points = this.state.points.filter((_, j) => j !== i);
        points.sort((a, b) => {
            if(a[0] === b[0]){
                if(a[0] === 0){
                    return a[1] - b[1];
                }
                if(a[1] === 1){
                    return b[1] - a[1];
                }
                return a[1] - b[1];
            }
            return a[0] - b[0];
        });
        this.setState({points});
        this.change(points);
        event.stopPropagation();
    }

    private handleLeave() {
        if(this.selected === undefined) {
            return;
        }

        document.addEventListener('mousemove', this.handleDrag, true);
        document.addEventListener('mouseup', this.handlePointUpdate, true);
    }

    private handleEnter() {
        document.removeEventListener('mousemove', this.handleDrag, true);
        document.removeEventListener('mouseup', this.handlePointUpdate, true);
    }

    private normalizePoint(point: Vec2) {
        const min = this.padding / 2;
        const maxX = this.width + min;
        const maxY = this.height + min;
        const normalizedX = (point[0] * (maxX - min)) + min;
        const normalizedY = (point[1] * (maxY - min)) + min;
        const reverseY = (this.height + this.padding) - normalizedY;
        const newPoint = Vec2.create(normalizedX, reverseY);
        return newPoint;
    }

    private unNormalizePoint(point: Vec2) {
        const min = this.padding / 2;
        const maxX = this.width + min;
        const maxY = this.height + min;
        const unNormalizedX = (point[0] - min) / (maxX - min);

        // we have to take into account that we reversed y when we first normalized it.
        const unNormalizedY = ((this.height + this.padding) - point[1] - min) / (maxY - min);

        return Vec2.create(unNormalizedX, unNormalizedY);
    }

    private refCallBack(element: any) {
        if(element){
            this.myRef = element;
        }
    }

    private renderPoints() {
        const points: any[] = [];
        let point: Vec2;
        for (let i = 0; i < this.state.points.length; i++){
            if(i !== 0 && i !== this.state.points.length - 1){
                point = this.normalizePoint(this.state.points[i]);
                points.push(<PointComponent
                    key={i}
                    id={i}
                    x={point[0]}
                    y={point[1]}
                    nX={this.state.points[i][0]}
                    nY={this.state.points[i][1]}
                    selected={false}
                    delete={this.deletePoint}
                    onmouseover={this.props.onHover}
                    onmousedown={this.handleMouseDown(i)}
                    onclick={this.handleClick(i)}
                />);
            }
        }
        return points;
    }

    private renderLines() {
        const points: Vec2[] = [];
        let lines = [];
        let min: number;
        let maxX: number;
        let maxY: number;
        let normalizedX: number;
        let normalizedY: number;
        let reverseY: number;

        for(const point of this.state.points){
            min = this.padding / 2;
            maxX = this.width + min;
            maxY = this.height + min;
            normalizedX = (point[0] * (maxX - min)) + min;
            normalizedY = (point[1] * (maxY - min)) + min;
            reverseY = this.height + this.padding - normalizedY;
            points.push(Vec2.create(normalizedX, reverseY));
        }

        const data = points;
        const size = data.length;

        for (let i = 0; i < size - 1;i++){
            const x1 = data[i][0];
            const y1 = data[i][1];
            const x2 = data[i + 1][0];
            const y2 = data[i + 1][1];

            lines.push(<line key={`lineOf${i}`} x1={x1} x2={x2} y1={y1} y2={y2} stroke="#cec9ba" strokeWidth="5"/>);
        }

        return lines;
    }
}