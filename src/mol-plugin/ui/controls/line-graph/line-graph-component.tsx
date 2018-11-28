/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Paul Luna <paulluna0215@gmail.com>
 */
import PointComponent from './point-component';

import * as React from 'react';
import { Vec2 } from 'mol-math/linear-algebra';

interface LineGraphComponentState {
    points: Vec2[],
    selected?: number,
    copyPoint: any,
    updatedX: number, 
    updatedY: number,
}

export default class LineGraphComponent extends React.Component<any, LineGraphComponentState> {
    private myRef:any;
    private height: number;
    private width: number;
    private padding: number;
    
    constructor(props: any) {
        super(props);
        this.myRef = React.createRef();
        this.state = {
            points:[
                Vec2.create(0, 0),
                Vec2.create(1, 0)
            ],
            selected: undefined,
            copyPoint: undefined,
            updatedX: 0,
            updatedY: 0,
        };
        this.height = 400;
        this.width = 600;
        this.padding = 70;
        
        for (const point of this.props.data){
            this.state.points.push(point);
        }
        
        this.state.points.sort((a, b) => { 
            if(a[0] === b[0]){
                if(a[0] === 0){
                    return a[1]-b[1];
                }
                if(a[1] === 1){
                    return b[1]-a[1];
                }
                return a[1]-b[1];
            }
            return a[0] - b[0];
        });

        this.handleDrag = this.handleDrag.bind(this);
        this.handleDoubleClick = this.handleDoubleClick.bind(this);
        this.refCallBack = this.refCallBack.bind(this);
        this.handlePointUpdate = this.handlePointUpdate.bind(this);
        this.change = this.change.bind(this);
    }

    public render() {
        const points = this.renderPoints();
        const ghostPoint = this.state.copyPoint;
        return ([
            <div key="LineGraph">                
                <svg
                    className="msp-canvas"
                    ref={this.refCallBack} 
                    viewBox={`0 0 ${this.width+this.padding} ${this.height+this.padding}`}
                    onMouseMove={this.handleDrag} 
                    onMouseUp={this.handlePointUpdate}
                    onDoubleClick={this.handleDoubleClick}>  
            
                    <g stroke="black" fill="black">
                        <Poly 
                            data={this.state.points} 
                            k={0.5}
                            height={this.height}
                            width={this.width}
                            padding={this.padding}/>
                        {points}
                        {ghostPoint}
                    </g>

                     <defs>
                        <linearGradient id="Gradient">
                            <stop offset="0%" stopColor="#d30000"/>
                            <stop offset="30%" stopColor="#ffff05"/>
                            <stop offset="50%" stopColor="#05ff05"/>
                            <stop offset="70%" stopColor="#05ffff"/>
                            <stop offset="100%" stopColor="#041ae0"/>
                        </linearGradient>
                    </defs>
                    
                </svg>
            </div>,
            <div key="modal" id="modal-root" />
        ]);
    }

    private change(points: Vec2[]){
        let copyPoints = points.slice();
        copyPoints.shift();
        copyPoints.pop();
        this.props.onChange(copyPoints);    
    }

    private handleMouseDown = (id:number) => (event: any) => {
        if(id === 0 || id === this.state.points.length-1){
            return;
        } 
        const copyPoint: Vec2 = this.normalizePoint(Vec2.create(this.state.points[id][0], this.state.points[id][1]));
        this.setState({
            selected: id,
            copyPoint: "ready",
            updatedX: copyPoint[0],
            updatedY: copyPoint[1],
        });

        event.preventDefault();
    }

    private handleDrag(event: any) {
        if(this.state.copyPoint === undefined){
            return
        }
        const pt = this.myRef.createSVGPoint();
        let updatedCopyPoint;
        const padding = this.padding/2;
        pt.x = event.clientX;
        pt.y = event.clientY;
        const svgP = pt.matrixTransform(this.myRef.getScreenCTM().inverse());

        if( svgP.x < (padding) || 
            svgP.x > (this.width+(padding)) || 
            svgP.y > (this.height+(padding)) || 
            svgP.y < (padding)) {
            return;
        }
        updatedCopyPoint = Vec2.create(svgP.x, svgP.y);
        this.setState({
            updatedX: updatedCopyPoint[0],
            updatedY: updatedCopyPoint[1],
        });
        const unNormalizePoint = this.unNormalizePoint(updatedCopyPoint);
        this.setState({
            copyPoint: <PointComponent 
                            selected={false}
                            key="copy" 
                            x={updatedCopyPoint[0]} 
                            y={updatedCopyPoint[1]} 
                            nX={unNormalizePoint[0]} 
                            nY={unNormalizePoint[1]}
                            delete={this.deletePoint}
                            onmouseover={this.props.onHover}/>
        });
        this.props.onDrag(unNormalizePoint);
        event.preventDefault()
    }

    private handlePointUpdate(event: any) {
        const selected = this.state.selected;
        if(selected === undefined || selected === 0 || selected === this.state.points.length-1) {
            this.setState({
                selected: undefined,
                copyPoint: undefined,
            });
            return
        }
        const updatedPoint = this.unNormalizePoint(Vec2.create(this.state.updatedX, this.state.updatedY));
        const points = this.state.points.filter((_,i) => i !== this.state.selected);
        points.push(updatedPoint);;
        points.sort((a, b) => { 
            if(a[0] === b[0]){
                if(a[0] === 0){
                    return a[1]-b[1];
                }
                if(a[1] === 1){
                    return b[1]-a[1];
                }
                return a[1]-b[1];
            }
            return a[0] - b[0];
        });
        this.setState({
            points,
            selected: undefined,
            copyPoint: undefined,
        });
        this.change(points);
        event.preventDefault();
    }

    private handleDoubleClick(event: any) {
        let newPoint;
        const pt = this.myRef.createSVGPoint();
        pt.x = event.clientX;
        pt.y = event.clientY;
        const svgP = pt.matrixTransform(this.myRef.getScreenCTM().inverse());
        const points = this.state.points;
        const padding = this.padding/2; 

        if( svgP.x < (padding) || 
            svgP.x > (this.width+(padding)) || 
            svgP.y > (this.height+(padding)) || 
            svgP.y < (this.padding/2)) {
            return;
        }
        newPoint = this.unNormalizePoint(Vec2.create(svgP.x, svgP.y));
        points.push(newPoint);
        points.sort((a, b) => { 
            if(a[0] === b[0]){
                if(a[0] === 0){
                    return a[1]-b[1];
                }
                if(a[1] === 1){
                    return b[1]-a[1];
                }
                return a[1]-b[1];
            }
            return a[0] - b[0];
        });
        this.setState({points})
        this.change(points);
        event.preventDefault();
    }
    private deletePoint = (i:number) => (event: any) => {
    if(i===0 || i===this.state.points.length-1){ return};
        const points = this.state.points.filter((_,j) => j !== i);
        points.sort((a, b) => { 
            if(a[0] === b[0]){
                if(a[0] === 0){
                    return a[1]-b[1];
                }
                if(a[1] === 1){
                    return b[1]-a[1];
                }
                return a[1]-b[1];
            }
            return a[0] - b[0];
        });
        this.setState({points});
        this.change(points);
        event.stopPropagation();
    }

    private normalizePoint(point: Vec2) {
        const min = this.padding/2;
        const maxX = this.width+min;
        const maxY = this.height+min; 
        const normalizedX = (point[0]*(maxX-min))+min; 
        const normalizedY = (point[1]*(maxY-min))+min;
        const reverseY = (this.height+this.padding)-normalizedY;
        const newPoint = Vec2.create(normalizedX, reverseY);
        return newPoint;
    }

    private unNormalizePoint(point: Vec2) {
        const min = this.padding/2;
        const maxX = this.width+min; 
        const maxY = this.height+min;
        const unNormalizedX = (point[0]-min)/(maxX-min);

        // we have to take into account that we reversed y when we first normalized it.
        const unNormalizedY = ((this.height+this.padding)-point[1]-min)/(maxY-min); 

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
            if(i != 0 && i != this.state.points.length-1){
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
                        onMouseDown={this.handleMouseDown(i)}
                    />);
            }
        }
        return points;
    }
}

function Poly(props: any) {

    const points: Vec2[] = [];
    let min:number;
    let maxX:number;
    let maxY: number;
    let normalizedX: number;
    let normalizedY: number;
    let reverseY: number;
    
    for(const point of props.data){
        min = parseInt(props.padding, 10)/2;
        maxX = parseInt(props.width, 10)+min;
        maxY = parseInt(props.height, 10)+min; 
        normalizedX = (point[0]*(maxX-min))+min; 
        normalizedY = (point[1]*(maxY-min))+min;
        reverseY = (props.height+props.padding)-normalizedY;
        points.push(Vec2.create(normalizedX, reverseY));
    }

    if (props.k == null) {props.k = 0.3};
    const data = points;
    const size = data.length;
    const last = size - 2;
    let path = "M" + [data[0][0], data[0][1]];

    for (let i=0; i<size-1;i++){
        const x0 = i ? data[i-1][0] : data[0][0];
        const y0 = i ? data[i-1][1] : data[0][1];

        const x1 = data[i][0];
        const y1 = data[i][1];

        const x2 = data[i+1][0];
        const y2 = data[i+1][1];

        const x3 = i !== last ? data[i+2][0] : x2;
        const y3 = i !== last ? data[i+2][1] : y2; 

        const cp1x = x1 + (x2 - x0)/6 * props.k;
        const cp1y = y1 + (y2 -y0)/6 * props.k;

        const cp2x = x2 - (x3 -x1)/6 * props.k;
        const cp2y = y2 - (y3 - y1)/6 * props.k;

        path += "C" + [cp1x, cp1y, cp2x, cp2y, x2, y2];
    }

    return <path d={path} strokeWidth="5" stroke="#cec9ba" fill="none"/>
}