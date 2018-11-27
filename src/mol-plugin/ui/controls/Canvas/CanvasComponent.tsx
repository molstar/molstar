/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Paul Luna <paulluna0215@gmail.com>
 */
import PointComponent from './PointComponent';

import * as React from 'react';
import { Vec2 } from 'mol-math/linear-algebra';

// interface Line {
//     x1: number;
//     x2: number;
//     y1: number;
//     y2: number;
// };

// interface Rect {
//     x: number;
//     y: number;
// }

// interface TextLabel {
//     x: number;
//     y: number;
// }

interface CanvasComponentState {
    points: Vec2[],
    selectedPoint: any,
    selected: number | undefined,
    copyPoint: any,
    updatedX: number, 
    updatedY: number,
}

export default class CanvasComponent extends React.Component<any, CanvasComponentState> {
    private myRef:any;
    private height: number;
    private width: number;
    private padding: number;
    // private topLine: Line;
    // private bottomLine: Line;
    // private rightLine: Line;
    // private leftLine: Line;
    // private normalizedXLabels: Point[];
    // private normalizedYLabels: Point[];
    // private normalizedRect: Rect;
    // private normalizedLabel: TextLabel;
    // private normalizedInfoPoint: Point;
    // private normalizedInfoLabel: Point;
    // private label: TextLabel = {x: 0, y: 0.1};
    // private infoPoint = {x: 0.95, y: 0.5};
    // private infoLabel = {x: 0.945, y: 0.2};
    // private rainbowRect: Rect = {x: 0.0, y: 0.5};
    // private xLabel: Point[] = [
    //     {x: 0, y: 0.65},
    //     {x: 0.07, y: 0.65},
    //     {x: 0.17, y: 0.65},
    //     {x: 0.27, y: 0.65},
    //     {x: 0.37, y: 0.65},
    //     {x: 0.47, y: 0.65},
    //     {x: 0.57, y: 0.65},
    //     {x: 0.67, y: 0.65},
    //     {x: 0.77, y: 0.65},
    //     {x: 0.87, y: 0.65},
    //     {x: 0.97, y: 0.65},
    // ];

    // private yLabel: Point[] = [
    //     {x: 0.45, y: 0},
    //     {x: 0.45, y: 0.1},
    //     {x: 0.45, y: 0.2},
    //     {x: 0.45, y: 0.3},
    //     {x: 0.45, y: 0.4},
    //     {x: 0.45, y: 0.5},
    //     {x: 0.45, y: 0.6},
    //     {x: 0.45, y: 0.7},
    //     {x: 0.45, y: 0.8},
    //     {x: 0.45, y: 0.9},
    //     {x: 0.45, y: 1},       
    // ];

    
    constructor(props: any) {
        super(props);
        this.myRef = React.createRef();
        this.state = {
            points:[
                Vec2.create(0, 0),
                Vec2.create(1, 0)
            ],
            selectedPoint: undefined,
            selected: undefined,
            copyPoint: undefined,
            updatedX: 0,
            updatedY: 0,
        };
        this.height = 400;
        this.width = 600;
        this.padding = 70;
        // this.normalizedXLabels = this.normalizeXLabel(this.xLabel);
        // this.normalizedYLabels = this.normalizeYLabel(this.yLabel);
        // this.topLine = {
        //     x1: this.padding/2,
        //     x2: this.width+this.padding/2,  
        //     y1: this.padding/2,
        //     y2: this.padding/2
        // };
        // this.rightLine = {
        //     x1: this.width+this.padding/2,
        //     x2: this.width+this.padding/2,
        //     y1: this.padding/2,
        //     y2: this.height+this.padding/2
        // };
        // this.bottomLine = {
        //     x1: this.padding/2,
        //     x2: this.width+this.padding/2,
        //     y1: this.height+this.padding/2,
        //     y2: this.height+this.padding/2
        // };
        // this.leftLine = {
        //     x1: this.padding/2,
        //     x2: this.padding/2, 
        //     y1: this.padding/2,
        //     y2: this.height+this.padding/2
        // };
        // this.normalizedLabel = this.normalizeGraphLabel(
        //     this.label, 
        //     this.padding/2, 
        //     this.height+this.padding/2, 
        //     this.width+this.padding/2, 
        //     this.height+this.padding
        // );
        // this.normalizedInfoPoint = this.normalizeGraphLabel(
        //     this.infoPoint, 
        //     this.padding/2, 
        //     this.height+this.padding/2, 
        //     this.width+this.padding/2, 
        //     this.height+this.padding
        // );
        // this.normalizedInfoLabel = this.normalizeGraphLabel(
        //     this.infoLabel, 
        //     this.padding/2, 
        //     this.height+this.padding/2, 
        //     this.width+this.padding/2, 
        //     this.height+this.padding
        // );
        // this.normalizedRect = this.normalizeGraphLabel(
        //     this.rainbowRect,
        //     this.padding/2,
        //     0,
        //     this.width+this.padding/2,
        //     this.padding/2
        // );
        
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
        this.handleClickCanvas = this.handleClickCanvas.bind(this);
        this.handleDoubleClick = this.handleDoubleClick.bind(this);
        this.refCallBack = this.refCallBack.bind(this);
        this.handlePointUpdate = this.handlePointUpdate.bind(this);
        this.change = this.change.bind(this);
    }

    public render() {
        const points = this.renderPoints();
        const ghostPoint = this.state.copyPoint;
        const selectedPoint = this.state.selectedPoint;
        return ([
            <div key="canvas" id="canvas">                
                <svg 
                    id="svg"
                    className="msp-canvas"
                    ref={this.refCallBack} 
                    viewBox={`0 0 ${this.width+this.padding} ${this.height+this.padding}`}
                    onMouseMove={this.handleDrag} 
                    onMouseUp={this.handlePointUpdate}
                    onClick={this.handleClickCanvas}
                    onDoubleClick={this.handleDoubleClick}>  
                    {/* <text x={this.normalizedLabel.x} y={this.normalizedLabel.y} fontSize="250%">YOUR LABEL HERE</text> */}
                    {/* <InfoComponent 
                        key="rootInfo" 
                        cx={this.normalizedInfoPoint.x} 
                        cy={this.normalizedInfoPoint.y} 
                        x={this.normalizedInfoLabel.x} 
                        y={this.normalizedInfoLabel.y} 
                        fill="white" 
                        stroke="black" 
                        strokeWidth="5" 
                    /> */}
                    <g stroke="black" fill="black">
                        <Poly 
                            data={this.state.points} 
                            k={0.5}
                            height={this.height}
                            width={this.width}
                            padding={this.padding}/>
                        {points}
                        {ghostPoint}
                        {selectedPoint}

                        {/* <line x1={this.topLine.x1} x2={this.topLine.x2} y1={this.topLine.y1} y2={this.topLine.y2} stroke="black" strokeWidth="5" strokeLinecap="square"/>
                        <line x1={this.bottomLine.x1} x2={this.bottomLine.x2} y1={this.bottomLine.y1} y2={this.bottomLine.y2} stroke="black" strokeWidth="5" strokeLinecap="square"/>
                        <line x1={this.rightLine.x1} x2={this.rightLine.x2} y1={this.rightLine.y1} y2={this.rightLine.y2} stroke="black" strokeWidth="5" strokeLinecap="square"/>
                        <line x1={this.leftLine.x1} x2={this.leftLine.x2} y1={this.leftLine.y1} y2={this.leftLine.y2} stroke="black" strokeWidth="5" strokeLinecap="square"/> */}
                    </g>
                    {/* <g className="x-labels">
                        <text x={this.normalizedXLabels[0].x} y={this.normalizedXLabels[0].y} fontSize="250%">0</text>
                        <text x={this.normalizedXLabels[1].x} y={this.normalizedXLabels[1].y} fontSize="250%">0.1</text>
                        <text x={this.normalizedXLabels[2].x} y={this.normalizedXLabels[2].y} fontSize="250%">0.2</text>
                        <text x={this.normalizedXLabels[3].x} y={this.normalizedXLabels[3].y} fontSize="250%">0.3</text>
                        <text x={this.normalizedXLabels[4].x} y={this.normalizedXLabels[4].y} fontSize="250%">0.4</text>
                        <text x={this.normalizedXLabels[5].x} y={this.normalizedXLabels[5].y} fontSize="250%">0.5</text>
                        <text x={this.normalizedXLabels[6].x} y={this.normalizedXLabels[6].y} fontSize="250%">0.6</text>
                        <text x={this.normalizedXLabels[7].x} y={this.normalizedXLabels[7].y} fontSize="250%">0.7</text>
                        <text x={this.normalizedXLabels[8].x} y={this.normalizedXLabels[8].y} fontSize="250%">0.8</text>
                        <text x={this.normalizedXLabels[9].x} y={this.normalizedXLabels[9].y} fontSize="250%">0.9</text>
                        <text x={this.normalizedXLabels[10].x} y={this.normalizedXLabels[10].y} fontSize="250%">1</text>
                    </g>
                    <g className="y-labels">
                        <text x={this.normalizedYLabels[0].x} y={this.normalizedYLabels[0].y} fontSize="250%">0</text>
                        <text x={this.normalizedYLabels[1].x} y={this.normalizedYLabels[1].y} fontSize="250%">0.1</text>
                        <text x={this.normalizedYLabels[2].x} y={this.normalizedYLabels[2].y} fontSize="250%">0.2</text>
                        <text x={this.normalizedYLabels[3].x} y={this.normalizedYLabels[3].y} fontSize="250%">0.3</text>
                        <text x={this.normalizedYLabels[4].x} y={this.normalizedYLabels[4].y} fontSize="250%">0.4</text>
                        <text x={this.normalizedYLabels[5].x} y={this.normalizedYLabels[5].y} fontSize="250%">0.5</text>
                        <text x={this.normalizedYLabels[6].x} y={this.normalizedYLabels[6].y} fontSize="250%">0.6</text>
                        <text x={this.normalizedYLabels[7].x} y={this.normalizedYLabels[7].y} fontSize="250%">0.7</text>
                        <text x={this.normalizedYLabels[8].x} y={this.normalizedYLabels[8].y} fontSize="250%">0.8</text>
                        <text x={this.normalizedYLabels[9].x} y={this.normalizedYLabels[9].y} fontSize="250%">0.9</text>
                        <text x={this.normalizedYLabels[9].x} y={this.normalizedYLabels[10].y} fontSize="250%">1</text>
                    </g> */}

                     <defs>
                        <linearGradient id="Gradient">
                            <stop offset="0%" stopColor="#d30000"/>
                            <stop offset="30%" stopColor="#ffff05"/>
                            <stop offset="50%" stopColor="#05ff05"/>
                            <stop offset="70%" stopColor="#05ffff"/>
                            <stop offset="100%" stopColor="#041ae0"/>
                        </linearGradient>
                    </defs>
 
                    {/* <rect id="rect1" x={this.normalizedRect.x} y={this.normalizedRect.y} rx="15" ry="15" width={this.width} height="30" fill="url(#Gradient)"/> */}
                    
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

    // private handleClickPoint = (id:number) => (event: any) => {
    //     const selectedPointX = this.state.points[id][0];
    //     const selectedPointY = this.state.points[id][1];
    //     const normalizePoint = this.normalizePoint(Vec2.create(selectedPointX, selectedPointY));
    //     const selectedPoint = <PointComponent
    //                         key={id} 
    //                         id={id}
    //                         x={normalizePoint[0]} 
    //                         y={normalizePoint[1]}
    //                         nX={selectedPointX}
    //                         nY={selectedPointY}
    //                         selected={true}
    //                         delete={this.deletePoint}
    //                         onMouseDown={this.handleMouseDown(id)}
    //                     />

    //     this.setState({selectedPoint});

    //     event.stopPropagation();
    // }
    
    private handleClickCanvas() {
        this.setState({selectedPoint: undefined});
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

    // private normalizeGraphLabel(point: Point, minX: number, minY: number, maxX: number, maxY: number) {
    //     const normalizedX = (point.x*(maxX-minX))+minX;
    //     const normalizedY = (point.y*(maxY-minY))+minY;
    //     const reverseY = (this.height+this.padding)-normalizedY;

    //     return {x: normalizedX, y: reverseY};
    // } 

    // private normalizeXLabel(points: Point[]) {
    //     const minX = this.padding/2;
    //     const minY = 0;
    //     const maxX = this.width+this.padding/2;
    //     const maxY = this.padding/2;
    //     const normalizedPoints: Point[] = [];
    //     for(const point of points){
    //         normalizedPoints.push(this.normalizeGraphLabel(point, minX, minY, maxX, maxY));
    //     }

    //     return normalizedPoints;
    // }

    // private normalizeYLabel(points: Point[]) {
    //     const minX = 0;
    //     const minY = this.padding/2;
    //     const maxX = this.padding/2;
    //     const maxY = this.height+this.padding/2;
    //     const normalizedPoints: Point[] = [];
    //     for(const point of points) {
    //         normalizedPoints.push(this.normalizeGraphLabel(point, minX, minY, maxX, maxY));
    //     }
    //     return normalizedPoints;
    // }

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
                        // onClick={this.handleClickPoint(i)}
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