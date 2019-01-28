/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Paul Luna <paulluna0215@gmail.com>
 */
import * as React from 'react';

import { Vec2 } from 'mol-math/linear-algebra';
import PointComponent from './point-component';

interface LineGraphComponentState {[x:string]:any,
    points: Vec2[],
    copyPoint: any,
    canSelectMultiple: boolean,
    x: string | number,
    y: string | number,
    selected: number[],
    isControlEnabled: {[name: string]: boolean};
}

export default class LineGraphComponent extends React.Component<any, LineGraphComponentState> {
    private myRef:any;
    private height: number;
    private width: number;
    private padding: number;
    private updatedX: number;
    private updatedY: number;
    private ghostPoints: {id: number, element: SVGElement}[];
    private leftMost:  SVGElement;
    private rightMost: SVGElement;
    private upperMost: SVGElement;
    private bottomMost: SVGElement;
    private gElement: SVGElement;
    private namespace: string;
    private userInput: {[name:string]: number} = {};
    private mouseStartPoint: Vec2;
    private mouseDown: boolean;
    private hasDragged: boolean;
    private data: any;
    

    constructor(props: any) {
        super(props);
        this.myRef = React.createRef();
        this.state = {
            points:[
                Vec2.create(0, 0),
                Vec2.create(1, 0)
            ],
            copyPoint: undefined,
            selected: [],
            canSelectMultiple: false,
            x: '',
            y: '',
            isControlEnabled: {'plus': false, 'minus': false},
        };
        this.height = 400;
        this.width = 600;
        this.padding = 70;

        this.ghostPoints = [];
        this.namespace = 'http://www.w3.org/2000/svg';
        this.userInput['x'] = -1;
        this.userInput['y'] = -1;
        this.mouseDown = false;
        this.hasDragged = false;        
        this.leftMost = document.createElementNS(this.namespace, 'circle') as SVGElement;
        this.leftMost.setAttribute('cx', '0');
        
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
        this.handleMultipleDrag = this.handleMultipleDrag.bind(this);
        this.handleDoubleClick = this.handleDoubleClick.bind(this);
        this.refCallBack = this.refCallBack.bind(this);
        this.change = this.change.bind(this);
        this.handleKeyUp = this.handleKeyUp.bind(this);
        this.handleLeave = this.handleLeave.bind(this);
        this.handleEnter = this.handleEnter.bind(this);
        this.handleSubmit = this.handleSubmit.bind(this);
        this.handleChange = this.handleChange.bind(this);
        this.addPoint = this.addPoint.bind(this);
        this.deletePoint = this.deletePoint.bind(this);
        this.handleCanvasClick = this.handleCanvasClick.bind(this); 
        this.addGhostPoint = this.addGhostPoint.bind(this);
        this.setControls = this.setControls.bind(this);
    }

    public render() {
        const points = this.renderPoints();
        const lines = this.renderLines();
        
        return ([
            <div key="LineGraph">                
                <svg
                    className="msp-canvas"
                    ref={this.refCallBack} 
                    viewBox={`0 0 ${this.width+this.padding} ${this.height+this.padding}`}
                    onMouseMove={this.state.canSelectMultiple? this.handleMultipleDrag : this.handleDrag} 
                    onClick={this.handleCanvasClick}
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
                <div key="line-graph-controls" className="line-graph-control">
                <ul style={{margin: "5px 50px", listStyle: "none"}}>
                    <li style={{display: "inline"}}><div style={{margin: "-2px 5px"}} title="Delete point" className={this.state.isControlEnabled['minus'] ? "control control-minus" : "control disabled-minus"} onClick={this.deletePoint}></div></li>
                    <li style={{display: "inline"}}><div style={{margin: "-2px 5px"}} title="Add point" className={this.state.isControlEnabled['plus'] ? "control control-plus" : "control disabled-plus"} onClick={this.handleSubmit}></div></li>
                    <li style={{display: "inline-block", margin: "auto 2px"}}><input min="0" max="1" step="0.1" type="number" placeholder="x" name="x" value={this.state.x} onChange={this.handleChange} required/></li>
                    <li style={{display: "inline-block", margin: "auto 2px"}}><input min="0" max="1" step="0.1" type="number" placeholder="y" name="y" value={this.state.y} onChange={this.handleChange} required/></li>
                </ul>
            </div>
            </div>,
        ]);
    }

    componentDidMount() {
        this.gElement = document.getElementsByClassName('ghost-points')[0] as SVGElement;
        document.addEventListener('keydown', this.handleKeyDown);
        document.addEventListener('keyup', this.handleKeyUp);
    }

    componentWillUnmount() {
        document.removeEventListener('keydown', this.handleKeyDown);
        document.removeEventListener('keyup', this.handleKeyUp);
    }

    private change(points: Vec2[]){
        let copyPoints = points.slice();
        copyPoints.shift();
        copyPoints.pop();
        this.props.onChange(copyPoints);    
    }

    private handleKeyDown = (event: any) => {
        if(event.keyCode == 16) {
            this.setControls(false, false);
            this.setState({
                canSelectMultiple: true, 
            });
        }
    }

    private handleKeyUp = (event: any) => {
        if(event.keyCode != 16) { return; }
        this.ghostPoints = [];
        this.gElement.innerHTML = '';
        this.setState({canSelectMultiple: false});
    }

    private handleCanvasClick() {
        if(this.state.canSelectMultiple) {
            return;
        }
        this.setControls(false, false);
        this.gElement.innerHTML = '';
        this.ghostPoints = [];
        this.leftMost.setAttribute('cx', '0');
        this.setState({
            selected: [],
            x: '',
            y: ''
        });
    }

    private handlePointClick = (event:any) => {
        event.stopPropagation();
        const id = parseInt(event.target.id);
        if(isNaN(id)) {
            for(let i=0; i<this.ghostPoints.length; i++) {
                if(this.ghostPoints[i].element == event.target) {
                    let ghostPointElements = this.gElement.getElementsByClassName('ghostPoint');
                    this.gElement.removeChild(ghostPointElements[i]);
                    this.ghostPoints = this.ghostPoints.filter((_, j) => j !== i);
                    return;
                }
            }
        }

        if(this.state.canSelectMultiple) {
            let size = this.ghostPoints.length;
            this.ghostPoints[size-1].element.setAttribute('style', 'display: visible');
            this.ghostPoints[size-1].element.addEventListener('mousedown', this.handleMouseDown); 
            this.ghostPoints[size-1].element.addEventListener('click', this.handlePointClick);
            this.gElement.appendChild(this.ghostPoints[size-1].element);
            return;
        } 
        
        this.setControls(true, true)
        this.ghostPoints[0].element.setAttribute('style', 'display: visible');
        this.ghostPoints[0].element.addEventListener('mousedown', this.handleMouseDown);
        this.setState({
            x: this.state.points[this.state.selected[0]][0],
            y: this.state.points[this.state.selected[0]][1]
        });

    }

    private handleMouseDown = (event: any) => {
        let selected;
        let size;
        const id = parseInt(event.target.id);
        const x = event.target.cx.animVal.value;
        const y = event.target.cy.animVal.value;

        if(id === 0 || id === this.state.points.length-1){
            return;
        }

        this.mouseDown = true;
        if(isNaN(id)) {
            this.mouseStartPoint = Vec2.create(x, y);
            return;
        }

        if (this.state.canSelectMultiple) {
            this.addGhostPoint(id, x, y);
            size = this.ghostPoints.length;
            this.gElement.appendChild(this.ghostPoints[size-1].element);
            this.mouseStartPoint = Vec2.create(x, y);
            selected = this.state.selected;
            selected.push(id);
            this.setState({selected: selected});
            return;
        }
        
        const copyPoint: Vec2 = Vec2.create(x, y);
        this.ghostPoints = [];
        this.gElement.innerHTML = '';
        this.addGhostPoint(id, x, y);
        this.gElement.appendChild(this.ghostPoints[0].element);
        this.updatedX = copyPoint[0];
        this.updatedY = copyPoint[1];
        this.setState({selected: [id]});

    }

    private handleDrag(event: any) {
        if(this.state.selected.length === 0 || !this.mouseDown){
            return
        }

        let updatedCopyPoint;
        let svgP;
        const pt = this.myRef.createSVGPoint();
        const padding = this.padding/2;
        this.setControls(false, false);
        pt.x = event.clientX;
        pt.y = event.clientY;
        svgP = pt.matrixTransform(this.myRef.getScreenCTM().inverse());
        
        if ((svgP.x < (padding) || svgP.x > (this.width+(padding))) && (svgP.y > (this.height+(padding)) || svgP.y < (padding))) {
            updatedCopyPoint = Vec2.create(this.updatedX, this.updatedY);
        }
        else if (svgP.x < padding) {
            updatedCopyPoint = Vec2.create(padding, svgP.y);
        }
        else if( svgP.x > (this.width+(padding))) {
            updatedCopyPoint = Vec2.create(this.width+padding, svgP.y);
        }
        else if (svgP.y > (this.height+(padding))) {
            updatedCopyPoint = Vec2.create(svgP.x, this.height+padding);
        }
        else if (svgP.y < (padding)) {
            updatedCopyPoint = Vec2.create(svgP.x, padding);
        } else {
            updatedCopyPoint = Vec2.create(svgP.x, svgP.y);
        }
        this.updatedX = updatedCopyPoint[0];
        this.updatedY = updatedCopyPoint[1];
        const unNormalizePoint = this.unNormalizePoint(updatedCopyPoint);
        this.ghostPoints[0].element.setAttribute('style', 'display: visible');
        this.ghostPoints[0].element.setAttribute('cx', `${updatedCopyPoint[0]}`);
        this.ghostPoints[0].element.setAttribute('cy', `${updatedCopyPoint[1]}`);
        this.props.onDrag(unNormalizePoint);
        this.hasDragged = true;
        this.setState({x: '', y: ''});
    }

    private handleMultipleDrag(event: any) {
        if(!this.mouseDown) {
            return;
        }

        const padding = this.padding/2;
        const pt = this.myRef.createSVGPoint(); // create a temp point for the cursor pointer
        let updatedGhostPoint: Vec2;
        let ghostPoint: Vec2;
        let selected = 0;

        pt.x = event.clientX;
        pt.y = event.clientY;
        const svgP = pt.matrixTransform(this.myRef.getScreenCTM().inverse());
        const directionalVector = Vec2.create(svgP.x-this.mouseStartPoint[0], svgP.y-this.mouseStartPoint[1]);
        const leftMostInt = parseInt(this.leftMost.getAttribute('cx') as string);
        const rightMostInt = parseInt(this.rightMost.getAttribute('cx') as string);
        const upperMostInt = parseInt(this.upperMost.getAttribute('cy') as string);
        const bottomMostInt = parseInt(this.bottomMost.getAttribute('cy') as string);
        if((directionalVector[0]+leftMostInt) <= padding) {
            directionalVector[0] = padding-leftMostInt;
        }
        
        if((directionalVector[0]+rightMostInt) >= this.width+padding) {
            directionalVector[0] = (this.width+padding)-rightMostInt;
        }

        if((directionalVector[1]+upperMostInt) <= padding) {
            directionalVector[1] = padding-upperMostInt;
        }

        if((directionalVector[1]+bottomMostInt) >= (this.height+padding)) {
            directionalVector[1] = (this.height+padding)-bottomMostInt;
        }

        for(let i = 0; i < this.ghostPoints.length; i++) {
            const tempX = this.ghostPoints[i].element.getAttribute('cx');
            const tempY = this.ghostPoints[i].element.getAttribute('cy');
            ghostPoint = Vec2.create(parseInt(tempX as string), parseInt(tempY as string));
            if(this.mouseStartPoint[0] == ghostPoint[0]) {selected = i;}
            updatedGhostPoint = Vec2.create(ghostPoint[0]+directionalVector[0], ghostPoint[1]+directionalVector[1]);
            this.ghostPoints[i].element.setAttribute('cx', `${updatedGhostPoint[0]}`);
            this.ghostPoints[i].element.setAttribute('cy', `${updatedGhostPoint[1]}`);
            this.ghostPoints[i].element.setAttribute('style', 'style: visible'); 
            this.sortOuterGhostPoints(this.ghostPoints[i].element);
        }

        const x = parseInt(this.ghostPoints[selected].element.getAttribute('cx') as string);
        const y = parseInt(this.ghostPoints[selected].element.getAttribute('cy') as string);
        this.mouseStartPoint = Vec2.create(x, y);
        this.hasDragged = true;
    }

    private handlePointUpdate = (event: any) => {
        const selected = this.state.selected;
        this.mouseDown = false;
        if ((this.state.canSelectMultiple && !this.hasDragged) || !this.hasDragged) { 
            return; 
        }
        if(selected.length === 0 || selected[0] === 0 || selected[selected.length-1] === this.state.points.length-1) {
            this.setState({
                copyPoint: undefined,
            });
            return;
        }
        let points = this.state.points;
        for(let i = 0; i < this.ghostPoints.length; i++) {
            const id = this.ghostPoints[i].id;
            const element = this.ghostPoints[i].element;
            const x = parseInt(element.getAttribute('cx') as string);
            const y = parseInt(element.getAttribute('cy') as string);
            const updatedPoint = this.unNormalizePoint(Vec2.create(x, y));
            points[id] = updatedPoint;
        }

        points.sort((a, b) => { 
            if(a[0] === b[0]){
                if(a[0] === 1){
                    return b[1]-a[1];
                }
                return a[1]-b[1];
            }
            return a[0] - b[0];
        });
        this.setState({
            points,
            selected: [],
        });
        this.change(points);
        this.gElement.innerHTML = '';
        this.ghostPoints.forEach(x => {
            x.element.removeEventListener('mousedown', this.handleMouseDown);
        });
        this.ghostPoints = [];
        this.hasDragged = false;
        this.leftMost.setAttribute('cx', '0');
        document.removeEventListener("mousemove", this.handleDrag, true);
        document.removeEventListener("mouseup", this.handlePointUpdate, true);
    }

    private handleDoubleClick(event: any) {
        let newPoint;
        const pt = this.myRef.createSVGPoint();
        pt.x = event.clientX;
        pt.y = event.clientY;
        const svgP = pt.matrixTransform(this.myRef.getScreenCTM().inverse());
        const padding = this.padding/2; 
        if( svgP.x < (padding) || 
            svgP.x > (this.width+(padding)) ||
            svgP.y > (this.height+(padding)) || 
            svgP.y < (this.padding/2)) {
            return;
        }
        newPoint = this.unNormalizePoint(Vec2.create(svgP.x, svgP.y));
        this.data = newPoint
        this.addPoint();
    }

    private handleSubmit(event: any) {
        const x = parseFloat(this.state.x as string);
        const y = parseFloat(this.state.y as string);
        const point = Vec2.create(x, y);
        this.userInput['x'] = -1;
        this.userInput['y'] = -1;
        this.setState({
            x: '',
            y: '',
        });

        this.data = point;
        if(this.state.selected.length != 0) {
            const originalX = this.state.points[this.state.selected[0]][0];
            const originalY = this.state.points[this.state.selected[0]][1];
            if(x < originalX || y < originalY) {
                this.deletePoint(event, this.addPoint);
            } else {
                this.addPoint(this.deletePoint);
            }
        } else {
            this.addPoint();
        }
        this.setControls(false, false);
    }

    private handleChange(event: any){
        this.userInput[event.target.name] = event.target.value;
        if(event.target.value === '') { this.userInput[event.target.name] = -1; }
        if(this.userInput['x'] > -1 && this.userInput['y'] > -1) {
            this.setControls(true, this.state.isControlEnabled['minus']);
            this.setState({[event.target.name]: event.target.value});
        } else {
            this.setState({[event.target.name]: event.target.value});
        }
    }

    private addPoint(callBack?: any) {
        const point = this.data;
        const points = this.state.points;
        points.push(point);
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
        
        this.change(points);
        this.setState({points}, callBack);
    }

    private deletePoint(event?: any, callBack?: any) {
        let points;
        const i = this.state.selected[0];
        if(i===0 || i===this.state.points.length-1){ return; }
        points = this.state.points.filter((_,j) => j !== i);
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
        this.gElement.innerHTML = '';
        this.setControls(false, false);
        this.change(points);

        this.setState({
            points,
            selected: [],
            x: '',
            y: ''
        }, callBack);

        if(event)
            event.stopPropagation();
    }

    private handleLeave() {
        if(this.state.selected.length === 0) {
            return;
        }

        if(this.state.canSelectMultiple){
            document.addEventListener('mousemove', this.handleMultipleDrag, true);
       } else {
           document.addEventListener('mousemove', this.handleDrag, true);
       }

        document.addEventListener('mouseup', this.handlePointUpdate, true);
    }

    private handleEnter() {
        if(this.state.canSelectMultiple) {
            document.removeEventListener('mousemove', this.handleMultipleDrag, true);
        } else {
            document.removeEventListener('mousemouse', this.handleDrag, true);
        }
        document.removeEventListener('mouseup', this.handlePointUpdate, true);
    }

    private addGhostPoint(id: number, x: number, y: number) {
        this.ghostPoints.push({id: id, element: document.createElementNS(this.namespace, 'circle') as SVGElement});
        const size = this.ghostPoints.length;
        this.ghostPoints[size-1].element.setAttribute('r', '10');
        this.ghostPoints[size-1].element.setAttribute('fill', 'orange');
        this.ghostPoints[size-1].element.setAttribute('cx', `${x}`);
        this.ghostPoints[size-1].element.setAttribute('cy', `${y}`);
        this.ghostPoints[size-1].element.setAttribute('style', 'display: none');
        this.ghostPoints[size-1].element.setAttribute('class', 'ghostPoint');
        this.sortOuterGhostPoints(this.ghostPoints[size-1].element);
    }

    private setControls(plus: boolean, minus: boolean) {
        let isControlEnabled = {...this.state.isControlEnabled};
        isControlEnabled['plus'] = plus;
        isControlEnabled['minus'] = minus;
        this.setState({isControlEnabled});
    }

    private sortOuterGhostPoints(element: SVGElement) {
        if(this.leftMost.getAttribute('cx') == '0') {
            this.leftMost = element;
            this.rightMost = element;
            this.upperMost = element;
            this.bottomMost = element;
        }

        const x = parseInt(element.getAttribute('cx') as string);
        const y = parseInt(element.getAttribute('cy') as string);

        if(x < parseInt(this.leftMost.getAttribute('cx') as string)) {
            this.leftMost = element;
        }
        
        if(x > parseInt(this.rightMost.getAttribute('cx') as string)) {
            this.rightMost = element;
        }

        if(y > parseInt(this.rightMost.getAttribute('cy') as string)) {
            this.bottomMost = element;
        }

        if(y < parseInt(this.upperMost.getAttribute('cy') as string)) {
            this.upperMost = element;
        }
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
                        onmouseover={this.props.onHover}
                        onmousedown={this.handleMouseDown}
                        onclick={this.handlePointClick}
                    />);
            }
        }
        return points;
    }

    private renderLines() {
        const points: Vec2[] = [];
        let lines = [];
        let min:number;
        let maxX:number;
        let maxY: number;
        let normalizedX: number;
        let normalizedY: number;
        let reverseY: number;

        for(const point of this.state.points){
            min = this.padding/2;
            maxX = this.width+min;
            maxY = this.height+min; 
            normalizedX = (point[0]*(maxX-min))+min; 
            normalizedY = (point[1]*(maxY-min))+min;
            reverseY = this.height+this.padding-normalizedY;
            points.push(Vec2.create(normalizedX, reverseY));
        }

        const data = points;
        const size = data.length;

        for (let i=0; i<size-1;i++){
            const x1 = data[i][0];
            const y1 = data[i][1];
            const x2 = data[i+1][0];
            const y2 = data[i+1][1];
            
            lines.push(<line key={`lineOf${i}`} x1={x1} x2={x2} y1={y1} y2={y2} stroke="#cec9ba" strokeWidth="5"/>)
        }
        
        return lines;
    }
}