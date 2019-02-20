/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Paul Luna <paulluna0215@gmail.com>
 */
import * as React from 'react';

import { Vec2 } from 'mol-math/linear-algebra';
import { getEventListener } from './line-graph-controls';
import PointComponent from './point-component';
import { Lines } from './lines';
import { normalize, unNormalize } from './normilization';

interface LineGraphComponentState {[x:string]:any,
    points: Vec2[],
    copyPoint: any,
    canSelectMultiple: boolean,
    x: string | number,
    y: string | number,
    selected: number[],
    isControlEnabled: {[name: string]: boolean};
    mouseDown: boolean,
    hasDragged: boolean,
    mouseStartPoint: Vec2,
    ghostPoints: {id: number, element: SVGElement}[];
    myRef: any;
    ghostPointsWrapper: any;
    padding: number;
    height: number;
    width: number;
    leftMost: SVGElement;
    rightMost: SVGElement;
    upperMost: SVGElement;
    bottomMost: SVGElement;
    displayPoint: Vec2;
    createdLasso: boolean;
}

export default class LineGraphComponent extends React.Component<any, LineGraphComponentState> {
    private namespace: string;
    private userInput: {[name:string]: number} = {};
    private data: any;
    

    constructor(props: any) {
        super(props);
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
            createdLasso: false,
            mouseDown: false,
            hasDragged: false,
            mouseStartPoint: Vec2.create(0, 0),
            ghostPoints: [],
            myRef: React.createRef,
            ghostPointsWrapper: '',
            padding: props.padding,
            height: props.height,
            width: props.width,
            displayPoint: Vec2.create(-1, -1),
            leftMost: document.createElementNS(this.namespace, 'circle') as SVGElement,
            rightMost: document.createElementNS(this.namespace, 'circle') as SVGElement,
            upperMost: document.createElementNS(this.namespace, 'circle') as SVGElement,
            bottomMost: document.createElementNS(this.namespace, 'circle') as SVGElement,
        };

        this.namespace = 'http://www.w3.org/2000/svg';
        this.userInput['x'] = -1;
        this.userInput['y'] = -1;
        this.state.leftMost.setAttribute('cx', '0');
        
        for (const point of this.props.data){
            this.state.points.push(point);
        }
        
        this.sortPoints(this.state.points);

        this.refCallBack = this.refCallBack.bind(this);
        this.change = this.change.bind(this);
        this.handleKeyUp = this.handleKeyUp.bind(this);
        this.handleLeave = this.handleLeave.bind(this);
        this.handleEnter = this.handleEnter.bind(this);
        this.handleSubmit = this.handleSubmit.bind(this);
        this.handleChange = this.handleChange.bind(this);
        this.addPoint = this.addPoint.bind(this);
        this.handlePointUpdate = this.handlePointUpdate.bind(this);
        this.deletePoint = this.deletePoint.bind(this);
        this.setControls = this.setControls.bind(this);
    }

    public render() {
        const points = this.renderPoints();
        
        return ([
            <div key="LineGraph">                
                <svg
                    className="msp-canvas"
                    ref={this.refCallBack} 
                    viewBox={`0 0 ${this.state.width+this.state.padding} ${this.state.height+this.state.padding}`}
                    onMouseMove={this.handleEvent}
                    onMouseDown={this.handleEvent}
                    onClick={this.handleEvent}
                    onMouseUp={this.handleEvent }
                    onMouseLeave={this.handleLeave}
                    onMouseEnter={this.handleEnter}
                    tabIndex={0}
                    onKeyDown={this.handleKeyDown}
                    onKeyUp={this.handleKeyUp}
                    onDoubleClick={this.handleEvent}>  
            
                    <g stroke="black" fill="black">
                        <Lines type='straight' points={this.state.points} width={this.state.width} height={this.state.height} padding={this.state.padding} />
                        {points}
                    </g>
                    <g className="ghost-points" stroke="black" fill="black">
                    </g>
                    <g className="rect-select" fillOpacity="0.25;">
                    </g>
                </svg>
                <div key="line-graph-controls" className="line-graph-control">
                <ul style={{margin: "5px 50px", listStyle: "none"}}>
                    <li style={{display: "inline"}}><div style={{margin: "-2px 5px"}} title="Delete point" className={this.state.isControlEnabled['minus'] ? "control control-minus" : "control disabled-minus"} onClick={this.deletePoint}></div></li>
                    <li style={{display: "inline-block", margin: "auto 2px"}}><input min="0" max="1" step="0.01" type="number" placeholder="x" name="x" value={this.state.x} onChange={this.handleChange} required/></li>
                    <li style={{display: "inline-block", margin: "auto 2px"}}><input min="0" max="1" step="0.01" type="number" placeholder="y" name="y" value={this.state.y} onChange={this.handleChange} required/></li>
                    <li style={{display: "inline"}}><div style={{margin: "-2px 5px"}} title="Add point" className={this.state.isControlEnabled['plus'] ? "control control-plus" : "control disabled-plus"} onClick={this.handleSubmit}></div></li>
                </ul>
            </div>
            </div>,
        ]);
    }

    
    componentDidMount() {
        this.setState({
            ghostPointsWrapper: document.getElementsByClassName('ghost-points')[0] as SVGElement,
            rectLassoWrapper: document.getElementsByClassName('rect-select')[0] as SVGElement
        });

        document.addEventListener('keydown', this.handleKeyDown);
        document.addEventListener('keyup', this.handleKeyUp);
    }

    componentDidUpdate() {
        this.userInput['x'] = this.state.x < '' ? -1 : this.state.x as number;
        this.userInput['y'] = this.state.y < '' ? -1 : this.state.y as number;
    }
    
    componentWillUnmount() {
        document.removeEventListener('keydown', this.handleKeyDown);
        document.removeEventListener('keyup', this.handleKeyUp);
    }

    private change(points: Vec2[]){
        let copyPoints = points.slice();
        this.sortPoints(copyPoints); 
        copyPoints.shift();
        copyPoints.pop();
        this.props.onChange(copyPoints);    
    }

    private handleKeyDown = (event: any) => {
        if(event.keyCode == 16) {
            this.setControls(false, false);
            this.userInput['x'] = -1;
            this.userInput['y'] = -1;
            this.setState({
                x: '',
                y: '',
                canSelectMultiple: true, 
            });
        }
    }

    private handleKeyUp = (event: any) => {
        if(event.keyCode != 16) { return; }
        this.state.ghostPointsWrapper.innerHTML = '';
        this.setState({canSelectMultiple: false, ghostPoints: [], mouseDown: false, selected: []});
    }

    private handleEvent = (event: any) => {
        if(!(event instanceof Event)) {event.persist()}
        event.stopPropagation();
        this.setState(function(prevState, props) { return getEventListener(event, prevState);}, () => {
                this.change(this.state.points);
                if(this.state.displayPoint[0] != -1 && this.state.displayPoint[1] != -1){
                    this.props.onDrag(this.state.displayPoint);
                }
            });
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
        if(event.target.value === '') { 
            this.userInput[event.target.name] = -1;
            this.setControls(false, this.state.isControlEnabled['minus']);
            return -1;
        }

        if(this.state.selected[0]) {
            this.handlePointUpdate();
            return 0;
        }

        if(this.userInput['x'] > -1 && this.userInput['y'] > -1) {
            this.setControls(true, this.state.isControlEnabled['minus']);
            this.setState({[event.target.name]: event.target.value});
            return 0;
        } else {
            this.setState({[event.target.name]: event.target.value});
            return 0; 
        }
    }

    private addPoint(callBack?: any) {
        const point = this.data;
        let points = this.state.points;
        points.push(point);
        points = this.sortPoints(points);
        this.change(points);
        this.setState({points}, callBack);
    }

    private handlePointUpdate() {
        let x: number = this.userInput['x'];
        let y: number = this.userInput['y'];
        if(x > 1 || y > 1) { return -1; }
        let id = this.state.selected[0];
        let unNormalizePoint;
        const updatePoint = normalize(this.state.height, this.state.width, Vec2.create(x, y), this.state.padding);
        
        this.state.ghostPoints[0].element.setAttribute('cx', `${updatePoint[0]}`);
        this.state.ghostPoints[0].element.setAttribute('cy', `${updatePoint[1]}`);
        x = parseInt(this.state.ghostPoints[0].element.getAttribute('cx') as string); // This, for some reason, was the only way to make
        y = parseInt(this.state.ghostPoints[0].element.getAttribute('cy') as string); // x and y be of the type 'number'
        unNormalizePoint = unNormalize(this.state.height, this.state.width, Vec2.create(x, y), this.state.padding);

        let points = this.state.points.filter((_, j) => j != id);
        points.push(unNormalizePoint);
        points = this.sortPoints(points);
        for(let i = 0; i < points.length; i++) {
            if(points[i][0] == unNormalizePoint[0] && points[i][1] == unNormalizePoint[1]) {
                id = i;
            }
        }
        this.change(points);
        this.setState({
            points,
            selected: [id], 
            'x': this.userInput['x'],
            'y': this.userInput['y']
        });
    }

    private deletePoint(event?: any, callBack?: any) {
        let points;
        const i = this.state.selected[0];
        if(i===0 || i===this.state.points.length-1){ return; }
        points = this.state.points.filter((_,j) => j !== i);
        this.state.ghostPointsWrapper.innerHTML = '';
        this.setControls(false, false);
        this.change(points);
        this.userInput['x'] = -1;
        this.userInput['y'] = -1;
        this.setState({
            selected: [],
            points,
            x: '',
            y: ''
        }, callBack);

        if(event)
            event.stopPropagation();
    }

    private handleLeave() {
        if(this.state.selected.length === 0) { return; }
        
        document.addEventListener('mousemove', this.handleEvent);
        document.addEventListener('mouseup', this.handleEvent);
    }

    private handleEnter() {
        document.removeEventListener('mousemove', this.handleEvent);
        document.removeEventListener('mouseup', this.handleEvent);
    }


    private setControls(plus: boolean, minus: boolean) {
        let isControlEnabled = {...this.state.isControlEnabled};
        isControlEnabled['plus'] = plus;
        isControlEnabled['minus'] = minus;
        this.setState({isControlEnabled});
    }

    private sortPoints(points: Vec2[]) {
        points.sort(function(a, b) { 
            if(a[0] === b[0]){
                if(a[0] === 0){
                    return a[1]-b[1];
                }
                if(a[0] === 1){
                    return 1;
                }
                return a[1]-b[1];
            }
            return a[0] - b[0];
        });  
        return points;
    }

    private refCallBack(element: any) {
        if(element){
            this.setState({myRef: element});
        }
    }

    private renderPoints() {
        const points: any[] = [];
        let point: Vec2;
        this.sortPoints(this.state.points);
        for (let i = 0; i < this.state.points.length; i++){
            if(i != 0 && i != this.state.points.length-1){
                point = normalize(this.state.height, this.state.width, this.state.points[i], this.state.padding);
                points.push(<PointComponent
                        key={i}
                        id={i}
                        x={point[0]} 
                        y={point[1]}
                        nX={this.state.points[i][0]}
                        nY={this.state.points[i][1]}
                        selected={i == this.state.selected[0] ? true : false}
                        onmouseover={this.props.onHover}
                        onmousedown={this.handleEvent}
                        onclick={this.handleEvent}
                    />);
            }
        }
        return points;
    }
}