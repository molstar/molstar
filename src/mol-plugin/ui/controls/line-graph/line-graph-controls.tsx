// import * as React from 'react';

import { Vec2 } from 'mol-math/linear-algebra';
import { unNormalize, normalize } from './normilization';

export function getEventListener(event:any, state:any) {
    switch(event.type) {
        case 'mousedown': return mouseDown(event, state);
        case 'mouseup': return mouseUp(event, state);
        case 'click': return click(event, state);
        case 'dblclick': return dbclick(event, state);
        case 'mousemove': return mousemove(event, state);
        default:
            console.warn(`There is no function for the ${event.type} event`);
    }
}

function mouseDown (event: any, state: any ) {
    if(!event) { return state; }
    const id = parseInt(event.target.id);
    state.mouseDown = true;
    if(event.target.tagName != 'circle') {
        let svgP;
        const pt = state.myRef.createSVGPoint();
        pt.x = event.clientX;
        pt.y = event.clientY;
        svgP = pt.matrixTransform(state.myRef.getScreenCTM().inverse());
        state.mouseStartPoint = Vec2.create(svgP.x, svgP.y);
        state.canSelectMultiple = false;
        return state;
    }
    let selected;
    let size;
    let copyPoint: Vec2;
    const x = event.target.cx.animVal.value;
    const y = event.target.cy.animVal.value;
    if(isNaN(id)) {
        state.isControlEnabled['plus'] = false;
        state.isControlEnabled['minus'] = false;
        state.x = '';
        state.y = '';
        state.mouseStartPoint = Vec2.create(x, y);
        return state;
    }
    state.mouseStartPoint = Vec2.create(x, y);
    if (state.canSelectMultiple) {
        addGhostPoint(state,id, x, y);
        size = state.ghostPoints.length;
        state.ghostPointsWrapper.appendChild(state.ghostPoints[size-1].element);
        selected = state.selected;
        selected.push(id);
        state.selected = selected;
        return state;
    }
    copyPoint = Vec2.create(x, y);
    state.ghostPoints = [];
    state.ghostPointsWrapper.innerHTML = '';
    addGhostPoint(state, id, x, y);
    state.ghostPointsWrapper.appendChild(state.ghostPoints[0].element);
    state.updatedX = copyPoint[0];
    state.updatedY = copyPoint[1];
    state.selected = [id];
    state.copyPoint = copyPoint;
    return state;
}

function mousemove(event:any, state:any) {
    if(!state.mouseDown) {
        return -1;
    }
    let updatedCopyPoint: Vec2;
    let svgP;
    const pt = state.myRef.createSVGPoint();
    state.isControlEnabled['plus'] = false;
    state.isControlEnabled['minus'] = false;
    state.x = '';
    state.y = '';
    pt.x = event.clientX;
    pt.y = event.clientY;
    svgP = pt.matrixTransform(state.myRef.getScreenCTM().inverse());
    state.hasDragged = true;
    if(state.selected.length == 0 && state.mouseDown) {
        svgP = Vec2.create(svgP.x, svgP.y);
        createRect(svgP, state);
        state.createdLasso = true;
        return state;
    }

        let selected = 0;
        let ghostPoint: Vec2;
        const directionalVector = Vec2.create(svgP.x-state.mouseStartPoint[0], svgP.y-state.mouseStartPoint[1]);

        const outerCoords = {
            leftMostInt: parseInt(state.leftMost.getAttribute('cx') as string),
            rightMostInt: parseInt(state.rightMost.getAttribute('cx') as string),
            upperMostInt: parseInt(state.upperMost.getAttribute('cy') as string),
            bottomMostInt: parseInt(state.bottomMost.getAttribute('cy') as string),
        }

        if((directionalVector[0]+outerCoords.leftMostInt) <= state.padding/2) {
            directionalVector[0] = (state.padding/2)-outerCoords.leftMostInt;
        }
        
        if((directionalVector[0]+outerCoords.rightMostInt) >= state.width+(state.padding/2)) {
            directionalVector[0] = (state.width+(state.padding/2))-outerCoords.rightMostInt;
        }

        if((directionalVector[1]+outerCoords.upperMostInt) <= (state.padding/2)) {
            directionalVector[1] = (state.padding/2)-outerCoords.upperMostInt;
        }

        if((directionalVector[1]+outerCoords.bottomMostInt) >= (state.height+(state.padding/2))) {
            directionalVector[1] = (state.height+(state.padding/2))-outerCoords.bottomMostInt;
        }

        for(let i = 0; i < state.ghostPoints.length; i++) {
            const tempX = state.ghostPoints[i].element.getAttribute('cx');
            const tempY = state.ghostPoints[i].element.getAttribute('cy');
            ghostPoint = Vec2.create(parseInt(tempX as string), parseInt(tempY as string));
            if(state.mouseStartPoint[0] == ghostPoint[0]) {selected = i;}
            updatedCopyPoint = Vec2.create(ghostPoint[0]+directionalVector[0], ghostPoint[1]+directionalVector[1]);
            state.ghostPoints[i].element.setAttribute('cx', `${updatedCopyPoint[0]}`);
            state.ghostPoints[i].element.setAttribute('cy', `${updatedCopyPoint[1]}`);
            state.ghostPoints[i].element.setAttribute('style', 'style: visible'); 
            sortOuterGhostPoints(state);
        }
        
        const x = parseInt(state.ghostPoints[selected].element.getAttribute('cx') as string);
        const y = parseInt(state.ghostPoints[selected].element.getAttribute('cy') as string);
        if(state.selected.length == 1) {
            const unNormalizePoint = unNormalize(state.height, state.width, Vec2.create(x, y), state.padding);
            state.displayPoint = unNormalizePoint;
        }
        state.mouseStartPoint = Vec2.create(x, y);
        
        return state;
}

function mouseUp(event: any, state: any) {
    if((state.canSelectMultiple && !state.hasDragged) || !state.hasDragged) { return -1}
    let i;
    let points = [...state.points];
    state.mouseDown = false;
    if(state.hasDragged && state.selected.length == 0) {
        let selected = [];
        let rect = document.getElementsByClassName('rectLasso')[0];
        let leftBound = parseInt(rect.getAttribute('x') as string);
        let upperBound = parseInt(rect.getAttribute('y') as string);
        let rightBound = leftBound + parseInt(rect.getAttribute('width') as string);
        let lowerBound = upperBound + parseInt(rect.getAttribute('height') as string);
        for(i = 1; i < points.length-1; i++){
            let size;
            let point = normalize(state.height, state.width, Vec2.create(points[i][0], points[i][1]), state.padding);
            
            if(point[0] >= leftBound && point[0] <= rightBound && point[1] >= upperBound && point[1] <= lowerBound) {
                selected.push(i);
                addGhostPoint(state, i, point[0], point[1]);
                size = state.ghostPoints.length;
                state.ghostPoints[size-1].element.setAttribute('style', 'display: visible');
                state.ghostPointsWrapper.appendChild(state.ghostPoints[size-1].element);
            }
        }

        state.rectLassoWrapper.innerHTML = '';
        state.hasDragged = false;
        if(selected.length > 0) { state.canSelectMultiple = true; }
        state.selected = [...selected];
        return state;
    }
    
    for(i = 0; i < state.ghostPoints.length; i++) {
        const id = state.ghostPoints[i].id;
        const element = state.ghostPoints[i].element;
        const x = parseInt(element.getAttribute('cx') as string);
        const y = parseInt(element.getAttribute('cy') as string);
        const updatedPoint = unNormalize(state.height, state.width, Vec2.create(x, y), state.padding);
        points[id] = updatedPoint;
    }

    state.ghostPointsWrapper.innerHTML = '';
    state.ghostPoints = [];
    state.hasDragged = false;
    state.selected = [];
    state.points = [...points];
    state.leftMost.setAttribute('cx', '0');
    state.canSelectMultiple = false;
    state.displayPoint = Vec2.create(-1, -1);
    return state;
}

function click(event: any, state: any) {
    state.mouseDown = false;
    if(event.target.tagName != 'circle') {
        if(state.canSeleMultiple) { return;}

        state.isControlEnabled['plus'] = false;
        state.isControlEnabled['minus'] = false;
        if(state.createdLasso) {
            state.createdLasso = false;
            return state;
        }
        state.ghostPointsWrapper.innerHTML = '';
        state.ghostPoints = [];
        state.hasDragged = false;
        state.leftMost.setAttribute('cx', '0');
        // this.userInput['x'] = -1;
        // this.userInput['y'] = -1;
        state.selected = [];
        state.x = '';
        state.y = '';
        return state;
    }
    const id = parseInt(event.target.id);
    if(isNaN(id)) {
        for(let i=0; i<state.ghostPoints.length; i++) {
            if(state.ghostPoints[i].element == event.target) {
                const ghostPointElements = state.ghostPointsWrapper.getElementsByClassName('ghostPoint');
                state.ghostPointsWrapper.removeChild(ghostPointElements[i]);
                state.ghostPoints = state.ghostPoints.filter((_: any, j:any) => j !== i);
                return state;
            }
        }
    }

    if(state.canSelectMultiple) {
        let size = state.ghostPoints.length;
        state.ghostPoints[size-1].element.setAttribute('style', 'display: visible');
        state.ghostPointsWrapper.appendChild(state.ghostPoints[size-1].element);
        sortOuterGhostPoints(state);
        return state;
    } 
    
    // this.setControls(false, true)
    state.isControlEnabled['plus'] = false;
    state.isControlEnabled['minus'] = true;
    state.ghostPoints[0].element.setAttribute('style', 'display: visible');
    // this.userInput['x'] = this.state.points[this.state.selected[0]][0] as number;
    // this.userInput['y'] = this.state.points[this.state.selected[0]][1] as number; 
    state.x = state.points[state.selected[0]][0];
    state.y = state.points[state.selected[0]][1];
    return state;
}

function dbclick(event: any, state: any) {
    let newPoint;
    let points;
    let svgP;
    const pt = state.myRef.createSVGPoint();
    pt.x = event.clientX;
    pt.y = event.clientY;
    svgP  = pt.matrixTransform(state.myRef.getScreenCTM().inverse());
    if(svgP.x < state.padding ||
        svgP.x > (state.width+state.padding) || 
        svgP.y > (state.height+state.padding) || 
        svgP.y < state.padding/2) { return -1; }
    newPoint = unNormalize(state.height, state.width, Vec2.create(svgP.x, svgP.y), state.padding);
    points = [...state.points, newPoint];
    state.points = [...points];
    return state;
}

// TODO: Move to helper file?
function addGhostPoint(state: any, id: number, x: number, y: number) {
    state.ghostPoints.push({id: id, element: document.createElementNS('http://www.w3.org/2000/svg', 'circle') as SVGElement});
    const size = state.ghostPoints.length;
    state.ghostPoints[size-1].element.setAttribute('r', '10');
    state.ghostPoints[size-1].element.setAttribute('fill', 'orange');
    state.ghostPoints[size-1].element.setAttribute('cx', `${x}`);
    state.ghostPoints[size-1].element.setAttribute('cy', `${y}`);
    state.ghostPoints[size-1].element.setAttribute('style', 'display: none');
    state.ghostPoints[size-1].element.setAttribute('class', 'ghostPoint');
    sortOuterGhostPoints(state);
}

// TODO: Move to sorting file?
export function sortOuterGhostPoints(state: any) {
    // Need to reset the values
    state.leftMost.setAttribute('cx', `${state.width+(state.padding/2)}`);
    state.rightMost.setAttribute('cx', `0`);
    state.upperMost.setAttribute('cy', `${state.height+(state.padding/2)}`);
    state.bottomMost.setAttribute('cy', `0`);
    
    for(let i=0; i < state.ghostPoints.length; i++) {
        let element = state.ghostPoints[i].element;
        const x = parseInt(element.getAttribute('cx') as string);
        const y = parseInt(element.getAttribute('cy') as string);
        
        if(x < parseInt(state.leftMost.getAttribute('cx') as string)) {
            state.leftMost.setAttribute('cx', `${element.getAttribute('cx')}`);
        }

        if(x > parseInt(state.rightMost.getAttribute('cx') as string)) {
            state.rightMost.setAttribute('cx', `${element.getAttribute('cx')}`);
        }

        if(y > parseInt(state.bottomMost.getAttribute('cy') as string)) {
            state.bottomMost.setAttribute('cy', `${element.getAttribute('cy')}`);
        }

        if(y < parseInt(state.upperMost.getAttribute('cy') as string)) {
            state.upperMost.setAttribute('cy', `${element.getAttribute('cy')}`);
        }
    }
}

function createRect(endPoint: Vec2, state: any) {
    const width = Math.abs(endPoint[0] - state.mouseStartPoint[0]);
    const height = Math.abs(endPoint[1] - state.mouseStartPoint[1]);
    let rectLasso = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
    state.rectLassoWrapper.innerHTML = '';
    if(endPoint[0] > state.mouseStartPoint[0] && endPoint[1] > state.mouseStartPoint[1]) {
        rectLasso.setAttribute('x', `${state.mouseStartPoint[0]}`);
        rectLasso.setAttribute('y', `${state.mouseStartPoint[1]}`);
    } else if( endPoint[0] < state.mouseStartPoint[0] && endPoint[1] > state.mouseStartPoint[1]) {
        rectLasso.setAttribute('x', `${state.mouseStartPoint[0]-width}`);
        rectLasso.setAttribute('y', `${state.mouseStartPoint[1]}`);
    } else if ( endPoint[0] > state.mouseStartPoint[0] && endPoint[1] < state.mouseStartPoint[1]){
        rectLasso.setAttribute('x', `${state.mouseStartPoint[0]}`);
        rectLasso.setAttribute('y', `${state.mouseStartPoint[1] - height}`);
    } else {
        rectLasso.setAttribute('x', `${state.mouseStartPoint[0] - width}`);
        rectLasso.setAttribute('y', `${state.mouseStartPoint[1] - height}`);
    }
    
    rectLasso.setAttribute('width', `${width}`);
    rectLasso.setAttribute('height', `${height}`);
    rectLasso.setAttribute('fill', '#4c66b2');
    rectLasso.setAttribute('fill-opacity', '0.25');
    rectLasso.setAttribute('stroke', 'blue');
    rectLasso.setAttribute('class', 'rectLasso');
    state.rectLassoWrapper.append(rectLasso);
}