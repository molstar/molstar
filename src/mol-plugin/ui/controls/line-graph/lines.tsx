import * as React from 'react';
import { Vec2 } from 'mol-math/linear-algebra';

export class Lines extends React.Component<any, any> {
    render(){
        const lines = this.renderLines();
        return lines;
    }

    renderLines() {
        const points: Vec2[] = [];
        let lines = [];
        let min:number;
        let maxX:number;
        let maxY: number;
        let normalizedX: number;
        let normalizedY: number;
        let reverseY: number;

        for(const point of this.props.points){
            min = this.props.padding/2;
            maxX = this.props.width+min;
            maxY = this.props.height+min; 
            normalizedX = (point[0]*(maxX-min))+min; 
            normalizedY = (point[1]*(maxY-min))+min;
            reverseY = this.props.height+this.props.padding-normalizedY;
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