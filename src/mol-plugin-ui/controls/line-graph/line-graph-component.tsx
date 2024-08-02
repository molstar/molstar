/**
 * Copyright (c) 2018-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Paul Luna <paulluna0215@gmail.com>
 * @author David Sehnal <david.sehnal@gmail.com>
 */
import { PointComponent } from './point-component';

import * as React from 'react';
import { Vec2 } from '../../../mol-math/linear-algebra';
import { Grid } from '../../../mol-model/volume';
import { arrayMax } from '../../../mol-util/array';
import { ControlGroup } from '../common';
import { CloseSvg } from '../icons';
import { ParamDefinition } from '../../../mol-util/param-definition';
import { Color } from '../../../mol-util/color';
import { CombinedColorControl } from '../color';
import { ColorNames, getRandomColor } from '../../../mol-util/color/names';
import { UUID } from '../../../mol-util';
import { ParamOnChange } from '../parameters';
import { ColorListRangesEntry } from '../../../mol-util/color/color';


export function cpsToColorListRangesEntry(controlPoints: ControlPoint[]): ColorListRangesEntry[] {
    const colors: ColorListRangesEntry[] = [];
    for (const controlPoint of controlPoints) {
        const { color, data, id } = controlPoint;
        colors.push([color, data.x, id]);
    };
    return colors;
};

export interface ControlPointData {
    x: number,
    alpha: number
}

export interface ControlPoint {
    id: UUID,
    color: Color,
    data: ControlPointData
    index: number
}

interface LineGraphComponentState {
    points: ControlPoint[],
    copyPoint: any,
    canSelectMultiple: boolean,
    showColorPicker: boolean,
    colored: boolean,
    clickedPointId?: UUID
}

const startEndPoints = [
    {
        data: {
            x: 0,
            alpha: 0
        },
        id: UUID.create22(),
        color: ColorNames.black,
        index: 0
    },
    {
        data: {
            x: 1,
            alpha: 0
        },
        id: UUID.create22(),
        color: ColorNames.black,
        index: 9
    }
];

function ColorPicker(props: any) {
    const isActive = props.isActive;
    const defaultColor = props.defaultColor;
    const color = props.color;
    const updateColor = props.updateColor;
    const toggleColorPicker = props.toggleColorPicker;
    return (isActive ? <div style={{ marginBottom: '6px', marginTop: 1 }} className='msp-accent-offset'>
        <ControlGroup header='Select Color' initialExpanded={true} hideExpander={true} hideOffset={true} onHeaderClick={toggleColorPicker}
            topRightIcon={CloseSvg} noTopMargin childrenClassName='msp-viewport-controls-panel-controls'>
            <CombinedColorControl param={defaultColor} value={color} onChange={updateColor} name='color' hideNameRow />
        </ControlGroup>
    </div> : null);
}

export class LineGraphComponent extends React.Component<any, LineGraphComponentState> {
    private myRef: any;
    private height: number;
    private width: number;
    private padding: number;
    private updatedX: number;
    private updatedY: number;
    private selectedPointId?: UUID;
    private ghostPoints: SVGElement[];
    private gElement: SVGElement;
    private namespace: string;
    constructor(props: any) {
        super(props);
        this.myRef = React.createRef();
        this.state = {
            points: startEndPoints.concat(this.props.data),
            copyPoint: undefined,
            canSelectMultiple: false,
            showColorPicker: false,
            colored: this.props.colored
        };
        this.height = 400;
        this.width = 600;
        this.padding = 70;
        this.selectedPointId = undefined;
        this.ghostPoints = [];
        this.namespace = 'http://www.w3.org/2000/svg';

        this.sortPoints(this.state.points);

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
    private getPoint(id: UUID) {
        const points = this.state.points;
        const point = points.find(p => p.id === id);
        if (!point) throw Error(`Point with id ${id} does not exist`);
        return point;
    }

    public render() {
        const points = this.renderPoints();
        const lines = this.renderLines();
        const histogram = this.renderHistogram();
        const color = this.state.clickedPointId ? this.getPoint(this.state.clickedPointId).color : void 0;
        const defaultColor = ParamDefinition.Color(Color(0x121212));
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
                    {/* renders points */}
                    <g stroke="black" fill="black">
                        {histogram}
                        {lines}
                        {points}
                    </g>
                    <g className="ghost-points" stroke="black" fill="black">
                    </g>
                </svg>
                <>
                    <ColorPicker isActive={this.state.showColorPicker} defaultColor={defaultColor} color={color} updateColor={this.updateColor} toggleColorPicker={this.toggleColorPicker}/>
                </>
            </div>,
            <div key="modal" id="modal-root" />
        ]);
    }

    toggleColorPicker = () => {
        this.setState({ showColorPicker: this.state.showColorPicker === true ? false : true });
    };

    componentDidMount() {
        this.gElement = document.getElementsByClassName('ghost-points')[0] as SVGElement;
        this.setState({
            colored: this.props.colored
        });
    }

    private change(points: ControlPoint[]) {
        const copyPoints = points.slice();
        copyPoints.shift();
        copyPoints.pop();
        this.props.onChange(copyPoints);
    }

    private updateColor: ParamOnChange = ({ value }: { value: Color }) => {
        if (!this.state.clickedPointId) throw Error('No point is selected');
        const point = this.getPoint(this.state.clickedPointId);
        if (!point) throw Error('Point should be selected');
        point.color = value;
        const points = this.state.points.map(p => p.id === this.state.clickedPointId ? point : p);
        this.setState({
            points: points
        });
        this.change(points);
    };

    private handleKeyDown = (event: any) => {
        // TODO: set canSelectMultiple = true
    };

    private handleKeyUp = (event: any) => {
        // TODO: SET canSelectMultiple = fasle
    };

    setColored() {
        this.setState({ colored: true });
    }

    private handleClick = (point: ControlPoint) => (event: any) => {
        this.setState({ clickedPointId: point.id });
        if (this.state.colored) {
            if (this.state.showColorPicker) {
                // TODO: what should happen there?
            } else {
                this.toggleColorPicker();
            }

        }
    };

    private sortPoints(points: ControlPoint[]) {
        points.sort((a, b) => {
            if (a.data.x === b.data.x) {
                if (a.data.x === 0) {
                    return a.data.alpha - b.data.alpha;
                }
                if (a.data.alpha === 1) {
                    return b.data.alpha - a.data.alpha;
                }
                return a.data.alpha - b.data.alpha;
            }
            return a.data.x - b.data.x;
        });
        return points;
    }

    private handleChangePoints(points: ControlPoint[]) {
        const pointsSorted = this.sortPoints(points);
        this.setState({ points: pointsSorted });
        this.change(points);
    }

    private handleMouseDown = (point: ControlPoint) => (event: any) => {
        const { id, index } = point;
        if (index === 0 || index === this.state.points.length - 1) {
            return;
        }
        if (this.state.canSelectMultiple) {
            return;
        }

        const copyPoint: Vec2 = this.normalizePoint(this.getPoint(id).data);
        this.ghostPoints.push(document.createElementNS(this.namespace, 'circle') as SVGElement);
        this.ghostPoints[0].setAttribute('r', '10');
        this.ghostPoints[0].setAttribute('fill', 'orange');
        this.ghostPoints[0].setAttribute('cx', `${copyPoint[0]}`);
        this.ghostPoints[0].setAttribute('cy', `${copyPoint[1]}`);
        this.ghostPoints[0].setAttribute('style', 'display: none');
        this.gElement.appendChild(this.ghostPoints[0]);
        this.updatedX = copyPoint[0];
        this.updatedY = copyPoint[1];
        this.selectedPointId = point.id;
    };

    private handleDrag(event: any) {
        if (this.selectedPointId === undefined) {
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
        } else if (svgP.x > (this.width + (padding))) {
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
    };

    private addPoint(point: ControlPoint) {
        this.state.points.push(point);
        this.handleChangePoints(this.state.points);
    }

    private replacePoint(point: ControlPoint) {
        // can get id from it
        let points = this.state.points.filter(p => p.id !== point.id);
        points.push(point);
        points = this.sortPoints(points);
        this.setState({
            points,
        });
        this.change(points);
    }

    private handlePointUpdate(event: any) {
        const selected = this.selectedPointId;
        if (this.state.canSelectMultiple) {
            return;
        }
        if (selected === undefined || this.getPoint(selected).index === 0 || this.getPoint(selected).index === this.state.points.length - 1) {
            this.setState({
                copyPoint: undefined,
            });
            return;
        }
        // "unselects" points
        this.selectedPointId = undefined;
        const pointData = this.getPoint(selected);
        const updatedPointData = this.unNormalizePoint(Vec2.create(this.updatedX, this.updatedY));
        const updatedPoint: ControlPoint = {
            data: updatedPointData,
            id: pointData.id,
            color: pointData.color,
            index: pointData.index
        };

        this.replacePoint(updatedPoint);
        this.gElement.innerHTML = '';
        this.ghostPoints = [];
        document.removeEventListener('mousemove', this.handleDrag, true);
        document.removeEventListener('mouseup', this.handlePointUpdate, true);
    }

    private handleDoubleClick(event: any) {
        const pt = this.myRef.createSVGPoint();
        pt.x = event.clientX;
        pt.y = event.clientY;
        const svgP = pt.matrixTransform(this.myRef.getScreenCTM().inverse());
        const points = this.state.points;
        const padding = this.padding / 2;

        if (svgP.x < (padding) ||
            svgP.x > (this.width + (padding)) ||
            svgP.y > (this.height + (padding)) ||
            svgP.y < (this.padding / 2)) {
            return;
        }
        const newPointData = this.unNormalizePoint(Vec2.create(svgP.x, svgP.y));
        // problem with index
        // should be the last one
        // find the last index and do + 1?
        const newPoint: ControlPoint = {
            data: newPointData,
            id: UUID.create22(),
            color: getRandomColor(),
            index: points.slice(-1)[0].index + 1
        };
        this.addPoint(newPoint);
    }

    private deletePoint = (id: UUID) => (event: any) => {
        const point = this.getPoint(id);
        if (point.index === 0 || point.index === this.state.points.length - 1) { return; }
        let points = this.state.points.filter(p => p.id !== point.id);
        points = this.sortPoints(points);
        this.setState({ points: points, clickedPointId: undefined, showColorPicker: false });
        this.change(points);
        console.log('Point', id, ' is deleted');
        event.stopPropagation();
    };

    private handleLeave() {
        if (this.selectedPointId === undefined) {
            return;
        }

        document.addEventListener('mousemove', this.handleDrag, true);
        document.addEventListener('mouseup', this.handlePointUpdate, true);
    }

    private handleEnter() {
        document.removeEventListener('mousemove', this.handleDrag, true);
        document.removeEventListener('mouseup', this.handlePointUpdate, true);
    }

    private normalizePoint(controlPointData: ControlPointData): Vec2 {
        const offset = this.padding / 2;
        const maxX = this.width + offset;
        const maxY = this.height + offset;
        const normalizedX = (controlPointData.x * (maxX - offset)) + offset;
        const normalizedY = (controlPointData.alpha * (maxY - offset)) + offset;
        const reverseY = (this.height + this.padding) - normalizedY;
        const newPoint = Vec2.create(normalizedX, reverseY);
        return newPoint;
    }

    private unNormalizePoint(vec2: Vec2): ControlPointData {
        // creates x and alpha from cartisian
        const min = this.padding / 2;
        const maxX = this.width + min;
        const maxY = this.height + min;
        const unNormalizedX = (vec2[0] - min) / (maxX - min);

        // we have to take into account that we reversed y when we first normalized it.
        const unNormalizedY = ((this.height + this.padding) - vec2[1] - min) / (maxY - min);

        return {
            x: unNormalizedX,
            alpha: unNormalizedY
        };
    }

    private refCallBack(element: any) {
        if (element) {
            this.myRef = element;
        }
    }

    private renderHistogram() {
        if (!this.props.volume) return null;
        const histogram = Grid.getHistogram(this.props.volume.grid, 40);
        const bars = [];
        const N = histogram.counts.length;
        const w = this.width / N;
        const offset = this.padding / 2;
        const max = arrayMax(histogram.counts) || 1;
        for (let i = 0; i < N; i++) {
            const x = this.width * i / (N - 1) + offset;
            const y1 = this.height + offset;
            const y2 = this.height * (1 - histogram.counts[i] / max) + offset;
            bars.push(<line key={`histogram${i}`} x1={x} x2={x} y1={y1} y2={y2} stroke="#ded9ca" strokeWidth={w} />);
        }
        return bars;
    }

    private renderPoints() {
        const points: any[] = [];
        let point: Vec2;
        debugger;
        for (let i = 0; i < this.state.points.length; i++) {
            if (i !== 0 && i !== this.state.points.length - 1) {
                const { data, color, id, index } = this.state.points[i];
                let finalColor = color;
                if (!this.state.colored) {
                    finalColor = ColorNames.black;
                };
                point = this.normalizePoint(data);
                points.push(<PointComponent
                    index={index}
                    key={id}
                    id={id}
                    x={point[0]}
                    y={point[1]}
                    nX={data.x}
                    nY={data.alpha}
                    selected={false}
                    // should be able to delete point
                    delete={this.deletePoint}
                    // onmouseover we provide props.onHover of line graph component
                    onmouseover={this.props.onHover}
                    onmousedown={this.handleMouseDown(this.state.points[i])}
                    onclick={this.handleClick(this.state.points[i])}
                    color={finalColor}
                />);
            }
        }
        return points;
    }

    private renderLines() {
        const points: Vec2[] = [];
        const lines = [];
        let maxX: number;
        let maxY: number;
        let normalizedX: number;
        let normalizedY: number;
        let reverseY: number;

        const o = this.padding / 2;
        for (const point of this.state.points) {
            maxX = this.width + o;
            maxY = this.height + this.padding;
            normalizedX = (point.data.x * (maxX - o)) + o;
            normalizedY = (point.data.alpha * (maxY - o)) + o;
            reverseY = this.height + this.padding - normalizedY;
            points.push(Vec2.create(normalizedX, reverseY));
        }

        const data = points;
        const size = data.length;

        for (let i = 0; i < size - 1; i++) {
            const x1 = data[i][0];
            const y1 = data[i][1];
            const x2 = data[i + 1][0];
            const y2 = data[i + 1][1];
            lines.push(<line key={`lineOf${i}`} x1={x1} x2={x2} y1={y1} y2={y2} stroke="#cec9ba" strokeWidth="5"/>);
        }

        return lines;
    }
}