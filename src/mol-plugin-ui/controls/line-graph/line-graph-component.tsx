/**
 * Copyright (c) 2018-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Paul Luna <paulluna0215@gmail.com>
 * @author David Sehnal <david.sehnal@gmail.com>
 */
import { PointComponent } from './point-component';

import * as React from 'react';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { Vec2 } from '../../../mol-math/linear-algebra';
import { Grid } from '../../../mol-model/volume';
import { arrayMax, arrayMin } from '../../../mol-util/array';
import { Button, ControlGroup } from '../common';
import { CloseSvg } from '../icons';
import { ParamDefinition } from '../../../mol-util/param-definition';
import { Color } from '../../../mol-util/color';
import { CombinedColorControl } from '../color';
import { ColorNames, getRandomColor } from '../../../mol-util/color/names';
import { UUID } from '../../../mol-util';
import { ParameterControls, ParamOnChange } from '../parameters';
import { ColorListRangesEntry } from '../../../mol-util/color/color';
import { generateGaussianControlPoints } from '../../../mol-geo/geometry/direct-volume/direct-volume';
import { capitalize } from '../../../mol-util/string';
// import { WaitingParameterControls } from '../../../extensions/volumes-and-segmentations/ui';
import { sleep } from '../../../mol-util/sleep';
import { useCallback, useEffect, useRef, useState } from 'react';

type ComponentParams<T extends React.Component<any, any, any> | ((props: any) => JSX.Element)> =
    T extends React.Component<infer P, any, any> ? P : T extends (props: infer P) => JSX.Element ? P : never;



function WaitingParameterControls<T extends PD.Params>({ values, onChangeValues, ...etc }: { values: PD.ValuesFor<T>, onChangeValues: (values: PD.ValuesFor<T>) => any } & ComponentParams<ParameterControls<T>>) {
    const [changing, currentValues, execute] = useAsyncChange(values);

    return <ParameterControls isDisabled={changing} values={currentValues} onChangeValues={newValue => execute(onChangeValues, newValue)} {...etc} />;
}


function useAsyncChange<T>(initialValue: T) {
    const [isExecuting, setIsExecuting] = useState(false);
    const [value, setValue] = useState(initialValue);
    const isMounted = useRef(false);

    useEffect(() => setValue(initialValue), [initialValue]);

    useEffect(() => {
        isMounted.current = true;
        return () => { isMounted.current = false; };
    }, []);

    const execute = useCallback(
        async (func: (val: T) => Promise<any>, val: T) => {
            setIsExecuting(true);
            setValue(val);
            try {
                await func(val);
            } catch (err) {
                if (isMounted.current) {
                    setValue(initialValue);
                }
                throw err;
            } finally {
                if (isMounted.current) {
                    setIsExecuting(false);
                }
            }
        },
        []
    );

    return [isExecuting, value, execute] as const;
}

class TFParamsWrapper extends React.Component<any> {
    // should accept params or rather use params from control
    state = {
        gaussianTFParamsValues: {
            gaussianCenter: 1.0,
            gaussianExtent: 1.0
        }
    };

    handleChange = (next: GaussianTFParamsValues) => {
        this.props.onChange('gaussian', next.gaussianExtent, next.gaussianCenter);

    };

    // check if need to hanlde state change
    render() {
        return <WaitingParameterControls params={GaussianTFParams} values={this.state.gaussianTFParamsValues} onChangeValues={async next => { await sleep(20); console.log(next); this.handleChange(next); }} />;

        // return (
        //     <Button onClick={this.handleChange}>{`Apply ${capitalize(this.props.kind)} Transfer Function`}</Button>
        // );
    }
}



class TFButton extends React.Component<any> {
    handleClick = () => {
        this.props.onClick(this.props.kind, this.props.sigmaMultiplierExtent, this.props.sigmaMultiplierCenter);
    };

    render() {
        return (
            <Button onClick={this.handleClick}>{`Apply ${capitalize(this.props.kind)} Transfer Function`}</Button>
        );
    }
}


export function areControlPointsColorsSame(newControlPoints: ControlPoint[], oldControlPoints: ControlPoint[]) {
    const newSorted = newControlPoints.sort((a, b) => a.index - b.index);
    const oldSorted = oldControlPoints.sort((a, b) => a.index - b.index);
    debugger;
    for (let i = 0, il = newSorted.length; i < il; ++i) {
        const newPoint = newSorted[i];
        const oldPoint = oldSorted[i];
        if (newPoint.color !== oldPoint.color) {
            return false;
        }
    }
};

export function controlPointsToColorListControlPointsEntry(controlPoints: ControlPoint[]): ColorListRangesEntry[] {
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

interface VolumeDescriptiveStatistics {
    min: number,
    max: number,
    mean: number,
    sigma: number
}

interface LineGraphComponentState {
    points: ControlPoint[],
    copyPoint: any,
    canSelectMultiple: boolean,
    showColorPicker: boolean,
    // colored: boolean,
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

export const GaussianTFParams = {
    gaussianCenter: PD.Numeric(0.2, { min: 0, max: 3, step: 0.05 }),
    gaussianExtent: PD.Numeric(0.2, { min: 0, max: 3, step: 0.05 })
};

export type GaussianTFParamsValues = PD.Values<typeof GaussianTFParams>;

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
    private descriptiveStatistics: VolumeDescriptiveStatistics;
    constructor(props: any) {
        super(props);
        this.myRef = React.createRef();
        this.state = {
            points: startEndPoints.concat(this.props.data),
            copyPoint: undefined,
            canSelectMultiple: false,
            showColorPicker: false,
            // colored: this.props.colored
        };
        this.height = 400;
        this.width = 600;
        this.padding = 70;
        this.selectedPointId = undefined;
        this.ghostPoints = [];
        this.namespace = 'http://www.w3.org/2000/svg';
        this.descriptiveStatistics = this.getDescriptiveStatistics();

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
        this.setPredefinedTransferFunction = this.setPredefinedTransferFunction.bind(this);
    }

    private getDescriptiveStatistics() {
        const s = this.props.volume.grid.stats;
        const d: VolumeDescriptiveStatistics = {
            min: s.min,
            max: s.max,
            mean: s.mean,
            sigma: s.sigma
        };
        return d;
    }

    private setPredefinedTransferFunction(type: 'gaussian', sigmaMultiplierExtent: number, sigmaMultiplierCenter: number) {
        const a = 0.2;
        const mean = this.descriptiveStatistics.mean;
        const min = this.descriptiveStatistics.min;
        const max = this.descriptiveStatistics.max;
        const sigma = this.descriptiveStatistics.sigma;
        const TFextent = max - min;
        // now it sort of works
        console.log(this.descriptiveStatistics);
        const b = sigmaMultiplierCenter * (mean + sigma) / TFextent;
        // fix this, should be just sigma,
        const c = sigmaMultiplierExtent * sigma / TFextent;
        // const l = (this.width * 2 * c / extent);
        debugger;
        switch (type) {
            case 'gaussian':
                const gaussianPoints: ControlPoint[] = generateGaussianControlPoints(a, b, c, TFextent);
                const currentPoints = this.state.points;
                gaussianPoints.push(currentPoints[currentPoints.length - 1]);
                gaussianPoints.unshift(currentPoints[0]);
                this.setState({
                    points: gaussianPoints
                });
                this.change(gaussianPoints);
                break;
            default: throw Error(`Transfer function type ${type} is not supported`);
        };
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
        const descriptiveStatisticsBars = this.renderDescriptiveStatisticsBars();
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
                        {descriptiveStatisticsBars}
                    </g>
                    <g className="ghost-points" stroke="black" fill="black">
                    </g>
                </svg>
                <>
                    <ColorPicker isActive={this.state.showColorPicker} defaultColor={defaultColor} color={color} updateColor={this.updateColor} toggleColorPicker={this.toggleColorPicker}/>
                </>
                <>
                    {/* can be select instead, then on select etc. */}
                    {/* change data in the UI somehow */}
                    {/* change gaussian on change of this */}
                    {/* make wrapper component similar to TFButton */}
                    {/* <WaitingParameterControls params={GaussianTFParams} values={GaussianTFParamsValues} onChangeValues={async next => { await sleep(20); console.log('stuff'); }} /> */}
                    <TFParamsWrapper onChange={this.setPredefinedTransferFunction}></TFParamsWrapper>
                    <TFButton onClick={this.setPredefinedTransferFunction} kind={'gaussian'} sigmaMultiplierExtent={0.25} sigmaMultiplierCenter={1.5}></TFButton>
                    {/* <Button onClick={() => this.setPredefinedTransferFunction('gaussian')}>Apply Gaussian Transfer Function</Button> */}
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
        // this.setState({
        //     colored: this.props.colored
        // });
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

    private handleClick = (point: ControlPoint) => (event: any) => {
        this.setState({ clickedPointId: point.id });
        if (this.state.showColorPicker) {
            // TODO: what should happen there?
        } else {
            this.toggleColorPicker();
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
        const data: number[] = this.props.volume.grid.cells.data;
        const bins = [];
        const increment = data.length / N;
        data.sort((a, b) => a - b);
        for (let i = 0; i < N; ++i) {
            // sort
            const chunk = data.slice(i * increment, (i + 1) * increment);
            bins.push(chunk);
        };
        debugger;
        for (let i = 0; i < N; i++) {
            const fromValue = arrayMin(bins[i]);
            const toValue = arrayMax(bins[i]);
            const x = this.width * i / (N - 1) + offset;
            const y1 = this.height + offset;
            const y2 = this.height * (1 - histogram.counts[i] / max) + offset;
            bars.push(<line key={`histogram${i}`} x1={x} x2={x} y1={y1} y2={y2} stroke="#ded9ca" strokeWidth={w}>
                <title>[{fromValue}; {toValue}]</title>
            </line>);
        }
        return bars;
    }

    private renderDescriptiveStatisticsBars() {
        if (!this.props.volume) return null;
        const offset = this.padding / 2;
        const mean = this.descriptiveStatistics.mean;
        const min = this.descriptiveStatistics.min;
        const max = this.descriptiveStatistics.max;
        const sigma = this.descriptiveStatistics.sigma;
        const extent = max - min;
        const x = this.width * (mean / extent);
        const w = offset / 5;
        const bars = [];
        const y1 = this.height + offset;
        const y2 = offset;
        const xPositive = this.width * ((mean + sigma) / extent);
        const xNegative = this.width * ((mean - sigma) / extent);
        bars.push(
            <line key={'meanBar'} x1={x} x2={x} y1={y1} y2={y2} stroke="#000000" strokeWidth={w}>
                <title>`Mean: ${mean}`</title>
            </line>);
        bars.push(<line key={'positiveSigmaBar'} x1={xPositive} x2={xPositive} y1={y1} y2={y2} stroke="#808080" strokeWidth={w} />);
        bars.push(<line key={'negativeSigmaBar'} x1={xNegative} x2={xNegative} y1={y1} y2={y2} stroke="#808080" strokeWidth={w} />);
        return bars;
    }

    private renderPoints() {
        const points: any[] = [];
        let point: Vec2;
        for (let i = 0; i < this.state.points.length; i++) {
            if (i !== 0 && i !== this.state.points.length - 1) {
                const { data, color, id, index } = this.state.points[i];
                const finalColor = color;
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