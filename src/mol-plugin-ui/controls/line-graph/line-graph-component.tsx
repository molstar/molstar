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
import { arrayMax } from '../../../mol-util/array';
import { Button, ControlGroup, ExpandGroup, IconButton, TextInput } from '../common';
import { CloseSvg, DeleteSvg, MinusBoxSvg, PlusBoxSvg } from '../icons';
import { ParamDefinition } from '../../../mol-util/param-definition';
import { Color } from '../../../mol-util/color';
import { CombinedColorControl } from '../color';
import { ColorNames, getRandomColor } from '../../../mol-util/color/names';
import { UUID } from '../../../mol-util';
import { ParameterControls, ParamOnChange } from '../parameters';
import { ColorListRangesEntry } from '../../../mol-util/color/color';
import { generateGaussianControlPoints } from '../../../mol-geo/geometry/direct-volume/direct-volume';
import { useCallback, useEffect, useRef, useState } from 'react';

type ComponentParams<T extends React.Component<any, any, any> | ((props: any) => JSX.Element)> =
    T extends React.Component<infer P, any, any> ? P : T extends (props: infer P) => JSX.Element ? P : never;

class PointButton extends React.Component<any> {
    onAbs = (v: number) => {
        // this changes that point in the state of linegraph component.points
        this.props.changeXValue(this.props.point.id, v);
        // this.props.value = v;
    };

    onAlpha = (v: number) => {
        this.props.changeAlphaValue(this.props.point.id, v);
    };

    handleClick = () => {
        this.props.onClick(this.props.point.id);
    };
    render() {
        debugger;
        const [absValue, relativeValue] = this.props.onExpandGroupOpen(this.props.point.data);
        debugger;
        console.log(absValue, relativeValue);
        const alpha = (this.props.point.data.alpha as number).toFixed(3);
        return (
            <div style={{ display: 'flex', marginBottom: 1 }} key={this.props.point.id}>
                <Button style={{ textAlign: 'start', textIndent: '20px' }}>Point</Button>
                <TextInput numeric
                    style={{ minWidth: 0 }} className='msp-form-control' onEnter={this.props.onEnter} blurOnEnter={true} blurOnEscape={true}
                    value={absValue} placeholder={'Some text'}
                    // not directly
                    // first onAbs, inside onAbs do call to this.props.changeXValue(this.props.point.id, value)
                    // and assign that value to one of the props (value)
                    isDisabled={false} onChange={(value) => { this.onAbs(value); }} />
                <TextInput numeric
                    style={{ minWidth: 0 }} className='msp-form-control' onEnter={this.props.onEnter} blurOnEnter={true} blurOnEscape={true}
                    value={alpha} placeholder={'Some text'}
                    isDisabled={false} onChange={(value) => { this.onAlpha(value as any); }} />

                <IconButton title={'Remove point'} svg={DeleteSvg} onClick={this.handleClick}></IconButton>
            </div>
        );
    }
}

function adjustTFParams(name: TFName, ds: VolumeDescriptiveStatistics) {
    switch (name) {
        // TODO: type for that
        case 'gaussian':
            const max = ds.max;
            const min = ds.min;
            // TODO: copy?
            GaussianTFParams.gaussianCenter.max = max;
            GaussianTFParams.gaussianCenter.min = min;
            GaussianTFParams.gaussianExtent.max = max / 2;
            GaussianTFParams.gaussianCenter.min = 0;
            return GaussianTFParams;
        case 'method2':
            return Method2TFParams;
        default: throw Error(`Transfer function ${name} is not supported`);
    }
}

class BaseLine extends React.Component<any> {
    // handleClick = () => {
    //     // this.props.onClick(this.props.kind, this.props.sigmaMultiplierExtent, this.props.sigmaMultiplierCenter);
    // };

    render() {
        const y = this.props.height;
        const x1 = this.props.offset;
        const x2 = this.props.offset + this.props.width;
        return (
            <>
                <line y1={y} y2={y} x1={x1} x2={x2} stroke='black' fill='black'></line>
            </>
        );
    }
}

class PointsPanel extends React.Component<any> {
    // handleClick = () => {
    //     // this.props.onClick(this.props.kind, this.props.sigmaMultiplierExtent, this.props.sigmaMultiplierCenter);
    // };

    render() {
        const points: ControlPoint[] = this.props.points;
        const realPoints = points.filter(p => p.isGhost !== true);
        debugger;
        // const realPoints = points.filter();
        // TODO: add ghost prop to points
        // const removeAllPointsButton = <Button style={{ position: 'absolute', top: 0, right: 0 }} onClick={this.props.removeAllPoints}>Remove All Points</Button>;
        const removeAllPointsButton = <IconButton style={{ position: 'absolute', top: 0, right: 0 }} svg={DeleteSvg} small onClick={() => {
            this.props.removeAllPoints();
        } }></IconButton>;
        const controlPointsButtons = realPoints.map(p => {
            return <PointButton key={p.id} point={p}
                onClick={this.props.onPointButtonClick}
                onExpandGroupOpen={this.props.onExpandGroupOpen}
                changeXValue={this.props.changeXValue}
                changeAlphaValue={this.props.changeAlphaValue}
            ></PointButton>;
        });
        return (
            <div style={{ position: 'relative' }}>
                <ExpandGroup header='Control Points Panel' initiallyExpanded={false}>
                    {controlPointsButtons}
                </ExpandGroup>
                {removeAllPointsButton}
                {/* TODO: flex? */}
                {/* <IconButton small onClick={() => {
                this.props.removeAllPoints();
            } }></IconButton> */}
            </div>
        );
    }
}

class TFMethodPanel extends React.Component<any> {
    // handleClick = () => {
    //     // this.props.onClick(this.props.kind, this.props.sigmaMultiplierExtent, this.props.sigmaMultiplierCenter);
    // };

    render() {
        return <ExpandGroup header={this.props.params.name} initiallyExpanded={true}>
            {/* render params */}
            <TFParamsWrapper onChange={this.props.onChange} params={this.props.params} descriptiveStatistics={this.props.descriptiveStatistics}></TFParamsWrapper>
        </ExpandGroup>;
    }
}

export type TFName = 'gaussian' | 'method2';

export interface TFMethod {
    name: TFName
    // TODO: find params
    params: TFParamsValues
}

class HelpersPanel extends React.Component<any> {
    // handleClick = () => {
    //     // this.props.onClick(this.props.kind, this.props.sigmaMultiplierExtent, this.props.sigmaMultiplierCenter);
    // };

    render() {
        const methods: TFMethod[] = this.props.methods;
        // const points: ControlPoint[] = this.props.points;
        // methods undefined
        const methodsUI = methods.map(p => {
            // should be expandgroup
            // pass on change method around
            return <TFMethodPanel key={p.name} onChange={this.props.onChange} params={p} descriptiveStatistics={this.props.descriptiveStatistics}></TFMethodPanel>;
        });
        return (
            <ExpandGroup header='Helpers' initiallyExpanded={false}>
                {methodsUI}
            </ExpandGroup>
        );
    }
}

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

// TODO: add explicit state type to it
class TFParamsWrapper extends React.Component<any> {
    state = { params: this.props.params };

    // TODO: check if the above is correct
    handleChange = (next: TFParamsValues) => {
        // on change should be generic as well
        this.props.onChange(next);
    };

    // rework for generic tf

    handleClick = () => {
        this.props.onChange(this.props.params);
        this.setState({ params: this.props.params });
    };

    render() {
        const adjustedParams = adjustTFParams(this.props.params.name, this.props.descriptiveStatistics);
        return (<ExpandGroup header='Transfer Function Settings' initiallyExpanded>
            {/* TODO: fix that as any */}
            <WaitingParameterControls params={adjustedParams} values={this.state.params} onChangeValues={async next => { this.handleChange(next as any); }} />
            <Button onClick={this.handleClick}>{`Apply ${this.props.params.name} Transfer Function`}</Button>
        </ExpandGroup>);
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
    absValue?: number
    relativeValue?: number
}

export interface ControlPoint {
    id: UUID,
    color: Color,
    data: ControlPointData
    index: number,
    isGhost?: boolean
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
    clickedPointIds?: UUID[],
    methodsParams: TFParamsValues[]
}

const startEndPoints: ControlPoint[] = [
    {
        data: {
            x: 0,
            alpha: 0
        },
        id: UUID.create22(),
        color: ColorNames.black,
        index: 0,
        isGhost: true
    },
    {
        data: {
            x: 1,
            alpha: 0
        },
        id: UUID.create22(),
        color: ColorNames.black,
        index: 9,
        isGhost: true
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
    gaussianCenter: PD.Numeric(0.2, { min: 0, max: 1, step: 0.01 }),
    gaussianExtent: PD.Numeric(0.2, { min: 0, max: 1, step: 0.01 }),
    gaussianHeight: PD.Numeric(0.2, { min: 0, max: 1, step: 0.01 }),
    // TODO: PD from object or something
    name: PD.Text('gaussian', { isHidden: true })
};

export const Method2TFParams = {
    param1: PD.Numeric(0.2, { min: 0, max: 1, step: 0.01 }),
    param2: PD.Numeric(0.2, { min: 0, max: 1, step: 0.01 }),
    // TODO: PD from object or something
    name: PD.Text('method2', { isHidden: true })
    // gaussianHeight: PD.Numeric(0.2, { min: 0, max: 1, step: 0.01 })
};

export type GaussianTFParamsValues = PD.Values<typeof GaussianTFParams>;

export type TFParamsValues = GaussianTFParamsValues | Method2ParamsValues;

export type Method2ParamsValues = PD.Values<typeof Method2TFParams>;

const DefaultGaussianParams: GaussianTFParamsValues = {
    gaussianCenter: 1.0,
    gaussianExtent: 0.25,
    gaussianHeight: 0.2,
    name: 'gaussian'
};

const DefaultMethod2Params: Method2ParamsValues = {
    name: 'method2',
    param1: 1.0,
    param2: 2.0
};

const DefaultTFParams: TFParamsValues[] = [
    DefaultGaussianParams,
    DefaultMethod2Params
];

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
            methodsParams: DefaultTFParams
            // colored: this.props.colored
        };
        this.height = 400;
        this.width = 600;
        this.padding = 70;
        this.selectedPointId = undefined;
        this.ghostPoints = [];
        this.namespace = 'http://www.w3.org/2000/svg';
        this.descriptiveStatistics = this.getDescriptiveStatistics();
        // TODO: do not show last and first point (shadow) in the list of point buttons
        this.sortPoints(this.state.points);
        this.removePoint = this.removePoint.bind(this);
        this.removeRightmostPoint = this.removeRightmostPoint.bind(this);
        this.createPoint = this.createPoint.bind(this);
        this.handleDrag = this.handleDrag.bind(this);
        this.handleMultipleDrag = this.handleMultipleDrag.bind(this);
        this.handleDoubleClick = this.handleDoubleClick.bind(this);
        this.refCallBack = this.refCallBack.bind(this);
        this.handlePointUpdate = this.handlePointUpdate.bind(this);
        this.change = this.change.bind(this);
        this.handleKeyUp = this.handleKeyUp.bind(this);
        this.handleLeave = this.handleLeave.bind(this);
        this.handleEnter = this.handleEnter.bind(this);
        this.setTF = this.setTF.bind(this);
        this.changeXValue = this.changeXValue.bind(this);
        this.changeAlphaValue = this.changeAlphaValue.bind(this);
        this.deleteAllPoints = this.deleteAllPoints.bind(this);
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

    private deleteAllPoints() {
        const points = this.state.points;
        const ghostPoints = [];
        for (const point of points) {
            if (point.isGhost === true) { ghostPoints.push(point); }
        };
        const ghostPointsSorted = this.sortPoints(ghostPoints);
        this.setState({ points: ghostPointsSorted, clickedPointIds: undefined, showColorPicker: false });
        this.change([]);
    }

    private _setMethod2TF(params: Method2ParamsValues) {
        // TODO: implement
    }

    private _setGaussianTF(params: GaussianTFParamsValues) {
        const yOffset = this.padding / 2 / this.height;
        const { name, gaussianExtent, gaussianCenter, gaussianHeight } = params;
        const a = gaussianHeight;
        const min = this.descriptiveStatistics.min;
        const max = this.descriptiveStatistics.max;
        // const sigma = this.descriptiveStatistics.sigma;
        const TFextent = max - min;
        // make it spread the points equally
        // console.log(this.descriptiveStatistics);
        const b = gaussianCenter / TFextent;
        // const b = gaussianCenter * (mean + sigma) / TFextent;
        const c = gaussianExtent / TFextent;
        // const l = (this.width * 2 * c / extent);
        // console.log(a, b, c);
        const gaussianPoints: ControlPoint[] = generateGaussianControlPoints(a, b, c, TFextent, yOffset, this.props.volume);
        const currentPoints = this.state.points;
        gaussianPoints.push(currentPoints[currentPoints.length - 1]);
        gaussianPoints.unshift(currentPoints[0]);
        this.setState({
            points: gaussianPoints
        });
        this.change(gaussianPoints);
    }

    changeAlphaValue(pointId: UUID, alpha: number) {
        const points = this.state.points;
        const modifiedPoints = points.map(p => {
            if (p.id === pointId) {
                const modifiedData: ControlPointData = {
                    x: p.data.x,
                    alpha: alpha
                };
                const modifiedP: ControlPoint = {
                    id: p.id,
                    color: p.color,
                    index: p.index,
                    data: modifiedData
                };
                return modifiedP;
            } else {
                return p;
            }
        });
        console.log(modifiedPoints);
        this.handleChangePoints(modifiedPoints);
    }

    // TODO: refactor to remove duplication with above function
    changeXValue(pointId: UUID, x: number) {
        const points = this.state.points;
        const modifiedPoints = points.map(p => {
            if (p.id === pointId) {
                const modifiedData: ControlPointData = {
                    x: x,
                    alpha: p.data.alpha
                };
                const modifiedP: ControlPoint = {
                    id: p.id,
                    color: p.color,
                    index: p.index,
                    data: modifiedData
                };
                return modifiedP;
            } else {
                return p;
            }
        });
        console.log(modifiedPoints);
        this.handleChangePoints(modifiedPoints);
    }

    private setTF(params: TFParamsValues) {
        const name = params.name;
        switch (name) {
            case 'gaussian':
                this._setGaussianTF(params as GaussianTFParamsValues);
                break;
            case 'method2':
                this._setMethod2TF(params as Method2ParamsValues);
                break;
            default: throw Error(`Transfer function ${name} is not supported`);
        };
    }

    private getPoint(id: UUID) {
        const points = this.state.points;
        const selectedPoints = points.find(p => p.id === id);
        if (!selectedPoints) throw Error(`Points with ids ${id} do not exist`);
        return selectedPoints;
    }

    private getPoints(id: UUID[]) {
        const points = this.state.points;
        const selectedPoints = points.filter(p => id.includes(p.id));
        if (!selectedPoints) throw Error(`Point with id ${id} does not exist`);
        return selectedPoints;
    }

    public render() {
        // TODO: fix keys somewhere here
        const points = this.renderPoints();
        // const baseline = this.renderBaseline();
        const lines = this.renderLines();
        const histogram = this.renderHistogram();
        const axes = this.renderAxes();
        const gridLines = this.renderGridLines();
        // const pointsButtons = this.renderPointsButtons();
        const descriptiveStatisticsBars = this.renderDescriptiveStatisticsBars();
        const firstPoint = this.state.clickedPointIds ? this.getPoint(this.state.clickedPointIds[0]) : void 0;
        const color = firstPoint ? firstPoint.color : void 0;
        const defaultColor = ParamDefinition.Color(Color(0x121212));
        // bind this
        const _createPoint = this.createPoint;
        const _removePoint = this.removeRightmostPoint;
        const plusIconButtonStyle = {
            position: ('absolute' as any),
            top: '100px',
            right: '30px',
            // transform: 'translate(-50%, -50%)'
        };
        // TODO: fix style
        const minusIconButtonStyle = {
            position: ('absolute' as any),
            top: '100px',
            right: '5px',
            // transform: 'translate(-50%, -50%)'
        };

        return ([
            <div key="LineGraph">
                <IconButton style={plusIconButtonStyle} small svg={PlusBoxSvg} onClick={function (e: React.MouseEvent<HTMLButtonElement>): void {
                    // throw new Error('Function not implemented.');
                    _createPoint();
                } } ></IconButton>
                <IconButton style={minusIconButtonStyle} small svg={MinusBoxSvg} onClick={function (e: React.MouseEvent<HTMLButtonElement>): void {
                    // throw new Error('Function not implemented.');
                    _removePoint();
                } } ></IconButton>
                <svg
                    className="msp-canvas"
                    ref={this.refCallBack}
                    // TODO: can e.g. change viewbox
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
                        {/* {baseline} */}
                        {histogram}
                        {lines}
                        {points}
                        {descriptiveStatisticsBars}
                        {axes}
                        {gridLines}
                    </g>
                    <g className="ghost-points" stroke="black" fill="black">
                    </g>
                </svg>
                {/* TODO: div with the same height as color picker */}
                <>{this.state.showColorPicker ? <ColorPicker isActive={this.state.showColorPicker} defaultColor={defaultColor} color={color} updateColor={this.updateColor} toggleColorPicker={this.toggleColorPicker}/> : <div>Select point to change color</div>}

                </>
                <>
                    {/* <TFParamsWrapper onChange={this.setTF} descriptiveStatistics={this.descriptiveStatistics}></TFParamsWrapper> */}
                    {/* <Button onClick={this.deleteAllPoints}>Remove All Points</Button> */}
                    <PointsPanel points={this.state.points}
                        removeAllPoints={this.deleteAllPoints}
                        onExpandGroupOpen={this.props.onExpandGroupOpen}
                        changeXValue={this.changeXValue}
                        changeAlphaValue={this.changeAlphaValue}
                        onPointButtonClick={this.removePoint}></PointsPanel>
                    {/* Connect this to existing code */}
                    {/* should render actual state of all methods
                    // add this to state then, params of all methods */}
                    {/* TODO: debug */}
                    <HelpersPanel onChange={this.setTF} methods={this.state.methodsParams} descriptiveStatistics={this.descriptiveStatistics}></HelpersPanel>
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
        const clickedPointIds = this.state.clickedPointIds;
        const currentPoints = this.state.points;
        if (!clickedPointIds || clickedPointIds.length === 0) throw Error('No point is selected');
        const clickedPoints = this.getPoints(clickedPointIds);
        if (!clickedPoints) throw Error('Point should be selected');
        for (const point of clickedPoints) {
            point.color = value;
        }
        // need to change just points in ps
        const clickedPointsIds = clickedPoints.map(p => p.id);
        const notClickedPoints = currentPoints.filter(p => !clickedPointsIds.includes(p.id));
        const newPoints = clickedPoints.concat(notClickedPoints);
        console.log('Points before color change', currentPoints);
        console.log('Points with color changed', newPoints);
        this.handleChangePoints(newPoints);
        debugger;
    };

    private handleKeyDown = (event: any) => {
        // TODO: set canSelectMultiple = true
        if (event.key === 'Shift') {
            this.setState({ canSelectMultiple: true });
            console.log('Shift is pressed');
            // TODO: allow to select another point
        } else {
            // console.log(`${event} is pressed`);
        }
    };

    private handleKeyUp = (event: any) => {
        // TODO: SET canSelectMultiple = fasle
        this.setState({ canSelectMultiple: false });
        if (event.shiftKey) {
            // do something
            console.log('Shift is released');
            // TODO: allow to select another point
        } else {
            // console.log(`${event} is released`);
        }
    };

    private handleClick = (point: ControlPoint) => (event: any) => {
        this.setState({ clickedPointIds: [point.id] });

        if (this.state.canSelectMultiple) {
            if (event.shiftKey) return;
            // TODO: function to execute on this
        }

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

        if ((svgP.x < (padding) || svgP.x > (this.width + (padding))) &&
        (svgP.y > (this.height + (0)) || svgP.y < (0))) {
            updatedCopyPoint = Vec2.create(this.updatedX, this.updatedY);
        } else if (svgP.x < padding) {
            // TODO: fix lines to start from true 0
            updatedCopyPoint = Vec2.create(padding, svgP.y);
        } else if (svgP.x > (this.width + (padding))) {
            updatedCopyPoint = Vec2.create(this.width + padding, svgP.y);
        } else if (svgP.y > (this.height + (0))) {
            // does not allow into such area
            // fix viewbox or?
            updatedCopyPoint = Vec2.create(svgP.x, this.height + 0);
        } else if (svgP.y < (0)) {
            updatedCopyPoint = Vec2.create(svgP.x, 0);
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
        // do not need to raise points, simply generate their original position
        // as x + offset
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

    private removePoint(id: UUID) {
        const point = this.getPoint(id);
        if (point.index === 0 || point.index === this.state.points.length - 1) { return; }
        let points = this.state.points.filter(p => p.id !== point.id);
        points = this.sortPoints(points);
        this.setState({ points: points, clickedPointIds: undefined, showColorPicker: false });
        this.change(points);
    }

    // TODO: fix movement of plus and minus sign buttons upon creation of points
    private removeRightmostPoint() {
        // remove based on data
        const points = this.state.points;
        const sortedPs = this.sortPoints(points);
        debugger;
        // rightmost is undefined
        // last is ghost
        // const rightmostP = sortedPs.slice(-1)[0];
        const rightmostP = sortedPs[sortedPs.length - 2];
        debugger;
        this.removePoint(rightmostP.id);
        debugger;
    };

    private createPoint() {
        const svgP = this.myRef.createSVGPoint();
        // TODO: resolve location, do something like this.width, height...
        svgP.x = this.width - 100;
        svgP.y = this.height - 100;

        // const svgP = pt.matrixTransform(this.myRef.getScreenCTM().inverse());
        const points = this.state.points;
        // const padding = this.padding / 2;
        // debugger;
        // // resolve this
        // const svtP = 0;
        // if (svgP.x < (padding) ||
        //     svgP.x > (this.width + (padding)) ||
        //     svgP.y > (this.height + (padding)) ||
        //     svgP.y < (this.padding / 2)) {
        //     debugger;
        //     return;
        // }
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

    private handleDoubleClick(event: any) {
        return;
        // const pt = this.myRef.createSVGPoint();
        // pt.x = event.clientX;
        // pt.y = event.clientY;
        // const svgP = pt.matrixTransform(this.myRef.getScreenCTM().inverse());
        // const points = this.state.points;
        // const padding = this.padding / 2;

        // if (svgP.x < (padding) ||
        //     svgP.x > (this.width + (padding)) ||
        //     svgP.y > (this.height + (padding)) ||
        //     svgP.y < (this.padding / 2)) {
        //     return;
        // }
        // const newPointData = this.unNormalizePoint(Vec2.create(svgP.x, svgP.y));
        // // problem with index
        // // should be the last one
        // // find the last index and do + 1?
        // const newPoint: ControlPoint = {
        //     data: newPointData,
        //     id: UUID.create22(),
        //     color: getRandomColor(),
        //     index: points.slice(-1)[0].index + 1
        // };
        // this.addPoint(newPoint);
    }

    private deletePoint = (id: UUID) => (event: any) => {
        this.removePoint(id);
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

    private renderGridLines() {
        const count = 5;
        const bars: any = [];
        const offset = this.padding / 2;
        const x1 = offset;
        const w = 0.5;
        // TODO: allow editing text field, may need some function similar to what is with RGB thing
        //
        // TODO: consider adding baseline height as attribute of LineGraphComponent
        const x2 = this.width + offset;
        // fix that
        const overallHeight = this.height;
        // TODO: limit the height of histogram bars
        // 1, 2, 3, should be 4
        // get histogram bars back
        for (let i = 0; i < count; ++i) {
            // adjust + - 1 etc.
            // since we have risen all elements, need to probably increase height or
            // decrease it by - offset everywhere where it is used
            const y = overallHeight * i / count;
            // const y = this.height + offset * 2 - i * this.height / count;
            bars.push(
                <line key={`${i / count}gridline`} x1={x1} x2={x2} y1={y} y2={y}
                    stroke="#A9A9A9" strokeWidth={w} strokeDasharray="15, 15">
                    <title>{(i / count).toFixed(1)}</title>
                </line>
            );
        }
        return bars;
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
            const fromValue = histogram.min + histogram.binWidth * i;
            const toValue = histogram.min + histogram.binWidth * (i + 1);
            const x = this.width * i / (N - 1) + offset;
            // may need to do +-offset on y1 or y2
            const y1 = this.height;// + 2 * offset;
            const y2 = this.height * (1 - histogram.counts[i] / max);// + 2 * offset;
            // console.log(y1, y2);
            bars.push(<line key={`histogram${i}`} x1={x} x2={x} y1={y1} y2={y2} stroke="#A9A9A9" strokeWidth={w}>
                <title>[{fromValue}; {toValue}]</title>
            </line>);
        }
        return bars;
    }

    private renderAxes() {
        if (!this.props.volume) return null;
        const offset = this.padding / 2;
        const mean = this.descriptiveStatistics.mean;
        const min = this.descriptiveStatistics.min;
        const max = this.descriptiveStatistics.max;
        // X: horizontal bar with arrow
        // need min max
        const x1HorizontalBar = 0;
        const x2HorizontalBar = this.width + offset;
        const y1HorizontalBar = this.height + offset;
        const y2HorizontalBar = this.height + offset;
        // vertical bar with arrow
        // x and y letters



        const x2VerticalBar = offset;
        const x1VerticalBar = offset;
        const y2VerticalBar = 25;
        const y1VerticalBar = this.height + offset;
        const w = offset / 10;
        const bars = [];
        // const hd = `M ${xh1} ${yh1} L ${xh2} ${yh2} L ${xh3} ${yh3} L ${xh4} ${yh4} L ${xh5} ${yh5} Z`;
        // const arrowheadd = `M ${xah1} ${yah1} L ${xah2} ${yah2} L ${xah3} ${yah3} L ${xah4} ${yah4} L ${xah5} ${yah5} Z`;
        bars.push(
            // TODO: color
            // TODO: may not work with <>/</>
            <>
                <defs>
                    <marker
                        id="head-horizontal"
                        viewBox="0 0 10 10"
                        refX="5"
                        refY="5"
                        markerWidth="6"
                        markerHeight="6"
                        orient="auto-start-reverse">
                        <path d="M 0 2 L 10 5 L 0 8 z" />
                    </marker>
                </defs>

                <line key={'horizontalAxis'} x1={x1HorizontalBar} x2={x2HorizontalBar}
                    y1={y1HorizontalBar} y2={y2HorizontalBar} stroke="#7d7f7c" strokeWidth={w} markerEnd="url(#head-horizontal)">
                </line>
            </>
        );
        // raise histogram above the base line somehow
        bars.push(
            <>
                <defs>
                    <marker
                        id="head-vertical"
                        viewBox="0 0 10 10"
                        refX="5"
                        refY="5"
                        markerWidth="6"
                        markerHeight="6"
                        orient="auto-start-reverse">
                        <path d="M 0 2 L 10 5 L 0 8 z" />
                        {/* <path d="M 2 10 L 5 0 L 8 10 z" /> */}
                    </marker>
                </defs>
                <line key={'verticalAxis'} x1={x1VerticalBar} x2={x2VerticalBar}
                    y1={y1VerticalBar} y2={y2VerticalBar} stroke="#7d7f7c" strokeWidth={w} markerEnd="url(#head-vertical)">
                    {/* <title>+ Sigma: {sigma}</title> */}
                </line>
            </>
        );
        return bars;
    }

    private makeXAxisLabel(value: number, x: number, y: number) {
        // TODO: make constants be dependant on precision and font size
        // TODO: font size
        // TODO: remove y dep
        return <text x={x - 20} y={y + 25}>{parseFloat(value.toFixed(5))} </text>;
    }

    private makeYAxisLabel(value: number, x: number, y: number) {
        // TODO: make constants be dependant on precision and font size
        // TODO: font size
        // TODO: remove x dep
        return <text x={x} y={y + 25} fontSize={25}>{parseFloat(value.toFixed(1))} </text>;
    }

    private renderDescriptiveStatisticsBars() {
        if (!this.props.volume) return null;
        const offset = this.padding / 2;
        const mean = this.descriptiveStatistics.mean;
        const min = this.descriptiveStatistics.min;
        const max = this.descriptiveStatistics.max;
        const sigma = this.descriptiveStatistics.sigma;
        const extent = max - min;
        const x = this.width * ((mean + Math.abs(min)) / extent) + offset;
        const w = offset / 10;
        const bars = [];
        const y1 = this.height + offset;// + offset * 2;
        const y2 = 0;// offset;
        const xPositive = this.width * ((mean + sigma + Math.abs(min)) / extent) + offset;
        const xNegative = this.width * ((mean - sigma + Math.abs(min)) / extent) + offset;
        bars.push(
            // raise points and lines
            //
            <>
                <line key={'meanBar'} x1={x} x2={x} y1={y1} y2={y2} stroke="#D3D3D3" strokeDasharray="5, 5" strokeWidth={w}>
                    <title>Mean: {mean}</title>
                </line>
                {/* TODO: center */}
                {this.makeXAxisLabel(mean, x, y1)}
                {/* <text x={x - 40} y={y1 + 25}>{parseFloat(mean.toFixed(7))} </text> */}
            </>);
        bars.push(<><line key={'positiveSigmaBar'} x1={xPositive} x2={xPositive} y1={y1} y2={y2} stroke="#D3D3D3" strokeWidth={w}>
            <title>+Sigma: {mean + sigma}</title>
        </line>
        {this.makeXAxisLabel(mean + sigma, xPositive, y1)}
        </>);
        bars.push(
            <><line key={'negativeSigmaBar'} x1={xNegative} x2={xNegative} y1={y1} y2={y2} stroke="#D3D3D3" strokeWidth={w}>
                <title>-Sigma: {-sigma}</title>
            </line>
            {this.makeXAxisLabel(mean - sigma, xNegative, y1)}
            </>);
        return bars;
    }

    // fix gridlines and labels
    // first fix points, allow them in area above and recalc the message or something
    // remove ghost points from the list
    private renderBaseline() {
        const offset = this.padding / 2;
        return <>
            <BaseLine height={this.height} offset={this.padding / 2} width={this.width} />
            {/* TODO: tune height need to be higher */}
            {this.makeYAxisLabel(0, offset / 2, this.height - offset)}
            {/*  */}
            {this.makeYAxisLabel(1, offset / 2, offset)}
            {/* TODO: 0.2 0.4 etc. */}
        </>;
    }

    // turn off handleclick too
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
                    // same for lines
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