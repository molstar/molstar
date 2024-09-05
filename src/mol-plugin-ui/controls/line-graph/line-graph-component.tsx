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
import { Grid, Volume } from '../../../mol-model/volume';
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
import { generateControlPoints, generateGaussianControlPoints } from '../../../mol-geo/geometry/direct-volume/direct-volume';
import { useCallback, useEffect, useRef, useState } from 'react';
import { LineGraphParams, startEndPoints } from './line-graph-params';
import { throwError } from 'rxjs';

type ComponentParams<T extends React.Component<any, any, any> | ((props: any) => JSX.Element)> =
    T extends React.Component<infer P, any, any> ? P : T extends (props: infer P) => JSX.Element ? P : never;

class PointButton extends React.Component<any> {
    onAbs = (v: number) => {
        this.props.changeXValue(this.props.point.id, v);
        // this.props.value = v;
    };

    onAlpha = (v: number) => {
        // set here based on value, that is fine, just truncate it then in visual
        this.props.changeAlphaValue(this.props.point.id, v);
    };

    handleClick = () => {
        this.props.onClick(this.props.point.id);
    };

    onPointIndexClick = (pointId: UUID) => {
        this.props.onPointIndexClick(pointId);
    };

    render() {
        const [absValue, relativeValue] = this.props.onExpandGroupOpen(this.props.point.data);
        const truncatedAbsValue = parseFloat((absValue as number).toFixed(3));
        // console.log(absValue, relativeValue);
        const alpha = (this.props.point.data.alpha as number).toFixed(3);
        const color = this.props.point.color as Color;
        return (
            <div style={{ display: 'flex', marginBottom: 1 }} key={this.props.point.id}>
                {/* TODO: fix indices */}
                <Button style={{ textAlign: 'start', textIndent: '20px' }}>{this.props.point.index}</Button>
                {/* TODO: pass to onColorSquareClick prop a function to invoke colorpicker */}
                {/* toggleColorPicker */}
                {<Button style={{ backgroundColor: Color.toStyle(color), minWidth: 32, width: 32 }}
                    onClick={() => { this.props.onColorSquareClick(this.props.point.id); } } />}
                {/* fix here number of points */}
                {/* so decouple from the values */}
                {/* just change it */}
                {/* got it need at list parsefloat */}
                <TextInput numeric
                    style={{ minWidth: 0 }} className='msp-form-control' onEnter={this.props.onEnter} blurOnEnter={true} blurOnEscape={true}
                    value={truncatedAbsValue} placeholder={'Some text'}
                    // fix that
                    // TODO: string to number?
                    isDisabled={false} onChange={(value) => { this.onAbs(value); }} />
                <TextInput numeric
                // ok set value to string repr parseFloat and toFixed, but store in state the true value instead
                // and set the value based on the state inside onAlpha onAbs
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
        case 'defaults':
            return DefaultParams;
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
        const realPoints = points.filter(p => p.isTerminal !== true);
        debugger;
        // const realPoints = points.filter();
        // TODO: add ghost prop to points
        // const removeAllPointsButton = <Button style={{ position: 'absolute', top: 0, right: 0 }} onClick={this.props.removeAllPoints}>Remove All Points</Button>;
        const iconStyle = { position: 'absolute', top: 0, right: 0, lineHeight: '24px', height: '24px', textAlign: 'right', width: '32px', paddingRight: '6px', background: 'none' };
        const removeAllPointsButton = <IconButton title='Remove All Points' style={iconStyle as any} svg={DeleteSvg} small onClick={() => {
            this.props.removeAllPoints();
        } }></IconButton>;
        const controlPointsButtons = realPoints.map(p => {
            return <PointButton key={p.id} point={p}
                onPointIndexClick={this.props.onPointIndexClick}
                onColorSquareClick={this.props.onColorSquareClick}
                onClick={this.props.onPointButtonClick}
                onExpandGroupOpen={this.props.onExpandGroupOpen}
                changeXValue={this.props.changeXValue}
                changeAlphaValue={this.props.changeAlphaValue}
            ></PointButton>;
        });
        return (
            // TODO: may need this className='msp-representation-entry'
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

export type TFName = 'gaussian' | 'method2' | 'defaults';

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
            <Button style={{ marginTop: 1 }} onClick={this.handleClick}>{`Apply ${this.props.params.name} Transfer Function`}</Button>
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
    // adjustedAlpha?: number
    absValue?: number
    relativeValue?: number
}

export interface ControlPoint {
    id: UUID,
    color: Color,
    data: ControlPointData
    index: number,
    isTerminal?: boolean
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
    methodsParams: TFParamsValues[],
}

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
    // TODO: PD from object or something, need list of values
    name: PD.Text('gaussian', { isHidden: true })
};

export const DefaultParams = {
    name: PD.Text('defaults', { isHidden: true })
};

export const Method2TFParams = {
    param1: PD.Numeric(0.2, { min: 0, max: 1, step: 0.01 }),
    param2: PD.Numeric(0.2, { min: 0, max: 1, step: 0.01 }),
    // TODO: PD from object or something
    name: PD.Text('method2', { isHidden: true })
    // gaussianHeight: PD.Numeric(0.2, { min: 0, max: 1, step: 0.01 })
};

export type GaussianTFParamsValues = PD.Values<typeof GaussianTFParams>;

export type TFParamsValues = GaussianTFParamsValues | Method2ParamsValues | DefaultParamsValues;
export type DefaultParamsValues = PD.Values<typeof DefaultParams>;
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

const DefaultDefaultParams: DefaultParamsValues = {
    name: 'defaults'
};

const DefaultTFParams: TFParamsValues[] = [
    DefaultGaussianParams,
    DefaultMethod2Params,
    DefaultDefaultParams
];

export interface LineGraphComponentProps {
    // TODO: better name
    data: ControlPoint[]
    // TODO: may need to have it as any
    volume: Volume | undefined
    // TODO: function types
    onChange: any
    onHover: any
    onDrag: any
    // TODO: why undefined?
    colored: boolean | undefined
    // TODO: better name
    onExpandGroupOpen: any
    onAbsValueToPointValue: any
};

// export interface LineGraphAttrs {
//     // TODO: correct types
//     myRef = React.createRef();
//     this.state = {
//         points: startEndPoints.concat(this.props.data),
//         copyPoint: undefined,
//         canSelectMultiple: false,
//         showColorPicker: false,
//         methodsParams: DefaultTFParams
//         // colored: this.props.colored
//     };
//     this.roof = LineGraphParams.roof;
//     this.roofUnnormalized = LineGraphParams.roofUnnormalized;
//     this.baseline = LineGraphParams.baseline;
//     this.baselineUnnormalized = LineGraphParams.baselineUnnormalized;
//     this.height = LineGraphParams.height;
//     this.width = LineGraphParams.width;
//     this.padding = LineGraphParams.padding;
//     // this.offsetY = LineGraphParams.offsetY;
//     // this.offsetY = 0;
//     this.selectedPointId = undefined;
//     this.ghostPoints = [];
//     this.namespace = 'http://www.w3.org/2000/svg';
//     this.descriptiveStatistics = this.getDescriptiveStatistics();

// }

export function LineGraphComponent(props: LineGraphComponentProps) {
    const [controlPoints, setControlPoints] = useState(startEndPoints.concat(props.data));
    const [copyPoint, setCopyPoint] = useState<Vec2 | undefined >(undefined);
    const [canSelectMultiple, setCanSelectMultiple] = useState(false);
    const [showColorPicker, setShowColorPicker] = useState(false);
    // TODO: correct type
    const [clickedPointIds, setClickedPointIds] = useState<UUID[]>([]);
    const [methodsParams, setMethodsParams] = useState(DefaultTFParams);

    // const attrs = useRef<LineGraphAttrs | undefined>();

    useEffect(() => {
        gElement.current = document.getElementsByClassName('ghost-points')[0] as SVGElement;
    });

    const myRef = useRef<React.RefObject<any> | undefined>(undefined);
    const height = useRef(LineGraphParams.height);
    const width = useRef(LineGraphParams.width);
    // TODO: revisit padding
    const padding = useRef(LineGraphParams.padding);
    // TODO: may have issues if should be undefined instead
    const updatedX = useRef(0);
    const updatedY = useRef(0);
    // TODO: compare with clickedPointIds
    const selectedPointId = useRef<UUID | undefined>(undefined);
    const ghostPoints = useRef<SVGElement[]>([]);
    const gElement = useRef<SVGElement | undefined>(undefined);
    const namespace = useRef('http://www.w3.org/2000/svg');
    const descriptiveStatistics = useRef(getDescriptiveStatistics());
    const baseline = useRef(LineGraphParams.baseline);
    const roof = useRef(LineGraphParams.roof);

    function getDescriptiveStatistics() {
        const v = props.volume;
        if (!v) throw Error('No volume is provided');
        const s = v.grid.stats;
        const d: VolumeDescriptiveStatistics = {
            min: s.min,
            max: s.max,
            mean: s.mean,
            sigma: s.sigma
        };
        return d;
    }

    function deleteAllPoints() {
        const terminalPoints = [];
        for (const point of controlPoints) {
            if (point.isTerminal === true) { terminalPoints.push(point); }
        };
        const terminalPointsSorted = sortPointsByXValues(terminalPoints);
        setControlPoints(terminalPointsSorted);
        setClickedPointIds([]);
        setShowColorPicker(false);
        // this.setState({ points: terminalPointsSorted, clickedPointIds: undefined, showColorPicker: false });
        change([]);
    }

    function highlightPoint(pointId: UUID) {
        const targetPoint = controlPoints.find(p => p.id === pointId);
        throw Error('Not implemented');
        // if (!targetPoint) throw Error('Cannot highlight inexisting point exist');
    }

    function _setMethod2TF(params: Method2ParamsValues) {
        // TODO: implement
    }

    function _setDefaultsTF(params: DefaultParamsValues) {
        const paddingUnnormalized = padding.current / height.current;
        // there two are the same values, reduce to a single one, e.g. could be be just padding, why not
        // should be e.g.
        const generatedPoints = generateControlPoints(paddingUnnormalized, ColorNames.black, undefined, props.volume);
        const currentPoints = controlPoints;
        generatedPoints.push(currentPoints[currentPoints.length - 1]);
        generatedPoints.unshift(currentPoints[0]);
        setControlPoints(generatedPoints);
        change(generatedPoints);
    }

    function change(points: ControlPoint[]) {
        const copyPoints = points.slice();
        copyPoints.shift();
        copyPoints.pop();
        props.onChange(copyPoints);
    }

    function _setGaussianTF(params: GaussianTFParamsValues) {
        const p = padding.current;
        const h = height.current;
        const paddingUnnormalized = (p / 2) / h;
        const { gaussianExtent, gaussianCenter, gaussianHeight } = params;
        const a = gaussianHeight;
        const ds = descriptiveStatistics.current;
        const min = ds.min;
        const max = ds.max;
        const TFextent = max - min;
        const b = gaussianCenter / TFextent;
        const c = gaussianExtent / TFextent;
        if (!props.volume) throw Error('No volume');
        const gaussianPoints: ControlPoint[] = generateGaussianControlPoints(a, b, c, TFextent, 0, props.volume, paddingUnnormalized);
        const currentPoints = controlPoints;
        gaussianPoints.push(currentPoints[currentPoints.length - 1]);
        gaussianPoints.unshift(currentPoints[0]);
        setControlPoints(gaussianPoints);
        change(gaussianPoints);
    }

    function changeAlphaValue(pointId: UUID, alpha: number) {
        const modifiedPoints = controlPoints.map(p => {
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
        handleChangePoints(modifiedPoints);
    }

    function handleChangePoints(points: ControlPoint[]) {
        const pointsSorted = sortPointsByXValues(points);
        setControlPoints(pointsSorted);
        change(pointsSorted);
    }

    function sortPointsByXValues(points: ControlPoint[]) {
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

    // TODO: refactor to remove duplication with changeAlphaValue function
    function changeXValue(pointId: UUID, absValue: number) {
        const x = props.onAbsValueToPointValue(absValue);
        const points = controlPoints;
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
        handleChangePoints(modifiedPoints);
    }

    function setTF(params: TFParamsValues) {
        const name = params.name as TFName;
        switch (name) {
            // TODO: separate into type
            case 'gaussian':
                _setGaussianTF(params as GaussianTFParamsValues);
                break;
            case 'method2':
                _setMethod2TF(params as Method2ParamsValues);
                break;
            case 'defaults':
                _setDefaultsTF(params as DefaultParamsValues);
                break;
            default: throw Error(`Transfer function ${name} is not supported`);
        };
    }

    function getPoint(id: UUID) {
        const selectedPoint = controlPoints.find(p => p.id === id);
        if (!selectedPoint) throw Error(`Point with id ${id} does not exist`);
        return selectedPoint;
    }

    function getPoints(ids: UUID[]) {
        const selectedPoints = controlPoints.filter(p => ids.includes(p.id));
        if (!selectedPoints) throw Error(`Points with ids ${ids} do not exist`);
        return selectedPoints;
    }

    function controlPointDataToSVGCoords(controlPointData: ControlPointData): Vec2 {
        const p = padding.current;
        const w = width.current;
        const h = height.current;
        const b = baseline.current;
        const r = roof.current;
        const offset = p / 2;
        const maxX = w + offset;
        const normalizedX = (controlPointData.x * (maxX - offset)) + offset;
        const space = _getYSpace();
        const yFromBottom = controlPointData.alpha * space;
        const fromTopToBaseline = h - b;
        const yFromTop = fromTopToBaseline - yFromBottom;
        const newPoint = Vec2.create(normalizedX, yFromTop);
        return newPoint;
    }

    const deletePoint = (id: UUID) => (event: any) => {
        removePoint(id);
        event.stopPropagation();
    };

    function removePoint(id: UUID) {
        const point = getPoint(id);
        if (point.isTerminal === true) { return; }
        let points = controlPoints.filter(p => p.id !== point.id);
        points = sortPointsByXValues(points);
        setControlPoints(points);
        // TODO: may not work, if it is the case - add undefined option to setState init
        setClickedPointIds([]);
        setShowColorPicker(false);
        change(points);
    }

    function handleLeave() {
        if (selectedPointId === undefined) {
            return;
        }

        document.addEventListener('mousemove', handleDrag, true);
        document.addEventListener('mouseup', handlePointUpdate, true);
    }
    function svgCoordsToPointData(vec2: Vec2): ControlPointData {
        const h = height.current;
        const b = baseline.current;
        const r = roof.current;
        const w = width.current;
        const p = padding.current;
        const space = _getYSpace();
        const min = p / 2;
        const maxX = w + min;
        const unNormalizedX = (vec2[0] - min) / (maxX - min);
        const cartesianDistanceBetweenPointAndBaseline = (h - b - vec2[1]);
        const dividedBySpace = cartesianDistanceBetweenPointAndBaseline / space;
        const unNormalizedY = dividedBySpace;
        return {
            x: unNormalizedX,
            alpha: unNormalizedY
        };
    }

    function handleDrag(event: any) {
        if (selectedPointId.current === undefined || !myRef.current) {
            return;
        }
        const w = width.current;
        const h = height.current;
        const b = baseline.current;
        const r = roof.current;
        const mr = (myRef.current as any);
        const pt = mr.createSVGPoint();
        let updatedCopyPoint;
        const offset = padding.current / 2;
        pt.x = event.clientX;
        pt.y = event.clientY;
        const svgP = pt.matrixTransform(mr.getScreenCTM().inverse());
        updatedCopyPoint = Vec2.create(svgP.x, svgP.y);
        if ((svgP.x < (offset) || svgP.x > (w + (offset))) &&
        (svgP.y > (h - b) || svgP.y < (0))) {
            updatedCopyPoint = Vec2.create(updatedX.current, updatedY.current);
        } else if (svgP.x < offset) {
            updatedCopyPoint = Vec2.create(offset, svgP.y);
        } else if (svgP.x > (w + (offset))) {
            updatedCopyPoint = Vec2.create(w + offset, svgP.y);
        } else if (svgP.y > (h - b)) {
            updatedCopyPoint = Vec2.create(svgP.x, h - b);
        } else if (svgP.y < (r)) {
            updatedCopyPoint = Vec2.create(svgP.x, r);
        } else {
            updatedCopyPoint = Vec2.create(svgP.x, svgP.y);
        }

        updatedX.current = updatedCopyPoint[0];
        updatedY.current = updatedCopyPoint[1];
        const unNormalizePoint = svgCoordsToPointData(updatedCopyPoint);
        // TODO: this.ghostPoints[0] is undefined when dragging a point towards right
        // border, why?
        ghostPoints.current[0].setAttribute('style', 'display: visible');
        ghostPoints.current[0].setAttribute('cx', `${updatedCopyPoint[0]}`);
        ghostPoints.current[0].setAttribute('cy', `${updatedCopyPoint[1]}`);

        props.onDrag(unNormalizePoint);
    }

    function replacePoint(point: ControlPoint) {
        let points = controlPoints.filter(p => p.id !== point.id);
        points.push(point);
        points = sortPointsByXValues(points);
        setControlPoints(points);
        change(points);
    }

    function handlePointUpdate(event: any) {
        const selected = selectedPointId.current;
        if (canSelectMultiple) {
            return;
        }
        if (selected === undefined || getPoint(selected).isTerminal === true) {
            setCopyPoint(undefined);
            return;
        }
        selectedPointId.current = undefined;
        const pointData = getPoint(selected);
        const updatedPointData = svgCoordsToPointData(Vec2.create(updatedX.current, updatedY.current));
        const updatedPoint: ControlPoint = {
            data: updatedPointData,
            id: pointData.id,
            color: pointData.color,
            index: pointData.index
        };

        replacePoint(updatedPoint);
        if (!gElement.current) throw Error('No gElement (SVGElement)');
        gElement.current.innerHTML = '';
        ghostPoints.current = [];
        document.removeEventListener('mousemove', handleDrag, true);
        document.removeEventListener('mouseup', handlePointUpdate, true);
    }

    const handleMouseDown = (point: ControlPoint) => (event: any) => {
        const { id, isTerminal } = point;
        if (isTerminal === true) {
            return;
        }
        if (canSelectMultiple) {
            return;
        }
        const gps = ghostPoints.current;
        const copyPoint: Vec2 = controlPointDataToSVGCoords(getPoint(id).data);
        gps.push(document.createElementNS(namespace.current, 'circle') as SVGElement);
        gps[0].setAttribute('r', '10');
        gps[0].setAttribute('fill', 'orange');
        gps[0].setAttribute('cx', `${copyPoint[0]}`);
        gps[0].setAttribute('cy', `${copyPoint[1]}`);
        gps[0].setAttribute('style', 'display: none');
        if (!gElement.current) throw Error('No gElement (SVGElement)');
        gElement.current.appendChild(gps[0]);
        updatedX.current = copyPoint[0];
        updatedY.current = copyPoint[1];
        selectedPointId.current = point.id;
    };

    const toggleColorPicker = () => {
        setShowColorPicker(showColorPicker === true ? false : true);
    };

    const handleClick = (point: ControlPoint) => (event: any) => {
        setClickedPointIds([point.id]);

        if (canSelectMultiple) {
            if (event.shiftKey) return;
            // TODO: function to execute on this
        }

        if (showColorPicker) {
            // TODO: what should happen there?
        } else {
            toggleColorPicker();
        }
    };

    function renderPoints() {
        const points: any[] = [];
        let point: Vec2;
        console.log(controlPoints);
        for (let i = 0; i < controlPoints.length; i++) {
            if (controlPoints[i].isTerminal !== true) {
            // if (i !== 0 && i !== controlPoints.length - 1) {
                const { data, color, id, index } = controlPoints[i];
                const finalColor = color;
                point = controlPointDataToSVGCoords(data);
                points.push(<PointComponent
                    index={index}
                    key={id}
                    id={id}
                    x={point[0]}
                    y={point[1]}
                    nX={data.x}
                    nY={data.alpha}
                    selected={false}
                    delete={deletePoint}
                    // onmouseover we provide props.onHover of line graph component
                    onmouseover={props.onHover}
                    onmousedown={handleMouseDown(controlPoints[i])}
                    onclick={handleClick(controlPoints[i])}
                    color={finalColor}
                />);
            }
        }
        return points;
    }

    function makeYAxisLabel(value: number, x: number, y: number) {
        // TODO: make constants be dependant on precision and font size
        // TODO: font size
        // TODO: remove x dep
        return <text x={x} y={y + 25} fontSize={25}>{parseFloat(value.toFixed(1))} </text>;
    }

    function renderBaseline() {
        const { h, b, r, w, p } = _getLineGraphAttributes();
        return <>
            <BaseLine height={h - b} offset={p / 2} width={w} />
            {makeYAxisLabel(0, b / 2, h)}
            {makeYAxisLabel(1, b / 2, b)}
        </>;
    }

    function _getYSpace() {
        const { h, b, r, w, p } = _getLineGraphAttributes();
        const space = h - b - r;
        return space;
    }

    function _getLineGraphAttributes() {
        const h = height.current;
        const b = baseline.current;
        const r = roof.current;
        const w = width.current;
        const p = padding.current;
        return { h, b, r, w, p };
    }

    function renderLines() {
        const points: Vec2[] = [];
        const lines = [];
        let maxX: number;
        let maxY: number;
        let normalizedX: number;
        let normalizedY: number;
        let reverseY: number;

        const { h, b, r, w, p } = _getLineGraphAttributes();

        const o = p / 2;
        for (const point of controlPoints) {
            maxX = w + o;
            maxY = h + p;
            normalizedX = (point.data.x * (maxX - o)) + o;
            const alpha = point.data.alpha;
            const space = _getYSpace();
            const pointHeightRealAboveBaseline = alpha * space;
            const pointHeightRealBelowRoof = space - pointHeightRealAboveBaseline;
            const pointHeightBelowSVGOrigin = pointHeightRealBelowRoof + r;
            points.push(Vec2.create(normalizedX, pointHeightBelowSVGOrigin));
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

    function renderHistogram() {
        if (!props.volume) return null;
        const { h, b, r, w, p } = _getLineGraphAttributes();
        const histogram = Grid.getHistogram(props.volume.grid, 40);
        const bars = [];
        const N = histogram.counts.length;
        const histogramW = w / N;
        const offset = p / 2;
        const max = arrayMax(histogram.counts) || 1;
        for (let i = 0; i < N; i++) {
            const fromValue = histogram.min + histogram.binWidth * i;
            const toValue = histogram.min + histogram.binWidth * (i + 1);
            const x = w * i / (N - 1) + offset;
            const y1 = h - b;
            const space = _getYSpace();
            const fraction = space * (1 - histogram.counts[i] / max);
            const y2 = fraction + r;
            bars.push(<line key={`histogram${i}`} x1={x} x2={x} y1={y1} y2={y2} stroke="#A9A9A9" strokeWidth={histogramW}>
                <title>[{fromValue}; {toValue}]</title>
            </line>);
        }
        return bars;
    }

    function renderAxes() {
        if (!props.volume) return null;
        const { h, b, r, w, p } = _getLineGraphAttributes();
        const offset = p / 2;
        const x1HorizontalBar = offset;
        const x2HorizontalBar = w + offset;
        const y1HorizontalBar = h + offset;
        const y2HorizontalBar = h + offset;
        const x2VerticalBar = offset;
        const x1VerticalBar = offset;
        const y2VerticalBar = 25;
        const y1VerticalBar = h + offset;
        const barW = offset / 10;
        const bars = [];
        bars.push(
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
                    y1={y1HorizontalBar} y2={y2HorizontalBar} stroke="#7d7f7c" strokeWidth={barW} markerEnd="url(#head-horizontal)">
                </line>
            </>
        );
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
                    y1={y1VerticalBar} y2={y2VerticalBar} stroke="#7d7f7c" strokeWidth={barW} markerEnd="url(#head-vertical)">
                    {/* <title>+ Sigma: {sigma}</title> */}
                </line>
            </>
        );
        return bars;
    }

    function makeXAxisLabel(value: number, x: number, y: number) {
        // TODO: make constants be dependant on precision and font size
        // TODO: font size
        // TODO: remove y dep
        return <text x={x - 20} y={y + 25}>{parseFloat(value.toFixed(5))} </text>;
    }

    function renderGridLines() {
        const { h, b, r, w, p } = _getLineGraphAttributes();
        const count = 5;
        const bars: any = [];
        const offset = p / 2;
        const x1 = offset;
        const lineW = 0.5;
        const x2 = w + offset;
        for (let i = 0; i < count; ++i) {
            const y = h * i / count;
            bars.push(
                <line key={`${1 - i / count}gridline`} x1={x1} x2={x2} y1={y} y2={y}
                    stroke="#A9A9A9" strokeWidth={lineW} strokeDasharray="15, 15">
                    <title>{(1 - i / count).toFixed(1)}</title>
                </line>
            );
        }
        return bars;
    }

    function renderDescriptiveStatisticsBars() {
        if (!props.volume) return null;
        const { h, b, r, w, p } = _getLineGraphAttributes();
        const offset = p / 2;
        const { mean, min, max, sigma } = descriptiveStatistics.current;
        const extent = max - min;
        const x = w * ((mean + Math.abs(min)) / extent) + offset;
        const barW = offset / 10;
        const bars = [];
        const y1 = h + offset;// + offset * 2;
        const y2 = 0;// offset;
        const xPositive = w * ((mean + sigma + Math.abs(min)) / extent) + offset;
        const xNegative = w * ((mean - sigma + Math.abs(min)) / extent) + offset;
        bars.push(
            <>
                <line key={'meanBar'} x1={x} x2={x} y1={y1} y2={y2} stroke="#D3D3D3" strokeDasharray="5, 5" strokeWidth={barW}>
                    <title>Mean: {mean}</title>
                </line>
                {makeXAxisLabel(mean, x, y1)}
            </>);
        bars.push(<>
            <line key={'positiveSigmaBar'} x1={xPositive} x2={xPositive} y1={y1} y2={y2} stroke="#D3D3D3" strokeWidth={barW}>
                <title>+Sigma: {mean + sigma}</title>
            </line>
            {makeXAxisLabel(mean + sigma, xPositive, y1)}
        </>);
        bars.push(
            <><line key={'negativeSigmaBar'} x1={xNegative} x2={xNegative} y1={y1} y2={y2} stroke="#D3D3D3" strokeWidth={barW}>
                <title>-Sigma: {-sigma}</title>
            </line>
            {makeXAxisLabel(mean - sigma, xNegative, y1)}
            </>);
        return bars;
    }

    function addPoint(point: ControlPoint) {
        controlPoints.push(point);
        handleChangePoints(controlPoints);
    }

    function createPoint() {
        const svgP = (myRef.current as any).createSVGPoint();
        // TODO: address constant 100, perhaps should be dependant on baseline, roof etc.
        const { h, b, r, w, p } = _getLineGraphAttributes();
        svgP.x = w - 100;
        svgP.y = h - 100;

        const points = controlPoints;
        const newPointData = svgCoordsToPointData(Vec2.create(svgP.x, svgP.y));
        const sorted = sortPointsByXValues(points);
        const realPoints = sorted.filter(p => p.isTerminal === undefined || false);
        const maxIndex = realPoints[realPoints.length - 1].index;
        const newPoint: ControlPoint = {
            data: newPointData,
            id: UUID.create22(),
            color: getRandomColor(),
            index: maxIndex + 1
        };
        addPoint(newPoint);
    }

    // TODO: fix movement of plus and minus sign buttons upon creation of points
    function removeRightmostPoint() {
        const points = controlPoints;
        const sortedPs = sortPointsByXValues(points);
        const rightmostP = sortedPs[sortedPs.length - 2];
        removePoint(rightmostP.id);
    };

    function refCallBack(element: any) {
        if (element) {
            myRef.current = element;
        }
    }

    function handleEnter() {
        document.removeEventListener('mousemove', handleDrag, true);
        document.removeEventListener('mouseup', handlePointUpdate, true);
    }

    const handleKeyDown = (event: any) => {
        // TODO: set canSelectMultiple = true
        if (event.key === 'Shift') {
            setCanSelectMultiple(true);
            // TODO: allow to select another point
        } else {
            // console.log(`${event} is pressed`);
        }
    };

    const handleKeyUp = (event: any) => {
        // TODO: SET canSelectMultiple = fasle
        setCanSelectMultiple(false);
        if (event.shiftKey) {
            // TODO: allow to select another point
        } else {
            // console.log(`${event} is released`);
        }
    };

    const updateColor: ParamOnChange = ({ value }: { value: Color }) => {
        const currentPoints = controlPoints;
        if (!clickedPointIds || clickedPointIds.length === 0) throw Error('No point is selected');
        const clickedPoints = getPoints(clickedPointIds);
        if (!clickedPoints) throw Error('Point should be selected');
        for (const point of clickedPoints) {
            point.color = value;
        }
        const clickedPointsIds = clickedPoints.map(p => p.id);
        const notClickedPoints = currentPoints.filter(p => !clickedPointsIds.includes(p.id));
        const newPoints = clickedPoints.concat(notClickedPoints);
        handleChangePoints(newPoints);
    };

    const onColorSquareClick = (pointId: UUID) => {
        // TODO: change color
        // use .change method
        // perhaps handle Change point or something
        setClickedPointIds([pointId]);
        toggleColorPicker();
    };

    const LineGraphRendered = () => {
        const points = renderPoints();
        const baseline = renderBaseline();
        const lines = renderLines();
        const histogram = renderHistogram();
        const axes = renderAxes();
        const gridLines = renderGridLines();
        const descriptiveStatisticsBars = renderDescriptiveStatisticsBars();
        const firstPoint = clickedPointIds.length > 0 ? getPoint(clickedPointIds[0]) : void 0;
        const color = firstPoint ? firstPoint.color : void 0;
        const defaultColor = ParamDefinition.Color(ColorNames.black);
        const _createPoint = createPoint;
        const _removePoint = removeRightmostPoint;
        const plusIconButtonStyle = {
            position: ('absolute' as any),
            top: '100px',
            right: '30px',
        };
        const minusIconButtonStyle = {
            position: ('absolute' as any),
            top: '100px',
            right: '5px',
        };

        const { h, b, r, w, p } = _getLineGraphAttributes();

        const rendered = ([
            <div key="LineGraph">
                <IconButton style={plusIconButtonStyle} small svg={PlusBoxSvg} onClick={function (e: React.MouseEvent<HTMLButtonElement>): void {
                    _createPoint();
                } } ></IconButton>
                <IconButton style={minusIconButtonStyle} small svg={MinusBoxSvg} onClick={function (e: React.MouseEvent<HTMLButtonElement>): void {
                    _removePoint();
                } } ></IconButton>
                <svg
                    className="msp-canvas"
                    ref={refCallBack}
                    viewBox={`0 0 ${w + p} ${h + p}`}
                    onMouseMove={handleDrag}
                    onMouseUp={handlePointUpdate}
                    onMouseLeave={handleLeave}
                    onMouseEnter={handleEnter}
                    tabIndex={0}
                    onKeyDown={handleKeyDown}
                    onKeyUp={handleKeyUp}
                    // onDoubleClick={handleDoubleClick}
                >
                    <g stroke="black" fill="black">
                        {baseline}
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
                <>{showColorPicker ? <ColorPicker isActive={showColorPicker} defaultColor={defaultColor}
                    color={color} updateColor={updateColor} toggleColorPicker={toggleColorPicker}/> :
                    // TODO: fix this div
                    <div>Select point to change color</div>}
                </>
                <>
                    {/* <TFParamsWrapper onChange={this.setTF} descriptiveStatistics={this.descriptiveStatistics}></TFParamsWrapper> */}
                    {/* <Button onClick={this.deleteAllPoints}>Remove All Points</Button> */}
                    <PointsPanel points={controlPoints}
                        onPointIndexClick={highlightPoint}
                        removeAllPoints={deleteAllPoints}
                        onExpandGroupOpen={props.onExpandGroupOpen}
                        changeXValue={changeXValue}
                        changeAlphaValue={changeAlphaValue}
                        onPointButtonClick={removePoint}
                        onColorSquareClick={onColorSquareClick}></PointsPanel>
                    <HelpersPanel onChange={setTF}
                        methods={methodsParams}
                        descriptiveStatistics={descriptiveStatistics}></HelpersPanel>
                </>
            </div>,
            <div key="modal" id="modal-root" />
        ]);

        return rendered;
    }

    return LineGraphRendered();
}

// export class LineGraphComponent extends React.Component<any, LineGraphComponentState> {
//     private myRef: any;
//     private height: number;
//     private width: number;
//     private padding: number;
//     private updatedX: number;
//     private updatedY: number;
//     private selectedPointId?: UUID;
//     private ghostPoints: SVGElement[];
//     private gElement: SVGElement;
//     private namespace: string;
//     private descriptiveStatistics: VolumeDescriptiveStatistics;
//     // private offsetY: number;
//     private baseline: number;
//     private baselineUnnormalized: number;
//     private roof: number;
//     private roofUnnormalized: number;
//     constructor(props: any) {
//         super(props);
//         this.myRef = React.createRef();
//         this.state = {
//             points: startEndPoints.concat(this.props.data),
//             copyPoint: undefined,
//             canSelectMultiple: false,
//             showColorPicker: false,
//             methodsParams: DefaultTFParams
//             // colored: this.props.colored
//         };
//         this.roof = LineGraphParams.roof;
//         this.roofUnnormalized = LineGraphParams.roofUnnormalized;
//         this.baseline = LineGraphParams.baseline;
//         this.baselineUnnormalized = LineGraphParams.baselineUnnormalized;
//         this.height = LineGraphParams.height;
//         this.width = LineGraphParams.width;
//         this.padding = LineGraphParams.padding;
//         // this.offsetY = LineGraphParams.offsetY;
//         // this.offsetY = 0;
//         this.selectedPointId = undefined;
//         this.ghostPoints = [];
//         this.namespace = 'http://www.w3.org/2000/svg';
//         this.descriptiveStatistics = this.getDescriptiveStatistics();
//         this.sortPointsByXValues(this.state.points);
//         this.highlightPoint = this.highlightPoint.bind(this);
//         this.removePoint = this.removePoint.bind(this);
//         this.removeRightmostPoint = this.removeRightmostPoint.bind(this);
//         this.createPoint = this.createPoint.bind(this);
//         this.handleDrag = this.handleDrag.bind(this);
//         this.handleMultipleDrag = this.handleMultipleDrag.bind(this);
//         this.handleDoubleClick = this.handleDoubleClick.bind(this);
//         this.refCallBack = this.refCallBack.bind(this);
//         this.handlePointUpdate = this.handlePointUpdate.bind(this);
//         this.change = this.change.bind(this);
//         this.handleKeyUp = this.handleKeyUp.bind(this);
//         this.handleLeave = this.handleLeave.bind(this);
//         this.handleEnter = this.handleEnter.bind(this);
//         this.setTF = this.setTF.bind(this);
//         this.changeXValue = this.changeXValue.bind(this);
//         this.changeAlphaValue = this.changeAlphaValue.bind(this);
//         this.deleteAllPoints = this.deleteAllPoints.bind(this);
//         this.toggleColorPicker = this.toggleColorPicker.bind(this);
//         this.onColorSquareClick = this.onColorSquareClick.bind(this);
//     }

//     private getDescriptiveStatistics() {
//         const s = this.props.volume.grid.stats;
//         const d: VolumeDescriptiveStatistics = {
//             min: s.min,
//             max: s.max,
//             mean: s.mean,
//             sigma: s.sigma
//         };
//         return d;
//     }

//     private deleteAllPoints() {
//         const points = this.state.points;
//         const ghostPoints = [];
//         for (const point of points) {
//             if (point.isTerminal === true) { ghostPoints.push(point); }
//         };
//         const ghostPointsSorted = this.sortPointsByXValues(ghostPoints);
//         this.setState({ points: ghostPointsSorted, clickedPointIds: undefined, showColorPicker: false });
//         this.change([]);
//     }

//     private highlightPoint(pointId: UUID) {
//         const points = this.state.points;
//         const targetPoint = points.find(p => p.id === pointId);
//         if (!targetPoint) throw Error('Cannot highlight inexisting point exist');
//         // TODO: highlight
//     }

//     private _setMethod2TF(params: Method2ParamsValues) {
//         // TODO: implement
//     }

//     private _setDefaultsTF(params: DefaultParamsValues) {
//         const paddingUnnormalized = this.padding / this.height;
//         // there two are the same values, reduce to a single one, e.g. could be be just padding, why not
//         // should be e.g.
//         const points = generateControlPoints(paddingUnnormalized, ColorNames.black, undefined, this.props.volume);
//         const currentPoints = this.state.points;
//         points.push(currentPoints[currentPoints.length - 1]);
//         points.unshift(currentPoints[0]);
//         this.setState({
//             points: points
//         });
//         this.change(points);
//         // TODO: set points to points similar to the below function
//     }

//     private _setGaussianTF(params: GaussianTFParamsValues) {
//         // TODO: simplify this
//         const offsetYNormalized = this.padding / 2;
//         const paddingNormalized = this.padding;
//         // will be 0
//         const offsetYUnnormalized = (offsetYNormalized - this.padding / 2) / this.height;
//         // will be ~0.09
//         const paddingUnnormalized = (paddingNormalized - this.padding / 2) / this.height;
//         // normalizedY = ([0;1] value * (this.height)) + offset
//         const { name, gaussianExtent, gaussianCenter, gaussianHeight } = params;
//         const a = gaussianHeight;
//         const min = this.descriptiveStatistics.min;
//         const max = this.descriptiveStatistics.max;
//         // const sigma = this.descriptiveStatistics.sigma;
//         const TFextent = max - min;
//         // make it spread the points equally
//         // console.log(this.descriptiveStatistics);
//         const b = gaussianCenter / TFextent;
//         // const b = gaussianCenter * (mean + sigma) / TFextent;
//         const c = gaussianExtent / TFextent;
//         // const l = (this.width * 2 * c / extent);
//         // console.log(a, b, c);
//         const gaussianPoints: ControlPoint[] = generateGaussianControlPoints(a, b, c, TFextent, offsetYUnnormalized, this.props.volume, paddingUnnormalized);
//         const currentPoints = this.state.points;
//         gaussianPoints.push(currentPoints[currentPoints.length - 1]);
//         gaussianPoints.unshift(currentPoints[0]);
//         this.setState({
//             points: gaussianPoints
//         });
//         // TODO: maybe add handleChange points instead?
//         this.change(gaussianPoints);
//     }

//     changeAlphaValue(pointId: UUID, alpha: number) {
//         const points = this.state.points;
//         const modifiedPoints = points.map(p => {
//             if (p.id === pointId) {
//                 const modifiedData: ControlPointData = {
//                     x: p.data.x,
//                     alpha: alpha
//                 };
//                 const modifiedP: ControlPoint = {
//                     id: p.id,
//                     color: p.color,
//                     index: p.index,
//                     data: modifiedData
//                 };
//                 return modifiedP;
//             } else {
//                 return p;
//             }
//         });
//         console.log(modifiedPoints);
//         this.handleChangePoints(modifiedPoints);
//     }

//     // TODO: refactor to remove duplication with above function
//     changeXValue(pointId: UUID, absValue: number) {
//         // got it, we have to convert abs value back to relative [0;1] value
//         // need value to point function
//         // now need to pass absValueToPointValue
//         // so that it reaches here;
//         const x = this.props.onAbsValueToPointValue(absValue);
//         const points = this.state.points;
//         const modifiedPoints = points.map(p => {
//             if (p.id === pointId) {
//                 const modifiedData: ControlPointData = {
//                     x: x,
//                     alpha: p.data.alpha
//                 };
//                 const modifiedP: ControlPoint = {
//                     id: p.id,
//                     color: p.color,
//                     index: p.index,
//                     data: modifiedData
//                 };
//                 return modifiedP;
//             } else {
//                 return p;
//             }
//         });
//         console.log(modifiedPoints);
//         this.handleChangePoints(modifiedPoints);
//     }

//     private setTF(params: TFParamsValues) {
//         const name = params.name;
//         switch (name) {
//             // TODO: separate into type
//             case 'gaussian':
//                 this._setGaussianTF(params as GaussianTFParamsValues);
//                 break;
//             case 'method2':
//                 this._setMethod2TF(params as Method2ParamsValues);
//                 break;
//             case 'defaults':
//                 this._setDefaultsTF(params as DefaultParamsValues);
//                 break;
//             default: throw Error(`Transfer function ${name} is not supported`);
//         };
//     }

//     private getPoint(id: UUID) {
//         const points = this.state.points;
//         const selectedPoints = points.find(p => p.id === id);
//         if (!selectedPoints) throw Error(`Points with ids ${id} do not exist`);
//         return selectedPoints;
//     }

//     private getPoints(id: UUID[]) {
//         const points = this.state.points;
//         const selectedPoints = points.filter(p => id.includes(p.id));
//         if (!selectedPoints) throw Error(`Point with id ${id} does not exist`);
//         return selectedPoints;
//     }

//     public render() {
//         const points = this.renderPoints();
//         const baseline = this.renderBaseline();
//         const lines = this.renderLines();
//         const histogram = this.renderHistogram();
//         const axes = this.renderAxes();
//         const gridLines = this.renderGridLines();
//         // const pointsButtons = this.renderPointsButtons();
//         const descriptiveStatisticsBars = this.renderDescriptiveStatisticsBars();
//         const firstPoint = this.state.clickedPointIds ? this.getPoint(this.state.clickedPointIds[0]) : void 0;
//         const color = firstPoint ? firstPoint.color : void 0;
//         const defaultColor = ParamDefinition.Color(Color(0x121212));
//         // bind this
//         const _createPoint = this.createPoint;
//         const _removePoint = this.removeRightmostPoint;
//         const plusIconButtonStyle = {
//             position: ('absolute' as any),
//             top: '100px',
//             right: '30px',
//             // transform: 'translate(-50%, -50%)'
//         };
//         // TODO: fix style
//         const minusIconButtonStyle = {
//             position: ('absolute' as any),
//             top: '100px',
//             right: '5px',
//             // transform: 'translate(-50%, -50%)'
//         };

//         return ([
//             <div key="LineGraph">
//                 <IconButton style={plusIconButtonStyle} small svg={PlusBoxSvg} onClick={function (e: React.MouseEvent<HTMLButtonElement>): void {
//                     // throw new Error('Function not implemented.');
//                     _createPoint();
//                 } } ></IconButton>
//                 <IconButton style={minusIconButtonStyle} small svg={MinusBoxSvg} onClick={function (e: React.MouseEvent<HTMLButtonElement>): void {
//                     // throw new Error('Function not implemented.');
//                     _removePoint();
//                 } } ></IconButton>
//                 <svg
//                     className="msp-canvas"
//                     ref={this.refCallBack}
//                     // TODO: can e.g. change viewbox
//                     viewBox={`0 0 ${this.width + this.padding} ${this.height + this.padding}`}
//                     onMouseMove={this.handleDrag}
//                     onMouseUp={this.handlePointUpdate}
//                     onMouseLeave={this.handleLeave}
//                     onMouseEnter={this.handleEnter}
//                     tabIndex={0}
//                     onKeyDown={this.handleKeyDown}
//                     onKeyUp={this.handleKeyUp}
//                     onDoubleClick={this.handleDoubleClick}>
//                     {/* renders points */}
//                     <g stroke="black" fill="black">
//                         {baseline}
//                         {histogram}
//                         {lines}
//                         {points}
//                         {descriptiveStatisticsBars}
//                         {axes}
//                         {gridLines}
//                     </g>
//                     <g className="ghost-points" stroke="black" fill="black">
//                     </g>
//                 </svg>
//                 {/* TODO: div with the same height as color picker */}
//                 <>{this.state.showColorPicker ? <ColorPicker isActive={this.state.showColorPicker} defaultColor={defaultColor} color={color} updateColor={this.updateColor} toggleColorPicker={this.toggleColorPicker}/> : <div>Select point to change color</div>}

//                 </>
//                 <>
//                     {/* <TFParamsWrapper onChange={this.setTF} descriptiveStatistics={this.descriptiveStatistics}></TFParamsWrapper> */}
//                     {/* <Button onClick={this.deleteAllPoints}>Remove All Points</Button> */}
//                     <PointsPanel points={this.state.points}
//                         onPointIndexClick={this.highlightPoint}
//                         removeAllPoints={this.deleteAllPoints}
//                         onExpandGroupOpen={this.props.onExpandGroupOpen}
//                         changeXValue={this.changeXValue}
//                         changeAlphaValue={this.changeAlphaValue}
//                         onPointButtonClick={this.removePoint}
//                         onColorSquareClick={this.onColorSquareClick}></PointsPanel>
//                     {/* Connect this to existing code */}
//                     {/* should render actual state of all methods
//                     // add this to state then, params of all methods */}
//                     {/* TODO: debug */}
//                     <HelpersPanel onChange={this.setTF} methods={this.state.methodsParams} descriptiveStatistics={this.descriptiveStatistics}></HelpersPanel>
//                 </>
//             </div>,
//             <div key="modal" id="modal-root" />
//         ]);
//     }

//     toggleColorPicker = () => {
//         this.setState({ showColorPicker: this.state.showColorPicker === true ? false : true });
//     };

//     private onColorSquareClick = (pointId: UUID) => {
//         this.setState({ clickedPointIds: [pointId] });
//         this.toggleColorPicker();
//     };

//     componentDidMount() {
//         this.gElement = document.getElementsByClassName('ghost-points')[0] as SVGElement;
//         // this.setState({
//         //     colored: this.props.colored
//         // });
//     }

//     private change(points: ControlPoint[]) {
//         const copyPoints = points.slice();
//         copyPoints.shift();
//         copyPoints.pop();
//         this.props.onChange(copyPoints);
//     }

//     private updateColor: ParamOnChange = ({ value }: { value: Color }) => {
//         const clickedPointIds = this.state.clickedPointIds;
//         const currentPoints = this.state.points;
//         if (!clickedPointIds || clickedPointIds.length === 0) throw Error('No point is selected');
//         const clickedPoints = this.getPoints(clickedPointIds);
//         if (!clickedPoints) throw Error('Point should be selected');
//         for (const point of clickedPoints) {
//             point.color = value;
//         }
//         // need to change just points in ps
//         const clickedPointsIds = clickedPoints.map(p => p.id);
//         const notClickedPoints = currentPoints.filter(p => !clickedPointsIds.includes(p.id));
//         const newPoints = clickedPoints.concat(notClickedPoints);
//         console.log('Points before color change', currentPoints);
//         console.log('Points with color changed', newPoints);
//         this.handleChangePoints(newPoints);
//         debugger;
//     };

//     private handleKeyDown = (event: any) => {
//         // TODO: set canSelectMultiple = true
//         if (event.key === 'Shift') {
//             this.setState({ canSelectMultiple: true });
//             console.log('Shift is pressed');
//             // TODO: allow to select another point
//         } else {
//             // console.log(`${event} is pressed`);
//         }
//     };

//     private handleKeyUp = (event: any) => {
//         // TODO: SET canSelectMultiple = fasle
//         this.setState({ canSelectMultiple: false });
//         if (event.shiftKey) {
//             // do something
//             console.log('Shift is released');
//             // TODO: allow to select another point
//         } else {
//             // console.log(`${event} is released`);
//         }
//     };

//     private handleClick = (point: ControlPoint) => (event: any) => {
//         this.setState({ clickedPointIds: [point.id] });

//         if (this.state.canSelectMultiple) {
//             if (event.shiftKey) return;
//             // TODO: function to execute on this
//         }

//         if (this.state.showColorPicker) {
//             // TODO: what should happen there?
//         } else {
//             this.toggleColorPicker();
//         }
//     };

//     private sortPointsByXValues(points: ControlPoint[]) {
//         points.sort((a, b) => {
//             if (a.data.x === b.data.x) {
//                 if (a.data.x === 0) {
//                     return a.data.alpha - b.data.alpha;
//                 }
//                 if (a.data.alpha === 1) {
//                     return b.data.alpha - a.data.alpha;
//                 }
//                 return a.data.alpha - b.data.alpha;
//             }
//             return a.data.x - b.data.x;
//         });
//         return points;
//     }

//     private sortPointsByIndices(points: ControlPoint[]) {
//         points.sort((a, b) => {
//             return a.index - b.index;
//             // if (a.data.x === b.data.x) {
//             //     if (a.data.x === 0) {
//             //         return a.data.alpha - b.data.alpha;
//             //     }
//             //     if (a.data.alpha === 1) {
//             //         return b.data.alpha - a.data.alpha;
//             //     }
//             //     return a.data.alpha - b.data.alpha;
//             // }
//             // return a.data.x - b.data.x;
//         });
//         return points;
//     }

//     private handleChangePoints(points: ControlPoint[]) {
//         const pointsSorted = this.sortPointsByXValues(points);
//         this.setState({ points: pointsSorted });
//         this.change(points);
//     }

//     private handleMouseDown = (point: ControlPoint) => (event: any) => {
//         const { id, isTerminal } = point;
//         if (isTerminal === true) {
//             return;
//         }
//         if (this.state.canSelectMultiple) {
//             return;
//         }

//         const copyPoint: Vec2 = this.controlPointDataToSVGCoords(this.getPoint(id).data);
//         this.ghostPoints.push(document.createElementNS(this.namespace, 'circle') as SVGElement);
//         this.ghostPoints[0].setAttribute('r', '10');
//         this.ghostPoints[0].setAttribute('fill', 'orange');
//         this.ghostPoints[0].setAttribute('cx', `${copyPoint[0]}`);
//         this.ghostPoints[0].setAttribute('cy', `${copyPoint[1]}`);
//         this.ghostPoints[0].setAttribute('style', 'display: none');
//         this.gElement.appendChild(this.ghostPoints[0]);
//         this.updatedX = copyPoint[0];
//         this.updatedY = copyPoint[1];
//         this.selectedPointId = point.id;
//     };

//     private handleDrag(event: any) {
//         if (this.selectedPointId === undefined) {
//             return;
//         }

//         // here use adjusted alpha
//         const pt = this.myRef.createSVGPoint();
//         let updatedCopyPoint;
//         const padding = this.padding / 2;
//         // this is true coords from the top left corner of svg
//         pt.x = event.clientX;
//         pt.y = event.clientY;
//         // still wrong poisition, shifted above on drag
//         const svgP = pt.matrixTransform(this.myRef.getScreenCTM().inverse());
//         // coords from top left corner of svg
//         updatedCopyPoint = Vec2.create(svgP.x, svgP.y);
//         console.log('created a svg point with coords', svgP.x, svgP.y);
//         if ((svgP.x < (padding) || svgP.x > (this.width + (padding))) &&
//         (svgP.y > (this.height - this.baseline) || svgP.y < (0))) {
//             debugger;
//             updatedCopyPoint = Vec2.create(this.updatedX, this.updatedY);
//         // right border is X > something
//         } else if (svgP.x < padding) {
//             debugger;
//             // TODO: fix lines to start from true 0
//             updatedCopyPoint = Vec2.create(padding, svgP.y);
//             // this probably
//         } else if (svgP.x > (this.width + (padding))) {
//             // this
//             console.log('Creating a point while crossing right border');
//             // should create on the same x that is fine, but on correct y
//             console.log(this.width + padding, svgP.y);
//             updatedCopyPoint = Vec2.create(this.width + padding, svgP.y);
//             // touch this
//         } else if (svgP.y > (this.height - this.baseline)) {
//             debugger;
//             // does not allow into such area
//             // fix viewbox or?
//             updatedCopyPoint = Vec2.create(svgP.x, this.height - this.baseline);
//         } else if (svgP.y < (this.roof)) {
//             debugger;
//             updatedCopyPoint = Vec2.create(svgP.x, this.roof);
//         } else {
//             // console.log('Point within the constraints, creating normally');
//             updatedCopyPoint = Vec2.create(svgP.x, svgP.y);
//         }

//         this.updatedX = updatedCopyPoint[0];
//         this.updatedY = updatedCopyPoint[1];
//         // TODO
//         const unNormalizePoint = this.svgCoordsToPointData(updatedCopyPoint);
//         // TODO: this.ghostPoints[0] is undefined when dragging a point towards right
//         // border, why?
//         this.ghostPoints[0].setAttribute('style', 'display: visible');
//         this.ghostPoints[0].setAttribute('cx', `${updatedCopyPoint[0]}`);
//         this.ghostPoints[0].setAttribute('cy', `${updatedCopyPoint[1]}`);


//         this.props.onDrag(unNormalizePoint);
//         // ok got it then point is "normalized" and thus got raised up or something
//     }

//     private handleMultipleDrag() {
//         // TODO
//     };

//     private addPoint(point: ControlPoint) {
//         this.state.points.push(point);
//         this.handleChangePoints(this.state.points);
//     }

//     private replacePoint(point: ControlPoint) {
//         // can get id from it
//         let points = this.state.points.filter(p => p.id !== point.id);
//         points.push(point);
//         points = this.sortPointsByXValues(points);
//         this.setState({
//             points,
//         });
//         this.change(points);
//     }

//     private handlePointUpdate(event: any) {
//         const selected = this.selectedPointId;
//         if (this.state.canSelectMultiple) {
//             return;
//         }
//         // or here
//         if (selected === undefined || this.getPoint(selected).isTerminal === true) {
//             this.setState({
//                 copyPoint: undefined,
//             });
//             return;
//         }
//         // "unselects" points
//         // do not need to raise points, simply generate their original position
//         // as x + offset
//         this.selectedPointId = undefined;
//         const pointData = this.getPoint(selected);
//         // I guess this one
//         const updatedPointData = this.svgCoordsToPointData(Vec2.create(this.updatedX, this.updatedY));
//         const updatedPoint: ControlPoint = {
//             data: updatedPointData,
//             id: pointData.id,
//             color: pointData.color,
//             index: pointData.index
//         };

//         this.replacePoint(updatedPoint);
//         this.gElement.innerHTML = '';
//         this.ghostPoints = [];
//         document.removeEventListener('mousemove', this.handleDrag, true);
//         document.removeEventListener('mouseup', this.handlePointUpdate, true);
//     }

//     private removePoint(id: UUID) {
//         const point = this.getPoint(id);
//         if (point.isTerminal === true) { return; }
//         let points = this.state.points.filter(p => p.id !== point.id);
//         points = this.sortPointsByXValues(points);
//         this.setState({ points: points, clickedPointIds: undefined, showColorPicker: false });
//         this.change(points);
//     }

//     // TODO: fix movement of plus and minus sign buttons upon creation of points
//     private removeRightmostPoint() {
//         // remove based on data
//         const points = this.state.points;
//         const sortedPs = this.sortPointsByXValues(points);
//         debugger;
//         // rightmost is undefined
//         // last is ghost
//         // const rightmostP = sortedPs.slice(-1)[0];
//         const rightmostP = sortedPs[sortedPs.length - 2];
//         debugger;
//         this.removePoint(rightmostP.id);
//         debugger;
//     };

//     private createPoint() {
//         const svgP = this.myRef.createSVGPoint();
//         // TODO: resolve location, do something like this.width, height...
//         svgP.x = this.width - 100;
//         svgP.y = this.height - 100;

//         // const svgP = pt.matrixTransform(this.myRef.getScreenCTM().inverse());
//         const points = this.state.points;
//         const newPointData = this.svgCoordsToPointData(Vec2.create(svgP.x, svgP.y));
//         const sorted = this.sortPointsByXValues(points);
//         const realPoints = sorted.filter(p => p.isTerminal === undefined || false);
//         const maxIndex = realPoints[realPoints.length - 1].index;
//         // const
//         const newPoint: ControlPoint = {
//             data: newPointData,
//             id: UUID.create22(),
//             color: getRandomColor(),
//             index: maxIndex + 1
//         };
//         this.addPoint(newPoint);
//     }

//     private handleDoubleClick(event: any) {
//         return;
//         // const pt = this.myRef.createSVGPoint();
//         // pt.x = event.clientX;
//         // pt.y = event.clientY;
//         // const svgP = pt.matrixTransform(this.myRef.getScreenCTM().inverse());
//         // const points = this.state.points;
//         // const padding = this.padding / 2;

//         // if (svgP.x < (padding) ||
//         //     svgP.x > (this.width + (padding)) ||
//         //     svgP.y > (this.height + (padding)) ||
//         //     svgP.y < (this.padding / 2)) {
//         //     return;
//         // }
//         // const newPointData = this.unNormalizePoint(Vec2.create(svgP.x, svgP.y));
//         // // problem with index
//         // // should be the last one
//         // // find the last index and do + 1?
//         // const newPoint: ControlPoint = {
//         //     data: newPointData,
//         //     id: UUID.create22(),
//         //     color: getRandomColor(),
//         //     index: points.slice(-1)[0].index + 1
//         // };
//         // this.addPoint(newPoint);
//     }

//     private deletePoint = (id: UUID) => (event: any) => {
//         this.removePoint(id);
//         console.log('Point', id, ' is deleted');
//         event.stopPropagation();
//     };

//     private handleLeave() {
//         if (this.selectedPointId === undefined) {
//             return;
//         }

//         document.addEventListener('mousemove', this.handleDrag, true);
//         document.addEventListener('mouseup', this.handlePointUpdate, true);
//     }

//     private handleEnter() {
//         document.removeEventListener('mousemove', this.handleDrag, true);
//         document.removeEventListener('mouseup', this.handlePointUpdate, true);
//     }

//     // this function prouces actual SVG coordinates (pixels from top left corner)
//     private controlPointDataToSVGCoords(controlPointData: ControlPointData): Vec2 {
//         console.log('controlPointDataToSVGCoords');
//         const offset = this.padding / 2;
//         const maxX = this.width + offset;
//         // fix this
//         const maxY = this.height + offset;
//         // normalizedY = ([0;1] value * (this.height)) + offset
//         const normalizedX = (controlPointData.x * (maxX - offset)) + offset;
//         // const normalizedY = (controlPointData.alpha * (maxY - offset)) + offset;
//         debugger;
//         // I guess should not be adjusted here or?
//         // if (controlPointData.adjustedAlpha === void 0) {
//         //     debugger;
//         //     console.log(controlPointData);
//         //     // TODO: add this adjusted alpha thing to copy point etc.
//         //     // fix lines
//         //     throw Error('Adjusted alpha should be provided');
//         // }
//         const space = this.height - this.baseline - this.roof;
//         // should be reversed or?
//         // before
//         // normalizedY = (point[1] * (maxY - o)) + o;
//         // reverseY = this.height + this.padding - normalizedY;

//         // Ok found a mistake
//         // alpha / adjusted alpha is proportion
//         // so we should stick to using just alpha
//         // adjustment is not needed anymore
//         const normalizedY = controlPointData.alpha * space;
//         // ok so here it is e.g. 0.2 alpha
//         // thefn space is e.g 350
//         // 0.2 * 350 = 70
//         // 70 is linear size from the baseline to the point
//         // nothing to do with the roof
//         // 70 so we need to place this 70
//         // such that it is calculated from the top corner
//         // i.e. get the distance from 0 0 to this 70
//         // to do this
//         const fromTopToBaseline = this.height - this.baseline;
//         const correctY = fromTopToBaseline - normalizedY;
//         // right so if this above is 200 and e.g. 300
//         // (based on e.g. 0.2 and 0.3 adjusted alpha)
//         // below it will become 350 - 200 = 150 and 350 - 300 = 50
//         // i.e. need to be reversed
//         // const y = (this.height - this.baseline) + normalizedY;
//         // const reverseY = this.hegight
//         // try adding removing baseline here or whatever
//         // const normalizedY = (controlPointData.adjustedAlpha * (maxY - offset)) + offset;
//         // const reverseY = (this.height + this.padding) - normalizedY;
//         // here we should provide something like
//         // baseline +
//         const newPoint = Vec2.create(normalizedX, correctY);
//         return newPoint;
//     }

//     // TODO: TODO: TODO: modify this too
//     private svgCoordsToPointData(vec2: Vec2): ControlPointData {
//         const space = this.height - this.baseline - this.roof;
//         const min = this.padding / 2;
//         const maxX = this.width + min;
//         const maxY = this.height + min;
//         const unNormalizedX = (vec2[0] - min) / (maxX - min);
//         const cartesianDistanceBetweenPointAndBaseline = (this.height - this.baseline - vec2[1]);
//         const dividedBySpace = cartesianDistanceBetweenPointAndBaseline / space;
//         const unNormalizedY = dividedBySpace;
//         return {
//             x: unNormalizedX,
//             alpha: unNormalizedY
//         };
//     }

//     private refCallBack(element: any) {
//         if (element) {
//             this.myRef = element;
//         }
//     }

//     private renderGridLines() {
//         const count = 5;
//         const bars: any = [];
//         const offset = this.padding / 2;
//         const x1 = offset;
//         const w = 0.5;
//         // TODO: allow editing text field, may need some function similar to what is with RGB thing
//         //
//         // TODO: consider adding baseline height as attribute of LineGraphComponent
//         const x2 = this.width + offset;
//         // fix that
//         const overallHeight = this.height;
//         // TODO: limit the height of histogram bars
//         // 1, 2, 3, should be 4
//         // get histogram bars back
//         for (let i = 0; i < count; ++i) {
//             // adjust + - 1 etc.
//             // since we have risen all elements, need to probably increase height or
//             // decrease it by - offset everywhere where it is used
//             const y = overallHeight * i / count;
//             // const y = this.height + offset * 2 - i * this.height / count;
//             bars.push(
//                 <line key={`${1 - i / count}gridline`} x1={x1} x2={x2} y1={y} y2={y}
//                     stroke="#A9A9A9" strokeWidth={w} strokeDasharray="15, 15">
//                     <title>{(1 - i / count).toFixed(1)}</title>
//                 </line>
//             );
//         }
//         return bars;
//     }

//     private renderHistogram() {
//         if (!this.props.volume) return null;
//         const histogram = Grid.getHistogram(this.props.volume.grid, 40);
//         const bars = [];
//         const N = histogram.counts.length;
//         const w = this.width / N;
//         const offset = this.padding / 2;
//         const max = arrayMax(histogram.counts) || 1;
//         for (let i = 0; i < N; i++) {
//             const fromValue = histogram.min + histogram.binWidth * i;
//             const toValue = histogram.min + histogram.binWidth * (i + 1);
//             const x = this.width * i / (N - 1) + offset;
//             // add roof
//             // MAYBE remove offset from here
//             // const y1 = this.height + offset - this.baseline;// + 2 * offset;
//             // lowest point
//             // for some reason, maybe this got divided by 2 somewhere
//             const y1 = this.height - this.baseline;
//             // highest point
//             // 1. get the distance between roofline and baseline
//             const space = this.height - this.baseline - this.roof;
//             // 2. get a fraction of that distance based on histogram counts etc.
//             // 2.5. do 1 - that fraction
//             const fraction = space * (1 - histogram.counts[i] / max);
//             // 3. add it to this.roof
//             const y2 = fraction + this.roof;
//             // you got y2
//             // const
//             // const y2 = (this.height - this.baseline + this.roof) * (1 - histogram.counts[i] / max);// + 2 * offset;
//             // console.log(y1, y2);
//             // console.log(this.height, this.roof, this.baseline);
//             bars.push(<line key={`histogram${i}`} x1={x} x2={x} y1={y1} y2={y2} stroke="#A9A9A9" strokeWidth={w}>
//                 <title>[{fromValue}; {toValue}]</title>
//             </line>);
//         }
//         return bars;
//     }

//     private renderAxes() {
//         if (!this.props.volume) return null;
//         const offset = this.padding / 2;
//         const mean = this.descriptiveStatistics.mean;
//         const min = this.descriptiveStatistics.min;
//         const max = this.descriptiveStatistics.max;
//         // X: horizontal bar with arrow
//         // need min max
//         const x1HorizontalBar = offset;
//         const x2HorizontalBar = this.width + offset;
//         const y1HorizontalBar = this.height + offset;
//         const y2HorizontalBar = this.height + offset;
//         // vertical bar with arrow
//         // x and y letters



//         const x2VerticalBar = offset;
//         const x1VerticalBar = offset;
//         const y2VerticalBar = 25;
//         const y1VerticalBar = this.height + offset;
//         const w = offset / 10;
//         const bars = [];
//         // const hd = `M ${xh1} ${yh1} L ${xh2} ${yh2} L ${xh3} ${yh3} L ${xh4} ${yh4} L ${xh5} ${yh5} Z`;
//         // const arrowheadd = `M ${xah1} ${yah1} L ${xah2} ${yah2} L ${xah3} ${yah3} L ${xah4} ${yah4} L ${xah5} ${yah5} Z`;
//         bars.push(
//             // TODO: color
//             // TODO: may not work with <>/</>
//             <>
//                 <defs>
//                     <marker
//                         id="head-horizontal"
//                         viewBox="0 0 10 10"
//                         refX="5"
//                         refY="5"
//                         markerWidth="6"
//                         markerHeight="6"
//                         orient="auto-start-reverse">
//                         <path d="M 0 2 L 10 5 L 0 8 z" />
//                     </marker>
//                 </defs>

//                 <line key={'horizontalAxis'} x1={x1HorizontalBar} x2={x2HorizontalBar}
//                     y1={y1HorizontalBar} y2={y2HorizontalBar} stroke="#7d7f7c" strokeWidth={w} markerEnd="url(#head-horizontal)">
//                 </line>
//             </>
//         );
//         // raise histogram above the base line somehow
//         bars.push(
//             <>
//                 <defs>
//                     <marker
//                         id="head-vertical"
//                         viewBox="0 0 10 10"
//                         refX="5"
//                         refY="5"
//                         markerWidth="6"
//                         markerHeight="6"
//                         orient="auto-start-reverse">
//                         <path d="M 0 2 L 10 5 L 0 8 z" />
//                         {/* <path d="M 2 10 L 5 0 L 8 10 z" /> */}
//                     </marker>
//                 </defs>
//                 <line key={'verticalAxis'} x1={x1VerticalBar} x2={x2VerticalBar}
//                     y1={y1VerticalBar} y2={y2VerticalBar} stroke="#7d7f7c" strokeWidth={w} markerEnd="url(#head-vertical)">
//                     {/* <title>+ Sigma: {sigma}</title> */}
//                 </line>
//             </>
//         );
//         return bars;
//     }

//     private makeXAxisLabel(value: number, x: number, y: number) {
//         // TODO: make constants be dependant on precision and font size
//         // TODO: font size
//         // TODO: remove y dep
//         return <text x={x - 20} y={y + 25}>{parseFloat(value.toFixed(5))} </text>;
//     }

//     private makeYAxisLabel(value: number, x: number, y: number) {
//         // TODO: make constants be dependant on precision and font size
//         // TODO: font size
//         // TODO: remove x dep
//         return <text x={x} y={y + 25} fontSize={25}>{parseFloat(value.toFixed(1))} </text>;
//     }

//     private renderDescriptiveStatisticsBars() {
//         if (!this.props.volume) return null;
//         const offset = this.padding / 2;
//         const mean = this.descriptiveStatistics.mean;
//         const min = this.descriptiveStatistics.min;
//         const max = this.descriptiveStatistics.max;
//         const sigma = this.descriptiveStatistics.sigma;
//         const extent = max - min;
//         const x = this.width * ((mean + Math.abs(min)) / extent) + offset;
//         const w = offset / 10;
//         const bars = [];
//         const y1 = this.height + offset;// + offset * 2;
//         const y2 = 0;// offset;
//         const xPositive = this.width * ((mean + sigma + Math.abs(min)) / extent) + offset;
//         const xNegative = this.width * ((mean - sigma + Math.abs(min)) / extent) + offset;
//         bars.push(
//             // raise points and lines
//             //
//             <>
//                 <line key={'meanBar'} x1={x} x2={x} y1={y1} y2={y2} stroke="#D3D3D3" strokeDasharray="5, 5" strokeWidth={w}>
//                     <title>Mean: {mean}</title>
//                 </line>
//                 {/* TODO: center */}
//                 {this.makeXAxisLabel(mean, x, y1)}
//                 {/* <text x={x - 40} y={y1 + 25}>{parseFloat(mean.toFixed(7))} </text> */}
//             </>);
//         bars.push(<><line key={'positiveSigmaBar'} x1={xPositive} x2={xPositive} y1={y1} y2={y2} stroke="#D3D3D3" strokeWidth={w}>
//             <title>+Sigma: {mean + sigma}</title>
//         </line>
//         {this.makeXAxisLabel(mean + sigma, xPositive, y1)}
//         </>);
//         bars.push(
//             <><line key={'negativeSigmaBar'} x1={xNegative} x2={xNegative} y1={y1} y2={y2} stroke="#D3D3D3" strokeWidth={w}>
//                 <title>-Sigma: {-sigma}</title>
//             </line>
//             {this.makeXAxisLabel(mean - sigma, xNegative, y1)}
//             </>);
//         return bars;
//     }

//     // fix gridlines and labels
//     // first fix points, allow them in area above and recalc the message or something
//     // remove ghost points from the list
//     private renderBaseline() {
//         // TODO: fix this function
//         const baseline = this.baseline;
//         return <>
//             <BaseLine height={this.height - baseline} offset={this.padding / 2} width={this.width} />
//             {/* TODO: tune height need to be higher */}
//             {/* or + offset */}
//             {this.makeYAxisLabel(0, baseline / 2, this.height)}
//             {/*  */}
//             {this.makeYAxisLabel(1, baseline / 2, baseline)}
//             {/* TODO: 0.2 0.4 etc. */}
//         </>;
//     }

//     // turn off handleclick too
//     private renderPoints() {
//         const points: any[] = [];
//         let point: Vec2;
//         for (let i = 0; i < this.state.points.length; i++) {
//             if (i !== 0 && i !== this.state.points.length - 1) {
//                 // render points such that they are at baseline
//                 const { data, color, id, index } = this.state.points[i];
//                 const finalColor = color;
//                 // TODO: here there is no adjusted alpha for some reason
//                 // here
//                 point = this.controlPointDataToSVGCoords(data);
//                 points.push(<PointComponent
//                     index={index}
//                     key={id}
//                     id={id}
//                     // so here you should provide correct positions
//                     // actual ones I mean
//                     // from the top left corner
//                     // of the SVG!
//                     // x={50}
//                     // y={50}
//                     // so we push here coords
//                     x={point[0]}
//                     y={point[1]}
//                     nX={data.x}
//                     nY={data.alpha}
//                     selected={false}
//                     // should be able to delete point
//                     delete={this.deletePoint}
//                     // onmouseover we provide props.onHover of line graph component
//                     onmouseover={this.props.onHover}
//                     onmousedown={this.handleMouseDown(this.state.points[i])}
//                     onclick={this.handleClick(this.state.points[i])}
//                     color={finalColor}
//                 />);
//             }
//         }
//         return points;
//     }

//     private renderLines() {
//         const points: Vec2[] = [];
//         const lines = [];
//         let maxX: number;
//         let maxY: number;
//         let normalizedX: number;
//         let normalizedY: number;
//         let reverseY: number;

//         const o = this.padding / 2;
//         for (const point of this.state.points) {
//             maxX = this.width + o;
//             // render lines is ~ fine
//             // so introduce here baseline and roof
//             // not padding I guess, or maybe padding + baseline and or roof
//             maxY = this.height + this.padding;
//             // normalize point data using a function
//             normalizedX = (point.data.x * (maxX - o)) + o;
//             // if (point.data.adjustedAlpha === void 0) throw Error('Adjusted alpha is not provided');
//             // this is [0;1] thing
//             const alpha = point.data.alpha;
//             const space = this.height - this.baseline - this.roof;
//             const pointHeightRealAboveBaseline = alpha * space;
//             const pointHeightRealBelowRoof = space - pointHeightRealAboveBaseline;
//             // TODO: may need to add padding
//             const pointHeightBelowSVGOrigin = pointHeightRealBelowRoof + this.roof;
//             // this gonna be fraction [0;1];
//             // const fraction
//             // const a = 0;
//             // normalizedY = point.data.alpha;
//             // normalizedY = ratio
//             // normalizedY = (ratio * (maxY - o)) + o;
//             // normalizedY = (point.data.alpha * (maxY - o)) + o;
//             // normalizedY = (point.data.alpha * (maxY - o)) + o;
//             // reverseY = this.height + this.padding - normalizedY;
//             points.push(Vec2.create(normalizedX, pointHeightBelowSVGOrigin));
//         }

//         const data = points;
//         const size = data.length;

//         for (let i = 0; i < size - 1; i++) {
//             const x1 = data[i][0];
//             const y1 = data[i][1];
//             const x2 = data[i + 1][0];
//             const y2 = data[i + 1][1];
//             lines.push(<line key={`lineOf${i}`} x1={x1} x2={x2} y1={y1} y2={y2} stroke="#cec9ba" strokeWidth="5"/>);
//         }

//         return lines;
//     }
// }