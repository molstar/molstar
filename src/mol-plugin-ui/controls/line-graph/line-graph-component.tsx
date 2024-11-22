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

type ComponentParams<T extends React.Component<any, any, any> | ((props: any) => JSX.Element)> =
    T extends React.Component<infer P, any, any> ? P : T extends (props: infer P) => JSX.Element ? P : never;

export interface PointButtonProps {
    changeXValue: any
    changeAlphaValue: any
    onClick: any
    onPointIndexClick: any
    getValueFromPoint: any
    onColorSquareClick: any
    onEnter: any
    point: ControlPoint
}

function PointButton(props: PointButtonProps) {
    const onAbs = (v: number) => {
        props.changeXValue(props.point.id, v);
    };

    const onAlpha = (v: number) => {
        props.changeAlphaValue(props.point.id, v);
    };

    const handleClick = () => {
        props.onClick(props.point.id);
    };

    const onPointIndexClick = (pointId: UUID) => {
        props.onPointIndexClick(pointId);
    };

    const render = () => {
        const [absValue, relativeValue] = props.getValueFromPoint(props.point.data);
        const truncatedAbsValue = parseFloat((absValue as number).toFixed(3));
        // console.log(absValue, relativeValue);
        const alpha = (props.point.data.alpha as number).toFixed(3);
        const color = props.point.color as Color;
        return (
            <div style={{ display: 'flex', marginBottom: 1 }} key={props.point.id}>
                <Button onClick={() => onPointIndexClick(props.point.id)} style={{ textAlign: 'start', textIndent: '20px' }}>{props.point.index}</Button>
                {<Button style={{ backgroundColor: Color.toStyle(color), minWidth: 32, width: 32 }}
                    onClick={() => { props.onColorSquareClick(props.point.id); } } />}
                <TextInput numeric
                    style={{ minWidth: 0 }} className='msp-form-control' onEnter={props.onEnter} blurOnEnter={true} blurOnEscape={true}
                    value={truncatedAbsValue} placeholder={'Some text'}
                    isDisabled={false} onChange={(value) => { onAbs(value); }} />
                <TextInput numeric
                // ok set value to string repr parseFloat and toFixed, but store in state the true value instead
                // and set the value based on the state inside onAlpha onAbs
                    style={{ minWidth: 0 }} className='msp-form-control' onEnter={props.onEnter} blurOnEnter={true} blurOnEscape={true}
                    value={alpha} placeholder={'Some text'}
                    isDisabled={false} onChange={(value) => { onAlpha(value as any); }} />

                <IconButton title={'Remove point'} svg={DeleteSvg} onClick={handleClick}></IconButton>
            </div>
        );
    };

    return render();
}

function adjustTFParams(name: TFName, ds: VolumeDescriptiveStatistics) {
    switch (name) {
        case 'gaussian':
            const max = ds.max;
            const min = ds.min;
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

interface BaselineProps {
    height: number
    offset: number
    width: number
}

function BaseLine(props: BaselineProps) {
    const render = () => {
        const y = props.height;
        const x1 = props.offset;
        const x2 = props.offset + props.width;
        return (
            <>
                <line y1={y} y2={y} x1={x1} x2={x2} stroke='black' fill='black'></line>
            </>
        );
    };
    return render();
}


interface PointPanelProps {
    points: ControlPoint[]
    removeAllPoints: any
    onPointIndexClick: any
    onColorSquareClick: any
    onPointButtonClick: any
    getValueFromPoint: any
    changeXValue: any
    changeAlphaValue: any
}

function PointsPanel(props: PointPanelProps) {
    const render = () => {
        const points: ControlPoint[] = props.points;
        const realPoints = points.filter(p => p.isTerminal !== true);
        const iconStyle = { position: 'absolute', top: 0, right: 0, lineHeight: '24px', height: '24px', textAlign: 'right', width: '32px', paddingRight: '6px', background: 'none' };
        const removeAllPointsButton = <IconButton title='Remove All Points' style={iconStyle as any} svg={DeleteSvg} small onClick={() => {
            props.removeAllPoints();
        } }></IconButton>;
        const controlPointsButtons = realPoints.map(p => {
            return <PointButton key={p.id} point={p}
                onPointIndexClick={props.onPointIndexClick}
                onColorSquareClick={props.onColorSquareClick}
                onClick={props.onPointButtonClick}
                getValueFromPoint={props.getValueFromPoint}
                changeXValue={props.changeXValue}
                changeAlphaValue={props.changeAlphaValue}
                onEnter={() => {}}
            ></PointButton>;
        });
        return (
            <div style={{ position: 'relative' }}>
                <ExpandGroup header='Control Points Panel' initiallyExpanded={false}>
                    {controlPointsButtons}
                </ExpandGroup>
                {removeAllPointsButton}
            </div>
        );
    };
    return render();
}

interface TFMethodPanelProps {
    params: TFParamsValues
    descriptiveStatistics: VolumeDescriptiveStatistics
    onChange: any
}

function TFMethodPanel(props: TFMethodPanelProps) {
    const render = () => {
        return <ExpandGroup header={props.params.name} initiallyExpanded={true}>
            <TFParamsWrapper onChange={props.onChange} params={props.params} descriptiveStatistics={props.descriptiveStatistics}></TFParamsWrapper>
        </ExpandGroup>;
    };
    return render();
}

export type TFName = 'gaussian' | 'method2' | 'defaults';

export interface TFMethod {
    name: TFName
    params: TFParamsValues
}

interface HelpersPanelProps {
    methodsParams: TFParamsValues[]
    onChange: any
    descriptiveStatistics: VolumeDescriptiveStatistics
}

function HelpersPanel(props: HelpersPanelProps) {
    const render = () => {
        const methodsParams: TFParamsValues[] = props.methodsParams;
        const methodsUI = methodsParams.map(p => {
            return <TFMethodPanel key={p.name} onChange={props.onChange} params={p} descriptiveStatistics={props.descriptiveStatistics}></TFMethodPanel>;
        });
        return (
            <ExpandGroup header='Helpers' initiallyExpanded={false}>
                {methodsUI}
            </ExpandGroup>
        );
    };

    return render();
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

interface TFParamsWrapperProps {
    params: TFParamsValues
    onChange: any
    descriptiveStatistics: VolumeDescriptiveStatistics
}

function TFParamsWrapper(props: TFParamsWrapperProps) {
    const [params, setParams] = useState(props.params);

    const handleChange = (next: TFParamsValues) => {
        props.onChange(next);
    };

    const handleClick = () => {
        props.onChange(props.params);
        setParams(props.params);
    };

    const render = () => {
        const adjustedParams = adjustTFParams((props.params.name as TFName), props.descriptiveStatistics);
        return (<ExpandGroup header='Transfer Function Settings' initiallyExpanded>
            <WaitingParameterControls params={adjustedParams} values={params} onChangeValues={async next => { handleChange(next as any); }} />
            <Button style={{ marginTop: 1 }} onClick={handleClick}>{`Apply ${props.params.name} Transfer Function`}</Button>
        </ExpandGroup>);
    };

    return render();

}

// class TFParamsWrapper extends React.Component<any> {
//     state = { params: this.props.params };

//     handleChange = (next: TFParamsValues) => {
//         this.props.onChange(next);
//     };

//     handleClick = () => {
//         this.props.onChange(this.props.params);
//         this.setState({ params: this.props.params });
//     };

//     render() {
//         const adjustedParams = adjustTFParams(this.props.params.name, this.props.descriptiveStatistics);
//         return (<ExpandGroup header='Transfer Function Settings' initiallyExpanded>
//             <WaitingParameterControls params={adjustedParams} values={this.state.params} onChangeValues={async next => { this.handleChange(next as any); }} />
//             <Button style={{ marginTop: 1 }} onClick={this.handleClick}>{`Apply ${this.props.params.name} Transfer Function`}</Button>
//         </ExpandGroup>);
//     }
// }

export function areControlPointsColorsSame(newControlPoints: ControlPoint[], oldControlPoints: ControlPoint[]) {
    const newSorted = newControlPoints.sort((a, b) => a.index - b.index);
    const oldSorted = oldControlPoints.sort((a, b) => a.index - b.index);
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

function ColorPicker(props: any) {
    const isActive = props.isActive;
    const defaultColor = props.defaultColor;
    const color = props.color;
    const updateColor = props.updateColor;
    const toggleColorPicker = props.toggleColorPicker;
    if (isActive) {
        return <div style={{ marginBottom: '6px', marginTop: 1 }} className='msp-accent-offset'>
            <ControlGroup header='Select Color' initialExpanded={true} hideExpander={true} hideOffset={true} onHeaderClick={toggleColorPicker}
                topRightIcon={CloseSvg} noTopMargin childrenClassName='msp-viewport-controls-panel-controls'>
                <CombinedColorControl param={defaultColor} value={color} onChange={updateColor} name='color' hideNameRow />
            </ControlGroup>
        </div>;
    } else {
        return <Button disabled >Select point to change color</Button>;
    }
}

export const GaussianTFParams = {
    gaussianCenter: PD.Numeric(0.2, { min: 0, max: 1, step: 0.01 }),
    gaussianExtent: PD.Numeric(0.2, { min: 0, max: 1, step: 0.01 }),
    gaussianHeight: PD.Numeric(0.2, { min: 0, max: 1, step: 0.01 }),
    name: PD.Text('gaussian', { isHidden: true })
};

export const DefaultParams = {
    name: PD.Text('defaults', { isHidden: true })
};

export const Method2TFParams = {
    param1: PD.Numeric(0.2, { min: 0, max: 1, step: 0.01 }),
    param2: PD.Numeric(0.2, { min: 0, max: 1, step: 0.01 }),
    name: PD.Text('method2', { isHidden: true })
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
    controlPoints: ControlPoint[]
    volume: Volume | undefined
    // TODO: function types
    onChange: any
    onHover: any
    onDrag: any
    colored: boolean | undefined
    getValueFromPoint: any
    onAbsValueToPointValue: any
    onNumberOfPointsChange: any
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

interface PointRef {
    ref: any
    id: UUID
}

export function LineGraphComponent(props: LineGraphComponentProps) {
    const [controlPoints, setControlPoints] = useState(startEndPoints.concat(props.controlPoints));
    const [copyPoint, setCopyPoint] = useState<Vec2 | undefined >(undefined);
    const [canSelectMultiple, setCanSelectMultiple] = useState(false);
    const [showColorPicker, setShowColorPicker] = useState(false);
    const [clickedPointIds, setClickedPointIds] = useState<UUID[]>([]);

    useEffect(() => {
        console.log(`Control points changed: number ${controlPoints.length}, ${controlPoints}`);
    }, [controlPoints]);

    useEffect(() => {
        gElement.current = document.getElementsByClassName('ghost-points')[0] as SVGElement;
    }, []);


    // https://stackoverflow.com/a/73796936/13136429
    const useRefs = () => {
        const refsByKey = useRef<Record<UUID, HTMLElement | null>>({});

        const setRef = (element: HTMLElement | null, key: UUID) => {
            refsByKey.current[key] = element;
        };

        return { refsByKey: refsByKey.current, setRef };
    };

    const { refsByKey, setRef } = useRefs();

    // const refs = Object.values(refsByKey).filter(Boolean);

    // const pointRefs = useRef<PointRef[]>([]);
    // useEffect(() => {
    //     const refs = controlPoints.map(a => {
    //         const r: PointRef = {
    //             id: a.id,
    //             ref: React.createRef()
    //         };
    //         return r;
    //     });
    //     // new refs each time
    //     pointRefs.current = refs;
    // }, [controlPoints]);

    // useEffect(() => {
    //     console.log(pointRefs.current)
    // }, [pointRefs]);

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
        console.log(`Highlighted point with ID ${pointId}`);
        const targetPoint = controlPoints.find(p => p.id === pointId);
        // throw Error('Not implemented');
        if (!targetPoint) throw Error('Cannot highlight inexisting point exist');
        // get ref
        // const ref = pointRefs.current.find(r => r.id === pointId);
        const ref = refsByKey[pointId];
        console.log('Highlighting ref', ref);
        if (!ref) throw Error(`No ref for point ID ${pointId}`);
        // hover
        const event = new MouseEvent('mouseover', {
            'view': window,
            'bubbles': true,
            'cancelable': true
        });
        // DISPATCHING THE EVENT, i.e., ACTUALLY HOVERING
        ref.dispatchEvent(event);
    }

    function _setMethod2TF(params: Method2ParamsValues) {
        // TODO: implement
    }

    function _setDefaultsTF(params: DefaultParamsValues) {
        const paddingUnnormalized = padding.current / height.current;
        // there two are the same values, reduce to a single one, e.g. could be be just padding, why not
        // should be e.g.
        const generatedPoints = generateControlPoints(ColorNames.black, undefined, props.volume);
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

    function handleChangePoints(points: ControlPoint[]) {
        const pointsSorted = sortPointsByXValues(points);
        const oldNumberOfPoints = controlPoints.filter(p => p.isTerminal !== true).length;
        const newNumberOfPoints = points.filter(p => p.isTerminal !== true).length;
        if (oldNumberOfPoints !== newNumberOfPoints) {
            props.onNumberOfPointsChange(newNumberOfPoints);
        }
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

    function _changePointDataInState(pointId: UUID, data: ControlPointData) {
        const modifiedPoints = controlPoints.map(p => {
            if (p.id === pointId) {

                const modifiedP: ControlPoint = {
                    id: p.id,
                    color: p.color,
                    index: p.index,
                    data: data
                };
                return modifiedP;
            } else {
                return p;
            }
        });
        handleChangePoints(modifiedPoints);
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
        setClickedPointIds([]);
        setShowColorPicker(false);
        handleChangePoints(points);
    }

    function handleLeave() {
        if (selectedPointId === undefined) {
            return;
        }
        debugger;
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
    // const _handleMouseDown = (point: ControlPoint) => {
    //     const { id, isTerminal } = point;
    //     if (isTerminal === true) {
    //         return;
    //     }
    //     if (canSelectMultiple) {
    //         return;
    //     }
    //     const gps = ghostPoints.current;
    //     const copyPoint: Vec2 = controlPointDataToSVGCoords(getPoint(id).data);
    //     gps.push(document.createElementNS(namespace.current, 'circle') as SVGElement);
    //     gps[0].setAttribute('r', '10');
    //     gps[0].setAttribute('fill', 'orange');
    //     gps[0].setAttribute('cx', `${copyPoint[0]}`);
    //     gps[0].setAttribute('cy', `${copyPoint[1]}`);
    //     gps[0].setAttribute('style', 'display: none');
    //     if (!gElement.current) throw Error('No gElement (SVGElement)');
    //     gElement.current.appendChild(gps[0]);
    //     updatedX.current = copyPoint[0];
    //     updatedY.current = copyPoint[1];
    //     selectedPointId.current = point.id;
    // };

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
        for (let i = 0; i < controlPoints.length; i++) {
            if (controlPoints[i].isTerminal !== true) {
                const { data, color, id, index } = controlPoints[i];
                const finalColor = color;
                point = controlPointDataToSVGCoords(data);
                points.push(<PointComponent
                    ref={(element: HTMLElement) => setRef(element, id)}
                    index={index}
                    key={id}
                    id={id}
                    x={point[0]}
                    y={point[1]}
                    nX={data.x}
                    nY={data.alpha}
                    // onmouseover we provide props.onHover of line graph component
                    onmouseover={props.onHover}
                    onmousedown={handleMouseDown(controlPoints[i])}
                    onclick={handleClick(controlPoints[i])}
                    color={finalColor}
                />);
                console.log(`Just pushed ${i}st/th control point`);
                debugger;
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
            {makeYAxisLabel(0, b / 2, h - b)}
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
        console.log('addPoint, points', controlPoints);
        // const points = controlPoints;
        // points.push(point);
        handleChangePoints([...controlPoints, point]);
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
        let maxIndex = 0;
        if (realPoints.length > 0) {
            maxIndex = realPoints[realPoints.length - 1].index;
        }
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
        console.log('removeRightmostPoint, points', controlPoints);
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
        console.log('Handle enter, points', controlPoints);
        document.removeEventListener('mousemove', handleDrag, true);
        document.removeEventListener('mouseup', handlePointUpdate, true);
    }

    const handleKeyDown = (event: any) => {
        console.log('handlekeydown, points', controlPoints);
        // TODO: set canSelectMultiple = true
        if (event.key === 'Shift') {
            setCanSelectMultiple(true);
            // TODO: allow to select another point
        } else {
            // console.log(`${event} is pressed`);
        }
    };

    const handleKeyUp = (event: any) => {
        console.log('handlekeyup, points', controlPoints);
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
        // remove point and add it again?
        const clickedPointsIds = clickedPoints.map(p => p.id);
        const notClickedPoints = currentPoints.filter(p => !clickedPointsIds.includes(p.id));
        handleChangePoints(notClickedPoints);
        const newPoints = clickedPoints.concat(notClickedPoints);
        handleChangePoints(newPoints);
        console.log('Color updated');
        props.onDrag(clickedPoints[0].data);
        // const shakedPoint = _shake();
        // _revertShake(shakedPoint);
        // _handleMouseDown(clickedPoints[0]);
    };

    const onColorSquareClick = (pointId: UUID) => {
        setClickedPointIds([pointId]);
        if (showColorPicker) {
            // TODO: what should happen there?
        } else {
            toggleColorPicker();
        }
    };
    // function _shake() {
    //     const points = controlPoints;
    //     const realPoints = points.filter(p => p.isTerminal !== true);
    //     const randomPoint = realPoints[Math.floor(Math.random() * realPoints.length)];
    //     randomPoint.data.x = randomPoint.data.x + 0.01;
    //     const otherPoints = points.filter(p => p.id !== randomPoint.id);
    //     const newPoints = [...otherPoints, randomPoint];
    //     setControlPoints(otherPoints);
    //     handleChangePoints(otherPoints);
    //     // setControlPoints(newPoints);
    //     // handleChangePoints(newPoints);
    //     console.log('Shaked point', randomPoint, randomPoint.data.x);
    //     return randomPoint;
    // }

    // function _revertShake(randomPoint: ControlPoint) {
    //     const points = controlPoints;
    //     // const realPoints = points.filter(p => p.isTerminal !== true);
    //     // const randomPoint = realPoints[Math.floor(Math.random() * realPoints.length)];
    //     // randomPoint.data.x = randomPoint.data.x - 0.01;
    //     const otherPoints = points.filter(p => p.id !== randomPoint.id);
    //     randomPoint.data.x = randomPoint.data.x - 0.01;
    //     const newPoints = [...otherPoints, randomPoint];
    //     console.log('Revert shaked point', randomPoint, randomPoint.data.x);
    //     setControlPoints(newPoints);
    //     handleChangePoints(newPoints);
    // }

    const LineGraphRendered = () => {
        // console.log(pointRefs);
        console.log(`Rendering ${controlPoints.length} points`, controlPoints);
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
        const plusIconButtonStyle = {
            position: ('absolute' as any),
            top: '100px',
            right: '30px',
            // backgroundColor: 'transparent'
            // border: '2px solid black',
        };
        const minusIconButtonStyle = {
            position: ('absolute' as any),
            top: '100px',
            right: '5px',
            // backgroundColor: 'transparent'
            // border: '2px solid black',
        };

        const { h, b, r, w, p } = _getLineGraphAttributes();

        const rendered = ([
            <div key="LineGraph">
                <IconButton style={plusIconButtonStyle} small svg={PlusBoxSvg} onClick={function (e: React.MouseEvent<HTMLButtonElement>): void {
                    createPoint();
                } } ></IconButton>
                <IconButton style={minusIconButtonStyle} small svg={MinusBoxSvg} onClick={function (e: React.MouseEvent<HTMLButtonElement>): void {
                    removeRightmostPoint();
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
                <ColorPicker isActive={showColorPicker} defaultColor={defaultColor}
                    color={color} updateColor={updateColor} toggleColorPicker={toggleColorPicker}/>
                <>
                    <PointsPanel points={controlPoints}
                        onPointIndexClick={highlightPoint}
                        removeAllPoints={deleteAllPoints}
                        getValueFromPoint={props.getValueFromPoint}
                        changeXValue={changeXValue}
                        changeAlphaValue={changeAlphaValue}
                        onPointButtonClick={removePoint}
                        onColorSquareClick={onColorSquareClick}></PointsPanel>
                    <HelpersPanel onChange={setTF}
                        methodsParams={DefaultTFParams}
                        descriptiveStatistics={descriptiveStatistics.current}></HelpersPanel>
                </>
            </div>,
            <div key="modal" id="modal-root" />
        ]);

        return rendered;
    };

    return LineGraphRendered();
}