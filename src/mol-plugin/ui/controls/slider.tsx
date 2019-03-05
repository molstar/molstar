/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react'
import { NumericInput } from './common';
import { noop } from 'mol-util';

export class Slider extends React.Component<{
    min: number,
    max: number,
    value: number,
    step?: number,
    onChange: (v: number) => void,
    disabled?: boolean,
    onEnter?: () => void
}, { isChanging: boolean, current: number }> {

    state = { isChanging: false, current: 0 }

    static getDerivedStateFromProps(props: { value: number }, state: { isChanging: boolean, current: number }) {
        if (state.isChanging || props.value === state.current) return null;
        return { current: props.value };
    }

    begin = () => {
        this.setState({ isChanging: true });
    }

    end = (v: number) => {
        this.setState({ isChanging: false });
        this.props.onChange(v);
    }

    updateCurrent = (current: number) => {
        this.setState({ current });
    }

    updateManually = (v: number) => {
        this.setState({ isChanging: true });

        let n = v;
        if (this.props.step === 1) n = Math.round(n);
        if (n < this.props.min) n = this.props.min;
        if (n > this.props.max) n = this.props.max;

        this.setState({ current: n, isChanging: true });
    }

    onManualBlur = () => {
        this.setState({ isChanging: false });
        this.props.onChange(this.state.current);
    }

    render() {
        let step = this.props.step;
        if (step === void 0) step = 1;
        return <div className='msp-slider'>
            <div>
                <SliderBase min={this.props.min} max={this.props.max} step={step} value={this.state.current} disabled={this.props.disabled}
                    onBeforeChange={this.begin}
                    onChange={this.updateCurrent as any} onAfterChange={this.end as any} />
            </div>
            <div>
                <NumericInput
                    value={this.state.current} blurOnEnter={true} onBlur={this.onManualBlur}
                    isDisabled={this.props.disabled} onChange={this.updateManually} />
            </div>
        </div>;
    }
}

export class Slider2 extends React.Component<{
    min: number,
    max: number,
    value: [number, number],
    step?: number,
    onChange: (v: [number, number]) => void,
    disabled?: boolean,
    onEnter?: () => void
}, { isChanging: boolean, current: [number, number] }> {

    state = { isChanging: false, current: [0, 1] as [number, number] }

    static getDerivedStateFromProps(props: { value: [number, number] }, state: { isChanging: boolean, current: [number, number] }) {
        if (state.isChanging || (props.value[0] === state.current[0]) && (props.value[1] === state.current[1])) return null;
        return { current: props.value };
    }

    begin = () => {
        this.setState({ isChanging: true });
    }

    end = (v: [number, number]) => {
        this.setState({ isChanging: false });
        this.props.onChange(v);
    }

    updateCurrent = (current: [number, number]) => {
        this.setState({ current });
    }

    updateMax = (v: number) => {
        let n = v;
        if (this.props.step === 1) n = Math.round(n);
        if (n < this.state.current[0]) n = this.state.current[0]
        else if (n < this.props.min) n = this.props.min;
        if (n > this.props.max) n = this.props.max;
        this.props.onChange([this.state.current[0], n]);
    }

    updateMin = (v: number) => {
        let n = v;
        if (this.props.step === 1) n = Math.round(n);
        if (n < this.props.min) n = this.props.min;
        if (n > this.state.current[1]) n = this.state.current[1];
        else if (n > this.props.max) n = this.props.max;
        this.props.onChange([n, this.state.current[1]]);
    }

    render() {
        let step = this.props.step;
        if (step === void 0) step = 1;
        return <div className='msp-slider2'>
            <div>
                <NumericInput
                    value={this.state.current[0]} onEnter={this.props.onEnter} blurOnEnter={true}
                    isDisabled={this.props.disabled} onChange={this.updateMin} />
            </div>
            <div>
                <SliderBase min={this.props.min} max={this.props.max} step={step} value={this.state.current} disabled={this.props.disabled}
                    onBeforeChange={this.begin} onChange={this.updateCurrent as any} onAfterChange={this.end as any} range={true} pushable={true} />
            </div>
            <div>
                <NumericInput
                    value={this.state.current[1]} onEnter={this.props.onEnter} blurOnEnter={true}
                    isDisabled={this.props.disabled} onChange={this.updateMax} />
            </div>
        </div>;
    }
}

/**
 * The following code was adapted from react-components/slider library.
 *
 * The MIT License (MIT)
 * Copyright (c) 2015-present Alipay.com, https://www.alipay.com/
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:

 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.

 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

function classNames(_classes: { [name: string]: boolean | number }) {
    let classes = [];
    let hasOwn = {}.hasOwnProperty;

    for (let i = 0; i < arguments.length; i++) {
        let arg = arguments[i];
        if (!arg) continue;

        let argType = typeof arg;

        if (argType === 'string' || argType === 'number') {
            classes.push(arg);
        } else if (Array.isArray(arg)) {
            classes.push(classNames.apply(null, arg));
        } else if (argType === 'object') {
            for (let key in arg) {
                if (hasOwn.call(arg, key) && arg[key]) {
                    classes.push(key);
                }
            }
        }
    }

    return classes.join(' ');
}

function isNotTouchEvent(e: TouchEvent) {
    return e.touches.length > 1 || (e.type.toLowerCase() === 'touchend' && e.touches.length > 0);
}

function getTouchPosition(vertical: boolean, e: TouchEvent) {
    return vertical ? e.touches[0].clientY : e.touches[0].pageX;
}

function getMousePosition(vertical: boolean, e: MouseEvent) {
    return vertical ? e.clientY : e.pageX;
}

function getHandleCenterPosition(vertical: boolean, handle: HTMLElement) {
    const coords = handle.getBoundingClientRect();
    return vertical ?
        coords.top + (coords.height * 0.5) :
        coords.left + (coords.width * 0.5);
}

function pauseEvent(e: Event) {
    e.stopPropagation();
    e.preventDefault();
}

export class Handle extends React.Component<Partial<HandleProps>, {}> {
    render() {
        const {
            className,
            tipFormatter,
            vertical,
            offset,
            value,
            index,
        } = this.props as HandleProps;

        const style = vertical ? { bottom: `${offset}%` } : { left: `${offset}%` };
        return (
            <div className={className} style={style} title={tipFormatter(value, index)}
            />
        );
    }
}

export interface SliderBaseProps {
    min: number,
    max: number,
    step?: number,
    defaultValue?: number | number[],
    value?: number | number[],
    marks?: any,
    included?: boolean,
    className?: string,
    prefixCls?: string,
    disabled?: boolean,
    children?: any,
    onBeforeChange?: (value: number | number[]) => void,
    onChange?: (value: number | number[]) => void,
    onAfterChange?: (value: number | number[]) => void,
    handle?: JSX.Element,
    tipFormatter?: (value: number, index: number) => any,
    dots?: boolean,
    range?: boolean | number,
    vertical?: boolean,
    allowCross?: boolean,
    pushable?: boolean | number,
}

export interface SliderBaseState {
    handle: number | null,
    recent: number,
    bounds: number[]
}

export class SliderBase extends React.Component<SliderBaseProps, SliderBaseState> {
    private sliderElement: HTMLElement | undefined = void 0;
    private handleElements: (HTMLElement | undefined)[] = [];

    constructor(props: SliderBaseProps) {
        super(props);

        const { range, min, max } = props;
        const initialValue = range ? Array.apply(null, Array(+range + 1)).map(() => min) : min;
        const defaultValue = ('defaultValue' in props ? props.defaultValue : initialValue);
        const value = (props.value !== undefined ? props.value : defaultValue);

        const bounds = (range ? value : [min, value]).map((v: number) => this.trimAlignValue(v));

        let recent;
        if (range && bounds[0] === bounds[bounds.length - 1] && bounds[0] === max) {
            recent = 0;
        } else {
            recent = bounds.length - 1;
        }

        this.state = {
            handle: null,
            recent,
            bounds,
        };
    }

    public static defaultProps: SliderBaseProps = {
        prefixCls: 'msp-slider-base',
        className: '',
        min: 0,
        max: 100,
        step: 1,
        marks: {},
        handle: <Handle className='' vertical={false} offset={0} tipFormatter={v => v} value={0} index={0} />,
        onBeforeChange: noop,
        onChange: noop,
        onAfterChange: noop,
        tipFormatter: (value, index) => value,
        included: true,
        disabled: false,
        dots: false,
        range: false,
        vertical: false,
        allowCross: true,
        pushable: false,
    };

    private dragOffset = 0;
    private startPosition = 0;
    private startValue = 0;
    private _getPointsCache: any = void 0;

    componentWillReceiveProps(nextProps: SliderBaseProps) {
        if (!('value' in nextProps || 'min' in nextProps || 'max' in nextProps)) return;

        const { bounds } = this.state;
        if (nextProps.range) {
            const value = nextProps.value || bounds;
            const nextBounds = (value as number[]).map((v: number) => this.trimAlignValue(v, nextProps));
            if (nextBounds.every((v: number, i: number) => v === bounds[i])) return;

            this.setState({ bounds: nextBounds } as SliderBaseState);
            if (bounds.some(v => this.isValueOutOfBounds(v, nextProps))) {
                this.props.onChange!(nextBounds);
            }
        } else {
            const value = nextProps.value !== undefined ? nextProps.value : bounds[1];
            const nextValue = this.trimAlignValue(value as number, nextProps);
            if (nextValue === bounds[1] && bounds[0] === nextProps.min) return;

            this.setState({ bounds: [nextProps.min, nextValue] } as SliderBaseState);
            if (this.isValueOutOfBounds(bounds[1], nextProps)) {
                this.props.onChange!(nextValue);
            }
        }
    }

    onChange(state: this['state']) {
        const props = this.props;
        const isNotControlled = !('value' in props);
        if (isNotControlled) {
            this.setState(state);
        } else if (state.handle !== undefined) {
            this.setState({ handle: state.handle } as SliderBaseState);
        }

        const data = { ...this.state, ...(state as any) };
        const changedValue = props.range ? data.bounds : data.bounds[1];
        props.onChange!(changedValue);
    }

    onMouseDown(e: MouseEvent) {
        if (e.button !== 0) { return; }

        let position = getMousePosition(this.props.vertical!, e);
        if (!this.isEventFromHandle(e)) {
            this.dragOffset = 0;
        } else {
            const handlePosition = getHandleCenterPosition(this.props.vertical!, e.target as HTMLElement);
            this.dragOffset = position - handlePosition;
            position = handlePosition;
        }
        this.onStart(position);
        this.addDocumentEvents('mouse');
        pauseEvent(e);
    }

    onMouseMove(e: MouseEvent) {
        const position = getMousePosition(this.props.vertical!, e);
        this.onMove(e, position - this.dragOffset);
    }

    onMove(e: MouseEvent | TouchEvent, position: number) {
        pauseEvent(e);
        const props = this.props;
        const state = this.state;

        let diffPosition = position - this.startPosition;
        diffPosition = this.props.vertical ? -diffPosition : diffPosition;
        const diffValue = diffPosition / this.getSliderLength() * (props.max - props.min);

        const value = this.trimAlignValue(this.startValue + diffValue);
        const oldValue = state.bounds[state.handle!];
        if (value === oldValue) return;

        const nextBounds = [...state.bounds];
        nextBounds[state.handle!] = value;
        let nextHandle = state.handle!;
        if (props.pushable !== false) {
            const originalValue = state.bounds[nextHandle];
            this.pushSurroundingHandles(nextBounds, nextHandle, originalValue);
        } else if (props.allowCross) {
            nextBounds.sort((a, b) => a - b);
            nextHandle = nextBounds.indexOf(value);
        }
        this.onChange({
            handle: nextHandle,
            bounds: nextBounds,
        } as SliderBaseState);
    }

    onStart(position: number) {
        const props = this.props;
        props.onBeforeChange!(this.getValue());

        const value = this.calcValueByPos(position);
        this.startValue = value;
        this.startPosition = position;

        const state = this.state;
        const { bounds } = state;

        let valueNeedChanging = 1;
        if (this.props.range) {
            let closestBound = 0;
            for (let i = 1; i < bounds.length - 1; ++i) {
                if (value > bounds[i]) { closestBound = i; }
            }
            if (Math.abs(bounds[closestBound + 1] - value) < Math.abs(bounds[closestBound] - value)) {
                closestBound = closestBound + 1;
            }
            valueNeedChanging = closestBound;

            const isAtTheSamePoint = (bounds[closestBound + 1] === bounds[closestBound]);
            if (isAtTheSamePoint) {
                valueNeedChanging = state.recent;
            }

            if (isAtTheSamePoint && (value !== bounds[closestBound + 1])) {
                valueNeedChanging = value < bounds[closestBound + 1] ? closestBound : closestBound + 1;
            }
        }

        this.setState({
            handle: valueNeedChanging,
            recent: valueNeedChanging,
        } as SliderBaseState);

        const oldValue = state.bounds[valueNeedChanging];
        if (value === oldValue) return;

        const nextBounds = [...state.bounds];
        nextBounds[valueNeedChanging] = value;
        this.onChange({ bounds: nextBounds } as SliderBaseState);
    }

    onTouchMove(e: TouchEvent) {
        if (isNotTouchEvent(e)) {
            this.end('touch');
            return;
        }

        const position = getTouchPosition(this.props.vertical!, e);
        this.onMove(e, position - this.dragOffset);
    }

    onTouchStart(e: TouchEvent) {
        if (isNotTouchEvent(e)) return;

        let position = getTouchPosition(this.props.vertical!, e);
        if (!this.isEventFromHandle(e)) {
            this.dragOffset = 0;
        } else {
            const handlePosition = getHandleCenterPosition(this.props.vertical!, e.target as HTMLElement);
            this.dragOffset = position - handlePosition;
            position = handlePosition;
        }
        this.onStart(position);
        this.addDocumentEvents('touch');
        pauseEvent(e);
    }

    /**
     * Returns an array of possible slider points, taking into account both
     * `marks` and `step`. The result is cached.
     */
    getPoints() {
        const { marks, step, min, max } = this.props;
        const cache = this._getPointsCache;
        if (!cache || cache.marks !== marks || cache.step !== step) {
            const pointsObject = { ...marks };
            if (step !== null) {
                for (let point = min; point <= max; point += step!) {
                    pointsObject[point] = point;
                }
            }
            const points = Object.keys(pointsObject).map(parseFloat);
            points.sort((a, b) => a - b);
            this._getPointsCache = { marks, step, points };
        }
        return this._getPointsCache.points;
    }

    getPrecision(step: number) {
        const stepString = step.toString();
        let precision = 0;
        if (stepString.indexOf('.') >= 0) {
            precision = stepString.length - stepString.indexOf('.') - 1;
        }
        return precision;
    }

    getSliderLength() {
        const slider = this.sliderElement;
        if (!slider) {
            return 0;
        }

        return this.props.vertical ? slider.clientHeight : slider.clientWidth;
    }

    getSliderStart() {
        const slider = this.sliderElement as HTMLElement;
        const rect = slider.getBoundingClientRect();

        return this.props.vertical ? rect.top : rect.left;
    }

    getValue(): number {
        const { bounds } = this.state;
        return (this.props.range ? bounds : bounds[1]) as number;
    }

    private eventHandlers = {
        'touchmove': (e: TouchEvent) => this.onTouchMove(e),
        'touchend': (e: TouchEvent) => this.end('touch'),
        'mousemove': (e: MouseEvent) => this.onMouseMove(e),
        'mouseup': (e: MouseEvent) => this.end('mouse'),
    }

    addDocumentEvents(type: 'touch' | 'mouse') {
        if (type === 'touch') {
            document.addEventListener('touchmove', this.eventHandlers.touchmove);
            document.addEventListener('touchend', this.eventHandlers.touchend);
        } else if (type === 'mouse') {
            document.addEventListener('mousemove', this.eventHandlers.mousemove);
            document.addEventListener('mouseup', this.eventHandlers.mouseup);
        }
    }

    calcOffset(value: number) {
        const { min, max } = this.props;
        const ratio = (value - min) / (max - min);
        return ratio * 100;
    }

    calcValue(offset: number) {
        const { vertical, min, max } = this.props;
        const ratio = Math.abs(offset / this.getSliderLength());
        const value = vertical ? (1 - ratio) * (max - min) + min : ratio * (max - min) + min;
        return value;
    }

    calcValueByPos(position: number) {
        const pixelOffset = position - this.getSliderStart();
        const nextValue = this.trimAlignValue(this.calcValue(pixelOffset));
        return nextValue;
    }

    end(type: 'mouse' | 'touch') {
        this.removeEvents(type);
        this.props.onAfterChange!(this.getValue());
        this.setState({ handle: null } as SliderBaseState);
    }

    isEventFromHandle(e: Event) {
        for (const h of this.handleElements) {
            if (h === e.target) return true;
        }
        return false;

        // return this.state.bounds.some((x, i) => e.target

        // (
        //     //this.handleElements[i] && e.target === ReactDOM.findDOMNode(this.handleElements[i])
        // ));
    }

    isValueOutOfBounds(value: number, props: SliderBaseProps) {
        return value < props.min || value > props.max;
    }

    pushHandle(bounds: number[], handle: number, direction: number, amount: number) {
        const originalValue = bounds[handle];
        let currentValue = bounds[handle];
        while (direction * (currentValue - originalValue) < amount) {
            if (!this.pushHandleOnePoint(bounds, handle, direction)) {
                // can't push handle enough to create the needed `amount` gap, so we
                // revert its position to the original value
                bounds[handle] = originalValue;
                return false;
            }
            currentValue = bounds[handle];
        }
        // the handle was pushed enough to create the needed `amount` gap
        return true;
    }

    pushHandleOnePoint(bounds: number[], handle: number, direction: number) {
        const points = this.getPoints();
        const pointIndex = points.indexOf(bounds[handle]);
        const nextPointIndex = pointIndex + direction;
        if (nextPointIndex >= points.length || nextPointIndex < 0) {
            // reached the minimum or maximum available point, can't push anymore
            return false;
        }
        const nextHandle = handle + direction;
        const nextValue = points[nextPointIndex];
        const { pushable: threshold } = this.props;
        const diffToNext = direction * (bounds[nextHandle] - nextValue);
        if (!this.pushHandle(bounds, nextHandle, direction, +threshold! - diffToNext)) {
            // couldn't push next handle, so we won't push this one either
            return false;
        }
        // push the handle
        bounds[handle] = nextValue;
        return true;
    }

    pushSurroundingHandles(bounds: number[], handle: number, originalValue: number) {
        const { pushable: threshold } = this.props;
        const value = bounds[handle];

        let direction = 0;
        if (bounds[handle + 1] - value < threshold!) {
            direction = +1;
        } else if (value - bounds[handle - 1] < threshold!) {
            direction = -1;
        }

        if (direction === 0) { return; }

        const nextHandle = handle + direction;
        const diffToNext = direction * (bounds[nextHandle] - value);
        if (!this.pushHandle(bounds, nextHandle, direction, +threshold! - diffToNext)) {
            // revert to original value if pushing is impossible
            bounds[handle] = originalValue;
        }
    }

    removeEvents(type: 'touch' | 'mouse') {
        if (type === 'touch') {
            document.removeEventListener('touchmove', this.eventHandlers.touchmove);
            document.removeEventListener('touchend', this.eventHandlers.touchend);
        } else if (type === 'mouse') {
            document.removeEventListener('mousemove', this.eventHandlers.mousemove);
            document.removeEventListener('mouseup', this.eventHandlers.mouseup);
        }
    }

    trimAlignValue(v: number, nextProps?: SliderBaseProps) {
        const { handle, bounds } = (this.state || {}) as this['state'];
        const { marks, step, min, max, allowCross } = { ...this.props, ...(nextProps || {}) } as SliderBaseProps;

        let val = v;
        if (val <= min) {
            val = min;
        }
        if (val >= max) {
            val = max;
        }
        /* eslint-disable eqeqeq */
        if (!allowCross && handle != null && handle > 0 && val <= bounds[handle - 1]) {
            val = bounds[handle - 1];
        }
        if (!allowCross && handle != null && handle < bounds.length - 1 && val >= bounds[handle + 1]) {
            val = bounds[handle + 1];
        }
        /* eslint-enable eqeqeq */

        const points = Object.keys(marks).map(parseFloat);
        if (step !== null) {
            const closestStep = (Math.round((val - min) / step!) * step!) + min;
            points.push(closestStep);
        }

        const diffs = points.map((point) => Math.abs(val - point));
        const closestPoint = points[diffs.indexOf(Math.min.apply(Math, diffs))];

        return step !== null ? parseFloat(closestPoint.toFixed(this.getPrecision(step!))) : closestPoint;
    }

    render() {
        const {
            handle,
            bounds,
        } = this.state;
        const {
            className,
            prefixCls,
            disabled,
            vertical,
            dots,
            included,
            range,
            step,
            marks,
            max, min,
            tipFormatter,
            children,
        } = this.props;

        const customHandle = this.props.handle;

        const offsets = bounds.map(v => this.calcOffset(v));

        const handleClassName = `${prefixCls}-handle`;

        const handlesClassNames = bounds.map((v, i) => classNames({
            [handleClassName]: true,
            [`${handleClassName}-${i + 1}`]: true,
            [`${handleClassName}-lower`]: i === 0,
            [`${handleClassName}-upper`]: i === bounds.length - 1,
        }));

        const isNoTip = (step === null) || (tipFormatter === null);

        const commonHandleProps = {
            prefixCls,
            noTip: isNoTip,
            tipFormatter,
            vertical,
        };

        this.handleElements = [];
        const handles = bounds.map((v, i) => React.cloneElement(customHandle!, {
            ...commonHandleProps,
            className: handlesClassNames[i],
            value: v,
            offset: offsets[i],
            dragging: handle === i,
            index: i,
            key: i,
            ref: (h: any) => this.handleElements.push(h)  // `handle-${i}`,
        }));
        if (!range) { handles.shift(); }

        const isIncluded = included || range;

        const tracks: JSX.Element[] = [];
        // for (let i = 1; i < bounds.length; ++i) {
        //     const trackClassName = classNames({
        //         [`${prefixCls}-track`]: true,
        //         [`${prefixCls}-track-${i}`]: true,
        //     });
        //     tracks.push(
        //         <Track className={trackClassName} vertical={vertical} included={isIncluded}
        //             offset={offsets[i - 1]} length={offsets[i] - offsets[i - 1]} key={i}
        //             />
        //     );
        // }

        const sliderClassName = classNames({
            [prefixCls!]: true,
            [`${prefixCls}-with-marks`]: Object.keys(marks).length,
            [`${prefixCls}-disabled`]: disabled!,
            [`${prefixCls}-vertical`]: this.props.vertical!,
            [className!]: !!className,
        });

        return (
            <div ref={e => this.sliderElement = e!} className={sliderClassName}
                onTouchStart={disabled ? noop : this.onTouchStart.bind(this)}
                onMouseDown={disabled ? noop : this.onMouseDown.bind(this)}
            >
                <div className={`${prefixCls}-rail`} />
                {tracks}
                <Steps prefixCls={prefixCls} vertical={vertical} marks={marks} dots={dots} step={step}
                    included={isIncluded} lowerBound={bounds[0]}
                    upperBound={bounds[bounds.length - 1]} max={max} min={min}
                />
                {handles}
                <Marks className={`${prefixCls}-mark`} vertical={vertical!} marks={marks}
                    included={isIncluded!} lowerBound={bounds[0]}
                    upperBound={bounds[bounds.length - 1]} max={max} min={min}
                />
                {children}
            </div>
        );
    }
}

export interface HandleProps {
    className: string,
    vertical: boolean,
    offset: number,
    tipFormatter: (v: number, index: number) => any,
    value: number,
    index: number,
}

interface MarksProps {
    className: string,
    vertical: boolean,
    marks: any,
    included: boolean | number,
    upperBound: number,
    lowerBound: number,
    max: number,
    min: number
}
const Marks = ({ className, vertical, marks, included, upperBound, lowerBound, max, min }: MarksProps) => {
    const marksKeys = Object.keys(marks);
    const marksCount = marksKeys.length;
    const unit = 100 / (marksCount - 1);
    const markWidth = unit * 0.9;

    const range = max - min;
    const elements = marksKeys.map(parseFloat).sort((a, b) => a - b).map((point) => {
        const isActived = (!included && point === upperBound) ||
            (included && point <= upperBound && point >= lowerBound);
        const markClassName = classNames({
            [`${className}-text`]: true,
            [`${className}-text-active`]: isActived,
        });

        const bottomStyle = {
            // height: markWidth + '%',
            marginBottom: '-50%',
            bottom: `${(point - min) / range * 100}%`,
        };

        const leftStyle = {
            width: `${markWidth}%`,
            marginLeft: `${-markWidth / 2}%`,
            left: `${(point - min) / range * 100}%`,
        };

        const style = vertical ? bottomStyle : leftStyle;

        const markPoint = marks[point];
        const markPointIsObject = typeof markPoint === 'object' && !React.isValidElement(markPoint);
        const markLabel = markPointIsObject ? markPoint.label : markPoint;
        const markStyle = markPointIsObject ? { ...style, ...markPoint.style } : style;
        return (<span className={markClassName} style={markStyle} key={point}>
            {markLabel}
        </span>);
    });

    return <div className={className}>{elements}</div>;
};

function calcPoints(vertical: boolean, marks: any, dots: boolean, step: number, min: number, max: number) {
    const points = Object.keys(marks).map(parseFloat);
    if (dots) {
        for (let i = min; i <= max; i = i + step) {
            if (points.indexOf(i) >= 0) continue;
            points.push(i);
        }
    }
    return points;
}

const Steps = ({ prefixCls, vertical, marks, dots, step, included,
    lowerBound, upperBound, max, min }: any) => {
    const range = max - min;
    const elements = calcPoints(vertical, marks, dots, step, min, max).map((point) => {
        const offset = `${Math.abs(point - min) / range * 100}%`;
        const style = vertical ? { bottom: offset } : { left: offset };

        const isActived = (!included && point === upperBound) ||
            (included && point <= upperBound && point >= lowerBound);
        const pointClassName = classNames({
            [`${prefixCls}-dot`]: true,
            [`${prefixCls}-dot-active`]: isActived,
        });

        return <span className={pointClassName} style={style} key={point} />;
    });

    return <div className={`${prefixCls}-step`}>{elements}</div>;
};

    // const Track = ({ className, included, vertical, offset, length }: any) => {
    //     const style: any = {
    //         visibility: included ? 'visible' : 'hidden'
    //     };
    //     if (vertical) {
    //         style.bottom = `${offset}%`;
    //         style.height = `${length}%`;
    //     } else {
    //         style.left = `${offset}%`;
    //         style.width = `${length}%`;
    //     }
    //     return <div className={className} style={style} />;
    // };