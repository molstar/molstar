/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Paul Luna <paulluna0215@gmail.com>
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */
import { PointComponent } from './point-component';

import * as React from 'react';
import { Vec2 } from '../../../mol-math/linear-algebra';
import { Grid } from '../../../mol-model/volume';
import { downsampleHistogram } from '../../../mol-math/histogram';
import { arrayMax } from '../../../mol-util/array';

interface LineGraphYAxisOptions {
    scale?: 'linear' | 'log',
    logMin?: number
}

const DefaultLogMin = 1e-3;

interface LineGraphComponentState {
    points: Vec2[],
    copyPoint: any,
    canSelectMultiple: boolean,
}

export class LineGraphComponent extends React.Component<any, LineGraphComponentState> {
    private myRef: any;
    private height: number;
    private width: number;
    private padding: number;
    private updatedX: number;
    private updatedY: number;
    private selected?: number[];
    private ghostPoints: SVGElement[];
    private gElement: SVGElement;
    private namespace: string;

    constructor(props: any) {
        super(props);
        this.myRef = React.createRef();
        this.state = {
            points: [
                Vec2.create(0, 0),
                Vec2.create(1, 0)
            ],
            copyPoint: undefined,
            canSelectMultiple: false,
        };
        this.height = 400;
        this.width = 600;
        this.padding = 70;
        this.selected = undefined;
        this.ghostPoints = [];
        this.namespace = 'http://www.w3.org/2000/svg';

        for (const point of this.props.data) {
            this.state.points.push(point);
        }

        this.state.points.sort((a, b) => {
            if (a[0] === b[0]) {
                if (a[0] === 0) {
                    return a[1] - b[1];
                }
                if (a[1] === 1) {
                    return b[1] - a[1];
                }
                return a[1] - b[1];
            }
            return a[0] - b[0];
        });

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

    public render() {
        const points = this.renderPoints();
        const lines = this.renderLines();
        const histogram = this.renderHistogram();
        const axes = this.renderAxes();

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

                    <g stroke="black" fill="black">
                        {histogram}
                        {axes}
                        {lines}
                        {points}
                    </g>
                    <g className="ghost-points" stroke="black" fill="black">
                    </g>
                </svg>
            </div>,
            <div key="modal" id="modal-root" />
        ]);
    }

    componentDidMount() {
        this.gElement = document.getElementsByClassName('ghost-points')[0] as SVGElement;
    }

    componentDidUpdate(prevProps: any) {
        // Re-sync internal points when the parent provides a new `data` array
        // (e.g. when the user picks a preset). The constructor seeds points
        // only once, so changes from outside need to be merged in here.
        if (prevProps.data !== this.props.data) {
            const points: Vec2[] = [Vec2.create(0, 0), Vec2.create(1, 0)];
            for (const p of this.props.data) points.push(p);
            points.sort((a, b) => {
                if (a[0] === b[0]) {
                    if (a[0] === 0) return a[1] - b[1];
                    if (a[1] === 1) return b[1] - a[1];
                    return a[1] - b[1];
                }
                return a[0] - b[0];
            });
            this.selected = undefined;
            this.setState({ points });
        }
    }

    private change(points: Vec2[]) {
        const copyPoints = points.slice();
        copyPoints.shift();
        copyPoints.pop();
        this.props.onChange(copyPoints);
    }

    private handleKeyDown = (event: any) => {
        // TODO: set canSelectMultiple = true
    };

    private handleKeyUp = (event: any) => {
        // TODO: SET canSelectMultiple = fasle
    };

    private handleClick = (id: number) => (event: any) => {
        // TODO: add point to selected array
    };

    private handleMouseDown = (id: number) => (event: any) => {
        if (id === 0 || id === this.state.points.length - 1) {
            return;
        }

        if (this.state.canSelectMultiple) {
            return;
        }

        const copyPoint: Vec2 = this.normalizePoint(Vec2.create(this.state.points[id][0], this.state.points[id][1]));
        this.ghostPoints.push(document.createElementNS(this.namespace, 'circle') as SVGElement);
        this.ghostPoints[0].setAttribute('r', '10');
        this.ghostPoints[0].setAttribute('fill', 'orange');
        this.ghostPoints[0].setAttribute('cx', `${copyPoint[0]}`);
        this.ghostPoints[0].setAttribute('cy', `${copyPoint[1]}`);
        this.ghostPoints[0].setAttribute('style', 'display: none');
        this.gElement.appendChild(this.ghostPoints[0]);
        this.updatedX = copyPoint[0];
        this.updatedY = copyPoint[1];
        this.selected = [id];
    };

    private handleDrag(event: any) {
        if (this.selected === undefined) {
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
    }

    private handlePointUpdate(event: any) {
        const selected = this.selected;
        if (this.state.canSelectMultiple) {
            return;
        }

        if (selected === undefined || selected[0] === 0 || selected[0] === this.state.points.length - 1) {
            this.setState({
                copyPoint: undefined,
            });
            return;
        }
        this.selected = undefined;

        const updatedPoint = this.unNormalizePoint(Vec2.create(this.updatedX, this.updatedY));
        const points = this.state.points.filter((_, i) => i !== selected[0]);
        points.push(updatedPoint);
        points.sort((a, b) => {
            if (a[0] === b[0]) {
                if (a[0] === 0) {
                    return a[1] - b[1];
                }
                if (a[1] === 1) {
                    return b[1] - a[1];
                }
                return a[1] - b[1];
            }
            return a[0] - b[0];
        });
        this.setState({
            points,
        });
        this.change(points);
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
        const newPoint = this.unNormalizePoint(Vec2.create(svgP.x, svgP.y));
        points.push(newPoint);
        points.sort((a, b) => {
            if (a[0] === b[0]) {
                if (a[0] === 0) {
                    return a[1] - b[1];
                }
                if (a[1] === 1) {
                    return b[1] - a[1];
                }
                return a[1] - b[1];
            }
            return a[0] - b[0];
        });
        this.setState({ points });
        this.change(points);
    }

    private deletePoint = (i: number) => (event: any) => {
        if (i === 0 || i === this.state.points.length - 1) { return; }
        const points = this.state.points.filter((_, j) => j !== i);
        points.sort((a, b) => {
            if (a[0] === b[0]) {
                if (a[0] === 0) {
                    return a[1] - b[1];
                }
                if (a[1] === 1) {
                    return b[1] - a[1];
                }
                return a[1] - b[1];
            }
            return a[0] - b[0];
        });
        this.setState({ points });
        this.change(points);
        event.stopPropagation();
    };

    private handleLeave() {
        if (this.selected === undefined) {
            return;
        }

        document.addEventListener('mousemove', this.handleDrag, true);
        document.addEventListener('mouseup', this.handlePointUpdate, true);
    }

    private handleEnter() {
        document.removeEventListener('mousemove', this.handleDrag, true);
        document.removeEventListener('mouseup', this.handlePointUpdate, true);
    }

    private get yAxisOptions(): LineGraphYAxisOptions {
        return (this.props.yAxis as LineGraphYAxisOptions | undefined) ?? {};
    }

    private get isLogY() {
        return this.yAxisOptions.scale === 'log';
    }

    private get logMin() {
        return this.yAxisOptions.logMin ?? DefaultLogMin;
    }

    /** Map a stored y in [0,1] to a display fraction in [0,1]. */
    private yToDisplay(y: number) {
        if (!this.isLogY) return y;
        const lo = this.logMin;
        if (y <= lo) return 0;
        if (y >= 1) return 1;
        const logLo = Math.log10(lo);
        return (Math.log10(y) - logLo) / (0 - logLo); // log10(1) = 0
    }

    /** Inverse of yToDisplay: display fraction in [0,1] to stored y in [0,1]. */
    private yFromDisplay(d: number) {
        if (!this.isLogY) return d;
        const lo = this.logMin;
        if (d <= 0) return 0;
        if (d >= 1) return 1;
        const logLo = Math.log10(lo);
        return Math.pow(10, logLo + d * (0 - logLo));
    }

    private normalizePoint(point: Vec2) {
        const offset = this.padding / 2;
        const maxX = this.width + offset;
        const maxY = this.height + offset;
        const normalizedX = (point[0] * (maxX - offset)) + offset;
        const normalizedY = (this.yToDisplay(point[1]) * (maxY - offset)) + offset;
        const reverseY = (this.height + this.padding) - normalizedY;
        const newPoint = Vec2.create(normalizedX, reverseY);
        return newPoint;
    }

    private unNormalizePoint(point: Vec2) {
        const min = this.padding / 2;
        const maxX = this.width + min;
        const maxY = this.height + min;
        const unNormalizedX = (point[0] - min) / (maxX - min);

        // we have to take into account that we reversed y when we first normalized it.
        const displayY = ((this.height + this.padding) - point[1] - min) / (maxY - min);
        const unNormalizedY = this.yFromDisplay(displayY);

        return Vec2.create(unNormalizedX, unNormalizedY);
    }

    private refCallBack(element: any) {
        if (element) {
            this.myRef = element;
        }
    }

    private renderHistogram() {
        if (!this.props.volume) return null;

        // Reuse the histogram computed for robust stats (same bin count) and
        // down-bin it for display, avoiding a second full scan of the volume.
        const fine = Grid.getHistogram(this.props.volume.grid, Grid.RobustStatsBinCount);
        const histogram = downsampleHistogram(fine, 40);
        const bars = [];
        const N = histogram.counts.length;
        const w = this.width / N;
        const offset = this.padding / 2;
        const max = arrayMax(histogram.counts) || 1;
        const isLog = this.isLogY;
        // For log scaling, compress counts in [1, max] to [0, 1]; counts of 0 map to 0.
        const logMax = isLog ? Math.log10(Math.max(2, max)) : 0;
        for (let i = 0; i < N; i++) {
            const x = this.width * (i + 0.5) / N + offset;
            const y1 = this.height + offset;
            const c = histogram.counts[i];
            let frac: number;
            if (isLog) {
                frac = c > 0 ? Math.log10(c) / logMax : 0;
                if (frac < 0) frac = 0;
                else if (frac > 1) frac = 1;
            } else {
                frac = c / max;
            }
            const y2 = this.height * (1 - frac) + offset;
            bars.push(<line key={`histogram${i}`} x1={x} x2={x} y1={y1} y2={y2} stroke="#ded9ca" strokeWidth={w} />);
        }
        return bars;
    }

    private renderAxes() {
        const offset = 0; // this.padding / 2;
        const elements: React.ReactNode[] = [];
        const ticks: { y: number, label: string }[] = this.isLogY
            ? [
                { y: this.logMin, label: this.logMin.toString() },
                { y: 0.01, label: '0.01' },
                { y: 0.1, label: '0.1' },
                { y: 1, label: '1' },
            ]
            : [
                { y: 0, label: '0' },
                { y: 0.5, label: '0.5' },
                { y: 1, label: '1' },
            ];
        for (let i = 0; i < ticks.length; i++) {
            const { y, label } = ticks[i];
            const py = this.normalizePoint(Vec2.create(0, y))[1];
            elements.push(
                <text key={`ticklbl${i}`} x={offset + 4} y={py - 4} textAnchor="start" fontSize="22" fill="#888" stroke="none">{label}</text>
            );
        }
        return elements;
    }

    private renderPoints() {
        const points: any[] = [];
        let point: Vec2;
        for (let i = 0; i < this.state.points.length; i++) {
            if (i !== 0 && i !== this.state.points.length - 1) {
                point = this.normalizePoint(this.state.points[i]);
                points.push(<PointComponent
                    key={i}
                    id={i}
                    x={point[0]}
                    y={point[1]}
                    nX={this.state.points[i][0]}
                    nY={this.state.points[i][1]}
                    selected={false}
                    delete={this.deletePoint}
                    onmouseover={this.props.onHover}
                    onmousedown={this.handleMouseDown(i)}
                    onclick={this.handleClick(i)}
                />);
            }
        }
        return points;
    }

    private renderLines() {
        const data: Vec2[] = [];
        for (const point of this.state.points) {
            data.push(this.normalizePoint(point));
        }

        const lines = [];
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