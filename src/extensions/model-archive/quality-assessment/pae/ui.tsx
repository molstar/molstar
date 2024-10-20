/**
 * Copyright (c) 2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Fragment, memo, useEffect, useRef } from 'react';
import { CollapsableControls, CollapsableState } from '../../../../mol-plugin-ui/base';
import { Button } from '../../../../mol-plugin-ui/controls/common';
import { PencilRulerSvg } from '../../../../mol-plugin-ui/controls/icons';
import { PluginContext } from '../../../../mol-plugin/context';
import { drawPAEPng, PAEDrawing } from './plot';
import { Model, ResidueIndex } from '../../../../mol-model/structure';
import { round } from '../../../../mol-util';
import { BehaviorSubject } from 'rxjs';
import { useBehavior } from '../../../../mol-plugin-ui/hooks/use-behavior';
import { clamp } from '../../../../mol-math/interpolate';
import { AtomicHierarchy } from '../../../../mol-model/structure/model/properties/atomic';

interface State {
    model?: Model;
}

export class PAEPlotUI extends CollapsableControls<{}, State> {
    protected defaultState(): State & CollapsableState {
        return {
            header: 'PAE',
            isCollapsed: false,
            brand: { accent: 'green', svg: PencilRulerSvg }
        };
    }

    private draw = () => {
        const model = this.plugin.managers.structure.hierarchy.current.models[0];
        this.setState({ model: model.cell?.obj?.data });
    };

    protected renderControls(): JSX.Element | null {
        return <>
            <Button onClick={this.draw}>Draw</Button>
            {this.state.model && <Plot plugin={this.plugin} model={this.state.model} />}
        </>;
    }
}

const PlotSize = 1000;
const PlotOffset = 120;

interface PlotInteractivityState {
    crosshairOffset?: [number, number];
    inside?: boolean;
    mouseDown?: boolean;
    boxStart?: [number, number];
    boxEnd?: [number, number];
}

const Plot = memo(({ plugin, model }: { plugin: PluginContext, model: Model }) => {
    const interactivityRect = useRef<SVGRectElement>();
    const _interactivity = useRef<BehaviorSubject<PlotInteractivityState>>();
    if (!_interactivity.current) _interactivity.current = new BehaviorSubject({});
    const interactity = _interactivity.current;

    useEffect(() => {
        const moveEvent = (ev: MouseEvent) => {
            const offset = getPlotMouseOffsetBase(interactivityRect.current!, ev.clientX, ev.clientY);
            interactity.next({ ...interactity.value, crosshairOffset: offset });
        };
        const mouseUpEvent = (ev: MouseEvent) => {
            if (!interactity.value.mouseDown) return;
            const offset = getPlotMouseOffsetBase(interactivityRect.current!, ev.clientX, ev.clientY);
            interactity.next({ ...interactity.value, mouseDown: false, boxEnd: offset });
        };
        window.addEventListener('mousemove', moveEvent);
        window.addEventListener('mouseup', mouseUpEvent);
        return () => {
            window.removeEventListener('mousemove', moveEvent);
            window.removeEventListener('mouseup', mouseUpEvent);
        };
    }, [model]);

    const drawing = drawPAEPng(model);
    if (!drawing) return <>Not available</>;


    const { info, colorRange, chains, png } = drawing;
    const nResidues = info.maxResidueIndex - info.minResidueIndex;

    const border = '#333';
    const line = '#000';

    const legendHeight = 80;
    const legendOffsetY = PlotOffset + PlotSize + 50;

    const viewBox = '0 0 1140 1270';

    return <div style={{ margin: '8px 8px 0 8px', position: 'relative' }}>
        <svg viewBox={viewBox} width='100%'>
            <image x={PlotOffset} y={PlotOffset} width={PlotSize} height={PlotSize} href={png} />
            <line x1={PlotOffset} x2={PlotOffset + PlotSize} y1={PlotOffset} y2={PlotOffset + PlotSize} style={{ stroke: line, strokeDasharray: '15,15' }} />
            <linearGradient id='legend-gradient' x1={0} x2={1} y1={0} y2={0}>
                <stop offset='0%' stopColor={colorRange[0]} />
                <stop offset='100%' stopColor={colorRange[1]} />
            </linearGradient>
            <rect x={PlotOffset} y={legendOffsetY} width={PlotSize} height={legendHeight} style={{ fill: 'url(#legend-gradient)', strokeWidth: 1, stroke: border }} />
            <text x={PlotOffset + 20} y={legendOffsetY + legendHeight - 22} style={{ fontSize: '45px', fill: 'white', fontWeight: 'bold' }}>{round(info.minMetric, 2)} Å</text>
            <text x={PlotOffset + PlotSize - 20} y={legendOffsetY + legendHeight - 22} style={{ fontSize: '45px', fill: 'black', fontWeight: 'bold' }} textAnchor='end'>{round(info.maxMetric, 2)} Å</text>

            <text x={PlotOffset + PlotSize / 2} y={50} className='msp-svg-text' style={{ fontSize: '45px', fontWeight: 'bold' }} textAnchor='middle'>Scored Residue</text>
            <text className='msp-svg-text' style={{ fontSize: '50px', fontWeight: 'bold' }} transform={`translate(50, ${PlotOffset + PlotSize / 2}) rotate(270)`} textAnchor='middle'>Aligned Residue</text>

            {chains.map(({ startOffset, endOffset, label }) => {
                const textOffset = PlotOffset + PlotSize * (startOffset + (endOffset - startOffset) / 2) / nResidues;
                const endLineOffset = PlotOffset + PlotSize * endOffset / nResidues;
                const startLineOffset = PlotOffset + PlotSize * startOffset / nResidues;

                const seq_id = model.atomicHierarchy.residues.label_seq_id;
                const startIndex = seq_id.value(info.minResidueIndex + startOffset);
                const endIndex = seq_id.value(info.minResidueIndex + endOffset - 1);

                return <Fragment key={startOffset}>
                    <text x={textOffset} y={PlotOffset - 15} className='msp-svg-text' style={{ fontSize: '40px' }} textAnchor='middle'>{label} {startIndex}-{endIndex}</text>
                    <text className='msp-svg-text' style={{ fontSize: '40px' }} transform={`translate(${PlotOffset - 15}, ${textOffset}) rotate(270)`} textAnchor='middle'>{label} {startIndex}-{endIndex}</text>
                    <line x1={startLineOffset} x2={startLineOffset} y1={PlotOffset - 20} y2={PlotOffset + PlotSize + 20} style={{ stroke: line, strokeDasharray: '15,15' }} />
                    <line x1={endLineOffset} x2={endLineOffset} y1={PlotOffset - 20} y2={PlotOffset + PlotSize + 20} style={{ stroke: line, strokeDasharray: '15,15' }} />
                    <line x1={PlotOffset - 20} x2={PlotOffset + PlotSize + 20} y1={startLineOffset} y2={startLineOffset} style={{ stroke: line, strokeDasharray: '15,15' }} />
                    <line x1={PlotOffset - 20} x2={PlotOffset + PlotSize + 20} y1={endLineOffset} y2={endLineOffset} style={{ stroke: line, strokeDasharray: '15,15' }} />
                </Fragment>;
            })}
        </svg>
        <svg viewBox={viewBox} style={{ position: 'absolute', inset: 0 }}>
            <rect x={PlotOffset} y={PlotOffset} width={PlotSize} height={PlotSize} style={{ fill: 'transparent', strokeWidth: 0, stroke: '#333', cursor: 'crosshair' }}
                ref={interactivityRect as any}
                onMouseMove={(ev) => {
                    interactity.next({ ...interactity.value, inside: true });
                    ev.currentTarget.style.stroke = 'black';
                    ev.currentTarget.style.strokeWidth = '4px';
                }}
                onMouseDown={(ev) => {
                    interactity.next({ ...interactity.value, mouseDown: true, boxStart: getPlotMouseOffset(ev) });
                }}
                onMouseLeave={(ev) => {
                    interactity.next({ ...interactity.value, inside: false, crosshairOffset: undefined });
                    ev.currentTarget.style.stroke = '#333';
                    ev.currentTarget.style.strokeWidth = '1px';
                }} />
            <PlotInteractivity drawing={drawing} interactity={interactity} />
        </svg>
    </div>;
}, (prev, next) => prev.model === next.model);

function PlotInteractivity({ drawing, interactity }: { drawing: PAEDrawing, interactity: BehaviorSubject<PlotInteractivityState> }) {
    const state = useBehavior(interactity);
    const { crosshairOffset, inside } = state;
    const box = getBox(state);
    const label = getCrosshairLabel(drawing, state);

    return <>
        {inside && crosshairOffset && <line x1={crosshairOffset[0] + PlotOffset} x2={crosshairOffset[0] + PlotOffset} y1={PlotOffset} y2={PlotOffset + PlotSize} style={{ pointerEvents: 'none', stroke: 'black', strokeDasharray: '5,5' }} />}
        {inside && crosshairOffset && <line x1={PlotOffset} x2={PlotOffset + PlotSize} y1={crosshairOffset[1] + PlotOffset} y2={crosshairOffset[1] + PlotOffset} style={{ pointerEvents: 'none', stroke: 'black', strokeDasharray: '5,5' }} />}
        {box && <rect x={PlotOffset + box[0]} y={PlotOffset + box[1]} width={box[2]} height={box[3]} style={{ stroke: '#333', fill: 'rgba(0, 0, 0, 0.1)', pointerEvents: 'none' }} />}
        {label && <text x={PlotOffset + crosshairOffset![0] + 20} y={PlotOffset + crosshairOffset![1] - 10} style={{ fontSize: '40px', fill: 'black', fontWeight: 'bold', pointerEvents: 'none', userSelect: 'none' }} textAnchor='start'>{label[0]}</text>}
        {label && <text x={PlotOffset + crosshairOffset![0] - 20} y={PlotOffset + crosshairOffset![1] + 50} style={{ fontSize: '40px', fill: 'black', fontWeight: 'bold', pointerEvents: 'none', userSelect: 'none' }} textAnchor='end'>{label[1]}</text>}
    </>;
}

function getCrosshairLabel(drawing: PAEDrawing, state: PlotInteractivityState) {
    if (!state.crosshairOffset || !state.inside) return;

    const x = clamp(state.crosshairOffset[0], 0, PlotSize);
    const y = clamp(state.crosshairOffset[1], 0, PlotSize);

    return [getPAEResidueLabel(drawing, x), getPAEResidueLabel(drawing, y)];
}

function getPAEResidueLabel(drawing: PAEDrawing, offset: number) {
    let rI = drawing!.info.minResidueIndex + Math.round(offset / PlotSize * (drawing!.info.maxResidueIndex - drawing!.info.minResidueIndex)) as ResidueIndex;
    rI = clamp(rI, drawing!.info.minResidueIndex, drawing!.info.maxResidueIndex) as ResidueIndex;
    const hierarchy = drawing!.model.atomicHierarchy;
    const asym_id = hierarchy.chains.label_asym_id;
    const seq_id = hierarchy.residues.label_seq_id;
    const comp_id = hierarchy.atoms.label_comp_id;

    return `${asym_id.value(AtomicHierarchy.residueChainIndex(hierarchy, rI))} ${seq_id.value(rI)} ${comp_id.value(AtomicHierarchy.residueFirstAtomIndex(hierarchy, rI))}`;
}

function getBox(state: PlotInteractivityState) {
    const start = state.boxStart;
    const end = state.mouseDown ? state.crosshairOffset : state.boxEnd;
    if (!start || !end) return undefined;

    const x = clamp(Math.min(start[0], end[0]), 0, PlotSize);
    const width = clamp(Math.max(start[0], end[0]), 0, PlotSize) - x;
    const y = clamp(Math.min(start[1], end[1]), 0, PlotSize);
    const height = clamp(Math.max(start[1], end[1]), 0, PlotSize) - y;

    if (width < 1 && height < 1) return undefined;

    return [x, y, width, height];
}

function getPlotMouseOffset(ev: React.MouseEvent<SVGRectElement, MouseEvent>) {
    return getPlotMouseOffsetBase(ev.currentTarget, ev.clientX, ev.clientY);
}

function getPlotMouseOffsetBase(target: HTMLElement | SVGRectElement, clientX: number, clientY: number) {
    const rect = target.getBoundingClientRect();
    const offsetX = PlotSize * (clientX - rect.left) / rect.width;
    const offsetY = PlotSize * (clientY - rect.top) / rect.height;
    return [offsetX, offsetY] as [number, number];
}