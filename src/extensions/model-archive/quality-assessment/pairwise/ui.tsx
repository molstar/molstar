/**
 * Copyright (c) 2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { CSSProperties, Fragment, memo, useEffect, useRef } from 'react';
import { BehaviorSubject, combineLatest } from 'rxjs';
import { clamp } from '../../../../mol-math/interpolate';
import { Model, ResidueIndex } from '../../../../mol-model/structure';
import { AtomicHierarchy } from '../../../../mol-model/structure/model/properties/atomic';
import { PluginStateObject } from '../../../../mol-plugin-state/objects';
import { CollapsableControls, CollapsableState } from '../../../../mol-plugin-ui/base';
import { ScatterPlotSvg } from '../../../../mol-plugin-ui/controls/icons';
import { ParameterControls } from '../../../../mol-plugin-ui/controls/parameters';
import { useBehavior } from '../../../../mol-plugin-ui/hooks/use-behavior';
import { PluginContext } from '../../../../mol-plugin/context';
import { StateTransform } from '../../../../mol-state';
import { round } from '../../../../mol-util';
import { ParamDefinition as PD } from '../../../../mol-util/param-definition';
import { QualityAssessment } from '../prop';
import { drawPairwiseMetricPNG, PAEDrawing } from './plot';

type State = ReturnType<typeof getPropsAndValues>

export class PairwiseMetricPlotUI extends CollapsableControls<{}, State> {
    protected defaultState(): State & CollapsableState {
        return {
            header: 'Pairwise Residue Score',
            isCollapsed: false,
            isHidden: true,
            brand: { accent: 'purple', svg: ScatterPlotSvg },
            params: { } as any,
            values: undefined as any,
            dataSources: [],
        };
    }

    componentDidMount() {
        this.subscribe(combineLatest([
            this.plugin.state.data.events.changed,
            this.plugin.behaviors.state.isAnimating
        ]), ([_, anim]) => {
            if (anim) return;
            const state = getPropsAndValues(this.plugin, this.state.values);
            this.setState({
                ...state,
                isHidden: state.params.data.options.length === 0 || state.params.model.options.length === 0
            });
        });
    }

    protected renderControls(): JSX.Element | null {
        const { params, values, dataSources } = this.state;
        return <>
            <ParameterControls params={params} values={values} onChangeValues={values => this.setState({ values })} />
            <PlotWrapper plugin={this.plugin} values={values} dataSources={dataSources} />
        </>;
    }
}

const PlotWrapper = memo(({ plugin, values, dataSources }: { plugin: PluginContext, values: State['values'], dataSources: State['dataSources'] }) => {
    const model: PluginStateObject.Molecule.Model | undefined = plugin.state.data.cells.get(values.model)?.obj;
    const src = dataSources.find(src => src.id === values.data);
    const cif: PluginStateObject.Format.Cif | undefined = plugin.state.data.cells.get(src?.dataRef!)?.obj;
    const block = cif?.data.blocks[src?.blockIndex!];

    if (!model || !block || !src) return <div className='msp-description'>Data not available</div>;

    const metric = QualityAssessment.pairwiseMetricFromModelArchiveCIF(model.data, block, src.metridId);
    if (!metric) return <div className='msp-description'>Data not available</div>;

    return <Plot plugin={plugin} model={model.data} pairwiseMetric={metric} />;

}, (prev, next) => prev.values.data === next.values.data && prev.values.model === next.values.model);

function getPropsAndValues(plugin: PluginContext, current?: { model?: string, data?: string }) {
    const models = plugin.managers.structure.hierarchy.current.models;
    const cifs = plugin.state.data.selectQ(q => q.root.subtree().ofType(PluginStateObject.Format.Cif));

    const dataSources: {
        id: string,
        label: string,
        metridId: number,
        dataRef: StateTransform.Ref,
        blockIndex: number,
    }[] = [];

    for (const cif of cifs) {
        if (!cif.obj?.data.blocks) continue;
        let blockIndex = 0;
        for (const block of cif.obj.data.blocks) {
            for (const pae of QualityAssessment.findModelArchiveCIFPAEMetrics(block)) {
                dataSources.push({
                    id: `${cif.transform.ref}:${blockIndex}:${pae.id}`,
                    metridId: pae.id,
                    label: `${block.header}: ${pae.name}`,
                    dataRef: cif.transform.ref,
                    blockIndex,
                });
            }
            blockIndex++;
        }
    }

    const params = {
        model: PD.Select(models[0]?.cell.transform.ref, models.map(m => [m.cell.transform.ref, m.cell.obj?.data.label!]), { isHidden: models.length <= 1 }),
        data: PD.Select(dataSources[0]?.id, dataSources.map(o => [o.id, o.label]), { isHidden: dataSources.length <= 1 })
    };

    const values = {
        model: params.model.options.find(o => o[0] === current?.model)?.[0] ?? params.model.options[0]?.[0],
        data: params.data.options.find(o => o[0] === current?.data)?.[0] ?? params.data.options[0]?.[0],
    };

    return { params, values, dataSources };
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

const Plot = memo(({ plugin, model, pairwiseMetric }: { plugin: PluginContext, model: Model, pairwiseMetric: QualityAssessment.Pairwise }) => {
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

    const drawing = drawPairwiseMetricPNG(model, pairwiseMetric);
    if (!drawing) return <>Not available</>;


    const { metric, colorRange, chains, png } = drawing;
    const nResidues = metric.residueRange[1] - metric.residueRange[0];

    const border = '#333';
    const line = '#000';

    const legendHeight = 80;
    const legendOffsetY = PlotOffset + PlotSize + 50;

    const viewBox = '0 0 1140 1270';

    return <div style={{ margin: '8px 8px 0 8px', position: 'relative' }}>
        <svg viewBox={viewBox} width='100%'>
            <image x={PlotOffset + 1} y={PlotOffset + 1} width={PlotSize - 1} height={PlotSize - 1} href={png} />
            <line x1={PlotOffset} x2={PlotOffset + PlotSize} y1={PlotOffset} y2={PlotOffset + PlotSize} style={{ stroke: line, strokeDasharray: '15,15' }} />
            <linearGradient id='legend-gradient' x1={0} x2={1} y1={0} y2={0}>
                <stop offset='0%' stopColor={colorRange[0]} />
                <stop offset='100%' stopColor={colorRange[1]} />
            </linearGradient>
            <rect x={PlotOffset} y={legendOffsetY} width={PlotSize} height={legendHeight} style={{ fill: 'url(#legend-gradient)', strokeWidth: 1, stroke: border }} />
            <text x={PlotOffset + 20} y={legendOffsetY + legendHeight - 22} style={{ fontSize: '45px', fill: 'white', fontWeight: 'bold' }}>{round(metric.valueRange[0], 2)} Å</text>
            <text x={PlotOffset + PlotSize - 20} y={legendOffsetY + legendHeight - 22} style={{ fontSize: '45px', fill: 'black', fontWeight: 'bold' }} textAnchor='end'>{round(metric.valueRange[1], 2)} Å</text>
            <text x={PlotOffset + PlotSize / 2} y={legendOffsetY + legendHeight - 22} style={{ fontSize: '45px', fill: 'black' }} textAnchor='middle'>{metric.name}</text>

            <text x={PlotOffset + PlotSize / 2} y={50} className='msp-svg-text' style={{ fontSize: '45px', fontWeight: 'bold' }} textAnchor='middle'>Scored Residue</text>
            <text className='msp-svg-text' style={{ fontSize: '50px', fontWeight: 'bold' }} transform={`translate(50, ${PlotOffset + PlotSize / 2}) rotate(270)`} textAnchor='middle'>Aligned Residue</text>

            {chains.map(({ startOffset, endOffset, label }) => {
                const textOffset = PlotOffset + PlotSize * (startOffset + (endOffset - startOffset) / 2) / nResidues;
                const endLineOffset = PlotOffset + PlotSize * endOffset / nResidues;
                const startLineOffset = PlotOffset + PlotSize * startOffset / nResidues;

                const seq_id = model.atomicHierarchy.residues.label_seq_id;
                const startIndex = seq_id.value(metric.residueRange[0] + startOffset);
                const endIndex = seq_id.value(metric.residueRange[0] + endOffset - 1);

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
            <rect x={PlotOffset} y={PlotOffset} width={PlotSize} height={PlotSize} style={{ fill: 'transparent', cursor: 'crosshair' }}
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
    const labelStyle: CSSProperties | undefined = label ? { fontSize: '40px', fill: 'black', fontWeight: 'bold', pointerEvents: 'none', userSelect: 'none' } : undefined;

    // TODO: position label depending crosshairOffset position to avoid clipping at edges
    return <>
        {inside && crosshairOffset && <line x1={crosshairOffset[0] + PlotOffset} x2={crosshairOffset[0] + PlotOffset} y1={PlotOffset} y2={PlotOffset + PlotSize} style={{ pointerEvents: 'none', stroke: 'black', strokeDasharray: '5,5' }} />}
        {inside && crosshairOffset && <line x1={PlotOffset} x2={PlotOffset + PlotSize} y1={crosshairOffset[1] + PlotOffset} y2={crosshairOffset[1] + PlotOffset} style={{ pointerEvents: 'none', stroke: 'black', strokeDasharray: '5,5' }} />}
        {box && <rect x={PlotOffset + box[0]} y={PlotOffset + box[1]} width={box[2]} height={box[3]} style={{ stroke: '#333', fill: 'rgba(0, 0, 0, 0.1)', pointerEvents: 'none' }} />}
        {label && <text x={PlotOffset + crosshairOffset![0] + 20} y={PlotOffset + crosshairOffset![1] - 10} style={labelStyle} textAnchor='start'>{label[0]}</text>}
        {label && <text x={PlotOffset + crosshairOffset![0] - 20} y={PlotOffset + crosshairOffset![1] + 50} style={labelStyle} textAnchor='end'>{label[1]}</text>}
        {label?.[2] && <text x={PlotOffset + crosshairOffset![0] - 20} y={PlotOffset + crosshairOffset![1] - 10} style={labelStyle} textAnchor='end'>{label[2]}</text>}
    </>;
}

function getCrosshairLabel(drawing: PAEDrawing, state: PlotInteractivityState) {
    if (!state.crosshairOffset || !state.inside) return;

    const rA = getResidueIndex(drawing, clamp(state.crosshairOffset[0], 0, PlotSize));
    const rB = getResidueIndex(drawing, clamp(state.crosshairOffset[1], 0, PlotSize));

    const value = drawing.metric.values[rA]?.[rB] ?? drawing.metric.values[rB]?.[rA];
    const valueLabel = typeof value === 'number' ? `${round(value, 2)} Å` : '';

    return [getResidueLabel(drawing, rA), getResidueLabel(drawing, rB), valueLabel];
}

function getResidueIndex(drawing: PAEDrawing, offset: number) {
    const rI = drawing.metric.residueRange[0] + Math.round(offset / PlotSize * (drawing.metric.residueRange[1] - drawing.metric.residueRange[0] + 1)) as ResidueIndex;
    return clamp(rI, drawing.metric.residueRange[0], drawing.metric.residueRange[1]) as ResidueIndex;
}

function getResidueLabel(drawing: PAEDrawing, rI: ResidueIndex) {
    const hierarchy = drawing.model.atomicHierarchy;
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