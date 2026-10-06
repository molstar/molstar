/**
 * Copyright (c) 2024-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */


import { createRoot } from 'react-dom/client';
import { Viewer } from '../../apps/viewer/app';
import { MAPairwiseScorePlot } from '../../extensions/model-archive/quality-assessment/pairwise/ui';
import { QualityAssessment } from '../../extensions/model-archive/quality-assessment/prop';
import { Model, ResidueIndex } from '../../mol-model/structure';
import './index.html';
import '../../mol-plugin-ui/skin/light.scss';

export class AlphaFoldPAEExample {
    viewer: Viewer;
    plotContainerId: string;


    async init(options: { pluginContainerId: string, plotContainerId: string }) {
        this.plotContainerId = options.plotContainerId;
        this.viewer = await Viewer.create(options.pluginContainerId, {
            layoutIsExpanded: false,
            layoutShowControls: false,
            layoutShowLeftPanel: false,
            layoutShowLog: false,
        });

        return this;
    }

    async load(afId: string) {
        const id = afId.trim().toUpperCase();

        const plotRoot = createRoot(document.getElementById(this.plotContainerId)!);
        plotRoot.render(<div>Loading...</div>);

        await this.viewer.plugin.clear();
        await this.viewer.loadAlphaFoldDb(id);

        try {
            // AlphaFold DB rotates the file version (the current records use
            // v6). Resolve the document URL through the public prediction API
            // instead of keeping a stale version in the example.
            const predictionReq = await fetch(`https://alphafold.ebi.ac.uk/api/prediction/${encodeURIComponent(id)}`);
            if (!predictionReq.ok) throw new Error(`AlphaFold DB prediction lookup failed (${predictionReq.status}).`);
            const predictions = await predictionReq.json();
            const paeUrl = Array.isArray(predictions) ? predictions.find((p: any) => typeof p?.paeDocUrl === 'string')?.paeDocUrl : undefined;
            if (!paeUrl) throw new Error('AlphaFold DB did not provide a predicted-aligned-error document URL.');

            const req = await fetch(paeUrl);
            if (!req.ok) throw new Error(`AlphaFold DB PAE download failed (${req.status}).`);
            const json = await req.json();

            const model = this.viewer.plugin.managers.structure.hierarchy.current.models[0]?.cell.obj?.data!;
            if (!model) throw new Error('AlphaFold DB model was not loaded.');
            const metric = pairwiseMetricFromAlphaFoldDbJson(model, json);
            if (!metric) throw new Error('AlphaFold DB returned an invalid predicted-aligned-error document.');

            plotRoot.render(
                <div className='msp-plugin' style={{ background: 'white' }}>
                    <MAPairwiseScorePlot plugin={this.viewer.plugin} pairwiseMetric={metric} model={model} />
                </div>
            );
        } catch (err) {
            plotRoot.render(<div>Error: {String(err)}</div>);
        }
    }
}

function pairwiseMetricFromAlphaFoldDbJson(model: Model, data: any): QualityAssessment.Pairwise | undefined {
    if (!Array.isArray(data) || !data[0]?.predicted_aligned_error) return undefined;

    const { residues, residueAtomSegments, atomSourceIndex } = model.atomicHierarchy;
    const sortedResidueIndices = new Array(residues._rowCount).fill(0).map((_, i) => i);
    sortedResidueIndices.sort((a, b) => {
        const idxA = atomSourceIndex.value(residueAtomSegments.offsets[a]);
        const idxB = atomSourceIndex.value(residueAtomSegments.offsets[b]);
        return idxA - idxB;
    });

    const metricData = data[0].predicted_aligned_error as number[][];

    const metric: QualityAssessment.Pairwise = {
        id: 0,
        name: 'AlphaFold DB PAE',
        residueRange: [0 as ResidueIndex, (residues._rowCount - 1) as ResidueIndex],
        valueRange: [0, data[0].max_predicted_aligned_error],
        values: {}
    };

    for (let i = 0; i < metricData.length; i++) {
        const rA = sortedResidueIndices[i];
        if (typeof rA !== 'number') continue;
        const row = metricData[i];
        const xs: any = (metric.values[rA as ResidueIndex] = {});
        for (let j = 0; j < row.length; j++) {
            const rB = sortedResidueIndices[j];
            if (typeof rB !== 'number') continue;
            xs[rB] = row[j];
        }
    }

    return metric;
}

(window as any).AlphaFoldPAEExample = new AlphaFoldPAEExample();
