/**
 * Copyright (c) 2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */


import { createRoot } from 'react-dom/client';
import { Viewer } from '../../apps/viewer/app';
import { MAPairwiseScorePlot } from '../../extensions/model-archive/quality-assessment/pairwise/ui';
import { QualityAssessment } from '../../extensions/model-archive/quality-assessment/prop';
import { Model, ResidueIndex } from '../../mol-model/structure';
import './index.html';
require('mol-plugin-ui/skin/light.scss');

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
            const req = await fetch(`https://alphafold.ebi.ac.uk/files/AF-${id}-F1-predicted_aligned_error_v4.json`);
            const json = await req.json();

            const model = this.viewer.plugin.managers.structure.hierarchy.current.models[0]?.cell.obj?.data!;
            const metric = pairwiseMetricFromAlphaFoldDbJson(model, json)!;

            createRoot(document.getElementById(this.plotContainerId)!).render(
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