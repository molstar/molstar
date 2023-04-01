/**
 * Copyright (c) 2020-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 */

import { SymmetryOperator } from '../../mol-math/geometry';
import { Mat4 } from '../../mol-math/linear-algebra';
import { SIFTSMapping } from '../../mol-model-props/sequence/sifts-mapping';
import { QueryContext, Structure, StructureElement, StructureProperties, StructureSelection } from '../../mol-model/structure';
import { alignAndSuperpose, superpose } from '../../mol-model/structure/structure/util/superposition';
import { alignAndSuperposeWithSIFTSMapping } from '../../mol-model/structure/structure/util/superposition-sifts-mapping';
import { StructureSelectionQueries } from '../../mol-plugin-state/helpers/structure-selection-query';
import { StructureSelectionHistoryEntry } from '../../mol-plugin-state/manager/structure/selection';
import { PluginStateObject } from '../../mol-plugin-state/objects';
import { StateTransforms } from '../../mol-plugin-state/transforms';
import { PluginCommands } from '../../mol-plugin/commands';
import { PluginConfig } from '../../mol-plugin/config';
import { StateObjectCell, StateObjectRef } from '../../mol-state';
import { elementLabel, structureElementStatsLabel } from '../../mol-theme/label';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { stripTags } from '../../mol-util/string';
import { CollapsableControls, PurePluginUIComponent } from '../base';
import { Button, IconButton, ToggleButton } from '../controls/common';
import { ArrowDownwardSvg, ArrowUpwardSvg, DeleteOutlinedSvg, HelpOutlineSvg, Icon, SuperposeAtomsSvg, SuperposeChainsSvg, SuperpositionSvg, TuneSvg } from '../controls/icons';
import { ParameterControls } from '../controls/parameters';
import { ToggleSelectionModeButton } from './selection';

export class StructureSuperpositionControls extends CollapsableControls {
    defaultState() {
        return {
            isCollapsed: false,
            header: 'Superposition',
            brand: { accent: 'gray' as const, svg: SuperpositionSvg },
            isHidden: true
        };
    }

    componentDidMount() {
        this.subscribe(this.plugin.managers.structure.hierarchy.behaviors.selection, sel => {
            this.setState({ isHidden: sel.structures.length < 2 });
        });
    }

    renderControls() {
        return <>
            <SuperpositionControls />
        </>;
    }
}

export const StructureSuperpositionParams = {
    alignSequences: PD.Boolean(true, { isEssential: true, description: 'For Chain-based 3D superposition, perform a sequence alignment and use the aligned residue pairs to guide the 3D superposition.' }),
    traceOnly: PD.Boolean(true, { description: 'For Chain- and Uniprot-based 3D superposition, base superposition only on CA (and equivalent) atoms.' })
};
const DefaultStructureSuperpositionOptions = PD.getDefaultValues(StructureSuperpositionParams);
export type StructureSuperpositionOptions = PD.ValuesFor<typeof StructureSuperpositionParams>

const SuperpositionTag = 'SuperpositionTransform';

type SuperpositionControlsState = {
    isBusy: boolean,
    action?: 'byChains' | 'byAtoms' | 'options',
    canUseDb?: boolean,
    options: StructureSuperpositionOptions
}

export interface LociEntry {
    loci: StructureElement.Loci,
    label: string,
    cell: StateObjectCell<PluginStateObject.Molecule.Structure>
};

interface AtomsLociEntry extends LociEntry {
    atoms: StructureSelectionHistoryEntry[]
};

export class SuperpositionControls extends PurePluginUIComponent<{ }, SuperpositionControlsState> {
    state: SuperpositionControlsState = {
        isBusy: false,
        canUseDb: false,
        action: undefined,
        options: DefaultStructureSuperpositionOptions
    };

    componentDidMount() {
        this.subscribe(this.selection.events.changed, () => {
            this.forceUpdate();
        });

        this.subscribe(this.selection.events.additionsHistoryUpdated, () => {
            this.forceUpdate();
        });

        this.subscribe(this.plugin.behaviors.state.isBusy, v => {
            this.setState({ isBusy: v });
        });

        this.subscribe(this.plugin.managers.structure.hierarchy.behaviors.selection, sel => {
            this.setState({ canUseDb: sel.structures.every(s => !!s.cell.obj?.data && s.cell.obj.data.models.some(m => SIFTSMapping.Provider.isApplicable(m))) });
        });
    }

    get selection() {
        return this.plugin.managers.structure.selection;
    }

    async transform(s: StateObjectRef<PluginStateObject.Molecule.Structure>, matrix: Mat4, coordinateSystem?: SymmetryOperator) {
        const r = StateObjectRef.resolveAndCheck(this.plugin.state.data, s);
        if (!r) return;
        const o = this.plugin.state.data.selectQ(q => q.byRef(r.transform.ref).subtree().withTransformer(StateTransforms.Model.TransformStructureConformation))[0];

        const transform = coordinateSystem && !Mat4.isIdentity(coordinateSystem.matrix)
            ? Mat4.mul(Mat4(), coordinateSystem.matrix, matrix)
            : matrix;

        const params = {
            transform: {
                name: 'matrix' as const,
                params: { data: transform, transpose: false }
            }
        };
        const b = o
            ? this.plugin.state.data.build().to(o).update(params)
            : this.plugin.state.data.build().to(s)
                .insert(StateTransforms.Model.TransformStructureConformation, params, { tags: SuperpositionTag });
        await this.plugin.runTask(this.plugin.state.data.updateTree(b));
    }

    private getRootStructure(s: Structure) {
        const parent = this.plugin.helpers.substructureParent.get(s)!;
        return this.plugin.state.data.selectQ(q => q.byValue(parent).rootOfType(PluginStateObject.Molecule.Structure))[0].obj?.data!;
    }

    superposeChains = async () => {
        const { query } = this.state.options.traceOnly ? StructureSelectionQueries.trace : StructureSelectionQueries.polymer;
        const entries = this.chainEntries;

        const locis = entries.map(e => {
            const s = StructureElement.Loci.toStructure(e.loci);
            const loci = StructureSelection.toLociWithSourceUnits(query(new QueryContext(s)));
            return StructureElement.Loci.remap(loci, this.getRootStructure(e.loci.structure));
        });

        const pivot = this.plugin.managers.structure.hierarchy.findStructure(locis[0]?.structure);
        const coordinateSystem = pivot?.transform?.cell.obj?.data.coordinateSystem;

        const transforms = this.state.options.alignSequences
            ? alignAndSuperpose(locis)
            : superpose(locis);

        const eA = entries[0];
        for (let i = 1, il = locis.length; i < il; ++i) {
            const eB = entries[i];
            const { bTransform, rmsd } = transforms[i - 1];
            await this.transform(eB.cell, bTransform, coordinateSystem);
            const labelA = stripTags(eA.label);
            const labelB = stripTags(eB.label);
            this.plugin.log.info(`Superposed [${labelA}] and [${labelB}] with RMSD ${rmsd.toFixed(2)}.`);
        }
        await this.cameraReset();
    };

    superposeAtoms = async () => {
        const entries = this.atomEntries;

        const atomLocis = entries.map(e => {
            return StructureElement.Loci.remap(e.loci, this.getRootStructure(e.loci.structure));
        });
        const transforms = superpose(atomLocis);

        const pivot = this.plugin.managers.structure.hierarchy.findStructure(atomLocis[0]?.structure);
        const coordinateSystem = pivot?.transform?.cell.obj?.data.coordinateSystem;

        const eA = entries[0];
        for (let i = 1, il = atomLocis.length; i < il; ++i) {
            const eB = entries[i];
            const { bTransform, rmsd } = transforms[i - 1];
            await this.transform(eB.cell, bTransform, coordinateSystem);
            const labelA = stripTags(eA.label);
            const labelB = stripTags(eB.label);
            const count = entries[i].atoms.length;
            this.plugin.log.info(`Superposed ${count} ${count === 1 ? 'atom' : 'atoms'} of [${labelA}] and [${labelB}] with RMSD ${rmsd.toFixed(2)}.`);
        }
        await this.cameraReset();
    };

    superposeDb = async () => {
        const input = this.plugin.managers.structure.hierarchy.behaviors.selection.value.structures;
        const traceOnly = this.state.options.traceOnly;

        const structures = input.map(s => s.cell.obj?.data!);
        const { entries, failedPairs, zeroOverlapPairs } = alignAndSuperposeWithSIFTSMapping(structures, { traceOnly });

        const coordinateSystem = input[0]?.transform?.cell.obj?.data.coordinateSystem;

        let rmsd = 0;

        for (const xform of entries) {
            await this.transform(input[xform.other].cell, xform.transform.bTransform, coordinateSystem);
            rmsd += xform.transform.rmsd;
        }

        rmsd /= Math.max(entries.length - 1, 1);

        const formatPairs = (pairs: [number, number][]) => {
            return `[${pairs.map(([i, j]) => `(${structures[i].models[0].entryId}, ${structures[j].models[0].entryId})`).join(', ')}]`;
        };

        if (zeroOverlapPairs.length) {
            this.plugin.log.warn(`Superposition: No UNIPROT mapping overlap between structures ${formatPairs(zeroOverlapPairs)}.`);
        }

        if (failedPairs.length) {
            this.plugin.log.error(`Superposition: Failed to superpose structures ${formatPairs(failedPairs)}.`);
        }

        if (entries.length) {
            this.plugin.log.info(`Superposed ${entries.length + 1} structures with avg. RMSD ${rmsd.toFixed(2)} Å.`);
            await this.cameraReset();
        }
    };

    async cameraReset() {
        await new Promise(res => requestAnimationFrame(res));
        PluginCommands.Camera.Reset(this.plugin);
    }

    toggleByChains = () => this.setState({ action: this.state.action === 'byChains' ? void 0 : 'byChains' });
    toggleByAtoms = () => this.setState({ action: this.state.action === 'byAtoms' ? void 0 : 'byAtoms' });
    toggleOptions = () => this.setState({ action: this.state.action === 'options' ? void 0 : 'options' });

    highlight(loci: StructureElement.Loci) {
        this.plugin.managers.interactivity.lociHighlights.highlightOnly({ loci }, false);
    }

    moveHistory(e: StructureSelectionHistoryEntry, direction: 'up' | 'down') {
        this.plugin.managers.structure.selection.modifyHistory(e, direction, void 0, true);
    }

    focusLoci(loci: StructureElement.Loci) {
        this.plugin.managers.camera.focusLoci(loci);
    }

    lociEntry(e: LociEntry, idx: number) {
        return <div className='msp-flex-row' key={idx}>
            <Button noOverflow title='Click to focus. Hover to highlight.' onClick={() => this.focusLoci(e.loci)} style={{ width: 'auto', textAlign: 'left' }} onMouseEnter={() => this.highlight(e.loci)} onMouseLeave={() => this.plugin.managers.interactivity.lociHighlights.clearHighlights()}>
                <span dangerouslySetInnerHTML={{ __html: e.label }} />
            </Button>
        </div>;
    }

    historyEntry(e: StructureSelectionHistoryEntry, idx: number) {
        const history = this.plugin.managers.structure.selection.additionsHistory;
        return <div className='msp-flex-row' key={e.id}>
            <Button noOverflow title='Click to focus. Hover to highlight.' onClick={() => this.focusLoci(e.loci)} style={{ width: 'auto', textAlign: 'left' }} onMouseEnter={() => this.highlight(e.loci)} onMouseLeave={() => this.plugin.managers.interactivity.lociHighlights.clearHighlights()}>
                {idx}. <span dangerouslySetInnerHTML={{ __html: e.label }} />
            </Button>
            {history.length > 1 && <IconButton svg={ArrowUpwardSvg} small={true} className='msp-form-control' onClick={() => this.moveHistory(e, 'up')} flex='20px' title={'Move up'} />}
            {history.length > 1 && <IconButton svg={ArrowDownwardSvg} small={true} className='msp-form-control' onClick={() => this.moveHistory(e, 'down')} flex='20px' title={'Move down'} />}
            <IconButton svg={DeleteOutlinedSvg} small={true} className='msp-form-control' onClick={() => this.plugin.managers.structure.selection.modifyHistory(e, 'remove')} flex title={'Remove'} />
        </div>;
    }

    atomsLociEntry(e: AtomsLociEntry, idx: number) {
        return <div key={idx}>
            <div className='msp-control-group-header'>
                <div className='msp-no-overflow' title={e.label}>{e.label}</div>
            </div>
            <div className='msp-control-offset'>
                {e.atoms.map((h, i) => this.historyEntry(h, i))}
            </div>
        </div>;
    }

    get chainEntries() {
        const location = StructureElement.Location.create();
        const entries: LociEntry[] = [];
        this.plugin.managers.structure.selection.entries.forEach(({ selection }, ref) => {
            const cell = StateObjectRef.resolveAndCheck(this.plugin.state.data, ref);
            if (!cell || StructureElement.Loci.isEmpty(selection)) return;

            // only single polymer chain selections
            const l = StructureElement.Loci.getFirstLocation(selection, location)!;
            if (selection.elements.length > 1 || StructureProperties.entity.type(l) !== 'polymer') return;

            const stats = StructureElement.Stats.ofLoci(selection);
            const counts = structureElementStatsLabel(stats, { countsOnly: true });
            const chain = elementLabel(l, { reverse: true, granularity: 'chain' }).split('|');
            const label = `${counts} | ${chain[0]} | ${chain[chain.length - 1]}`;
            entries.push({ loci: selection, label, cell });
        });
        return entries;
    }

    get atomEntries() {
        const structureEntries = new Map<Structure, StructureSelectionHistoryEntry[]>();
        const history = this.plugin.managers.structure.selection.additionsHistory;

        for (let i = 0, il = history.length; i < il; ++i) {
            const e = history[i];
            if (StructureElement.Loci.size(e.loci) !== 1) continue;

            const k = e.loci.structure;
            if (structureEntries.has(k)) structureEntries.get(k)!.push(e);
            else structureEntries.set(k, [e]);
        }

        const entries: AtomsLociEntry[] = [];
        structureEntries.forEach((atoms, structure) => {
            const cell = this.plugin.helpers.substructureParent.get(structure)!;

            const elements: StructureElement.Loci['elements'][0][] = [];
            for (let i = 0, il = atoms.length; i < il; ++i) {
                // note, we don't do loci union here to keep order of selected atoms
                // for atom pairing during superposition
                elements.push(atoms[i].loci.elements[0]);
            }

            const loci = StructureElement.Loci(atoms[0].loci.structure, elements);
            const label = loci.structure.label.split(' | ')[0];
            entries.push({ loci, label, cell, atoms });
        });
        return entries;
    }

    toggleHint() {
        const shouldShowToggleHint = this.plugin.config.get(PluginConfig.Viewport.ShowSelectionMode);
        return shouldShowToggleHint ? (<>{' '}(toggle <ToggleSelectionModeButton inline /> mode)</>) : null;
    }

    addByChains() {
        const entries = this.chainEntries;
        return <>
            {entries.length > 0 && <div className='msp-control-offset'>
                {entries.map((e, i) => this.lociEntry(e, i))}
            </div>}
            {entries.length < 2 && <div className='msp-control-offset msp-help-text'>
                <div className='msp-help-description'><Icon svg={HelpOutlineSvg} inline />Add 2 or more selections{this.toggleHint()} from separate structures. Selections must be limited to single polymer chains or residues therein.</div>
            </div>}
            {entries.length > 1 && <Button title='Superpose structures by selected chains.' className='msp-btn-commit msp-btn-commit-on' onClick={this.superposeChains} style={{ marginTop: '1px' }}>
                Superpose
            </Button>}
        </>;
    }

    addByAtoms() {
        const entries = this.atomEntries;
        return <>
            {entries.length > 0 && <div className='msp-control-offset'>
                {entries.map((e, i) => this.atomsLociEntry(e, i))}
            </div>}
            {entries.length < 2 && <div className='msp-control-offset msp-help-text'>
                <div className='msp-help-description'><Icon svg={HelpOutlineSvg} inline />Add 1 or more selections{this.toggleHint()} from
                separate structures. Selections must be limited to single atoms.</div>
            </div>}
            {entries.length > 1 && <Button title='Superpose structures by selected atoms.' className='msp-btn-commit msp-btn-commit-on' onClick={this.superposeAtoms} style={{ marginTop: '1px' }}>
                Superpose
            </Button>}
        </>;
    }

    superposeByDbMapping() {
        return <>
            <Button icon={SuperposeChainsSvg} title='Superpose structures using intersection of residues from SIFTS UNIPROT mapping.' className='msp-btn msp-btn-block' onClick={this.superposeDb} style={{ marginTop: '1px' }} disabled={this.state.isBusy}>
                Uniprot
            </Button>
        </>;
    }

    private setOptions = (values: StructureSuperpositionOptions) => {
        this.setState({ options: values });
    };

    render() {
        return <>
            <div className='msp-flex-row'>
                <ToggleButton icon={SuperposeChainsSvg} label='Chains' toggle={this.toggleByChains} isSelected={this.state.action === 'byChains'} disabled={this.state.isBusy} />
                <ToggleButton icon={SuperposeAtomsSvg} label='Atoms' toggle={this.toggleByAtoms} isSelected={this.state.action === 'byAtoms'} disabled={this.state.isBusy} />
                {this.state.canUseDb && this.superposeByDbMapping()}
                <ToggleButton icon={TuneSvg} label='' title='Options' toggle={this.toggleOptions} isSelected={this.state.action === 'options'} disabled={this.state.isBusy} style={{ flex: '0 0 40px', padding: 0 }} />
            </div>
            {this.state.action === 'byChains' && this.addByChains()}
            {this.state.action === 'byAtoms' && this.addByAtoms()}
            {this.state.action === 'options' && <div className='msp-control-offset'>
                <ParameterControls params={StructureSuperpositionParams} values={this.state.options} onChangeValues={this.setOptions} isDisabled={this.state.isBusy} />
            </div>}
        </>;
    }
}