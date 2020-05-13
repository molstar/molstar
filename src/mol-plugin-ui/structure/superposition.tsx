/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react';
import { CollapsableControls, PurePluginUIComponent } from '../base';
import { Icon, ArrowUpwardSvg, ArrowDownwardSvg, DeleteOutlinedSvg, HelpOutlineSvg, TuneSvg, SuperposeAtomsSvg, SuperposeChainsSvg, SuperpositionSvg } from '../controls/icons';
import { Button, ToggleButton, IconButton } from '../controls/common';
import { StructureElement, StructureSelection, QueryContext, Structure, StructureProperties } from '../../mol-model/structure';
import { Mat4 } from '../../mol-math/linear-algebra';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { StateObjectRef, StateObjectCell, StateSelection } from '../../mol-state';
import { StateTransforms } from '../../mol-plugin-state/transforms';
import { PluginStateObject } from '../../mol-plugin-state/objects';
import { alignAndSuperpose, superpose } from '../../mol-model/structure/structure/util/superposition';
import { StructureSelectionQueries } from '../../mol-plugin-state/helpers/structure-selection-query';
import { structureElementStatsLabel, elementLabel } from '../../mol-theme/label';
import { ParameterControls } from '../controls/parameters';
import { stripTags } from '../../mol-util/string';
import { StructureSelectionHistoryEntry } from '../../mol-plugin-state/manager/structure/selection';
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
    alignSequences: PD.Boolean(true, { isEssential: true, description: 'Perform a sequence alignment and use the aligned residue pairs to guide the 3D superposition.' }),
};
const DefaultStructureSuperpositionOptions = PD.getDefaultValues(StructureSuperpositionParams);
export type StructureSuperpositionOptions = PD.ValuesFor<typeof StructureSuperpositionParams>

const SuperpositionTag = 'SuperpositionTransform';

type SuperpositionControlsState = {
    isBusy: boolean,
    action?: 'byChains' | 'byAtoms' | 'options',
    options: StructureSuperpositionOptions
}

interface LociEntry {
    loci: StructureElement.Loci,
    label: string,
    cell: StateObjectCell<PluginStateObject.Molecule.Structure>
};

interface AtomsLociEntry extends LociEntry {
    atoms: StructureSelectionHistoryEntry[]
};

export class SuperpositionControls extends PurePluginUIComponent<{}, SuperpositionControlsState> {
    state: SuperpositionControlsState = {
        isBusy: false,
        action: undefined,
        options: DefaultStructureSuperpositionOptions
    }

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
    }

    get selection() {
        return this.plugin.managers.structure.selection;
    }

    async transform(s: StateObjectRef<PluginStateObject.Molecule.Structure>, matrix: Mat4) {
        const r = StateObjectRef.resolveAndCheck(this.plugin.state.data, s);
        if (!r) return;
        // TODO should find any TransformStructureConformation decorator instance
        const o = StateSelection.findTagInSubtree(this.plugin.state.data.tree, r.transform.ref, SuperpositionTag);

        const params = {
            transform: {
                name: 'matrix' as const,
                params: { data: matrix, transpose: false }
            }
        };
        // TODO add .insertOrUpdate to StateBuilder?
        const b = o
            ? this.plugin.state.data.build().to(o).update(params)
            : this.plugin.state.data.build().to(s)
                .insert(StateTransforms.Model.TransformStructureConformation, params, { tags: SuperpositionTag });
        await this.plugin.runTask(this.plugin.state.data.updateTree(b));
    }

    superposeChains = async () => {
        const { query } = StructureSelectionQueries.trace;
        const entries = this.chainEntries;

        const traceLocis = entries.map((e, i) => {
            const s = StructureElement.Loci.toStructure(e.loci);
            const loci = StructureSelection.toLociWithSourceUnits(query(new QueryContext(s)));
            return StructureElement.Loci.remap(loci, i === 0
                ? this.plugin.helpers.substructureParent.get(e.loci.structure.root)!.obj!.data
                : loci.structure.root
            );
        });

        const transforms = this.state.options.alignSequences
            ? alignAndSuperpose(traceLocis)
            : superpose(traceLocis);

        const eA = entries[0];
        for (let i = 1, il = traceLocis.length; i < il; ++i) {
            const eB = entries[i];
            const { bTransform, rmsd } = transforms[i - 1];
            await this.transform(eB.cell, bTransform);
            const labelA = stripTags(eA.label);
            const labelB = stripTags(eB.label);
            this.plugin.log.info(`Superposed [${labelA}] and [${labelB}] with RMSD ${rmsd.toFixed(2)}.`);
        }
    }

    superposeAtoms = async () => {
        const entries = this.atomEntries;

        const atomLocis = entries.map((e, i) => {
            return StructureElement.Loci.remap(e.loci, i === 0
                ? this.plugin.helpers.substructureParent.get(e.loci.structure.root)!.obj!.data
                : e.loci.structure.root
            );
        });
        const transforms = superpose(atomLocis);

        const eA = entries[0];
        for (let i = 1, il = atomLocis.length; i < il; ++i) {
            const eB = entries[i];
            const { bTransform, rmsd } = transforms[i - 1];
            await this.transform(eB.cell, bTransform);
            const labelA = stripTags(eA.label);
            const labelB = stripTags(eB.label);
            const count = entries[i].atoms.length;
            this.plugin.log.info(`Superposed ${count} ${count === 1 ? 'atom' : 'atoms'} of [${labelA}] and [${labelB}] with RMSD ${rmsd.toFixed(2)}.`);
        }
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
            <Button noOverflow title='Click to focus. Hover to highlight.' onClick={() => this.focusLoci(e.loci)} style={{ width: 'auto', textAlign: 'left' }} onMouseEnter={() => this.highlight(e.loci)} onMouseLeave={this.plugin.managers.interactivity.lociHighlights.clearHighlights}>
                <span dangerouslySetInnerHTML={{ __html: e.label }} />
            </Button>
        </div>;
    }

    historyEntry(e: StructureSelectionHistoryEntry, idx: number) {
        const history = this.plugin.managers.structure.selection.additionsHistory;
        return <div className='msp-flex-row' key={e.id}>
            <Button noOverflow title='Click to focus. Hover to highlight.' onClick={() => this.focusLoci(e.loci)} style={{ width: 'auto', textAlign: 'left' }} onMouseEnter={() => this.highlight(e.loci)} onMouseLeave={this.plugin.managers.interactivity.lociHighlights.clearHighlights}>
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

    addByChains() {
        const entries = this.chainEntries;
        return <>
            {entries.length > 0 && <div className='msp-control-offset'>
                {entries.map((e, i) => this.lociEntry(e, i))}
            </div>}
            {entries.length < 2 && <div className='msp-control-offset msp-help-text'>
                <div className='msp-help-description'><Icon svg={HelpOutlineSvg} inline />Add 2 or more selections (toggle <ToggleSelectionModeButton inline /> mode) from separate structures. Selections must be limited to single polymer chains or residues therein.</div>
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
                <div className='msp-help-description'><Icon svg={HelpOutlineSvg} inline />Add 1 or more selections (toggle <ToggleSelectionModeButton inline /> mode) from
                separate structures. Selections must be limited to single atoms.</div>
            </div>}
            {entries.length > 1 && <Button title='Superpose structures by selected atoms.' className='msp-btn-commit msp-btn-commit-on' onClick={this.superposeAtoms} style={{ marginTop: '1px' }}>
                Superpose
            </Button>}
        </>;
    }

    private setOptions = (values: StructureSuperpositionOptions) => {
        this.setState({ options: values });
    }

    render() {
        return <>
            <div className='msp-flex-row'>
                <ToggleButton icon={SuperposeChainsSvg} label='By Chains' toggle={this.toggleByChains} isSelected={this.state.action === 'byChains'} disabled={this.state.isBusy} />
                <ToggleButton icon={SuperposeAtomsSvg} label='By Atoms' toggle={this.toggleByAtoms} isSelected={this.state.action === 'byAtoms'} disabled={this.state.isBusy} />
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