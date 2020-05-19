/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react';
import { Loci } from '../../mol-model/loci';
import { StructureElement } from '../../mol-model/structure';
import { StructureMeasurementCell, StructureMeasurementOptions, StructureMeasurementParams } from '../../mol-plugin-state/manager/structure/measurement';
import { StructureSelectionHistoryEntry } from '../../mol-plugin-state/manager/structure/selection';
import { PluginStateObject } from '../../mol-plugin-state/objects';
import { PluginCommands } from '../../mol-plugin/commands';
import { angleLabel, dihedralLabel, distanceLabel, lociLabel } from '../../mol-theme/label';
import { FiniteArray } from '../../mol-util/type-helpers';
import { CollapsableControls, PurePluginUIComponent } from '../base';
import { ActionMenu } from '../controls/action-menu';
import { Button, ExpandGroup, IconButton, ToggleButton } from '../controls/common';
import { Icon, PencilRulerSvg, SetSvg, ArrowUpwardSvg, ArrowDownwardSvg, DeleteOutlinedSvg, HelpOutlineSvg, AddSvg, TuneSvg, VisibilityOffOutlinedSvg, VisibilityOutlinedSvg, MoreHorizSvg } from '../controls/icons';
import { ParameterControls } from '../controls/parameters';
import { UpdateTransformControl } from '../state/update-transform';
import { ToggleSelectionModeButton } from './selection';

// TODO details, options (e.g. change text for labels)

export class StructureMeasurementsControls extends CollapsableControls {
    defaultState() {
        return {
            isCollapsed: false,
            header: 'Measurements',
            brand: { accent: 'gray' as const, svg: PencilRulerSvg }
        };
    }

    renderControls() {
        return <>
            <MeasurementControls />
            <MeasurementList />
        </>;
    }
}

export class MeasurementList extends PurePluginUIComponent {
    componentDidMount() {
        this.subscribe(this.plugin.managers.structure.measurement.behaviors.state, () => {
            this.forceUpdate();
        });
    }

    renderGroup(cells: ReadonlyArray<StructureMeasurementCell>, header: string) {
        const group: JSX.Element[] = [];
        for (const cell of cells) {
            if (cell.obj) group.push(<MeasurementEntry key={cell.obj.id} cell={cell} />);
        }
        return group.length ? <ExpandGroup header={header} initiallyExpanded={true}>{group}</ExpandGroup> : null;
    }

    render() {
        const measurements = this.plugin.managers.structure.measurement.state;

        return <div style={{ marginTop: '6px' }}>
            {this.renderGroup(measurements.labels, 'Labels')}
            {this.renderGroup(measurements.distances, 'Distances')}
            {this.renderGroup(measurements.angles, 'Angles')}
            {this.renderGroup(measurements.dihedrals, 'Dihedrals')}
        </div>;
    }
}

export class MeasurementControls extends PurePluginUIComponent<{}, { isBusy: boolean, action?: 'add' | 'options' }> {
    state = { isBusy: false, action: void 0 as 'add' | 'options' | undefined }

    componentDidMount() {
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

    measureDistance = () => {
        const loci = this.plugin.managers.structure.selection.additionsHistory;
        this.plugin.managers.structure.measurement.addDistance(loci[0].loci, loci[1].loci);
    }

    measureAngle = () => {
        const loci = this.plugin.managers.structure.selection.additionsHistory;
        this.plugin.managers.structure.measurement.addAngle(loci[0].loci, loci[1].loci, loci[2].loci);
    }

    measureDihedral = () => {
        const loci = this.plugin.managers.structure.selection.additionsHistory;
        this.plugin.managers.structure.measurement.addDihedral(loci[0].loci, loci[1].loci, loci[2].loci, loci[3].loci);
    }

    addLabel = () => {
        const loci = this.plugin.managers.structure.selection.additionsHistory;
        this.plugin.managers.structure.measurement.addLabel(loci[0].loci);
    }

    get actions(): ActionMenu.Items {
        const history = this.selection.additionsHistory;
        const ret: ActionMenu.Item[] = [
            { kind: 'item', label: `Label ${history.length === 0 ? ' (1 selection required)' : ' (1st selection)'}`, value: this.addLabel, disabled: history.length === 0 },
            { kind: 'item', label: `Distance ${history.length < 2 ? ' (2 selections required)' : ' (top 2 selections)'}`, value: this.measureDistance, disabled: history.length < 2 },
            { kind: 'item', label: `Angle ${history.length < 3 ? ' (3 selections required)' : ' (top 3 selections)'}`, value: this.measureAngle, disabled: history.length < 3 },
            { kind: 'item', label: `Dihedral ${history.length < 4 ? ' (4 selections required)' : ' (top 4 selections)'}`, value: this.measureDihedral, disabled: history.length < 4 },
        ];
        return ret;
    }

    selectAction: ActionMenu.OnSelect = item => {
        this.toggleAdd();
        if (!item) return;
        (item?.value as any)();
    }

    toggleAdd = () => this.setState({ action: this.state.action === 'add' ? void 0 : 'add' });
    toggleOptions = () => this.setState({ action: this.state.action === 'options' ? void 0 : 'options'  });

    highlight(loci: StructureElement.Loci) {
        this.plugin.managers.interactivity.lociHighlights.highlightOnly({ loci }, false);
    }

    moveHistory(e: StructureSelectionHistoryEntry, direction: 'up' | 'down') {
        this.plugin.managers.structure.selection.modifyHistory(e, direction, 4);
    }

    focusLoci(loci: StructureElement.Loci) {
        this.plugin.managers.camera.focusLoci(loci);
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

    add() {
        const history = this.plugin.managers.structure.selection.additionsHistory;

        const entries: JSX.Element[] = [];
        for (let i = 0, _i = Math.min(history.length, 4); i < _i; i++) {
            entries.push(this.historyEntry(history[i], i + 1));
        }

        return <>
            <ActionMenu items={this.actions} onSelect={this.selectAction} />
            {entries.length > 0 && <div className='msp-control-offset'>
                {entries}
            </div>}
            {entries.length === 0 && <div className='msp-control-offset msp-help-text'>
                <div className='msp-help-description'><Icon svg={HelpOutlineSvg} inline />Add one or more selections (toggle <ToggleSelectionModeButton inline /> mode)</div>
            </div>}
        </>;
    }

    render() {
        return <>
            <div className='msp-flex-row'>
                <ToggleButton icon={AddSvg} label='Add' toggle={this.toggleAdd} isSelected={this.state.action === 'add'} disabled={this.state.isBusy} className='msp-btn-apply-simple' />
                <ToggleButton icon={TuneSvg} label='' title='Options' toggle={this.toggleOptions} isSelected={this.state.action === 'options'} disabled={this.state.isBusy} style={{ flex: '0 0 40px', padding: 0 }} />
            </div>
            {this.state.action === 'add' && this.add()}
            {this.state.action === 'options' && <MeasurementsOptions />}
        </>;
    }
}

class MeasurementsOptions extends PurePluginUIComponent<{}, { isDisabled: boolean }> {
    state = { isDisabled: false }

    componentDidMount() {
        this.subscribe(this.plugin.managers.structure.measurement.behaviors.state, () => {
            this.forceUpdate();
        });

        this.subscribe(this.plugin.behaviors.state.isBusy, v => {
            this.setState({ isDisabled: v });
        });
    }

    changed = (options: StructureMeasurementOptions) => {
        this.plugin.managers.structure.measurement.setOptions(options);
    }

    render() {
        const measurements = this.plugin.managers.structure.measurement.state;

        return <div className='msp-control-offset'>
            <ParameterControls params={StructureMeasurementParams} values={measurements.options} onChangeValues={this.changed} isDisabled={this.state.isDisabled} />
        </div>;
    }
}

class MeasurementEntry extends PurePluginUIComponent<{ cell: StructureMeasurementCell }, { showUpdate: boolean }> {
    state = { showUpdate: false }

    componentDidMount() {
        this.subscribe(this.plugin.state.events.cell.stateUpdated, e => {
            this.forceUpdate();
        });
    }

    get selections() {
        return this.props.cell.obj?.data.source as PluginStateObject.Molecule.Structure.Selections | undefined;
    }

    delete = () => {
        PluginCommands.State.RemoveObject(this.plugin, { state: this.props.cell.parent!, ref: this.props.cell.transform.parent, removeParentGhosts: true });
    };

    toggleVisibility = (e: React.MouseEvent<HTMLElement>) => {
        e.preventDefault();
        PluginCommands.State.ToggleVisibility(this.plugin, { state: this.props.cell.parent!, ref: this.props.cell.transform.parent });
        e.currentTarget.blur();
    }

    highlight = () => {
        const selections = this.selections;
        if (!selections) return;

        this.plugin.managers.interactivity.lociHighlights.clearHighlights();
        for (const d of selections.data) {
            this.plugin.managers.interactivity.lociHighlights.highlight({ loci: d.loci }, false);
        }
        this.plugin.managers.interactivity.lociHighlights.highlight({ loci: this.props.cell.obj?.data.repr.getLoci()! }, false);
    }

    clearHighlight = () => {
        this.plugin.managers.interactivity.lociHighlights.clearHighlights();
    }

    toggleUpdate = () => this.setState({ showUpdate: !this.state.showUpdate });

    focus = () => {
        const selections = this.selections;
        if (!selections) return;

        const sphere = Loci.getBundleBoundingSphere(toLociBundle(selections.data));
        if (sphere) {
            this.plugin.managers.camera.focusSphere(sphere);
        }
    }

    get label() {
        const selections = this.selections;
        switch (selections?.data.length) {
            case 1: return lociLabel(selections.data[0].loci, { condensed: true });
            case 2: return distanceLabel(toLociBundle(selections.data), { condensed: true, unitLabel: this.plugin.managers.structure.measurement.state.options.distanceUnitLabel });
            case 3: return angleLabel(toLociBundle(selections.data), { condensed: true });
            case 4: return dihedralLabel(toLociBundle(selections.data), { condensed: true });
            default: return '';
        }
    }

    get actions(): ActionMenu.Items {
        this.props.cell.sourceRef;
        return [ActionMenu.Item('Select This', () => this.plugin.managers.structure.selection.fromSelections(this.props.cell.sourceRef!), { icon: SetSvg })];
    }

    selectAction: ActionMenu.OnSelect = item => {
        if (!item) return;
        this.setState({ showUpdate: false });
        (item?.value as any)();
    }

    render() {
        const { cell } = this.props;
        const { obj } = cell;
        if (!obj) return null;

        return <>
            <div className='msp-flex-row' key={obj.id} onMouseEnter={this.highlight} onMouseLeave={this.clearHighlight}>
                <button className='msp-form-control msp-control-button-label msp-no-overflow' title='Click to focus. Hover to highlight.' onClick={this.focus} style={{ width: 'auto', textAlign: 'left' }}>
                    <span dangerouslySetInnerHTML={{ __html: this.label }} />
                </button>
                <IconButton svg={cell.state.isHidden ? VisibilityOffOutlinedSvg : VisibilityOutlinedSvg} toggleState={false} small className='msp-form-control' onClick={this.toggleVisibility} flex title={cell.state.isHidden ? 'Show' : 'Hide'} />
                <IconButton svg={DeleteOutlinedSvg} small className='msp-form-control' onClick={this.delete} flex title='Delete' toggleState={false} />
                <IconButton svg={MoreHorizSvg} className='msp-form-control' onClick={this.toggleUpdate} flex title='Actions' toggleState={this.state.showUpdate} />
            </div>
            {this.state.showUpdate && cell.parent && <>
                <div className='msp-accent-offset'>
                    <ActionMenu items={this.actions} onSelect={this.selectAction} noOffset />
                    <ExpandGroup header='Options' noOffset>
                        <UpdateTransformControl state={cell.parent} transform={cell.transform} customHeader='none' autoHideApply />
                    </ExpandGroup>
                </div>
            </>}
        </>;
    }
}

function toLociBundle(data: FiniteArray<{ loci: Loci }, any>): { loci: FiniteArray<Loci, any> } {
    return { loci: (data.map(d => d.loci) as unknown as FiniteArray<Loci, any>) };
}