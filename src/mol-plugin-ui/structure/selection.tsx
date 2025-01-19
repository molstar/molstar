/**
 * Copyright (c) 2019-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Jason Pattle <jpattle.exscientia.co.uk>
 * @author Ludovic Autin <ludovic.autin@gmail.com>
 * @author Ventura Rivera <venturaxrivera@gmail.com>
 */

import * as React from 'react';
import { Structure } from '../../mol-model/structure/structure/structure';
import { getElementQueries, getNonStandardResidueQueries, getPolymerAndBranchedEntityQueries, StructureSelectionQueries, StructureSelectionQuery } from '../../mol-plugin-state/helpers/structure-selection-query';
import { InteractivityManager } from '../../mol-plugin-state/manager/interactivity';
import { StructureComponentManager } from '../../mol-plugin-state/manager/structure/component';
import { StructureComponentRef, StructureRef } from '../../mol-plugin-state/manager/structure/hierarchy-state';
import { StructureSelectionModifier } from '../../mol-plugin-state/manager/structure/selection';
import { PluginConfig } from '../../mol-plugin/config';
import { PluginContext } from '../../mol-plugin/context';
import { compileIdListSelection } from '../../mol-script/util/id-list';
import { memoizeLatest } from '../../mol-util/memoize';
import { ParamDefinition } from '../../mol-util/param-definition';
import { capitalize, stripTags } from '../../mol-util/string';
import { PluginUIComponent, PurePluginUIComponent } from '../base';
import { ActionMenu } from '../controls/action-menu';
import { Button, ControlGroup, IconButton, ToggleButton } from '../controls/common';
import { BrushSvg, CancelOutlinedSvg, CloseSvg, CubeOutlineSvg, HelpOutlineSvg, Icon, IntersectSvg, RemoveSvg, RestoreSvg, SelectionModeSvg, SetSvg, SubtractSvg, UnionSvg } from '../controls/icons';
import { ParameterControls, ParamOnChange, PureSelectControl } from '../controls/parameters';
import { HelpGroup, HelpText, ViewportHelpContent } from '../viewport/help';
import { AddComponentControls } from './components';


export class ToggleSelectionModeButton extends PurePluginUIComponent<{ inline?: boolean }> {
    componentDidMount() {
        this.subscribe(this.plugin.events.canvas3d.settingsUpdated, () => this.forceUpdate());
        this.subscribe(this.plugin.layout.events.updated, () => this.forceUpdate());
        this.subscribe(this.plugin.behaviors.interaction.selectionMode, () => this.forceUpdate());
    }

    _toggleSelMode = () => {
        this.plugin.selectionMode = !this.plugin.selectionMode;
    };

    render() {
        const style = this.props.inline
            ? { background: 'transparent', width: 'auto', height: 'auto', lineHeight: 'unset' }
            : { background: 'transparent' };
        return <IconButton svg={SelectionModeSvg} onClick={this._toggleSelMode} title={'Toggle Selection Mode'} style={style} toggleState={this.plugin.selectionMode} />;
    }
}

const StructureSelectionParams = {
    granularity: InteractivityManager.Params.granularity,
};

type SelectionHelperType = 'residue-list'

interface StructureSelectionActionsControlsState {
    isEmpty: boolean,
    isBusy: boolean,
    canUndo: boolean,

    action?: StructureSelectionModifier | 'theme' | 'add-component' | 'help',
    helper?: SelectionHelperType,

    structureSelectionParams?: typeof StructureSelectionParams,
}

const ActionHeader = new Map<StructureSelectionModifier, string>([
    ['add', 'Add/Union Selection'],
    ['remove', 'Remove/Subtract Selection'],
    ['intersect', 'Intersect Selection'],
    ['set', 'Set Selection']
] as const);

export class StructureSelectionActionsControls extends PluginUIComponent<{}, StructureSelectionActionsControlsState> {
    state = {
        action: void 0 as StructureSelectionActionsControlsState['action'],
        helper: void 0 as StructureSelectionActionsControlsState['helper'],

        isEmpty: true,
        isBusy: false,
        canUndo: false,

        structureSelectionParams: StructureSelectionParams,
    };

    componentDidMount() {
        this.subscribe(this.plugin.managers.structure.hierarchy.behaviors.selection, c => {
            const isEmpty = c.hierarchy.structures.length === 0;
            if (this.state.isEmpty !== isEmpty) {
                this.setState({ isEmpty });
            }
            // trigger elementQueries and nonStandardResidueQueries recalculation
            this.queriesVersion = -1;
            this.forceUpdate();
        });

        this.subscribe(this.plugin.behaviors.state.isBusy, v => {
            this.setState({ isBusy: v, action: void 0 });
        });

        this.subscribe(this.plugin.managers.interactivity.events.propsUpdated, () => {
            this.forceUpdate();
        });

        this.subscribe(this.plugin.state.data.events.historyUpdated, ({ state }) => {
            this.setState({ canUndo: state.canUndo });
        });

        // Update structureSelectionParams state if there are custom-defined granularityOptions
        const granularityOptions = this.plugin.spec.components?.selectionTools?.granularityOptions;
        if (granularityOptions) {
            const granularitySet = new Set((granularityOptions));
            const structureSelectionParams = {
                ...StructureSelectionParams,
                granularity: {
                    ...StructureSelectionParams.granularity,
                    options: StructureSelectionParams.granularity.options.filter(([firstItem]) => granularitySet.has(firstItem)),
                },
            };
            this.setState({ structureSelectionParams: structureSelectionParams });
        }
    }

    get isDisabled() {
        return this.state.isBusy || this.state.isEmpty;
    }

    set = (modifier: StructureSelectionModifier, selectionQuery: StructureSelectionQuery) => {
        this.plugin.managers.structure.selection.fromSelectionQuery(modifier, selectionQuery, false);
    };

    selectQuery: ActionMenu.OnSelect = (item, e) => {
        if (!item || !this.state.action) {
            this.setState({ action: void 0 });
            return;
        }
        const q = this.state.action! as StructureSelectionModifier;
        if (e?.shiftKey) {
            this.set(q, item.value as StructureSelectionQuery);
        } else {
            this.setState({ action: void 0 }, () => {
                this.set(q, item.value as StructureSelectionQuery);
            });
        }
    };

    selectHelper: ActionMenu.OnSelect = (item, e) => {
        console.log(item);
        if (!item || !this.state.action) {
            this.setState({ action: void 0, helper: void 0 });
            return;
        }
        this.setState({ helper: (item.value as { kind: SelectionHelperType }).kind });
    };

    get structures() {
        const structures: Structure[] = [];
        for (const s of this.plugin.managers.structure.hierarchy.selection.structures) {
            const structure = s.cell.obj?.data;
            if (structure) structures.push(structure);
        }
        return structures;
    }

    private queriesItems: ActionMenu.Items[] = [];
    private queriesVersion = -1;
    get queries() {
        const { registry } = this.plugin.query.structure;
        if (registry.version !== this.queriesVersion) {
            const structures = this.structures;
            const queries = [
                ...registry.list,
                ...getPolymerAndBranchedEntityQueries(structures),
                ...getNonStandardResidueQueries(structures),
                ...getElementQueries(structures)
            ].sort((a, b) => b.priority - a.priority);
            this.queriesItems = ActionMenu.createItems(queries, {
                filter: q => q !== StructureSelectionQueries.current && !q.isHidden,
                label: q => q.label,
                category: q => q.category,
                description: q => q.description
            });
            this.queriesVersion = registry.version;
        }
        return this.queriesItems;
    }

    private helpersItems?: ActionMenu.Items[] = void 0;
    get helpers() {
        if (this.helpersItems) return this.helpersItems;
        // TODO: this is an initial implementation of the helper UI
        //       the plan is to add support to input queries in different languages
        //       after this has been implemented in mol-script
        const helpers = [
            { kind: 'residue-list' as SelectionHelperType, category: 'Helpers', label: 'Atom/Residue Identifier List', description: 'Create a selection from a list of atom/residue ranges.' }
        ];
        this.helpersItems = ActionMenu.createItems(helpers, {
            label: q => q.label,
            category: q => q.category,
            description: q => q.description
        });
        return this.helpersItems;
    }

    private showAction(q: StructureSelectionActionsControlsState['action']) {
        return () => this.setState({ action: this.state.action === q ? void 0 : q, helper: void 0 });
    }

    toggleAdd = this.showAction('add');
    toggleRemove = this.showAction('remove');
    toggleIntersect = this.showAction('intersect');
    toggleSet = this.showAction('set');
    toggleTheme = this.showAction('theme');
    toggleAddComponent = this.showAction('add-component');
    toggleHelp = this.showAction('help');

    setGranuality: ParamOnChange = ({ value }) => {
        this.plugin.managers.interactivity.setProps({ granularity: value });
    };

    turnOff = () => this.plugin.selectionMode = false;

    undo = () => {
        const task = this.plugin.state.data.undo();
        if (task) this.plugin.runTask(task);
    };

    subtract = () => {
        const sel = this.plugin.managers.structure.hierarchy.getStructuresWithSelection();
        const components: StructureComponentRef[] = [];
        for (const s of sel) components.push(...s.components);
        if (components.length === 0) return;
        this.plugin.managers.structure.component.modifyByCurrentSelection(components, 'subtract');
    };

    render() {
        const granularity = this.plugin.managers.interactivity.props.granularity;
        const hide = this.plugin.spec.components?.selectionTools?.hide;
        const undoTitle = this.state.canUndo
            ? `Undo ${this.plugin.state.data.latestUndoLabel}`
            : 'Some mistakes of the past can be undone.';

        let children: React.ReactNode | undefined = void 0;

        if (this.state.action && !this.state.helper) {
            children = <>
                {(this.state.action && this.state.action !== 'theme' && this.state.action !== 'add-component' && this.state.action !== 'help') && <div className='msp-selection-viewport-controls-actions'>
                    <ActionMenu header={ActionHeader.get(this.state.action as StructureSelectionModifier)} title='Click to close.' items={this.queries} onSelect={this.selectQuery} noOffset />
                    <ActionMenu items={this.helpers} onSelect={this.selectHelper} noOffset />
                </div>}
                {this.state.action === 'theme' && <div className='msp-selection-viewport-controls-actions'>
                    <ControlGroup header='Theme' title='Click to close.' initialExpanded={true} hideExpander={true} hideOffset={true} onHeaderClick={this.toggleTheme} topRightIcon={CloseSvg}>
                        <ApplyThemeControls onApply={this.toggleTheme} />
                    </ControlGroup>
                </div>}
                {this.state.action === 'add-component' && <div className='msp-selection-viewport-controls-actions'>
                    <ControlGroup header='Add Component' title='Click to close.' initialExpanded={true} hideExpander={true} hideOffset={true} onHeaderClick={this.toggleAddComponent} topRightIcon={CloseSvg}>
                        <AddComponentControls onApply={this.toggleAddComponent} forSelection />
                    </ControlGroup>
                </div>}
                {this.state.action === 'help' && <div className='msp-selection-viewport-controls-actions'>
                    <ControlGroup header='Help' title='Click to close.' initialExpanded={true} hideExpander={true} hideOffset={true} onHeaderClick={this.toggleHelp} topRightIcon={CloseSvg} maxHeight='300px'>
                        <HelpGroup header='Selection Operations'>
                            <HelpText>Use <Icon svg={UnionSvg} inline /> <Icon svg={SubtractSvg} inline /> <Icon svg={IntersectSvg} inline /> <Icon svg={SetSvg} inline /> to modify the selection.</HelpText>
                        </HelpGroup>
                        <HelpGroup header='Representation Operations'>
                            <HelpText>Use <Icon svg={BrushSvg} inline /> <Icon svg={CubeOutlineSvg} inline /> <Icon svg={RemoveSvg} inline /> <Icon svg={RestoreSvg} inline /> to color, create components, remove from components, or undo actions.</HelpText>
                        </HelpGroup>
                        <ViewportHelpContent selectOnly={true} />
                    </ControlGroup>
                </div>}
            </>;
        } else if (ActionHeader.has(this.state.action as any) && this.state.helper === 'residue-list') {
            const close = () => this.setState({ action: void 0, helper: void 0 });
            children = <div className='msp-selection-viewport-controls-actions'>
                <ControlGroup header='Atom/Residue Identifier List' title='Click to close.' initialExpanded={true} hideExpander={true} hideOffset={true} onHeaderClick={close} topRightIcon={CloseSvg}>
                    <ResidueListSelectionHelper modifier={this.state.action as any} plugin={this.plugin} close={close} />
                </ControlGroup>
            </div>;
        }

        return <>
            <div className='msp-flex-row' style={{ background: 'none' }}>
                {(!hide?.granularity) && <PureSelectControl title={`Picking Level for selecting and highlighting`} param={this.state.structureSelectionParams.granularity} name='granularity' value={granularity} onChange={this.setGranuality} isDisabled={this.isDisabled} />}
                {(!hide?.union) && <ToggleButton icon={UnionSvg} title={`${ActionHeader.get('add')}. Hold shift key to keep menu open.`} toggle={this.toggleAdd} isSelected={this.state.action === 'add'} disabled={this.isDisabled} />}
                {(!hide?.subtract) && <ToggleButton icon={SubtractSvg} title={`${ActionHeader.get('remove')}. Hold shift key to keep menu open.`} toggle={this.toggleRemove} isSelected={this.state.action === 'remove'} disabled={this.isDisabled} />}
                {(!hide?.intersect) && <ToggleButton icon={IntersectSvg} title={`${ActionHeader.get('intersect')}. Hold shift key to keep menu open.`} toggle={this.toggleIntersect} isSelected={this.state.action === 'intersect'} disabled={this.isDisabled} />}
                {(!hide?.set) && <ToggleButton icon={SetSvg} title={`${ActionHeader.get('set')}. Hold shift key to keep menu open.`} toggle={this.toggleSet} isSelected={this.state.action === 'set'} disabled={this.isDisabled} />}

                {(!hide?.theme) && <ToggleButton icon={BrushSvg} title='Apply Theme to Selection' toggle={this.toggleTheme} isSelected={this.state.action === 'theme'} disabled={this.isDisabled} style={{ marginLeft: '10px' }} />}
                {(!hide?.componentAdd) && <ToggleButton icon={CubeOutlineSvg} title='Create Component of Selection with Representation' toggle={this.toggleAddComponent} isSelected={this.state.action === 'add-component'} disabled={this.isDisabled} />}
                {(!hide?.componentRemove) && <IconButton svg={RemoveSvg} title='Remove/subtract Selection from all Components' onClick={this.subtract} disabled={this.isDisabled} />}
                {(!hide?.undo) && <IconButton svg={RestoreSvg} onClick={this.undo} disabled={!this.state.canUndo || this.isDisabled} title={undoTitle} />}

                {(!hide?.help) && <ToggleButton icon={HelpOutlineSvg} title='Show/hide help' toggle={this.toggleHelp} style={{ marginLeft: '10px' }} isSelected={this.state.action === 'help'} />}
                {((!hide?.cancel) && this.plugin.config.get(PluginConfig.Viewport.ShowSelectionMode)) && (<IconButton svg={CancelOutlinedSvg} title='Turn selection mode off' onClick={this.turnOff} />)}
            </div>
            {children}
        </>;
    }
}

export class StructureSelectionStatsControls extends PluginUIComponent<{ hideOnEmpty?: boolean }, { isEmpty: boolean, isBusy: boolean }> {
    state = {
        isEmpty: true,
        isBusy: false
    };

    componentDidMount() {
        this.subscribe(this.plugin.managers.structure.selection.events.changed, () => {
            this.forceUpdate();
        });

        this.subscribe(this.plugin.managers.structure.hierarchy.behaviors.selection, c => {
            const isEmpty = c.structures.length === 0;
            if (this.state.isEmpty !== isEmpty) {
                this.setState({ isEmpty });
            }
        });

        this.subscribe(this.plugin.behaviors.state.isBusy, v => {
            this.setState({ isBusy: v });
        });
    }

    get isDisabled() {
        return this.state.isBusy || this.state.isEmpty;
    }

    get stats() {
        const stats = this.plugin.managers.structure.selection.stats;
        if (stats.structureCount === 0 || stats.elementCount === 0) {
            return 'Nothing Selected';
        } else {
            return `${stripTags(stats.label)} Selected`;
        }
    }

    clear = () => this.plugin.managers.interactivity.lociSelects.deselectAll();

    focus = () => {
        if (this.plugin.managers.structure.selection.stats.elementCount === 0) return;
        const { sphere } = this.plugin.managers.structure.selection.getBoundary();
        this.plugin.managers.camera.focusSphere(sphere);
    };

    highlight = (e: React.MouseEvent<HTMLElement>) => {
        this.plugin.managers.interactivity.lociHighlights.clearHighlights();
        this.plugin.managers.structure.selection.entries.forEach(e => {
            this.plugin.managers.interactivity.lociHighlights.highlight({ loci: e.selection }, false);
        });
    };

    clearHighlight = () => {
        this.plugin.managers.interactivity.lociHighlights.clearHighlights();
    };

    render() {
        const stats = this.plugin.managers.structure.selection.stats;
        const empty = stats.structureCount === 0 || stats.elementCount === 0;

        if (empty && this.props.hideOnEmpty) return null;

        return <>
            <div className='msp-flex-row'>
                <Button noOverflow onClick={this.focus} title='Click to Focus Selection' disabled={empty} onMouseEnter={this.highlight} onMouseLeave={this.clearHighlight}
                    style={{ textAlignLast: !empty ? 'left' : void 0 }}>
                    {this.stats}
                </Button>
                {!empty && <IconButton svg={CancelOutlinedSvg} onClick={this.clear} title='Clear' className='msp-form-control' flex />}
            </div>
        </>;
    }
}

interface ApplyThemeControlsState {
    values: StructureComponentManager.ThemeParams
}

interface ApplyThemeControlsProps {
    onApply?: () => void
}

class ApplyThemeControls extends PurePluginUIComponent<ApplyThemeControlsProps, ApplyThemeControlsState> {
    _params = memoizeLatest((pivot: StructureRef | undefined) => StructureComponentManager.getThemeParams(this.plugin, pivot));
    get params() { return this._params(this.plugin.managers.structure.component.pivotStructure); }

    state = { values: ParamDefinition.getDefaultValues(this.params) };

    apply = () => {
        this.plugin.managers.structure.component.applyTheme(this.state.values, this.plugin.managers.structure.hierarchy.current.structures);
        this.props.onApply?.();
    };

    paramsChanged = (values: any) => this.setState({ values });

    render() {
        return <>
            <ParameterControls params={this.params} values={this.state.values} onChangeValues={this.paramsChanged} />
            <Button icon={BrushSvg} className='msp-btn-commit msp-btn-commit-on' onClick={this.apply} style={{ marginTop: '1px' }}>
                Apply Theme
            </Button>
        </>;
    }
}

const ResidueListIdTypeParams = {
    idType: ParamDefinition.Select<'auth' | 'label' | 'atom-id' | 'element-symbol'>('auth', ParamDefinition.arrayToOptions(['auth', 'label', 'atom-id', 'element-symbol'])),
    identifiers: ParamDefinition.Text('', { description: 'A comma separated list of atom identifiers (e.g. 10, 15-25), element symbols (e.g. N, C or 20-200) or residue ranges in given chain (e.g. A 10-15, B 25, C 30:i)' })
};

const DefaultResidueListIdTypeParams = ParamDefinition.getDefaultValues(ResidueListIdTypeParams);

function ResidueListSelectionHelper({ modifier, plugin, close }: { modifier: StructureSelectionModifier, plugin: PluginContext, close: () => void }) {
    const [state, setState] = React.useState(DefaultResidueListIdTypeParams);

    const apply = () => {
        if (state.identifiers.trim().length === 0) return;

        try {
            close();
            const query = compileIdListSelection(state.identifiers, state.idType);
            plugin.managers.structure.selection.fromCompiledQuery(modifier, query, false);
        } catch (e) {
            console.error(e);
            plugin.log.error('Failed to create selection');
        }
    };

    return <>
        <ParameterControls params={ResidueListIdTypeParams} values={state} onChangeValues={setState} onEnter={apply} />
        <Button className='msp-btn-commit msp-btn-commit-on' disabled={state.identifiers.trim().length === 0} onClick={apply} style={{ marginTop: '1px' }}>
            {capitalize(modifier)} Selection
        </Button>
    </>;
}