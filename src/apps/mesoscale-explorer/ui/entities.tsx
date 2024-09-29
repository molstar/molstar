/**
 * Copyright (c) 2022-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginReactContext, PluginUIComponent } from '../../../mol-plugin-ui/base';
import { Button, ControlGroup, IconButton } from '../../../mol-plugin-ui/controls/common';
import { ArrowDropDownSvg, ArrowRightSvg, CloseSvg, VisibilityOffOutlinedSvg, VisibilityOutlinedSvg, ContentCutSvg, BrushSvg, SearchSvg, TooltipTextSvg, TooltipTextOutlineSvg, PlusBoxSvg, MinusBoxSvg } from '../../../mol-plugin-ui/controls/icons';
import { PluginCommands } from '../../../mol-plugin/commands';
import { State, StateObjectCell, StateSelection, StateTransformer } from '../../../mol-state';
import { ParameterControls, ParameterMappingControl, ParamOnChange, SelectControl } from '../../../mol-plugin-ui/controls/parameters';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { Clip } from '../../../mol-util/clip';
import { StructureRepresentation3D } from '../../../mol-plugin-state/transforms/representation';
import { Color } from '../../../mol-util/color';
import { CombinedColorControl } from '../../../mol-plugin-ui/controls/color';
import { MarkerAction } from '../../../mol-util/marker-action';
import { EveryLoci, Loci } from '../../../mol-model/loci';
import { deepEqual } from '../../../mol-util';
import { ColorValueParam, ColorParams, ColorProps, DimLightness, LightnessParams, LodParams, MesoscaleGroup, MesoscaleGroupProps, OpacityParams, SimpleClipParams, SimpleClipProps, createClipMapping, getClipObjects, getDistinctGroupColors, RootParams, MesoscaleState, getRoots, getAllGroups, getAllLeafGroups, getFilteredEntities, getAllFilteredEntities, getGroups, getEntities, getAllEntities, getEntityLabel, updateColors, getGraphicsModeProps, GraphicsMode, MesoscaleStateParams, setGraphicsCanvas3DProps, PatternParams, expandAllGroups, EmissiveParams, IllustrativeParams, getCellDescription, getEntityDescription, getEveryEntity } from '../data/state';
import React, { useState } from 'react';
import { MesoscaleExplorerState } from '../app';
import { StructureElement } from '../../../mol-model/structure/structure/element';
import { PluginStateObject as PSO } from '../../../mol-plugin-state/objects';
import { Structure } from '../../../mol-model/structure';
import { PluginContext } from '../../../mol-plugin/context';
import { Sphere3D } from '../../../mol-math/geometry';
import { MesoFocusLoci } from '../behavior/camera';
import Markdown from 'react-markdown';
import { combineLatest } from 'rxjs';

function centerLoci(plugin: PluginContext, loci: Loci, durationMs = 250) {
    const { canvas3d } = plugin;
    if (!canvas3d) return;

    const sphere = Loci.getBoundingSphere(loci) || Sphere3D();
    const snapshot = canvas3d.camera.getCenter(sphere.center);
    canvas3d.requestCameraReset({ durationMs, snapshot });
}

export class ModelInfo extends PluginUIComponent<{}, { isDisabled: boolean }> {
    state = {
        isDisabled: false,
    };

    componentDidMount() {
        this.subscribe(this.plugin.state.data.behaviors.isUpdating, v => {
            this.setState({ isDisabled: v });
        });

        this.subscribe(this.plugin.state.events.cell.stateUpdated, e => {
            if (!this.state.isDisabled && MesoscaleState.has(this.plugin) && MesoscaleState.ref(this.plugin) === e.ref) {
                this.forceUpdate();
            }
        });
    }

    get info() {
        if (!MesoscaleState.has(this.plugin)) return;

        const state = MesoscaleState.get(this.plugin);
        if (!state.description && !state.link) return;

        return {
            selectionDescription: state.focusInfo,
            description: state.description,
            link: state.link,
        };
    }

    render() {
        const info = this.info;
        return info && <>
            <div id='modelinfo' className='msp-help-text'>
                <div>{info.description}</div>
                <div><a href={info.link} target='_blank'>Source</a></div>
            </div>
        </>;
    }
}

const SelectionStyleParam = PD.Select('color+outline', PD.objectToOptions({
    'color+outline': 'Color & Outline',
    'color': 'Color',
    'outline': 'Outline'
} as const));
type SelectionStyle = typeof SelectionStyleParam['defaultValue']

export class SelectionInfo extends PluginUIComponent<{}, { isDisabled: boolean }> {
    state = {
        isDisabled: false,
    };

    componentDidMount() {
        this.subscribe(this.plugin.state.data.behaviors.isUpdating, v => {
            this.setState({ isDisabled: v });
        });

        this.subscribe(this.plugin.managers.structure.selection.events.changed, e => {
            if (!this.state.isDisabled) {
                this.forceUpdate();
            }
        });
    }

    get info() {
        const infos: { label: string, key: string, description?: string }[] = [];
        this.plugin.managers.structure.selection.entries.forEach((e, k) => {
            if (StructureElement.Loci.is(e.selection) && !StructureElement.Loci.isEmpty(e.selection)) {
                const cell = this.plugin.helpers.substructureParent.get(e.selection.structure);
                const { entities } = e.selection.structure.model;
                const description = entities.data.pdbx_description.value(0)[0] || 'model';
                infos.push({
                    description: description,
                    label: cell?.obj?.label || 'Unknown',
                    key: k,
                });
            }
        });
        return infos;
    }

    find(label: string) {
        MesoscaleState.set(this.plugin, { filter: `"${label}"` });
        if (label) expandAllGroups(this.plugin);
    };

    remove(key: string) {
        const e = this.plugin.managers.structure.selection.entries.get(key);
        if (!e) return;

        const loci = Structure.toStructureElementLoci(e.selection.structure);
        this.plugin.managers.interactivity.lociSelects.deselect({ loci }, false);
    }

    center(key: string) {
        const e = this.plugin.managers.structure.selection.entries.get(key);
        if (!e) return;

        const loci = Structure.toStructureElementLoci(e.selection.structure);
        centerLoci(this.plugin, loci);
        const cell = this.plugin.helpers.substructureParent.get(loci.structure);
        const d = getCellDescription(cell!); // '### ' + cell?.obj?.label + '\n\n' + cell?.obj?.description;
        MesoscaleState.set(this.plugin, { focusInfo: `${d}` });
    }

    get selection() {
        const info = this.info;
        const help_selection = <><div>Use <i>ctrl+left</i> to select entities, either on the 3D canvas or in the tree below</div><div>Use <i>shift+left</i> to select individual chain on the 3D canvas</div></>;
        if (!info.length) return <>
            <div id='seleinfo' className='msp-help-text'>
                {help_selection}
            </div>
        </>;
        return <>
            <div id='seleinfo'>
                {info.map((entry, index) => {
                    const label = <Button className={`msp-btn-tree-label`} noOverflow disabled={this.state.isDisabled}
                        onClick={() => this.center(entry.key)}
                    >
                        <span title={entry.label}>
                            {entry.label}
                        </span>
                    </Button>;
                    const find = <IconButton svg={SearchSvg} toggleState={false} disabled={this.state.isDisabled} small onClick={() => this.find(entry.label)} />;
                    const remove = <IconButton svg={CloseSvg} toggleState={false} disabled={this.state.isDisabled} onClick={() => this.remove(entry.key)} />;
                    return <>
                        <div key={index} className={`msp-flex-row`} style={{ margin: `1px 5px 1px ${1 * 10 + 5}px` }}>
                            {label}
                            {find}
                            {remove}
                        </div>
                    </>;
                })};
            </div>
        </>;
    }

    get style() {
        const p = this.plugin.canvas3d?.props;
        if (!p) return;

        if (p.renderer.dimStrength === 1 && p.marking.enabled) return 'color+outline';
        if (p.renderer.dimStrength === 1) return 'color';
        if (p.marking.enabled) return 'outline';
    }

    setStyle(value: SelectionStyle) {
        if (value.includes('color') && value.includes('outline')) {
            this.plugin.canvas3d?.setProps({
                renderer: {
                    dimStrength: 1,
                },
                marking: {
                    enabled: true
                }
            });
        } else if (value.includes('color')) {
            this.plugin.canvas3d?.setProps({
                renderer: {
                    dimStrength: 1,
                },
                marking: {
                    enabled: false
                }
            });
        } else if (value.includes('outline')) {
            this.plugin.canvas3d?.setProps({
                renderer: {
                    dimStrength: 0,
                    selectStrength: 0.3,
                },
                marking: {
                    enabled: true
                }
            });
        } else {
            this.plugin.canvas3d?.setProps({
                renderer: {
                    dimStrength: 0,
                    selectStrength: 0,
                },
                marking: {
                    enabled: false
                }
            });
        }

        this.forceUpdate();
    }

    renderStyle() {
        const style = this.style || '';
        return <div id='selestyle' style={{ margin: '5px', marginBottom: '10px' }}>
            <SelectControl name={'Style'} param={SelectionStyleParam} value={style} onChange={(e) => { this.setStyle(e.value); }} />
        </div>;
    }

    render() {
        return <>
            {this.renderStyle()}
            {this.selection}
        </>;
    }
}

export function MesoMarkdownAnchor({ href, children, element }: { href?: string, children?: any, element?: any }) {
    const plugin = React.useContext(PluginReactContext);
    if (!href) return element;
    // Decode the href to handle encoded spaces and other characters
    const decodedHref = href ? decodeURIComponent(href) : '';
    const handleHover = (e: React.MouseEvent<HTMLAnchorElement, MouseEvent>) => {
        e.preventDefault();
        if (decodedHref.startsWith('i')) {
            plugin.canvas3d?.mark({ loci: EveryLoci }, MarkerAction.RemoveHighlight);
            const query_names = decodedHref.substring(1).split(',');
            for (const query_name of query_names) {
                const entities = getEveryEntity(plugin, query_name);
                for (const r of entities) {
                    const repr = r.obj?.data.repr;
                    if (repr) {
                        plugin.canvas3d?.mark({ repr, loci: EveryLoci }, MarkerAction.Highlight);
                    }
                }
            }
        } else if (decodedHref.startsWith('g')) {
            plugin.canvas3d?.mark({ loci: EveryLoci }, MarkerAction.RemoveHighlight);
            const qindex = decodedHref.indexOf('.');
            const query = decodedHref.substring(1, qindex) + ':';
            const query_names = decodedHref.substring(qindex + 1).split(',');
            for (const query_name of query_names) {
                const e = getAllEntities(plugin, query + query_name);
                for (const r of e) {
                    const repr = r.obj?.data.repr;
                    if (repr) {
                        plugin.canvas3d?.mark({ repr, loci: EveryLoci }, MarkerAction.Highlight);
                    }
                }
            }
        }
    };
    const handleLeave = (e: React.MouseEvent<HTMLAnchorElement, MouseEvent>) => {
        // Implement your hover off logic here
        // Example: Perform an action if the href starts with 'h'
        if (decodedHref.startsWith('i') || decodedHref.startsWith('g')) {
            // Example hover off action
            e.preventDefault();
            plugin.canvas3d?.mark({ loci: EveryLoci }, MarkerAction.RemoveHighlight);
        }
    };
    const handleClick = (e: React.MouseEvent<HTMLAnchorElement, MouseEvent>) => {
        e.preventDefault();
        if (href.startsWith('#')) {
            plugin.managers.snapshot.applyKey(decodedHref.substring(1));
        } else if (decodedHref.startsWith('i')) {
            plugin.managers.interactivity.lociSelects.deselectAll();
            plugin.canvas3d?.mark({ loci: EveryLoci }, MarkerAction.RemoveHighlight);
            const query_names = decodedHref.substring(1).split(',');
            for (const query_name of query_names) {
                const entities = getFilteredEntities(plugin, '', query_name);
                for (const r of entities) {
                    const repr = r.obj?.data.repr;
                    if (repr) {
                        plugin.canvas3d?.mark({ repr, loci: EveryLoci }, MarkerAction.Highlight);
                    }
                    const cell = r as StateObjectCell<PSO.Molecule.Structure.Representation3D | PSO.Shape.Representation3D> | undefined;
                    if (!(cell?.obj?.data.sourceData instanceof Structure)) {
                        return;
                    }
                    const loci = Structure.toStructureElementLoci(cell.obj.data.sourceData);
                    plugin.managers.interactivity.lociSelects.toggle({ loci }, false);
                }
            }
        } else if (decodedHref.startsWith('g')) {
            plugin.managers.interactivity.lociSelects.deselectAll();
            plugin.canvas3d?.mark({ loci: EveryLoci }, MarkerAction.RemoveHighlight);
            const qindex = decodedHref.indexOf('.');
            const query = decodedHref.substring(1, qindex) + ':';
            const query_names = decodedHref.substring(qindex + 1).split(',');
            for (const query_name of query_names) {
                const entities = getAllEntities(plugin, query + query_name);
                for (const r of entities) {
                    const repr = r.obj?.data.repr;
                    if (repr) {
                        plugin.canvas3d?.mark({ repr, loci: EveryLoci }, MarkerAction.Highlight);
                    }
                    const cell = r as StateObjectCell<PSO.Molecule.Structure.Representation3D | PSO.Shape.Representation3D> | undefined;
                    if (!(cell?.obj?.data.sourceData instanceof Structure)) return;
                    const loci = Structure.toStructureElementLoci(cell.obj.data.sourceData);
                    plugin.managers.interactivity.lociSelects.toggle({ loci }, false);
                }
            }
        } else {
            // open the link in a new tab
            window.open(decodedHref, '_blank');
        }
    };

    if (decodedHref[0] === '#') {
        return <a href={decodedHref[0]} onMouseOver={handleHover} onClick={handleClick}>{children}</a>;
    }
    if (decodedHref[0] === 'i' || decodedHref[0] === 'g') {
        return <a href={decodedHref[0]} onMouseLeave={handleLeave} onMouseOver={handleHover} onClick={handleClick}>{children}</a>;
    }
    if (decodedHref[0] === 'h') {
        return <a href={decodedHref[0]} onClick={handleClick} rel='noopener noreferrer'>{children}</a>;
    }
    return element;
}

export function MesoViewportSnapshotDescription() {
    let tSize = 14;
    const plugin = React.useContext(PluginReactContext);
    if (MesoscaleState.has(plugin)) {
        const state = MesoscaleState.get(plugin);
        tSize = state.textSizeDescription;
    }
    const [_, setV] = React.useState(0);
    const [isShown, setIsShown] = useState(true);
    const [textSize, setTextSize] = useState(tSize);
    const toggleVisibility = () => {
        setIsShown(!isShown);
    };

    const increaseTextSize = () => {
        setTextSize(prevSize => Math.min(prevSize + 2, 50)); // Increase the text size by 2px, but not above 50px
    };

    const decreaseTextSize = () => {
        setTextSize(prevSize => Math.max(prevSize - 2, 2)); // Decrease the text size by 2px, but not below 2px
    };


    React.useEffect(() => {
        const sub = plugin.managers.snapshot.events.changed.subscribe(() => setV(v => v + 1));
        return () => sub.unsubscribe();
    }, [plugin]);

    const current = plugin.managers.snapshot.state.current;
    if (!current) return null;

    const e = plugin.managers.snapshot.getEntry(current)!;
    if (!e?.description?.trim()) return null;
    if (MesoscaleState.has(plugin)) {
        MesoscaleState.set(plugin, { textSizeDescription: textSize });
    }
    const showInfo = <IconButton svg={isShown ? TooltipTextSvg : TooltipTextOutlineSvg} flex='20px' onClick={toggleVisibility} title={isShown ? 'Hide Description' : 'Show Description'}/>;
    const increasePoliceSize = <IconButton svg={PlusBoxSvg} flex='20px' onClick={increaseTextSize} title='Bigger Text' />;
    const decreasePoliceSize = <IconButton svg={MinusBoxSvg} flex='20px' onClick={decreaseTextSize} title='Smaller Text' />;
    return (
        <>
            <div id='snapinfoctrl' className="msp-state-snapshot-viewport-controls" style={{ marginRight: '30px' }}>
                {showInfo}{increasePoliceSize}{decreasePoliceSize}
            </div>
            <div id='snapinfo' className={`msp-snapshot-description-me ${isShown ? 'shown' : 'hidden'}`} style={{ fontSize: `${textSize}px` }}>
                {<Markdown skipHtml={false} components={{ a: MesoMarkdownAnchor }}>{e.description}</Markdown>}
            </div>
        </>
    );
}

export class FocusInfo extends PluginUIComponent<{}, { isDisabled: boolean }> {
    componentDidMount() {
        this.subscribe(combineLatest([
            this.plugin.state.data.behaviors.isUpdating,
            this.plugin.managers.structure.selection.events.changed
        ]), ([isUpdating]) => {
            if (!isUpdating) this.forceUpdate();
        });
    }

    get info() {
        let focusInfo = '';
        if (MesoscaleState.has(this.plugin)) {
            const state = MesoscaleState.get(this.plugin);
            if (state.focusInfo) focusInfo = state.focusInfo;
        }
        return focusInfo;
    }

    render() {
        const focusInfo = this.info;
        const description = (focusInfo !== '') ? <Markdown skipHtml components={{ a: MesoMarkdownAnchor }}>{focusInfo}</Markdown> : '';
        return <>
            <div id='focusinfo' className='msp-help-text'>
                {description}
            </div>
        </>;
    }
}


export class EntityControls extends PluginUIComponent<{}, { isDisabled: boolean }> {
    filterRef = React.createRef<HTMLInputElement>();
    prevFilter = '';
    filterFocus = false;

    state = {
        isDisabled: false,
    };

    componentDidMount() {
        this.subscribe(this.plugin.state.events.object.created, e => {
            this.forceUpdate();
        });

        this.subscribe(this.plugin.state.events.object.removed, e => {
            this.forceUpdate();
        });

        this.subscribe(this.plugin.state.data.behaviors.isUpdating, v => {
            this.setState({ isDisabled: v });
        });

        this.subscribe(this.plugin.state.events.cell.stateUpdated, e => {
            if (!this.state.isDisabled && this.roots.some(r => e.cell === r) || (MesoscaleState.has(this.plugin) && MesoscaleState.ref(this.plugin) === e.ref)) {
                this.forceUpdate();
            }
        });
    }

    componentDidUpdate(): void {
        const filter = this.filter;
        if (this.filterFocus) {
            this.filterRef.current?.focus();
            this.prevFilter = filter;
        }
    }

    get roots() {
        return getRoots(this.plugin);
    }

    setGroupBy = (value: number) => {
        this.roots.forEach((c, i) => {
            if (c.state.isHidden && value === i || !c.state.isHidden && value !== i) {
                PluginCommands.State.ToggleVisibility(this.plugin, { state: c.parent!, ref: c.transform.ref });
            }
        });
    };

    get groupBy() {
        const roots = this.roots;
        for (let i = 0, il = roots.length; i < il; ++i) {
            if (!roots[i].state.isHidden) return i;
        }
        return 0;
    }

    setFilter = (value: string) => {
        this.filterFocus = true;
        const filter = value.trim().replace(/\s+/gi, ' ');
        MesoscaleState.set(this.plugin, { filter });
        if (filter) expandAllGroups(this.plugin);
    };

    get filter() {
        return MesoscaleState.has(this.plugin) ? MesoscaleState.get(this.plugin).filter : '';
    }

    setGraphics = (graphics: GraphicsMode) => {
        MesoscaleState.set(this.plugin, { graphics });
        (this.plugin.customState as MesoscaleExplorerState).graphicsMode = graphics;

        if (graphics === 'custom') return;

        const update = this.plugin.state.data.build();

        const { lodLevels, approximate, alphaThickness } = getGraphicsModeProps(graphics);

        for (const r of getAllEntities(this.plugin)) {
            update.to(r).update(old => {
                if (old.type) {
                    old.type.params.lodLevels = lodLevels;
                    old.type.params.approximate = approximate;
                    old.type.params.alphaThickness = alphaThickness;
                }
            });
        }

        for (const g of getAllGroups(this.plugin)) {
            update.to(g).update(old => {
                old.lod.lodLevels = lodLevels;
                old.lod.approximate = approximate;
            });
        }

        update.commit();

        setGraphicsCanvas3DProps(this.plugin, graphics);
    };

    get graphics() {
        const customState = this.plugin.customState as MesoscaleExplorerState;
        return MesoscaleState.has(this.plugin) ? MesoscaleState.get(this.plugin).graphics : customState.graphicsMode;
    }

    renderGraphics() {
        const graphics = this.graphics;
        return <div id='graphicsquality' style={{ margin: '5px', marginBottom: '10px' }}>
            <SelectControl name={'Graphics'} param={MesoscaleStateParams.graphics} value={`${graphics}`} onChange={(e) => { this.setGraphics(e.value); }} />
        </div>;
    }

    render() {
        const roots = this.roots;
        if (roots.length === 0 || !MesoscaleState.has(this.plugin)) {
            return <>
                {this.renderGraphics()}
            </>;
        }

        const disabled = this.state.isDisabled;
        const groupBy = this.groupBy;

        const options: [string, string][] = [];
        roots.forEach((c, i) => {
            options.push([`${i}`, c.obj!.label]);
        });
        const groupParam = PD.Select(options[0][0], options);
        const root = roots.length === 1 ? roots[0] : roots[groupBy];

        const filter = this.filter;

        return <>
            {this.renderGraphics()}
            <div id='searchtree' className={`msp-flex-row msp-control-row`} style={{ margin: '5px', marginBottom: '10px' }}>
                <input type='text' ref={this.filterRef}
                    value={filter}
                    placeholder='Search'
                    onChange={e => this.setFilter(e.target.value)}
                    disabled={disabled}
                    onBlur={() => this.filterFocus = false}
                />
                <IconButton svg={CloseSvg} toggleState={false} disabled={disabled} onClick={() => this.setFilter('')} />
            </div>
            {options.length > 1 && <div id='grouptree' style={{ margin: '5px', marginBottom: '10px' }}>
                <SelectControl name={'Group By'} param={groupParam} value={`${groupBy}`} onChange={(e) => { this.setGroupBy(parseInt(e.value)); }} />
            </div>}
            <div id='tree' style={{ position: 'relative', overflowY: 'auto', borderBottom: '1px  solid #000', maxHeight: '600px' }}>
                <GroupNode filter={filter} cell={root} depth={0} />
            </div>
        </>;
    }
}

class Node<P extends {}, S extends { isDisabled: boolean }> extends PluginUIComponent<P & { cell: StateObjectCell, depth: number }, S> {

    is(e: State.ObjectEvent) {
        return e.ref === this.ref && e.state === this.props.cell.parent;
    }

    get ref() {
        return this.props.cell.transform.ref;
    }

    get cell() {
        return this.props.cell;
    }

    get roots() {
        return getRoots(this.plugin);
    }

    componentDidMount() {
        this.subscribe(this.plugin.state.data.behaviors.isUpdating, v => {
            this.setState({ isDisabled: v });
        });

        this.subscribe(this.plugin.state.events.cell.stateUpdated, e => {
            if (!this.state.isDisabled && this.is(e)) {
                this.forceUpdate();
            }
        });
    }
}

export class GroupNode extends Node<{ filter: string }, { isCollapsed: boolean, action?: 'color' | 'clip' | 'root', isDisabled: boolean }> {
    state = {
        isCollapsed: !!this.props.cell.state.isCollapsed,
        action: undefined,
        isDisabled: false,
    };

    toggleExpanded = (e: React.MouseEvent<HTMLElement>) => {
        PluginCommands.State.ToggleExpanded(this.plugin, { state: this.cell.parent!, ref: this.ref });
    };

    toggleColor = (e?: React.MouseEvent<HTMLButtonElement, MouseEvent>) => {
        this.setState({ action: this.state.action === 'color' ? undefined : 'color' });
    };

    toggleClip = () => {
        this.setState({ action: this.state.action === 'clip' ? undefined : 'clip' });
    };

    toggleRoot = () => {
        this.setState({ action: this.state.action === 'root' ? undefined : 'root' });
    };

    showInfo = (e: React.MouseEvent<HTMLElement>) => {
        e.preventDefault();
        const d = getCellDescription(this.cell); // '### ' + this.cell?.obj?.label + '\n\n' + this.cell?.obj?.description;
        MesoscaleState.set(this.plugin, { focusInfo: `${d}` });
    };

    highlight = (e: React.MouseEvent<HTMLElement>) => {
        e.preventDefault();
        this.plugin.canvas3d?.mark({ loci: EveryLoci }, MarkerAction.RemoveHighlight);
        for (const r of this.allFilteredEntities) {
            const repr = r.obj?.data.repr;
            if (repr) {
                this.plugin.canvas3d?.mark({ repr, loci: EveryLoci }, MarkerAction.Highlight);
            }
        }
    };

    clearHighlight = (e: React.MouseEvent<HTMLElement>) => {
        e.preventDefault();
        this.plugin.canvas3d?.mark({ loci: EveryLoci }, MarkerAction.RemoveHighlight);
        e.currentTarget.blur();
    };

    get groups() {
        return getGroups(this.plugin, this.cell.params?.values.tag);
    }

    get allGroups() {
        const allGroups = getAllGroups(this.plugin, this.cell.params?.values.tag);
        allGroups.push(this.cell);
        return allGroups;
    }

    get entities() {
        return getEntities(this.plugin, this.cell.params?.values.tag);
    }

    get filteredEntities() {
        return getFilteredEntities(this.plugin, this.cell.params?.values.tag, this.props.filter);
    }

    get allEntities() {
        return getAllEntities(this.plugin, this.cell.params?.values.tag);
    }

    get allFilteredEntities() {
        return getAllFilteredEntities(this.plugin, this.cell.params?.values.tag, this.props.filter);
    }

    toggleVisible = (e: React.MouseEvent<HTMLElement>) => {
        PluginCommands.State.ToggleVisibility(this.plugin, { state: this.cell.parent!, ref: this.ref });
        const isHidden = this.cell.state.isHidden;

        for (const r of this.allFilteredEntities) {
            this.plugin.state.data.updateCellState(r.transform.ref, { isHidden });
        }

        this.plugin.build().to(this.ref).update(old => {
            old.hidden = isHidden;
        }).commit();
    };

    updateColor = (values: ColorProps) => {
        const update = this.plugin.state.data.build();
        const { value, illustrative, type, lightness, alpha, emissive } = values;

        const entities = this.filteredEntities;

        let groupColors: Color[] = [];

        if (type === 'generate') {
            groupColors = getDistinctGroupColors(entities.length, value, values.variability, values.shift);
        }

        for (let i = 0; i < entities.length; ++i) {
            const c = type === 'generate' ? groupColors[i] : value;
            update.to(entities[i]).update(old => {
                if (old.type) {
                    if (illustrative) {
                        old.colorTheme = { name: 'illustrative', params: { style: { name: 'uniform', params: { value: c, lightness } } } };
                    } else {
                        old.colorTheme = { name: 'uniform', params: { value: c, lightness } };
                    }
                    old.type.params.alpha = alpha;
                    old.type.params.xrayShaded = alpha < 1 ? 'inverted' : false;
                    old.type.params.emissive = emissive;
                } else {
                    old.coloring.params.color = c;
                    old.coloring.params.lightness = lightness;
                    old.alpha = alpha;
                    old.xrayShaded = alpha < 1 ? true : false;
                    old.emissive = emissive;
                }
            });
        }

        update.to(this.ref).update(old => {
            old.color = values;
        });

        for (const r of this.roots) {
            update.to(r).update(old => {
                old.color.type = 'custom';
            });
        }

        update.commit();
    };

    updateRoot = async (values: PD.Values) => {
        await updateColors(this.plugin, values, this.cell.params?.values.tag, this.props.filter);

        const update = this.plugin.state.data.build();

        for (const r of this.roots) {
            if (r !== this.cell) {
                update.to(r).update(old => {
                    old.color.type = 'custom';
                });
                const others = getAllLeafGroups(this.plugin, r.params?.values.tag);
                for (const o of others) {
                    update.to(o).update(old => {
                        old.color.type = 'custom';
                    });
                }
            }
        }

        update.to(this.ref).update(old => {
            old.color = values;
        });

        update.commit();
    };

    updateClip = (values: PD.Values) => {
        const update = this.plugin.state.data.build();
        const clipObjects = getClipObjects(values as SimpleClipProps, this.plugin.canvas3d!.boundingSphere);

        for (const r of this.allFilteredEntities) {
            update.to(r).update(old => {
                if (old.type) {
                    old.type.params.clip.objects = clipObjects;
                } else {
                    old.clip.objects = clipObjects;
                }
            });
        }

        for (const g of this.allGroups) {
            update.to(g).update(old => {
                old.clip = values;
            });
        }

        update.commit();
    };

    updateLod = (values: PD.Values) => {
        MesoscaleState.set(this.plugin, { graphics: 'custom' });
        (this.plugin.customState as MesoscaleExplorerState).graphicsMode = 'custom';

        const update = this.plugin.state.data.build();

        for (const r of this.allFilteredEntities) {
            update.to(r).update(old => {
                if (old.type) {
                    old.type.params.lodLevels = values.lodLevels;
                    old.type.params.cellSize = values.cellSize;
                    old.type.params.batchSize = values.batchSize;
                    old.type.params.approximate = values.approximate;
                }
            });
        }

        for (const g of this.allGroups) {
            update.to(g).update(old => {
                old.lod = values;
            });
        }

        update.commit();
    };

    update = (props: MesoscaleGroupProps) => {
        this.plugin.state.data.build().to(this.ref).update(props);
    };

    renderColor() {
        const color = this.cell.params?.values.color;
        if (this.cell.params?.values.color.type === 'uniform') {
            const style = {
                backgroundColor: Color.toStyle(color.value),
                minWidth: 32,
                width: 32,
                borderRight: `6px solid ${Color.toStyle(Color.lighten(color.value, color.lightness))}`
            };
            return <Button style={style} onClick={this.toggleColor} />;
        } else if (this.cell.params?.values.color.type === 'generate') {
            const style = {
                minWidth: 32,
                width: 32,
                borderRight: `6px solid ${Color.toStyle(Color.lighten(color.value, color.lightness))}`
            };
            return <IconButton style={style} svg={BrushSvg} toggleState={false} small onClick={this.toggleColor} />;
        } else {
            return <IconButton svg={BrushSvg} toggleState={false} small onClick={this.toggleColor} />;
        }
    }

    render() {
        if (this.allFilteredEntities.length === 0) return;

        const state = this.cell.state;
        const disabled = false;
        const groupLabel = this.cell.obj!.label;
        const depth = this.props.depth;
        const colorValue = this.cell.params?.values.color;
        const rootValue = this.cell.params?.values.color;
        const clipValue = this.cell.params?.values.clip;
        const lodValue = this.cell.params?.values.lod;
        const isRoot = this.cell.params?.values.root;

        const groups = this.groups;
        const entities = this.entities;

        const label = <Button className={`msp-btn-tree-label`} noOverflow disabled={disabled}
            onMouseEnter={this.highlight}
            onMouseLeave={this.clearHighlight}
            onClick={this.showInfo}
        >
            <span title={groupLabel}>{groupLabel}</span>
        </Button>;

        const expand = <IconButton svg={state.isCollapsed ? ArrowRightSvg : ArrowDropDownSvg} flex='20px' disabled={disabled} onClick={this.toggleExpanded} transparent className='msp-no-hover-outline' style={{ visibility: groups.length > 0 || entities.length > 0 ? 'visible' : 'hidden' }} />;
        const color = (entities.length > 0 && !isRoot) && this.renderColor();
        const root = (isRoot && this.allGroups.length > 1) && <IconButton svg={BrushSvg} toggleState={false} disabled={disabled} small onClick={this.toggleRoot} />;
        const clip = <IconButton svg={ContentCutSvg} toggleState={false} disabled={disabled} small onClick={this.toggleClip} />;
        const visibility = <IconButton svg={state.isHidden ? VisibilityOffOutlinedSvg : VisibilityOutlinedSvg} toggleState={false} disabled={disabled} small onClick={this.toggleVisible} />;

        return <>
            <div className={`msp-flex-row`} style={{ margin: `1px 5px 1px ${depth * 10 + 5}px` }}>
                {expand}
                {label}
                {root || color}
                {clip}
                {visibility}
            </div>
            {this.state.action === 'color' && <div style={{ marginRight: 5 }} className='msp-accent-offset'>
                <ControlGroup header='Color' initialExpanded={true} hideExpander={true} hideOffset={true} onHeaderClick={this.toggleColor}
                    topRightIcon={CloseSvg} noTopMargin childrenClassName='msp-viewport-controls-panel-controls'>
                    <ParameterControls params={ColorParams} values={colorValue} onChangeValues={this.updateColor} />
                </ControlGroup>
            </div>}
            {this.state.action === 'clip' && <div style={{ marginRight: 5 }} className='msp-accent-offset'>
                <ControlGroup header='Clip' initialExpanded={true} hideExpander={true} hideOffset={true} onHeaderClick={this.toggleClip}
                    topRightIcon={CloseSvg} noTopMargin childrenClassName='msp-viewport-controls-panel-controls'>
                    <ParameterControls params={SimpleClipParams} values={clipValue} onChangeValues={this.updateClip} />
                    <ParameterControls params={LodParams} values={lodValue} onChangeValues={this.updateLod} />
                </ControlGroup>
            </div>}
            {this.state.action === 'root' && <div style={{ marginRight: 5 }} className='msp-accent-offset'>
                <ControlGroup header='Color' initialExpanded={true} hideExpander={true} hideOffset={true} onHeaderClick={this.toggleRoot}
                    topRightIcon={CloseSvg} noTopMargin childrenClassName='msp-viewport-controls-panel-controls'>
                    <ParameterControls params={RootParams} values={rootValue} onChangeValues={this.updateRoot} />
                </ControlGroup>
            </div>}
            {(!state.isCollapsed) && <>
                {groups.map(c => {
                    return <GroupNode filter={this.props.filter} cell={c} depth={depth + 1} key={c.transform.ref} />;
                })}
                {this.filteredEntities.map(c => {
                    return <EntityNode cell={c} depth={depth + 1} key={c.transform.ref} />;
                })}
            </>}
        </>;
    }
}

export class EntityNode extends Node<{}, { action?: 'color' | 'clip', isDisabled: boolean }> {
    state = {
        action: undefined,
        isDisabled: false,
    };

    clipMapping = createClipMapping(this);

    get groups() {
        return this.plugin.state.data.select(StateSelection.Generators.ofTransformer(MesoscaleGroup)
            .filter(c => !!this.cell.transform.tags?.includes(c.params?.values.tag)));
    }

    toggleVisible = (e: React.MouseEvent<HTMLElement>) => {
        e.preventDefault();
        PluginCommands.State.ToggleVisibility(this.plugin, { state: this.props.cell.parent!, ref: this.ref });
        e.currentTarget.blur();
    };

    toggleColor = (e?: React.MouseEvent<HTMLButtonElement, MouseEvent>) => {
        if (e?.ctrlKey) {
            this.updateLightness({ lightness: this.lightnessValue?.lightness ? 0 : DimLightness });
            e.preventDefault();
        } else {
            this.setState({ action: this.state.action === 'color' ? undefined : 'color' });
        }
    };

    toggleClip = () => {
        this.setState({ action: this.state.action === 'clip' ? undefined : 'clip' });
    };

    highlight = (e: React.MouseEvent<HTMLElement>) => {
        e.preventDefault();
        this.plugin.canvas3d?.mark({ loci: EveryLoci }, MarkerAction.RemoveHighlight);
        const repr = this.cell?.obj?.data.repr;
        if (repr) {
            this.plugin.canvas3d?.mark({ repr, loci: EveryLoci }, MarkerAction.Highlight);
        }
        e.currentTarget.blur();
    };

    clearHighlight = (e: React.MouseEvent<HTMLElement>) => {
        e.preventDefault();
        this.plugin.canvas3d?.mark({ loci: EveryLoci }, MarkerAction.RemoveHighlight);
        e.currentTarget.blur();
    };

    toggleSelect = (e: React.MouseEvent<HTMLElement>) => {
        e.preventDefault();
        const cell = this.cell as StateObjectCell<PSO.Molecule.Structure.Representation3D | PSO.Shape.Representation3D> | undefined;
        if (!(cell?.obj?.data.sourceData instanceof Structure)) return;

        const loci = Structure.toStructureElementLoci(cell.obj.data.sourceData);
        this.plugin.managers.interactivity.lociSelects.toggle({ loci }, false);
    };

    center = (e: React.MouseEvent<HTMLElement>) => {
        e.preventDefault();
        const cell = this.cell as StateObjectCell<PSO.Molecule.Structure.Representation3D | PSO.Shape.Representation3D> | undefined;
        if (!(cell?.obj?.data.sourceData instanceof Structure)) return;

        const loci = Structure.toStructureElementLoci(cell.obj.data.sourceData);
        centerLoci(this.plugin, loci);
    };

    handleClick = (e: React.MouseEvent<HTMLElement>) => {
        if (e.ctrlKey) {
            this.toggleSelect(e);
        } else {
            const d = getEntityDescription(this.plugin, this.cell);
            if (this.cell?.obj?.data.sourceData.state.models.length !== 0) {
                const repr = this.cell?.obj?.data.repr;
                if (repr) {
                    // for fiber need to think how to handle.
                    const aloci = repr.getAllLoci()[0];
                    const locis = Loci.normalize(aloci, 'chainInstances') as StructureElement.Loci;
                    const nChain = aloci.structure.state.unitSymmetryGroups.length;
                    let index = MesoscaleState.get(this.plugin).index + 1;
                    if (index * nChain >= locis.elements.length) index = 0;
                    const elems = locis.elements.slice(index * nChain, ((index + 1) * nChain)); // end index is not included
                    const loci = StructureElement.Loci(aloci.structure, elems);
                    const sphere = Loci.getBoundingSphere(loci) || Sphere3D();
                    const state = this.plugin.state.behaviors;
                    const selections = state.select(StateSelection.Generators.ofTransformer(MesoFocusLoci));
                    const params = selections.length === 1 ? selections[0].obj?.data.params : undefined;
                    if (!params.centerOnly) {
                        this.plugin.managers.camera.focusSphere(sphere, params);
                    } else {
                        const snapshot = this.plugin.canvas3d?.camera.getCenter(sphere.center);
                        this.plugin.canvas3d?.requestCameraReset({ durationMs: params.durationMs, snapshot });
                    }
                    MesoscaleState.set(this.plugin, { index: index, focusInfo: `${d}` });
                }
            } else {
                this.center(e);
                MesoscaleState.set(this.plugin, { focusInfo: `${d}` });
            }
        }
    };

    get colorValue(): Color | undefined {
        if (this.cell.transform.params?.colorTheme?.params.value) {
            return this.cell.transform.params?.colorTheme?.params.value;
        } else if (this.cell.transform.params?.colorTheme?.name === 'illustrative') {
            return this.cell.transform.params?.colorTheme?.params.style.params.value;
        } else {
            return this.cell.transform.params?.colorTheme?.params.value ?? this.cell.transform.params?.coloring?.params.color;
        }
    }

    get illustrativeValue(): { illustrative: boolean } | undefined {
        return {
            illustrative: (this.cell.transform.params?.colorTheme?.name === 'illustrative')
        };
    }

    get lightnessValue(): { lightness: number } | undefined {
        if (this.cell.transform.params?.colorTheme?.name === 'illustrative') {
            return {
                lightness: this.cell.transform.params?.colorTheme?.params.style.params.lightness ?? 0
            };
        } else {
            return {
                lightness: this.cell.transform.params?.colorTheme?.params.lightness ?? this.cell.transform.params?.coloring?.params.lightness ?? 0
            };
        }
    }

    get opacityValue(): { alpha: number } | undefined {
        return {
            alpha: this.cell.transform.params?.type?.params.alpha ?? this.cell.transform.params?.alpha ?? 1
        };
    }

    get emissiveValue(): { emissive: number } | undefined {
        return {
            emissive: this.cell.transform.params?.type?.params.emissive ?? this.cell.transform.params?.emissive ?? 0
        };
    }

    get clipValue(): Clip.Props | undefined {
        return this.cell.transform.params.type?.params.clip ?? this.cell.transform.params.clip;
    }

    get lodValue(): PD.Values<typeof LodParams> | undefined {
        const p = this.cell.transform.params?.type?.params;
        if (!p) return;
        return {
            lodLevels: p.lodLevels,
            cellSize: p.cellSize,
            batchSize: p.batchSize,
            approximate: p.approximate,
        };
    }

    get patternValue(): { amplitude: number, frequency: number } | undefined {
        const p = this.cell.transform.params;
        if (p.type) return;
        return {
            amplitude: p.bumpAmplitude,
            frequency: p.bumpFrequency * 10,
        };
    }

    updateColor: ParamOnChange = ({ value }) => {
        const update = this.plugin.state.data.build();
        for (const g of this.groups) {
            update.to(g.transform.ref).update(old => {
                old.color.type = 'custom';
            });
        }
        for (const r of this.roots) {
            update.to(r).update(old => {
                old.color.type = 'custom';
            });
        }
        update.to(this.ref).update(old => {
            if (old.colorTheme) {
                if (old.colorTheme.name === 'illustrative') {
                    old.colorTheme.params.style.params.value = value;
                } else {
                    old.colorTheme.params.value = value;
                }
            } else if (old.coloring) {
                old.coloring.params.color = value;
            }
        });
        update.commit();
    };

    updateIllustrative = (values: PD.Values) => {
        return this.plugin.build().to(this.ref).update(old => {
            if (old.colorTheme) {
                if (old.colorTheme.name !== 'illustrative' && values.illustrative) {
                    old.colorTheme = { name: 'illustrative', params: { style: { name: 'uniform', params: { value: old.colorTheme.params.value, lightness: old.colorTheme.params.lightness } } } };
                } else if (old.colorTheme.name === 'illustrative' && !values.illustrative) {
                    old.colorTheme = { name: 'uniform', params: { value: old.colorTheme.params.style.params.value, lightness: old.colorTheme.params.style.params.lightness } };
                }
            }
        }).commit();
    };

    updateLightness = (values: PD.Values) => {
        return this.plugin.build().to(this.ref).update(old => {
            if (old.colorTheme) {
                if (old.colorTheme.name === 'illustrative') {
                    old.colorTheme.params.style.params.lightness = values.lightness;
                } else {
                    old.colorTheme.params.lightness = values.lightness;
                }
            } else if (old.coloring) {
                old.coloring.params.lightness = values.lightness;
            }
        }).commit();
    };

    updateOpacity = (values: PD.Values) => {
        return this.plugin.build().to(this.ref).update(old => {
            if (old.type) {
                old.type.params.alpha = values.alpha;
                old.type.params.xrayShaded = values.alpha < 1 ? 'inverted' : false;
            } else {
                old.alpha = values.alpha;
                old.xrayShaded = values.alpha < 1 ? true : false;
            }
        }).commit();
    };

    updateEmissive = (values: PD.Values) => {
        return this.plugin.build().to(this.ref).update(old => {
            if (old.type) {
                old.type.params.emissive = values.emissive;
            } else {
                old.emissive = values.emissive;
            }
        }).commit();
    };

    updateClip = (props: Clip.Props) => {
        const params = this.cell.transform.params;
        const clip = params.type ? params.type.params.clip : params.clip;
        if (!PD.areEqual(Clip.Params, clip, props)) {
            this.plugin.build().to(this.ref).update(old => {
                if (old.type) {
                    old.type.params.clip = props;
                } else {
                    old.clip = props;
                }
            }).commit();
        }
    };

    updateLod = (values: PD.Values) => {
        const params = this.cell.transform.params as StateTransformer.Params<StructureRepresentation3D>;
        if (!params.type) return;

        MesoscaleState.set(this.plugin, { graphics: 'custom' });
        (this.plugin.customState as MesoscaleExplorerState).graphicsMode = 'custom';

        if (!deepEqual(params.type.params.lodLevels, values.lodLevels) || params.type.params.cellSize !== values.cellSize || params.type.params.batchSize !== values.batchSize || params.type.params.approximate !== values.approximate) {
            this.plugin.build().to(this.ref).update(old => {
                old.type.params.lodLevels = values.lodLevels;
                old.type.params.cellSize = values.cellSize;
                old.type.params.batchSize = values.batchSize;
                old.type.params.approximate = values.approximate;
            }).commit();
        }
    };

    updatePattern = (values: PD.Values) => {
        return this.plugin.build().to(this.ref).update(old => {
            if (!old.type) {
                old.bumpAmplitude = values.amplitude;
                old.bumpFrequency = values.frequency / 10;
            }
        }).commit();
    };

    render() {
        const cellState = this.cell.state;
        const disabled = this.cell.status !== 'error' && this.cell.status !== 'ok';
        const depth = this.props.depth;
        const colorValue = this.colorValue;
        const lightnessValue = this.lightnessValue;
        const illustrativeValue = this.illustrativeValue;
        const opacityValue = this.opacityValue;
        const emissiveValue = this.emissiveValue;
        const lodValue = this.lodValue;
        const patternValue = this.patternValue;

        const l = getEntityLabel(this.plugin, this.cell);
        const label = <Button className={`msp-btn-tree-label msp-type-class-${this.cell.obj!.type.typeClass}`} noOverflow disabled={disabled}
            onClick={this.handleClick}
            onMouseEnter={this.highlight}
            onMouseLeave={this.clearHighlight}
        >
            <span title={l}>{l}</span>
        </Button>;

        const color = colorValue !== undefined && <Button style={{ backgroundColor: Color.toStyle(colorValue), minWidth: 32, width: 32, borderRight: `6px solid ${Color.toStyle(Color.lighten(colorValue, lightnessValue?.lightness || 0))}` }} onClick={this.toggleColor} />;
        const clip = <IconButton svg={ContentCutSvg} toggleState={false} disabled={disabled} small onClick={this.toggleClip} />;
        const visibility = <IconButton svg={cellState.isHidden ? VisibilityOffOutlinedSvg : VisibilityOutlinedSvg} toggleState={false} disabled={disabled} small onClick={this.toggleVisible} />;

        return <>
            <div className={`msp-flex-row`} style={{ margin: `1px 5px 1px ${depth * 10 + 5}px` }}>
                {label}
                {color}
                {clip}
                {visibility}
            </div>
            {this.state.action === 'color' && colorValue !== void 0 && <div style={{ marginRight: 5 }} className='msp-accent-offset'>
                <ControlGroup header='Color' initialExpanded={true} hideExpander={true} hideOffset={true} onHeaderClick={this.toggleColor}
                    topRightIcon={CloseSvg} noTopMargin childrenClassName='msp-viewport-controls-panel-controls'>
                    <CombinedColorControl param={ColorValueParam} value={colorValue ?? Color(0xFFFFFF)} onChange={this.updateColor} name='color' hideNameRow />
                    <ParameterControls params={IllustrativeParams} values={illustrativeValue} onChangeValues={this.updateIllustrative} />
                    <ParameterControls params={LightnessParams} values={lightnessValue} onChangeValues={this.updateLightness} />
                    <ParameterControls params={OpacityParams} values={opacityValue} onChangeValues={this.updateOpacity} />
                    <ParameterControls params={EmissiveParams} values={emissiveValue} onChangeValues={this.updateEmissive} />
                    {patternValue && <ParameterControls params={PatternParams} values={patternValue} onChangeValues={this.updatePattern} />}
                </ControlGroup>
            </div>}
            {this.state.action === 'clip' && <div style={{ marginRight: 5 }} className='msp-accent-offset'>
                <ControlGroup header='Clip' initialExpanded={true} hideExpander={true} hideOffset={true} onHeaderClick={this.toggleClip}
                    topRightIcon={CloseSvg} noTopMargin childrenClassName='msp-viewport-controls-panel-controls'>
                    <ParameterMappingControl mapping={this.clipMapping} />
                    {lodValue && <ParameterControls params={LodParams} values={lodValue} onChangeValues={this.updateLod} />}
                </ControlGroup>
            </div>}
        </>;
    }
}


