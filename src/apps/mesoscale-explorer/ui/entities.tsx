/**
 * Copyright (c) 2022-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginUIComponent } from '../../../mol-plugin-ui/base';
import { Button, ControlGroup, IconButton } from '../../../mol-plugin-ui/controls/common';
import { ArrowDropDownSvg, ArrowRightSvg, CloseSvg, VisibilityOffOutlinedSvg, VisibilityOutlinedSvg, ContentCutSvg, BrushSvg } from '../../../mol-plugin-ui/controls/icons';
import { PluginCommands } from '../../../mol-plugin/commands';
import { State, StateObjectCell, StateSelection, StateTransformer } from '../../../mol-state';
import { debounceTime, filter } from 'rxjs/operators';
import { escapeRegExp } from '../../../mol-util/string';
import { StructureElement, StructureProperties } from '../../../mol-model/structure';
import { ParameterControls, ParameterMappingControl, ParamOnChange, SelectControl } from '../../../mol-plugin-ui/controls/parameters';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { Clip } from '../../../mol-util/clip';
import { StructureRepresentation3D } from '../../../mol-plugin-state/transforms/representation';
import { Color } from '../../../mol-util/color';
import { CombinedColorControl } from '../../../mol-plugin-ui/controls/color';
import { MarkerAction } from '../../../mol-util/marker-action';
import { EveryLoci } from '../../../mol-model/loci';
import { deepEqual } from '../../../mol-util';
import { ColorValueParam, ColorParams, ColorProps, DimLightness, LightnessParams, LodParams, MesoscaleGroup, MesoscaleGroupObject, MesoscaleGroupProps, OpacityParams, SimpleClipParams, SimpleClipProps, createClipMapping, getClipProps, getDistinctBaseColors, getDistinctGroupColors, RootParams } from '../data/state';
import { PluginContext } from '../../../mol-plugin/context';

export class EntityControls extends PluginUIComponent<{}, { filter: string, isDisabled: boolean, groupBy: number }> {
    state = {
        filter: '',
        isDisabled: false,
        groupBy: 0,
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
    }

    get roots() {
        return getRoots(this.plugin);
    }

    render() {
        const disabled = this.state.isDisabled;
        const roots = this.roots;
        if (roots.length === 0) return;

        const options = roots.map((c, i) => [`${i}`, c.obj!.label] as [string, string]);
        const groupParam = PD.Select(options[0][0], options);
        const root = roots.length === 1 ? roots[0] : roots[this.state.groupBy];

        return <>
            <div className={`msp-flex-row msp-control-row`} style={{ margin: '5px', marginBottom: '10px' }}>
                <input type='text'
                    value={this.state.filter}
                    placeholder='Search'
                    onChange={e => this.setState({ filter: e.target.value.trim().replace(/\s+/gi, ' ') })}
                    disabled={disabled}
                />
                <IconButton svg={CloseSvg} toggleState={false} disabled={disabled} onClick={() => this.setState({ filter: '' })} />
            </div>
            <div style={{ margin: '5px', marginBottom: '10px' }}>
                <SelectControl name={'Group By'} param={groupParam} value={`${this.state.groupBy}`} onChange={(e) => { this.setState({ groupBy: parseInt(e.value) }); }} />
            </div>
            <GroupNode filter={this.state.filter} cell={root} depth={0} />
        </>;
    }
}

class Node<P extends {} = {}, S = {}, SS = {}> extends PluginUIComponent<P & { cell: StateObjectCell, depth: number }, S, SS> {
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
        this.subscribe(this.plugin.state.events.cell.stateUpdated.pipe(filter(e => this.is(e)), debounceTime(33)), e => {
            this.forceUpdate();
        });
    }
}

function getRoots(plugin: PluginContext) {
    return plugin.state.data.select(StateSelection.Generators.rootsOfType(MesoscaleGroupObject));
}

function getGroups(plugin: PluginContext, tag: string) {
    return plugin.state.data.select(StateSelection.Generators.ofTransformer(MesoscaleGroup).withTag(tag));
}

function _getAllGroups(plugin: PluginContext, tag: string, list: StateObjectCell[]) {
    const groups = getGroups(plugin, tag);
    list.push(...groups);
    for (const g of groups) {
        _getAllGroups(plugin, g.params?.values.tag, list);
    }
    return list;
}

function getAllGroups(plugin: PluginContext, tag: string) {
    return _getAllGroups(plugin, tag, []);
}

function getAllLeafGroups(plugin: PluginContext, tag: string) {
    const allGroups = getAllGroups(plugin, tag);
    allGroups.sort((a, b) => a.params?.values.index - b.params?.values.index);
    return allGroups.filter(g => {
        return plugin.state.data.select(StateSelection.Generators.ofTransformer(StructureRepresentation3D).withTag(g.params?.values.tag)).length > 0;
    });
}

function getEntities(plugin: PluginContext, tag: string) {
    return plugin.state.data.select(StateSelection.Generators.ofTransformer(StructureRepresentation3D).withTag(tag)).filter(c => c.obj!.data.sourceData.elementCount > 0);
}

function _getAllEntities(plugin: PluginContext, tag: string, list: StateObjectCell[]) {
    list.push(...getEntities(plugin, tag));
    for (const g of getGroups(plugin, tag)) {
        _getAllEntities(plugin, g.params?.values.tag, list);
    }
    return list;
}

function getAllEntities(plugin: PluginContext, tag: string) {
    return _getAllEntities(plugin, tag, []);
}

function getEntityLabel(cell: StateObjectCell) {
    try {
        const s = cell!.obj!.data.sourceData;
        const l = StructureElement.Location.create(s, s.units[0], s.units[0].elements[0]);
        const d = StructureProperties.entity.pdbx_description(l)[0] || 'model';
        return d.split('.').at(-1) || 'Entity';
    } catch (e) {
        console.log(e, cell!.obj!.data);
    }
    return 'Entity';
}

export class GroupNode extends Node<{ filter: string }, { isCollapsed: boolean, action?: 'color' | 'clip' | 'root' }> {
    state = {
        isCollapsed: !!this.props.cell.state.isCollapsed,
        action: undefined,
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

    highlight = (e: React.MouseEvent<HTMLElement>) => {
        this.plugin.canvas3d?.setProps({
            renderer: { dimStrength: e?.shiftKey ? 1.0 : 0.0 },
            marking: { enabled: !e?.shiftKey },
        });
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
        const reFilter = new RegExp(escapeRegExp(this.props.filter), 'gi');
        return this.entities.filter(c => getEntityLabel(c).match(reFilter) !== null);
    }

    get allEntities() {
        return getAllEntities(this.plugin, this.cell.params?.values.tag);
    }

    get allFilteredEntities() {
        const reFilter = new RegExp(escapeRegExp(this.props.filter), 'gi');
        return this.allEntities.filter(c => getEntityLabel(c).match(reFilter) !== null);
    }

    toggleVisible = (e: React.MouseEvent<HTMLElement>) => {
        PluginCommands.State.ToggleVisibility(this.plugin, { state: this.cell.parent!, ref: this.ref });
        const isHidden = this.cell.state.isHidden;

        for (const r of this.allEntities) {
            this.plugin.state.data.updateCellState(r.transform.ref, { isHidden });
        }

        this.plugin.build().to(this.ref).update(old => {
            old.hidden = isHidden;
        }).commit();
    };

    updateColor = (values: ColorProps) => {
        const update = this.plugin.state.data.build();
        const { value, type, lightness, alpha } = values;

        if (type === 'custom') {
            update.to(this.ref).update(old => {
                old.color.type = type;
            });
        } else {
            const entities = this.entities;

            let groupColors: Color[] = [];

            if (type === 'generate') {
                groupColors = getDistinctGroupColors(entities.length, value);
            }

            for (let i = 0; i < entities.length; ++i) {
                const c = type === 'generate' ? groupColors[i] : value;
                update.to(entities[i]).update(old => {
                    old.colorTheme.params.value = c;
                    old.colorTheme.params.lightness = lightness;
                    old.type.params.alpha = alpha;
                    old.type.params.xrayShaded = alpha < 1;
                });
            }

            update.to(this.ref).update(old => {
                old.color = values;
            });
        }

        for (const r of this.roots) {
            update.to(r).update(old => {
                old.color.type = 'custom';
            });
        }

        update.commit();
    };

    updateRoot = (values: PD.Values) => {
        const update = this.plugin.state.data.build();

        const { type } = values;
        const groups = getAllLeafGroups(this.plugin, this.cell.params?.values.tag);
        const baseColors = getDistinctBaseColors(groups.length);

        for (let i = 0; i < groups.length; ++i) {
            const g = groups[i];
            const entities = getEntities(this.plugin, g.params?.values.tag);
            let groupColors: Color[] = [];

            if (type === 'generate') {
                groupColors = getDistinctGroupColors(entities.length, baseColors[i]);
            }

            for (let j = 0; j < entities.length; ++j) {
                const c = type === 'generate' ? groupColors[j] : baseColors[i];
                update.to(entities[j]).update(old => {
                    old.colorTheme.params.value = c;
                });
            }

            update.to(g.transform.ref).update(old => {
                old.color.type = type;
                old.color.value = baseColors[i];
            });
        }

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
            old.color.type = type;
        });

        update.commit();
    };

    updateClip = (values: PD.Values) => {
        const update = this.plugin.state.data.build();
        const clip = getClipProps(values as SimpleClipProps, this.plugin.canvas3d!.boundingSphere);

        for (const r of this.allEntities) {
            update.to(r).update(old => {
                old.type.params.clip = clip;
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
        const update = this.plugin.state.data.build();

        for (const r of this.allEntities) {
            update.to(r).update(old => {
                old.type.params.lodLevels = values.lodLevels;
                old.type.params.cellSize = values.cellSize;
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
        const rootValue = { type: this.cell.params?.values.color.type };
        const clipValue = this.cell.params?.values.clip;
        const lodValue = this.cell.params?.values.lod;
        const isRoot = this.cell.params?.values.root;

        const groups = this.groups;
        const entities = this.entities;

        const label = <Button className={`msp-btn-tree-label`} noOverflow disabled={disabled}>
            <span>{groupLabel}</span>
        </Button>;

        const expand = <IconButton svg={state.isCollapsed ? ArrowRightSvg : ArrowDropDownSvg} flex='20px' disabled={disabled} onClick={this.toggleExpanded} transparent className='msp-no-hover-outline' style={{ visibility: groups.length > 0 || entities.length > 0 ? 'visible' : 'hidden' }} />;
        const color = (entities.length > 0 && !isRoot) && this.renderColor();
        const root = (isRoot && this.allGroups.length > 1) && <IconButton svg={BrushSvg} toggleState={false} disabled={disabled} small onClick={this.toggleRoot} />;
        const clip = <IconButton svg={ContentCutSvg} toggleState={false} disabled={disabled} small onClick={this.toggleClip} />;
        const visibility = <IconButton svg={state.isHidden ? VisibilityOffOutlinedSvg : VisibilityOutlinedSvg} toggleState={false} disabled={disabled} small onClick={this.toggleVisible} />;

        return <>
            <div className={`msp-flex-row`} style={{ margin: `1px 5px 1px ${depth * 10 + 5}px` }} onMouseEnter={this.highlight} onMouseLeave={this.clearHighlight}>
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
            {(!state.isCollapsed || this.props.filter) && <>
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

export class EntityNode extends Node<{}, { action?: 'color' | 'clip' }> {
    state = {
        action: undefined,
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
        this.plugin.canvas3d?.setProps({
            renderer: { dimStrength: e?.shiftKey ? 1.0 : 0.0 },
            marking: { enabled: !e?.shiftKey },
        });
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

    get colorValue(): Color | undefined {
        const hasValue = this.cell.transform.params?.colorTheme.params.value !== undefined;
        if (!hasValue) return;
        return this.cell.transform.params?.colorTheme.params.value;
    }

    get lightnessValue(): { lightness: number } | undefined {
        const hasLightness = this.cell.transform.params?.colorTheme.params.value !== undefined;
        if (!hasLightness) return;
        return {
            lightness: this.cell.transform.params?.colorTheme.params.lightness
        };
    }

    get opacityValue(): { alpha: number } | undefined {
        const hasOpacity = this.cell.transform.params?.type.params.alpha !== undefined;
        if (!hasOpacity) return;
        return {
            alpha: this.cell.transform.params?.type.params.alpha
        };
    }

    get clipValue(): Clip.Props | undefined {
        return this.cell.transform.params.type.params.clip;
    }

    get lodValue(): PD.Values<typeof LodParams> | undefined {
        const hasLod = this.cell.transform.params?.type.params.lodLevels !== undefined && this.cell.transform.params?.type.params.cellSize !== undefined;
        if (!hasLod) return;
        return {
            lodLevels: this.cell.transform.params.type.params.lodLevels,
            cellSize: this.cell.transform.params.type.params.cellSize,
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
            old.colorTheme.params.value = value;
        });
        update.commit();
    };

    updateLightness = (values: PD.Values) => {
        return this.plugin.build().to(this.ref).update(old => {
            old.colorTheme.params.lightness = values.lightness;
        }).commit();
    };

    updateOpacity = (values: PD.Values) => {
        return this.plugin.build().to(this.ref).update(old => {
            old.type.params.alpha = values.alpha;
            old.type.params.xrayShaded = values.alpha < 1;
        }).commit();
    };

    updateClip = (props: Clip.Props) => {
        const params = this.cell.transform.params as StateTransformer.Params<StructureRepresentation3D>;
        if (!PD.areEqual(Clip.Params, params.type.params.clip, props)) {
            this.plugin.build().to(this.ref).update(old => {
                old.type.params.clip = props;
            }).commit();
        }
    };

    updateLod = (values: PD.Values) => {
        const t = this.cell?.transform;
        if (!t) return;

        const params = t.params as StateTransformer.Params<StructureRepresentation3D>;
        if (!deepEqual(params.type.params.lodLevels, values.lodLevels) || params.type.params.cellSize !== values.cellSize) {
            this.plugin.build().to(t.ref).update(old => {
                old.type.params.lodLevels = values.lodLevels;
                old.type.params.cellSize = values.cellSize;
            }).commit();
        }
    };

    render() {
        const cellState = this.cell.state;
        const disabled = this.cell.status !== 'error' && this.cell.status !== 'ok';
        const depth = this.props.depth;
        const colorValue = this.colorValue;
        const lightnessValue = this.lightnessValue;
        const opacityValue = this.opacityValue;
        const lodValue = this.lodValue;

        const label = <Button className={`msp-btn-tree-label msp-type-class-${this.cell.obj!.type.typeClass}`} noOverflow disabled={disabled}>
            <span>{getEntityLabel(this.cell)}</span>
        </Button>;

        const color = colorValue !== undefined && <Button style={{ backgroundColor: Color.toStyle(colorValue), minWidth: 32, width: 32, borderRight: `6px solid ${Color.toStyle(Color.lighten(colorValue, lightnessValue?.lightness || 0))}` }} onClick={this.toggleColor} />;
        const clip = <IconButton svg={ContentCutSvg} toggleState={false} disabled={disabled} small onClick={this.toggleClip} />;
        const visibility = <IconButton svg={cellState.isHidden ? VisibilityOffOutlinedSvg : VisibilityOutlinedSvg} toggleState={false} disabled={disabled} small onClick={this.toggleVisible} />;

        return <>
            <div className={`msp-flex-row`} style={{ margin: `1px 5px 1px ${depth * 10 + 5}px` }} onMouseEnter={this.highlight} onMouseLeave={this.clearHighlight}>
                {label}
                {color}
                {clip}
                {visibility}
            </div>
            {this.state.action === 'color' && colorValue !== void 0 && <div style={{ marginRight: 5 }} className='msp-accent-offset'>
                <ControlGroup header='Color' initialExpanded={true} hideExpander={true} hideOffset={true} onHeaderClick={this.toggleColor}
                    topRightIcon={CloseSvg} noTopMargin childrenClassName='msp-viewport-controls-panel-controls'>
                    <CombinedColorControl param={ColorValueParam} value={colorValue ?? Color(0xFFFFFF)} onChange={this.updateColor} name='color' hideNameRow />
                    <ParameterControls params={LightnessParams} values={lightnessValue} onChangeValues={this.updateLightness} />
                    <ParameterControls params={OpacityParams} values={opacityValue} onChangeValues={this.updateOpacity} />
                </ControlGroup>
            </div>}
            {this.state.action === 'clip' && <div style={{ marginRight: 5 }} className='msp-accent-offset'>
                <ControlGroup header='Clip' initialExpanded={true} hideExpander={true} hideOffset={true} onHeaderClick={this.toggleClip}
                    topRightIcon={CloseSvg} noTopMargin childrenClassName='msp-viewport-controls-panel-controls'>
                    <ParameterMappingControl mapping={this.clipMapping} />
                    <ParameterControls params={LodParams} values={lodValue} onChangeValues={this.updateLod} />
                </ControlGroup>
            </div>}
        </>;
    }
}

