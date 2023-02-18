/**
 * Copyright (c) 2022-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginStateObject as PSO } from '../../../mol-plugin-state/objects';
import { PluginUIComponent } from '../../../mol-plugin-ui/base';
import { Button, ControlGroup, IconButton } from '../../../mol-plugin-ui/controls/common';
import { ArrowDropDownSvg, ArrowRightSvg, CloseSvg, VisibilityOffOutlinedSvg, VisibilityOutlinedSvg, ContentCutSvg } from '../../../mol-plugin-ui/controls/icons';
import { PluginCommands } from '../../../mol-plugin/commands';
import { State, StateObjectCell, StateSelection, StateTransformer } from '../../../mol-state';
import { debounceTime, filter } from 'rxjs/operators';
import { escapeRegExp, stringToWords } from '../../../mol-util/string';
import { StructureElement, StructureProperties } from '../../../mol-model/structure';
import { ParameterControls, ParameterMappingControl, ParamOnChange } from '../../../mol-plugin-ui/controls/parameters';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { Clip } from '../../../mol-util/clip';
import { StructureRepresentation3D } from '../../../mol-plugin-state/transforms/representation';
import { Color } from '../../../mol-util/color';
import { CombinedColorControl } from '../../../mol-plugin-ui/controls/color';
import { MarkerAction } from '../../../mol-util/marker-action';
import { EveryLoci } from '../../../mol-model/loci';
import { PluginContext } from '../../../mol-plugin/context';
import { Hcl } from '../../../mol-util/color/spaces/hcl';
import { distinctColors } from '../../../mol-util/color/distinct';
import { ParamMapping } from '../../../mol-util/param-mapping';
import { PluginUIContext } from '../../../mol-plugin-ui/context';
import { Vec3 } from '../../../mol-math/linear-algebra';
import { Spheres } from '../../../mol-geo/geometry/spheres/spheres';
import { deepEqual } from '../../../mol-util';

export class EntityControls extends PluginUIComponent<{}, { filter: string, isDisabled: boolean }> {
    state = {
        filter: '',
        isDisabled: false,
    };

    componentDidMount() {
        this.subscribe(this.plugin.state.events.object.created, e => {
            if (PSO.Molecule.Structure.is(e.obj)) this.forceUpdate();
        });

        this.subscribe(this.plugin.state.events.object.removed, e => {
            if (PSO.Molecule.Structure.is(e.obj)) this.forceUpdate();
        });

        this.subscribe(this.plugin.state.data.behaviors.isUpdating, v => {
            this.setState({ isDisabled: v });
        });
    }

    render() {
        const disabled = this.state.isDisabled;
        const roots = this.plugin.state.data.select(StateSelection.Generators.rootsOfType(PSO.Group).withTag('Entity'));

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
            {roots.map(s => {
                return <GroupNode filter={this.state.filter} cell={s} depth={0} key={s.transform.ref} />;
            })}
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

    componentDidMount() {
        this.subscribe(this.plugin.state.events.cell.stateUpdated.pipe(filter(e => this.is(e)), debounceTime(33)), e => {
            this.forceUpdate();
        });
    }

    toggleVisible = (e: React.MouseEvent<HTMLElement>) => {
        e.preventDefault();
        PluginCommands.State.ToggleVisibility(this.plugin, { state: this.props.cell.parent!, ref: this.ref });
        e.currentTarget.blur();
    };
}

function getGroups(plugin: PluginContext, cell: StateObjectCell, ref: string) {
    const children = cell.parent!.tree.children.get(ref);
    return children.toArray().filter(r => {
        const c = plugin.state.data.cells.get(r);
        if (c?.obj?.type !== PSO.Group.type) return false;
        if (!c.transform.tags?.includes('Entity')) return false;
        return true;
    }).map(r => plugin.state.data.cells.get(r)!);
}

function getEntities(plugin: PluginContext, cell: StateObjectCell, ref: string) {
    const children = cell.parent!.tree.children.get(ref);
    return children.toArray().filter(r => {
        const c = plugin.state.data.cells.get(r);
        if (c?.obj?.type !== PSO.Molecule.Structure.type) return false;
        if (c.obj.data.elementCount === 0) return false;
        if (!c.transform.tags?.includes('Entity')) return false;
        return true;
    }).map(r => plugin.state.data.cells.get(r)!);
}

function getAllEntities(plugin: PluginContext, cell: StateObjectCell, ref: string, list: StateObjectCell[]) {
    list.push(...getEntities(plugin, cell, ref));
    for (const g of getGroups(plugin, cell, ref)) {
        getAllEntities(plugin, cell, g.transform.ref, list);
    }
    return list;
}

const ColorParam = PD.Color(Color(0xFFFFFF));

const LightnessParams = {
    lightness: PD.Numeric(0, { min: -6, max: 6, step: 0.1 }),
};
const DimLightness = 6;

const OpacityParams = {
    alpha: PD.Numeric(1, { min: 0, max: 1, step: 0.01 }),
};

const LodParams = {
    lodLevels: Spheres.Params.lodLevels,
    cellSize: Spheres.Params.cellSize,
};

const SimpleClipParams = {
    type: PD.Select('plane', PD.objectToOptions(Clip.Type, t => stringToWords(t))),
    invert: PD.Boolean(false, { hideIf: g => g.type === 'none' }),
    position: PD.Group({
        x: PD.Numeric(0, { min: -100, max: 100, step: 1 }, { immediateUpdate: true }),
        y: PD.Numeric(0, { min: -100, max: 100, step: 1 }, { immediateUpdate: true }),
        z: PD.Numeric(0, { min: -100, max: 100, step: 1 }, { immediateUpdate: true }),
    }, { hideIf: g => g.type === 'none', isExpanded: true }),
    rotation: PD.Group({
        axis: PD.Vec3(Vec3.create(1, 0, 0)),
        angle: PD.Numeric(0, { min: -180, max: 180, step: 1 }, { immediateUpdate: true }),
    }, { hideIf: g => g.type === 'none', isExpanded: true }),
    scale: PD.Group({
        x: PD.Numeric(100, { min: 0, max: 100, step: 1 }, { immediateUpdate: true }),
        y: PD.Numeric(100, { min: 0, max: 100, step: 1 }, { immediateUpdate: true }),
        z: PD.Numeric(100, { min: 0, max: 100, step: 1 }, { immediateUpdate: true }),
    }, { hideIf: g => g.type === 'none', isExpanded: true }),
};
type SimpleClipParams = typeof SimpleClipParams

function createClipMapping(node: GroupNode | EntityNode) {
    return ParamMapping({
        params: SimpleClipParams,
        target: (ctx: PluginUIContext) => {
            return node.clipValue;
        }
    })({
        values(props, ctx) {
            if (!props || props.objects.length === 0) {
                return {
                    type: 'none',
                    invert: false,
                    position: { x: 0, y: 0, z: 0 },
                    rotation: { axis: Vec3.create(1, 0, 0), angle: 0 },
                    scale: { x: 100, y: 100, z: 100 },
                };
            }

            const { center, radius } = node.plugin.canvas3d!.boundingSphere;
            const { invert, position, scale, rotation, type } = props.objects[0];

            const p = Vec3.clone(position);
            Vec3.sub(p, p, center);
            Vec3.scale(p, p, 100 / radius);
            Vec3.round(p, p);

            const s = Vec3.clone(scale);
            Vec3.scale(s, s, 100 / radius / 2);
            Vec3.round(s, s);

            return {
                type,
                invert,
                position: { x: p[0], y: p[1], z: p[2] },
                rotation,
                scale: { x: s[0], y: s[1], z: s[2] },
            };
        },
        update: (s, props) => {
            if (!props) return;

            const { center, radius } = node.plugin.canvas3d!.boundingSphere;

            const position = Vec3.clone(center);
            Vec3.add(position, position, Vec3.create(
                radius * s.position.x / 100,
                radius * s.position.y / 100,
                radius * s.position.z / 100
            ));

            const scale = Vec3.create(s.scale.x, s.scale.y, s.scale.z);
            Vec3.scale(scale, scale, 2 * radius / 100);

            props.variant = 'instance';
            props.objects.length = 1;
            props.objects[0] = {
                type: s.type,
                invert: s.invert,
                position,
                scale,
                rotation: s.rotation
            };
        },
        apply: async (props, ctx) => {
            if (props) await node.updateClip(props);
        }
    });
}

export class GroupNode extends Node<{ filter: string }, { isCollapsed: boolean, action?: 'color' | 'clip' }> {
    state = {
        isCollapsed: !!this.props.cell.state.isCollapsed,
        action: undefined,
    };

    clipMapping = createClipMapping(this);

    toggleExpanded = (e: React.MouseEvent<HTMLElement>) => {
        e.preventDefault();
        PluginCommands.State.ToggleExpanded(this.plugin, { state: this.props.cell.parent!, ref: this.ref });
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
        for (const e of this.allEntities) {
            const repr = this.plugin.state.data.select(StateSelection.Generators.ofTransformer(StructureRepresentation3D, e.transform.ref))[0]?.obj?.data.repr;
            if (repr) {
                this.plugin.canvas3d?.mark({ repr, loci: EveryLoci }, MarkerAction.Highlight);
            }
        }
        e.currentTarget.blur();
    };

    clearHighlight = (e: React.MouseEvent<HTMLElement>) => {
        e.preventDefault();
        this.plugin.canvas3d?.mark({ loci: EveryLoci }, MarkerAction.RemoveHighlight);
        e.currentTarget.blur();
    };

    get groups() {
        return getGroups(this.plugin, this.props.cell, this.ref);
    }

    get entities() {
        const reFilter = new RegExp(escapeRegExp(this.props.filter), 'gi');
        return getEntities(this.plugin, this.props.cell, this.ref)
            .filter(c => c.obj!.label.match(reFilter) !== null);
    }

    get hasEntities() {
        return getEntities(this.plugin, this.props.cell, this.ref).length !== 0;
    }

    get allEntities() {
        return getAllEntities(this.plugin, this.props.cell, this.ref, []);
    }

    get pivotEntity(): StateObjectCell | undefined {
        const entities = getEntities(this.plugin, this.props.cell, this.ref);
        if (entities.length) {
            return entities[Math.floor(entities.length / 2)];
        } else {
            return getAllEntities(this.plugin, this.props.cell, this.ref, [])[0];
        }
    }

    get groupLabel() {
        return this.props.cell.obj!.label;
    }

    get pivotRepr(): StateObjectCell | undefined {
        const entity = this.pivotEntity;
        return entity && this.plugin.state.data.select(StateSelection.Generators.ofTransformer(StructureRepresentation3D, entity.transform.ref))[0];
    }

    get colorValue(): Color | undefined {
        const repr = this.pivotRepr;
        const hasValue = repr?.transform.params?.colorTheme.params.value !== undefined;
        if (!hasValue) return;
        return repr.transform.params?.colorTheme.params.value;
    }

    get lightnessValue(): { lightness: number } | undefined {
        const repr = this.pivotRepr;
        const hasLightness = repr?.transform.params?.colorTheme.params.lightness !== undefined;
        if (!hasLightness) return;
        return {
            lightness: repr.transform.params?.colorTheme.params.lightness
        };
    }

    get opacityValue(): { alpha: number } | undefined {
        const repr = this.pivotRepr;
        const hasOpacity = repr?.transform.params?.type.params.alpha !== undefined;
        if (!hasOpacity) return;
        return {
            alpha: repr.transform.params?.type.params.alpha
        };
    }

    get clipValue(): Clip.Props | undefined {
        return this.pivotRepr?.transform.params.type.params.clip;
    }

    get lodValue(): PD.Values<typeof LodParams> | undefined {
        const repr = this.pivotRepr;
        const hasLod = repr?.transform.params?.type.params.lodLevels !== undefined && repr?.transform.params?.type.params.cellSize !== undefined;
        if (!hasLod) return;
        return {
            lodLevels: repr.transform.params.type.params.lodLevels,
            cellSize: repr.transform.params.type.params.cellSize,
        };
    }

    updateColor: ParamOnChange = ({ value }) => {
        const update = this.plugin.state.data.build();
        const reprCells = this.plugin.state.data.select(StateSelection.Generators.ofTransformer(StructureRepresentation3D, this.ref));

        const hcl = Hcl.fromColor(Hcl(), value);
        const hue = value === 0
            ? [1, 360] as [number, number]
            : [Math.max(1, hcl[0] - 35), Math.min(360, hcl[0] + 35)] as [number, number];
        const groupColors = distinctColors(reprCells.length, {
            hue,
            chroma: [30, 80],
            luminance: [15, 85],
            clusteringStepCount: 50,
            minSampleCount: 800,
        });
        if (value !== 0) {
            groupColors[Math.floor(groupColors.length / 2)] = value;
        }

        for (let i = 0; i < reprCells.length; ++i) {
            const cell = reprCells[i];
            const params = cell.transform.params as StateTransformer.Params<StructureRepresentation3D>;
            if (params.colorTheme.params.value !== value) {
                update.to(cell).update(old => {
                    old.colorTheme.params.value = groupColors[i];
                });
            }
        }
        update.commit();
    };

    updateLightness = (values: PD.Values) => {
        const update = this.plugin.state.data.build();
        const reprCells = this.plugin.state.data.select(StateSelection.Generators.ofTransformer(StructureRepresentation3D, this.ref));

        for (let i = 0; i < reprCells.length; ++i) {
            const cell = reprCells[i];
            const params = cell.transform.params as StateTransformer.Params<StructureRepresentation3D>;
            if (params.colorTheme.params.lightness !== values.lightness) {
                update.to(cell).update(old => {
                    old.colorTheme.params.lightness = values.lightness;
                });
            }
        }
        update.commit();
    };

    updateOpacity = (values: PD.Values) => {
        const update = this.plugin.state.data.build();
        const reprCells = this.plugin.state.data.select(StateSelection.Generators.ofTransformer(StructureRepresentation3D, this.ref));

        for (let i = 0; i < reprCells.length; ++i) {
            const cell = reprCells[i];
            const params = cell.transform.params as StateTransformer.Params<StructureRepresentation3D>;
            if (params.type.params.alpha !== values.alpha) {
                update.to(cell).update(old => {
                    old.type.params.alpha = values.alpha;
                    old.type.params.xrayShaded = values.alpha < 1;
                });
            }
        }
        update.commit();
    };

    updateClip = async (props: Clip.Props) => {
        const update = this.plugin.state.data.build();
        this.plugin.state.data.select(StateSelection.Generators.ofTransformer(StructureRepresentation3D, this.ref)).forEach(cell => {

            const params = cell.transform.params as StateTransformer.Params<StructureRepresentation3D>;
            if (!PD.areEqual(Clip.Params, params.type.params.clip, props)) {
                update.to(cell).update(old => {
                    old.type.params.clip = props;
                });
            }
        });
        await update.commit();
    };

    updateLod = (values: PD.Values) => {
        const update = this.plugin.state.data.build();
        const reprCells = this.plugin.state.data.select(StateSelection.Generators.ofTransformer(StructureRepresentation3D, this.ref));

        for (let i = 0; i < reprCells.length; ++i) {
            const cell = reprCells[i];
            const params = cell.transform.params as StateTransformer.Params<StructureRepresentation3D>;
            if (!deepEqual(params.type.params.lodLevels, values.lodLevels) || params.type.params.cellSize !== values.cellSize) {
                update.to(cell).update(old => {
                    old.type.params.lodLevels = values.lodLevels;
                    old.type.params.cellSize = values.cellSize;
                });
            }
        }
        update.commit();
    };

    render() {
        const { cell } = this.props;
        const { obj } = cell;
        if (!obj) return null;

        const cellState = cell.state;
        const disabled = cell.status !== 'error' && cell.status !== 'ok';
        const depth = this.props.depth;
        const colorValue = this.colorValue;
        const lightnessValue = this.lightnessValue;
        const opacityValue = this.opacityValue;
        const lodValue = this.lodValue;
        const hasEntities = this.hasEntities;

        const groups = this.groups;
        const entities = this.entities;

        const label = <Button className={`msp-btn-tree-label msp-type-class-${obj.type.typeClass}`} noOverflow disabled={disabled}>
            <span>{this.groupLabel}</span>
        </Button>;

        const expand = <IconButton svg={cellState.isCollapsed ? ArrowRightSvg : ArrowDropDownSvg} flex='20px' disabled={disabled} onClick={this.toggleExpanded} transparent className='msp-no-hover-outline' style={{ visibility: groups.length > 0 || entities.length > 0 ? 'visible' : 'hidden' }} />;
        const color = hasEntities && colorValue !== undefined && <Button style={{ backgroundColor: Color.toStyle(colorValue), minWidth: 32, width: 32, borderRight: `6px solid ${Color.toStyle(Color.lighten(colorValue, lightnessValue?.lightness || 0))}` }} onClick={this.toggleColor} />;
        const clip = <IconButton svg={ContentCutSvg} toggleState={false} disabled={disabled} small onClick={this.toggleClip} />;
        const visibility = <IconButton svg={cellState.isHidden ? VisibilityOffOutlinedSvg : VisibilityOutlinedSvg} toggleState={false} disabled={disabled} small onClick={this.toggleVisible} />;

        return <>
            <div className={`msp-flex-row`} style={{ margin: `1px 5px 1px ${depth * 10 + 5}px` }} onMouseEnter={this.highlight} onMouseLeave={this.clearHighlight}>
                {expand}
                {label}
                {color}
                {clip}
                {visibility}
            </div>
            {this.state.action === 'color' && colorValue !== void 0 && <div style={{ marginRight: 5 }} className='msp-accent-offset'>
                <ControlGroup header='Color' initialExpanded={true} hideExpander={true} hideOffset={true} onHeaderClick={this.toggleColor}
                    topRightIcon={CloseSvg} noTopMargin childrenClassName='msp-viewport-controls-panel-controls'>
                    <CombinedColorControl param={ColorParam} value={colorValue ?? Color(0xFFFFFF)} onChange={this.updateColor} name='color' hideNameRow />
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
            {!cellState.isCollapsed && <>
                {groups.map(s => {
                    return <GroupNode filter={this.props.filter} cell={s} depth={depth + 1} key={s.transform.ref} />;
                })}
                {entities.map(s => {
                    return <EntityNode cell={s} depth={depth + 1} key={s.transform.ref} />;
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
        const repr = this.repr?.obj?.data.repr;
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

    get entityLabel() {
        try {
            const s = this.props.cell.obj!.data;
            const l = StructureElement.Location.create(s, s.units[0], s.units[0].elements[0]);
            const d = StructureProperties.entity.pdbx_description(l)[0] || 'model';
            return d.split('.').at(-1);
        } catch (e) {
            console.log(e, this.props.cell.obj!.data);
            return 'Entity';
        }
    }

    get repr(): StateObjectCell | undefined {
        return this.plugin.state.data.select(StateSelection.Generators.ofTransformer(StructureRepresentation3D, this.ref))[0];
    }

    get colorValue(): Color | undefined {
        const repr = this.repr;
        const hasValue = repr?.transform.params?.colorTheme.params.value !== undefined;
        if (!hasValue) return;
        return repr.transform.params?.colorTheme.params.value;
    }

    get lightnessValue(): { lightness: number } | undefined {
        const repr = this.repr;
        const hasLightness = repr?.transform.params?.colorTheme.params.value !== undefined;
        if (!hasLightness) return;
        return {
            lightness: repr.transform.params?.colorTheme.params.lightness
        };
    }

    get opacityValue(): { alpha: number } | undefined {
        const repr = this.repr;
        const hasOpacity = repr?.transform.params?.type.params.alpha !== undefined;
        if (!hasOpacity) return;
        return {
            alpha: repr.transform.params?.type.params.alpha
        };
    }

    get clipValue(): Clip.Props | undefined {
        return this.repr?.transform.params.type.params.clip;
    }

    get lodValue(): PD.Values<typeof LodParams> | undefined {
        const repr = this.repr;
        const hasLod = repr?.transform.params?.type.params.lodLevels !== undefined && repr?.transform.params?.type.params.cellSize !== undefined;
        if (!hasLod) return;
        return {
            lodLevels: repr.transform.params.type.params.lodLevels,
            cellSize: repr.transform.params.type.params.cellSize,
        };
    }

    updateColor: ParamOnChange = ({ value }) => {
        const t = this.repr?.transform;
        if (!t) return;

        return this.plugin.build().to(t.ref).update(old => {
            old.colorTheme.params.value = value;
        }).commit();
    };

    updateLightness = (values: PD.Values) => {
        const t = this.repr?.transform;
        if (!t) return;

        return this.plugin.build().to(t.ref).update(old => {
            old.colorTheme.params.lightness = values.lightness;
        }).commit();
    };

    updateOpacity = (values: PD.Values) => {
        const t = this.repr?.transform;
        if (!t) return;

        return this.plugin.build().to(t.ref).update(old => {
            old.type.params.alpha = values.alpha;
            old.type.params.xrayShaded = values.alpha < 1;
        }).commit();
    };

    updateClip = (props: Clip.Props) => {
        const t = this.repr?.transform;
        if (!t) return;

        const params = t.params as StateTransformer.Params<StructureRepresentation3D>;
        if (!PD.areEqual(Clip.Params, params.type.params.clip, props)) {
            this.plugin.build().to(t.ref).update(old => {
                old.type.params.clip = props;
            }).commit();
        }
    };

    updateLod = (values: PD.Values) => {
        const t = this.repr?.transform;
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
        const { cell } = this.props;
        const { obj } = cell;
        if (!obj) return null;

        const cellState = cell.state;
        const disabled = cell.status !== 'error' && cell.status !== 'ok';
        const depth = this.props.depth;
        const colorValue = this.colorValue;
        const lightnessValue = this.lightnessValue;
        const opacityValue = this.opacityValue;
        const lodValue = this.lodValue;

        const label = <Button className={`msp-btn-tree-label msp-type-class-${obj.type.typeClass}`} noOverflow disabled={disabled}>
            <span>{this.entityLabel}</span>
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
                    <CombinedColorControl param={ColorParam} value={colorValue ?? Color(0xFFFFFF)} onChange={this.updateColor} name='color' hideNameRow />
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


