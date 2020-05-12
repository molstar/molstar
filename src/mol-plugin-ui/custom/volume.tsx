/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginUIComponent } from '../base';
import { StateTransformParameters } from '../state/common';
import * as React from 'react';
import { VolumeStreaming } from '../../mol-plugin/behavior/dynamic/volume-streaming/behavior';
import { ExpandableControlRow, IconButton } from '../controls/common';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ParameterControls, ParamOnChange } from '../controls/parameters';
import { Slider } from '../controls/slider';
import { Volume, Grid } from '../../mol-model/volume';
import { Vec3 } from '../../mol-math/linear-algebra';
import { ColorNames } from '../../mol-util/color/names';
import { toPrecision } from '../../mol-util/number';
import { StateSelection, StateObjectCell } from '../../mol-state';
import { setSubtreeVisibility } from '../../mol-plugin/behavior/static/state';
import { VisibilityOutlinedSvg, VisibilityOffOutlinedSvg } from '../controls/icons';

const ChannelParams = {
    color: PD.Color(ColorNames.black, { description: 'Display color of the volume.' }),
    wireframe: PD.Boolean(false, { description: 'Control display of the volume as a wireframe.' }),
    opacity: PD.Numeric(0.3, { min: 0, max: 1, step: 0.01 }, { description: 'Opacity of the volume.' })
};
type ChannelParams = PD.Values<typeof ChannelParams>

const Bounds = new Map<VolumeStreaming.ChannelType, [number, number]>([
    ['em', [-5, 5]],
    ['2fo-fc', [0, 3]],
    ['fo-fc(+ve)', [1, 5]],
    ['fo-fc(-ve)', [-5, -1]],
]);

class Channel extends PluginUIComponent<{
    label: string,
    name: VolumeStreaming.ChannelType,
    channels: { [k: string]: VolumeStreaming.ChannelParams },
    isRelative: boolean,
    params: StateTransformParameters.Props,
    stats: Grid['stats'],
    changeIso: (name: string, value: number, isRelative: boolean) => void,
    changeParams: (name: string, param: string, value: any) => void,
    bCell: StateObjectCell,
    isDisabled?: boolean,
    isUnbounded?: boolean
}> {
    private ref = StateSelection.findTagInSubtree(this.plugin.state.data.tree, this.props.bCell!.transform.ref, this.props.name);

    componentDidUpdate() {
        this.ref = StateSelection.findTagInSubtree(this.plugin.state.data.tree, this.props.bCell!.transform.ref, this.props.name);
    }

    componentDidMount() {
        this.subscribe(this.plugin.state.data.events.cell.stateUpdated, e => {
            if (this.ref === e.ref) this.forceUpdate();
        });
    }

    getVisible = () => {
        const state = this.plugin.state.data;
        const ref = this.ref;
        if (!ref) return false;
        return !state.cells.get(ref)!.state.isHidden;
    };

    toggleVisible = () => {
        const state = this.plugin.state.data;
        const ref = this.ref;
        if (!ref) return;
        setSubtreeVisibility(state, ref, !state.cells.get(ref)!.state.isHidden);
    };

    render() {
        const props = this.props;
        const { isRelative, stats } = props;
        const channel = props.channels[props.name]!;

        const { min, max, mean, sigma } = stats;
        const value = Math.round(100 * (channel.isoValue.kind === 'relative' ? channel.isoValue.relativeValue : channel.isoValue.absoluteValue)) / 100;
        let relMin = (min - mean) / sigma;
        let relMax = (max - mean) / sigma;

        if (!this.props.isUnbounded) {
            const bounds = Bounds.get(this.props.name)!;
            if (this.props.name === 'em') {
                relMin = Math.max(bounds[0], relMin);
                relMax = Math.min(bounds[1], relMax);
            } else {
                relMin = bounds[0];
                relMax = bounds[1];
            }
        }

        const vMin = mean + sigma * relMin, vMax = mean + sigma * relMax;

        const step = toPrecision(isRelative ? Math.round(((vMax - vMin) / sigma)) / 100 : sigma / 100, 2);
        const ctrlMin = isRelative ? relMin : vMin;
        const ctrlMax = isRelative ? relMax : vMax;

        return <ExpandableControlRow
            label={props.label + (props.isRelative ? ' \u03C3' : '')}
            colorStripe={channel.color}
            pivot={<div className='msp-volume-channel-inline-controls'>
                <Slider value={value} min={ctrlMin} max={ctrlMax} step={step}
                    onChange={v => props.changeIso(props.name, v, isRelative)} disabled={props.params.isDisabled} onEnter={props.params.events.onEnter} />
                <IconButton svg={this.getVisible() ? VisibilityOutlinedSvg : VisibilityOffOutlinedSvg} onClick={this.toggleVisible} toggleState={false} disabled={props.params.isDisabled} />
            </div>}
            controls={<ParameterControls onChange={({ name, value }) => props.changeParams(props.name, name, value)} params={ChannelParams} values={channel} onEnter={props.params.events.onEnter} isDisabled={props.params.isDisabled} />}
        />;
    }
}

export class VolumeStreamingCustomControls extends PluginUIComponent<StateTransformParameters.Props> {
    private areInitial(params: any) {
        return PD.areEqual(this.props.info.params, params, this.props.info.initialValues);
    }

    private newParams(params: VolumeStreaming.Params) {
        this.props.events.onChange(params, this.areInitial(params));
    }

    changeIso = (name: string, value: number, isRelative: boolean) => {
        const old = this.props.params as VolumeStreaming.Params;
        this.newParams({
            ...old,
            entry: {
                name: old.entry.name,
                params: {
                    ...old.entry.params,
                    channels: {
                        ...old.entry.params.channels,
                        [name]: {
                            ...(old.entry.params.channels as any)[name],
                            isoValue: isRelative ? Volume.IsoValue.relative(value) : Volume.IsoValue.absolute(value)
                        }
                    }
                }
            }
        });
    };

    changeParams = (name: string, param: string, value: any) => {
        const old = this.props.params;
        this.newParams({
            ...old,
            entry: {
                name: old.entry.name,
                params: {
                    ...old.entry.params,
                    channels: {
                        ...old.entry.params.channels,
                        [name]: {
                            ...(old.entry.params.channels as any)[name],
                            [param]: value
                        }
                    }
                }
            }
        });
    };

    convert(channel: any, stats: Grid['stats'], isRelative: boolean) {
        return {
            ...channel, isoValue: isRelative
                ? Volume.IsoValue.toRelative(channel.isoValue, stats)
                : Volume.IsoValue.toAbsolute(channel.isoValue, stats)
        };
    }

    changeOption: ParamOnChange = ({ name, value }) => {
        const old = this.props.params as VolumeStreaming.Params;

        if (name === 'entry') {
            this.newParams({
                ...old,
                entry: {
                    name: value,
                    params: old.entry.params,
                }
            });
        } else {
            const b = (this.props.b as VolumeStreaming).data;
            const isEM = b.info.kind === 'em';

            const isRelative = value.params.isRelative;

            const sampling = b.info.header.sampling[0];
            const oldChannels = old.entry.params.channels as any;

            const oldView = old.entry.params.view.name === value.name
                ? old.entry.params.view.params
                : (((this.props.info.params as VolumeStreaming.ParamDefinition)
                    .entry.map(old.entry.name) as PD.Group<VolumeStreaming.EntryParamDefinition>)
                    .params as VolumeStreaming.EntryParamDefinition)
                    .view.map(value.name).defaultValue;

            const viewParams = { ...oldView };
            if (value.name === 'selection-box') {
                viewParams.radius = value.params.radius;
            } else if (value.name === 'box') {
                viewParams.bottomLeft = value.params.bottomLeft;
                viewParams.topRight = value.params.topRight;
            }
            viewParams.isUnbounded = !!value.params.isUnbounded;

            this.newParams({
                ...old,
                entry: {
                    name: old.entry.name,
                    params: {
                        ...old.entry.params,
                        view: {
                            name: value.name,
                            params: viewParams
                        },
                        detailLevel: value.params.detailLevel,
                        channels: isEM
                            ? { em: this.convert(oldChannels.em, sampling.valuesInfo[0], isRelative) }
                            : {
                                '2fo-fc': this.convert(oldChannels['2fo-fc'], sampling.valuesInfo[0], isRelative),
                                'fo-fc(+ve)': this.convert(oldChannels['fo-fc(+ve)'], sampling.valuesInfo[1], isRelative),
                                'fo-fc(-ve)': this.convert(oldChannels['fo-fc(-ve)'], sampling.valuesInfo[1], isRelative)
                            }
                    }
                }
            });
        }
    };

    render() {
        if (!this.props.b) return null;
        const b = (this.props.b as VolumeStreaming).data;

        const isEM = b.info.kind === 'em';
        const pivot = isEM ? 'em' : '2fo-fc';

        const params = this.props.params as VolumeStreaming.Params;
        const detailLevel = ((this.props.info.params as VolumeStreaming.ParamDefinition)
            .entry.map(params.entry.name) as PD.Group<VolumeStreaming.EntryParamDefinition>).params.detailLevel;
        const isRelative = ((params.entry.params.channels as any)[pivot].isoValue as Volume.IsoValue).kind === 'relative';

        const sampling = b.info.header.sampling[0];

        const isRelativeParam = PD.Boolean(isRelative, { description: 'Use normalized or absolute isocontour scale.', label: 'Normalized' });

        const isUnbounded = !!(params.entry.params.view.params as any).isUnbounded;
        const isUnboundedParam = PD.Boolean(isUnbounded, { description: 'Show full/limited range of iso-values for more fine-grained control.', label: 'Unbounded' });

        const isOff = params.entry.params.view.name === 'off';
        // TODO: factor common things out, cache
        const OptionsParams = {
            entry: PD.Select(params.entry.name, b.data.entries.map(info => [info.dataId, info.dataId] as [string, string]), { isHidden: isOff, description: 'Which entry with volume data to display.' }),
            view: PD.MappedStatic(params.entry.params.view.name, {
                'off': PD.Group({
                    isRelative: PD.Boolean(isRelative, { isHidden: true }),
                    isUnbounded: PD.Boolean(isUnbounded, { isHidden: true }),
                }, { description: 'Display off.' }),
                'box': PD.Group({
                    bottomLeft: PD.Vec3(Vec3.zero()),
                    topRight: PD.Vec3(Vec3.zero()),
                    detailLevel,
                    isRelative: isRelativeParam,
                    isUnbounded: isUnboundedParam,
                }, { description: 'Static box defined by cartesian coords.' }),
                'selection-box': PD.Group({
                    radius: PD.Numeric(5, { min: 0, max: 50, step: 0.5 }, { description: 'Radius in \u212B within which the volume is shown.' }),
                    detailLevel,
                    isRelative: isRelativeParam,
                    isUnbounded: isUnboundedParam,
                }, { description: 'Box around focused element.' }),
                'cell': PD.Group({
                    detailLevel,
                    isRelative: isRelativeParam,
                    isUnbounded: isUnboundedParam,
                }, { description: 'Box around the structure\'s bounding box.' }),
                // 'auto': PD.Group({  }), // TODO based on camera distance/active selection/whatever, show whole structure or slice.
            }, { options: VolumeStreaming.ViewTypeOptions, description: 'Controls what of the volume is displayed. "Off" hides the volume alltogether. "Bounded box" shows the volume inside the given box. "Around Focus" shows the volume around the element/atom last interacted with. "Whole Structure" shows the volume for the whole structure.' })
        };
        const options = {
            entry: params.entry.name,
            view: {
                name: params.entry.params.view.name,
                params: {
                    detailLevel: params.entry.params.detailLevel,
                    radius: (params.entry.params.view.params as any).radius,
                    bottomLeft: (params.entry.params.view.params as any).bottomLeft,
                    topRight: (params.entry.params.view.params as any).topRight,
                    isRelative,
                    isUnbounded
                }
            }
        };

        if (isOff) {
            return <ParameterControls onChange={this.changeOption} params={OptionsParams} values={options} onEnter={this.props.events.onEnter} isDisabled={this.props.isDisabled} />;
        }

        return <>
            {!isEM && <Channel label='2Fo-Fc' name='2fo-fc' bCell={this.props.bCell!} channels={params.entry.params.channels} changeIso={this.changeIso} changeParams={this.changeParams} isRelative={isRelative} params={this.props} stats={sampling.valuesInfo[0]} isUnbounded={isUnbounded} />}
            {!isEM && <Channel label='Fo-Fc(+ve)' name='fo-fc(+ve)' bCell={this.props.bCell!} channels={params.entry.params.channels} changeIso={this.changeIso} changeParams={this.changeParams} isRelative={isRelative} params={this.props} stats={sampling.valuesInfo[1]} isUnbounded={isUnbounded} />}
            {!isEM && <Channel label='Fo-Fc(-ve)' name='fo-fc(-ve)' bCell={this.props.bCell!} channels={params.entry.params.channels} changeIso={this.changeIso} changeParams={this.changeParams} isRelative={isRelative} params={this.props} stats={sampling.valuesInfo[1]} isUnbounded={isUnbounded} />}
            {isEM && <Channel label='EM' name='em' bCell={this.props.bCell!} channels={params.entry.params.channels} changeIso={this.changeIso} changeParams={this.changeParams} isRelative={isRelative} params={this.props} stats={sampling.valuesInfo[0]} isUnbounded={isUnbounded} />}

            <ParameterControls onChange={this.changeOption} params={OptionsParams} values={options} onEnter={this.props.events.onEnter} isDisabled={this.props.isDisabled} />
        </>;
    }
}