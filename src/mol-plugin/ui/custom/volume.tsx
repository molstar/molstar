/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginUIComponent } from '../base';
import { StateTransformParameters } from '../state/common';
import * as React from 'react';
import { VolumeStreaming } from '../../../mol-plugin/behavior/dynamic/volume-streaming/behavior';
import { ExpandableGroup } from '../controls/common';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { ParameterControls, ParamOnChange } from '../controls/parameters';
import { Slider } from '../controls/slider';
import { VolumeIsoValue, VolumeData } from '../../../mol-model/volume';
import { Vec3 } from '../../../mol-math/linear-algebra';
import { ColorNames } from '../../../mol-util/color/names';

const ChannelParams = {
    color: PD.Color(ColorNames.black, { description: 'Display color of the volume.' }),
    wireframe: PD.Boolean(false, { description: 'Control display of the volume as a wireframe.' }),
    opacity: PD.Numeric(0.3, { min: 0, max: 1, step: 0.01 }, { description: 'Opacity of the volume.' })
};
type ChannelParams = PD.Values<typeof ChannelParams>

function Channel(props: {
    label: string,
    name: VolumeStreaming.ChannelType,
    channels: { [k: string]: VolumeStreaming.ChannelParams },
    isRelative: boolean,
    params: StateTransformParameters.Props,
    stats: VolumeData['dataStats'],
    changeIso: (name: string, value: number, isRelative: boolean) => void
    changeParams: (name: string, param: string, value: any) => void
}) {
    const { isRelative, stats } = props;
    const channel = props.channels[props.name]!;

    const { min, max, mean, sigma } = stats;
    const value = channel.isoValue.kind === 'relative' ? channel.isoValue.relativeValue : channel.isoValue.absoluteValue;
    const relMin = (min - mean) / sigma;
    const relMax = (max - mean) / sigma;

    return <ExpandableGroup
        label={props.label + (props.isRelative ? ' \u03C3' : '')}
        colorStripe={channel.color}
        pivot={<Slider value={value} min={isRelative ? relMin : min} max={isRelative ? relMax : max}
            step={isRelative ? sigma / 100 : Math.round(((max - min) / sigma)) / 100}
            onChange={v => props.changeIso(props.name, v, isRelative)} disabled={props.params.isDisabled} onEnter={props.params.events.onEnter} />}
        controls={<ParameterControls onChange={({ name, value }) => props.changeParams(props.name, name, value)} params={ChannelParams} values={channel} onEnter={props.params.events.onEnter} />}
    />;
}

export class VolumeStreamingCustomControls extends PluginUIComponent<StateTransformParameters.Props> {

    private areInitial(params: any) {
        return PD.areEqual(this.props.info.params, params, this.props.info.initialValues);
    }

    private newParams(params: VolumeStreaming.Params) {
        this.props.events.onChange(params, this.areInitial(params));
    }

    changeIso = (name: string, value: number, isRelative: boolean) => {
        const old = this.props.params as VolumeStreaming.Params
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
                            isoValue: isRelative ? VolumeIsoValue.relative(value) : VolumeIsoValue.absolute(value)
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

    convert(channel: any, stats: VolumeData['dataStats'], isRelative: boolean) {
        return { ...channel, isoValue: isRelative
            ? VolumeIsoValue.toRelative(channel.isoValue, stats)
            : VolumeIsoValue.toAbsolute(channel.isoValue, stats) }
    }

    changeOption: ParamOnChange = ({ name, value }) => {
        const old = this.props.params as VolumeStreaming.Params

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
            .entry.map(params.entry.name) as PD.Group<VolumeStreaming.EntryParamDefinition>).params.detailLevel
        const isRelative = ((params.entry.params.channels as any)[pivot].isoValue as VolumeIsoValue).kind === 'relative';

        const sampling = b.info.header.sampling[0];

        // TODO: factor common things out
        const OptionsParams = {
            entry: PD.Select(params.entry.name, b.data.entries.map(info => [info.dataId, info.dataId] as [string, string]), { description: 'Which entry with volume data to display.' }),
            view: PD.MappedStatic(params.entry.params.view.name, {
                'off': PD.Group({}, { description: 'Display off.' }),
                'box': PD.Group({
                    bottomLeft: PD.Vec3(Vec3.zero()),
                    topRight: PD.Vec3(Vec3.zero()),
                    detailLevel,
                    isRelative: PD.Boolean(isRelative, { description: 'Use relative or absolute iso values.' })
                }, { description: 'Static box defined by cartesian coords.' }),
                'selection-box': PD.Group({
                    radius: PD.Numeric(5, { min: 0, max: 50, step: 0.5 }, { description: 'Radius in \u212B within which the volume is shown.' }),
                    detailLevel,
                    isRelative: PD.Boolean(isRelative, { description: 'Use relative or absolute iso values.' })
                }, { description: 'Box around last-interacted element.' }),
                'cell': PD.Group({
                    detailLevel,
                    isRelative: PD.Boolean(isRelative, { description: 'Use relative or absolute iso values.' })
                }, { description: 'Box around the structure\'s bounding box.' }),
                // 'auto': PD.Group({  }), // TODO based on camera distance/active selection/whatever, show whole structure or slice.
            }, { options: [['off', 'Off'], ['box', 'Bounded Box'], ['selection-box', 'Surroundings'], ['cell', 'Whole Structure']], description: 'Controls what of the volume is displayed. "Off" hides the volume alltogether. "Bounded box" shows the volume inside the given box. "Arround Interaction" shows the volume around the element/atom last interacted with. "Whole Structure" shows the volume for the whole structure.' })
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
                    isRelative
                }
            }
        };

        return <>
            {!isEM && <Channel label='2Fo-Fc' name='2fo-fc' channels={params.entry.params.channels} changeIso={this.changeIso} changeParams={this.changeParams} isRelative={isRelative} params={this.props} stats={sampling.valuesInfo[0]} />}
            {!isEM && <Channel label='Fo-Fc(+ve)' name='fo-fc(+ve)' channels={params.entry.params.channels} changeIso={this.changeIso} changeParams={this.changeParams} isRelative={isRelative} params={this.props} stats={sampling.valuesInfo[1]} />}
            {!isEM && <Channel label='Fo-Fc(-ve)' name='fo-fc(-ve)' channels={params.entry.params.channels} changeIso={this.changeIso} changeParams={this.changeParams} isRelative={isRelative} params={this.props} stats={sampling.valuesInfo[1]} />}
            {isEM && <Channel label='EM' name='em' channels={params.entry.params.channels} changeIso={this.changeIso} changeParams={this.changeParams} isRelative={isRelative} params={this.props} stats={sampling.valuesInfo[0]} />}

            <ParameterControls onChange={this.changeOption} params={OptionsParams} values={options} onEnter={this.props.events.onEnter} />
        </>
    }
}