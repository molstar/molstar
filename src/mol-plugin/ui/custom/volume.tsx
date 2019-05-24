/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginUIComponent } from '../base';
import { StateTransformParameters } from '../state/common';
import * as React from 'react';
import { VolumeStreaming } from 'mol-plugin/behavior/dynamic/volume-streaming/behavior';
import { ExpandableGroup } from '../controls/common';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { ParameterControls, ParamOnChange } from '../controls/parameters';
import { Slider } from '../controls/slider';
import { VolumeIsoValue, VolumeData } from 'mol-model/volume';
import { Vec3 } from 'mol-math/linear-algebra';

const ChannelParams = {
    color: PD.Color(0 as any),
    wireframe: PD.Boolean(false),
    opacity: PD.Numeric(0.3, { min: 0, max: 1, step: 0.01 })
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
        const old = this.props.params;
        this.newParams({
            ...old,
            channels: {
                ...old.channels,
                [name]: {
                    ...old.channels[name],
                    isoValue: isRelative ? VolumeIsoValue.relative(value) : VolumeIsoValue.absolute(value)
                }
            }
        });
    };

    changeParams = (name: string, param: string, value: any) => {
        const old = this.props.params;
        this.newParams({
            ...old,
            channels: {
                ...old.channels,
                [name]: {
                    ...old.channels[name],
                    [param]: value
                }
            }
        });
    };

    convert(channel: any, stats: VolumeData['dataStats'], isRelative: boolean) {
        return { ...channel, isoValue: isRelative
            ? VolumeIsoValue.toRelative(channel.isoValue, stats)
            : VolumeIsoValue.toAbsolute(channel.isoValue, stats) }
    }

    changeOption: ParamOnChange = ({ value }) => {
        const b = (this.props.b as VolumeStreaming).data;
        const isEM = b.info.kind === 'em';

        const isRelative = value.params.isRelative;
        const sampling = b.info.header.sampling[0];
        const old = this.props.params as VolumeStreaming.Params, oldChannels = old.channels as any;

        const oldView = old.view.name === value.name
            ? old.view.params
            : (this.props.info.params as VolumeStreaming.ParamDefinition).view.map(value.name).defaultValue;

        const viewParams = { ...oldView };
        if (value.name === 'selection-box') {
            viewParams.radius = value.params.radius;
        } else if (value.name === 'box') {
            viewParams.bottomLeft = value.params.bottomLeft;
            viewParams.topRight = value.params.topRight;
        }

        this.newParams({
            ...old,
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
        });
    };

    render() {
        if (!this.props.b) return null;

        const b = (this.props.b as VolumeStreaming).data;
        const isEM = b.info.kind === 'em';
        const pivot = isEM ? 'em' : '2fo-fc';

        const params = this.props.params as VolumeStreaming.Params;
        const isRelative = ((params.channels as any)[pivot].isoValue as VolumeIsoValue).kind === 'relative';

        const sampling = b.info.header.sampling[0];

        // TODO: factor common things out
        const OptionsParams = {
            view: PD.MappedStatic(params.view.name, {
                'box': PD.Group({
                    bottomLeft: PD.Vec3(Vec3.zero()),
                    topRight: PD.Vec3(Vec3.zero()),
                    detailLevel: this.props.info.params.detailLevel,
                    isRelative: PD.Boolean(isRelative, { description: 'Use relative or absolute iso values.' })
                }, { description: 'Static box defined by cartesian coords.' }),
                'selection-box': PD.Group({
                    radius: PD.Numeric(5, { min: 0, max: 50, step: 0.5 }),
                    detailLevel: this.props.info.params.detailLevel,
                    isRelative: PD.Boolean(isRelative, { description: 'Use relative or absolute iso values.' })
                }, { description: 'Box around last-interacted element.' }),
                'cell': PD.Group({
                    detailLevel: this.props.info.params.detailLevel,
                    isRelative: PD.Boolean(isRelative, { description: 'Use relative or absolute iso values.' })
                }, { description: 'Box around the structure\'s bounding box.' }),
                // 'auto': PD.Group({  }), // based on camera distance/active selection/whatever, show whole structure or slice.
            }, { options: [['box', 'Bounded Box'], ['selection-box', 'Selection'], ['cell', 'Whole Structure']] })
        };
        const options = {
            view: {
                name: params.view.name,
                params: {
                    detailLevel: params.detailLevel,
                    radius: (params.view.params as any).radius,
                    bottomLeft: (params.view.params as any).bottomLeft,
                    topRight: (params.view.params as any).topRight,
                    isRelative
                }
            }
        };

        return <>
            {!isEM && <Channel label='2Fo-Fc' name='2fo-fc' channels={params.channels} changeIso={this.changeIso} changeParams={this.changeParams} isRelative={isRelative} params={this.props} stats={sampling.valuesInfo[0]} />}
            {!isEM && <Channel label='Fo-Fc(+ve)' name='fo-fc(+ve)' channels={params.channels} changeIso={this.changeIso} changeParams={this.changeParams} isRelative={isRelative} params={this.props} stats={sampling.valuesInfo[1]} />}
            {!isEM && <Channel label='Fo-Fc(-ve)' name='fo-fc(-ve)' channels={params.channels} changeIso={this.changeIso} changeParams={this.changeParams} isRelative={isRelative} params={this.props} stats={sampling.valuesInfo[1]} />}
            {isEM && <Channel label='EM' name='em' channels={params.channels} changeIso={this.changeIso} changeParams={this.changeParams} isRelative={isRelative} params={this.props} stats={sampling.valuesInfo[0]} />}

            <ParameterControls onChange={this.changeOption} params={OptionsParams} values={options} onEnter={this.props.events.onEnter} />
        </>
    }
}