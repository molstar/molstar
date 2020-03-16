/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { produce } from 'immer';
import * as React from 'react';
import { Canvas3DParams, Canvas3DProps } from '../../mol-canvas3d/canvas3d';
import { PluginCommands } from '../../mol-plugin/commands';
import { PluginContext } from '../../mol-plugin/context';
import { StateTransform } from '../../mol-state';
import { Color } from '../../mol-util/color';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ParamMapping } from '../../mol-util/param-mapping';
import { Mutable } from '../../mol-util/type-helpers';
import { PluginUIComponent } from '../base';
import { ParameterMappingControl } from '../controls/parameters';
import { ViewportHelpContent } from './help';

export class SimpleSettingsControl extends PluginUIComponent {
    componentDidMount() {
        this.subscribe(this.plugin.events.canvas3d.settingsUpdated, () => this.forceUpdate());

        this.plugin.canvas3d!.camera.stateChanged.subscribe(state => {
            if (state.radiusMax !== undefined || state.radius !== undefined) {
                this.forceUpdate()
            }
        })
    }

    render() {
        if (!this.plugin.canvas3d) return null;
        return <>
            <ParameterMappingControl mapping={SimpleSettingsMapping} />
            <ViewportHelpContent />
        </>
    }
}

const SimpleSettingsParams = {
    spin: PD.Group({
        spin: Canvas3DParams.trackball.params.spin,
        speed: Canvas3DParams.trackball.params.spinSpeed
    }, { pivot: 'spin' }),
    camera: Canvas3DParams.cameraMode,
    background: PD.Group({
        color: PD.Color(Color(0xFCFBF9), { label: 'Background', description: 'Custom background color' }),
        transparent: PD.Boolean(false)
    }, { pivot: 'color' }),
    lighting: PD.Group({
        renderStyle: Canvas3DParams.renderer.params.style,
        occlusion: Canvas3DParams.postprocessing.params.occlusion,
        outline: Canvas3DParams.postprocessing.params.outline,
        fog: Canvas3DParams.cameraFog,
    }, { pivot: 'renderStyle' }),
    clipping: Canvas3DParams.cameraClipping,
    layout: PD.MultiSelect<'sequence' | 'log' | 'left'>([], [['sequence', 'Sequence'], ['log', 'Log'], ['left', 'Left Panel']] as const),
};

type SimpleSettingsParams = typeof SimpleSettingsParams
const SimpleSettingsMapping = ParamMapping({
    params: SimpleSettingsParams,
    target(ctx: PluginContext) {
        const layout: SimpleSettingsParams['layout']['defaultValue'] = [];
        if (ctx.layout.state.regionState.top !== 'hidden') layout.push('sequence');
        if (ctx.layout.state.regionState.bottom !== 'hidden') layout.push('log');
        if (ctx.layout.state.regionState.left !== 'hidden') layout.push('left');
        return { canvas: ctx.canvas3d?.props!, layout };
    }
})({
    values(props, ctx) {
        const { canvas } = props;
        const renderer = canvas.renderer;

        return {
            layout: props.layout,
            spin: { spin: !!canvas.trackball.spin, speed: canvas.trackball.spinSpeed },
            camera: canvas.cameraMode,
            background: {
                color: renderer.backgroundColor,
                transparent: canvas.transparentBackground
            },
            lighting: {
                renderStyle: renderer.style,
                occlusion: canvas.postprocessing.occlusion,
                outline: canvas.postprocessing.outline,
                fog: canvas.cameraFog
            },
            clipping: canvas.cameraClipping
        };
    },
    update(s, props) {
        const canvas = props.canvas as Mutable<Canvas3DProps>;
        canvas.trackball.spin = s.spin.spin;
        canvas.trackball.spinSpeed = s.spin.speed;
        canvas.cameraMode = s.camera;
        canvas.transparentBackground = s.background.transparent;
        canvas.renderer.backgroundColor = s.background.color;
        canvas.renderer.style = s.lighting.renderStyle
        canvas.postprocessing.occlusion = s.lighting.occlusion;
        canvas.postprocessing.outline = s.lighting.outline;
        canvas.cameraFog = s.lighting.fog
        canvas.cameraClipping = s.clipping

        props.layout = s.layout;
    },
    async apply(props, ctx) {
        await PluginCommands.Canvas3D.SetSettings(ctx, { settings: props.canvas });

        const hideLeft = props.layout.indexOf('left') < 0;
        const state = produce(ctx.layout.state, s => {
            s.regionState.top = props.layout.indexOf('sequence') >= 0 ? 'full' : 'hidden';
            s.regionState.bottom = props.layout.indexOf('log') >= 0 ? 'full' : 'hidden';
            s.regionState.left = hideLeft ? 'hidden' : ctx.behaviors.layout.leftPanelTabName.value === 'none' ? 'collapsed' : 'full';
        });
        await PluginCommands.Layout.Update(ctx, { state });

        if (hideLeft) {
            PluginCommands.State.SetCurrentObject(ctx, { state: ctx.state.data, ref: StateTransform.RootRef });
        }
    }
})