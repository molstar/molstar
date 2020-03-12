/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react';
import { Canvas3DParams, Canvas3DProps } from '../../mol-canvas3d/canvas3d';
import { PluginCommands } from '../../mol-plugin/commands';
import { ColorNames } from '../../mol-util/color/names';
import { ParameterMappingControl } from '../controls/parameters';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { PluginUIComponent } from '../base';
import { Color } from '../../mol-util/color';
import { ParamMapping } from '../../mol-util/param-mapping';
import { PluginContext } from '../../mol-plugin/context';
import { Mutable } from '../../mol-util/type-helpers';
import { produce } from 'immer';
import { StateTransform } from '../../mol-state';

export class SimpleSettingsControl extends PluginUIComponent {
    componentDidMount() {
        this.subscribe(this.plugin.events.canvas3d.settingsUpdated, () => this.forceUpdate());
    }

    render() {
        if (!this.plugin.canvas3d) return null;
        return <ParameterMappingControl mapping={SimpleSettingsMapping} />
    }
}

const SimpleSettingsParams = {
    layout: PD.MultiSelect<'sequence' | 'log' | 'left'>([], [['sequence', 'Sequence'], ['log', 'Log'], ['left', 'Left Panel']] as const),
    spin: Canvas3DParams.trackball.params.spin,
    camera: Canvas3DParams.cameraMode,
    background: PD.MappedStatic('opaque', {
        'transparent': PD.EmptyGroup(),
        'opaque': PD.Group({ color: PD.Color(Color(0xFCFBF9), { description: 'Custom background color' }) }, { isFlat: true })
    }, { description: 'Background of the 3D canvas' }),
    renderStyle: PD.Select('glossy', PD.arrayToOptions(['flat', 'matte', 'glossy', 'metallic']), { description: 'Style in which the 3D scene is rendered' }),
    occlusion: PD.Boolean(false, { description: 'Darken occluded crevices with the ambient occlusion effect' }),
    outline: PD.Boolean(false, { description: 'Draw outline around 3D objects' }),
    fog: PD.Boolean(false, { description: 'Show fog in the distance' }),
    clipFar: PD.Boolean(true, { description: 'Clip scene in the distance' }),
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

        let renderStyle: SimpleSettingsParams['renderStyle']['defaultValue'] = 'custom' as any;
        if (renderer) {
            if (renderer.lightIntensity === 0 && renderer.ambientIntensity === 1 && renderer.roughness === 0.4 && renderer.metalness === 0) {
                renderStyle = 'flat'
            } else if (renderer.lightIntensity === 0.6 && renderer.ambientIntensity === 0.4) {
                if (renderer.roughness === 1 && renderer.metalness === 0) {
                    renderStyle = 'matte'
                } else if (renderer.roughness === 0.4 && renderer.metalness === 0) {
                    renderStyle = 'glossy'
                } else if (renderer.roughness === 0.6 && renderer.metalness === 0.4) {
                    renderStyle = 'metallic'
                }
            }
        }

        return {
            layout: props.layout,
            spin: !!canvas.trackball.spin,
            camera: canvas.cameraMode,
            background:  (renderer.backgroundColor === ColorNames.white && canvas.transparentBackground) 
                ? { name: 'transparent', params: { } }
                : { name: 'opaque', params: { color: renderer.backgroundColor } },
            renderStyle,
            occlusion: canvas.postprocessing.occlusionEnable,
            outline: canvas.postprocessing.outlineEnable,
            fog: ctx.canvas3d ? canvas.cameraFog > 1 : false,
            clipFar: canvas.cameraClipFar
        };
    },
    update(s, props) {
        const canvas = props.canvas as Mutable<Canvas3DProps>;
        canvas.trackball.spin = s.spin;
        canvas.cameraMode = s.camera;
        canvas.transparentBackground = s.background.name === 'transparent';
        canvas.renderer.backgroundColor = s.background.name === 'transparent' ? ColorNames.white : s.background.params.color;
        switch (s.renderStyle) {
            case 'flat': Object.assign(canvas.renderer, { lightIntensity: 0, ambientIntensity: 1, roughness: 0.4, metalness: 0 }); break;
            case 'matte':  Object.assign(canvas.renderer, { lightIntensity: 0.6, ambientIntensity: 0.4, roughness: 1, metalness: 0 }); break;
            case 'glossy':  Object.assign(canvas.renderer, { lightIntensity: 0.6, ambientIntensity: 0.4, roughness: 0.4, metalness: 0 }); break;
            case 'metallic':  Object.assign(canvas.renderer, { lightIntensity: 0.6, ambientIntensity: 0.4, roughness: 0.6, metalness: 0.4 }); break;
        }
        canvas.postprocessing.occlusionEnable = s.occlusion;
        if (s.occlusion) { 
            canvas.postprocessing.occlusionBias = 0.5;
            canvas.postprocessing.occlusionRadius = 64;
        }
        canvas.postprocessing.outlineEnable = s.outline;
        canvas.cameraFog = s.fog ? 50 : 0;
        canvas.cameraClipFar = s.clipFar;

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
            PluginCommands.State.SetCurrentObject(ctx, { state: ctx.state.dataState, ref: StateTransform.RootRef });
        }
    }
})