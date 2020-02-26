/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react';
import { Canvas3DParams } from '../../mol-canvas3d/canvas3d';
import { PluginCommands } from '../../mol-plugin/command';
import { ColorNames } from '../../mol-util/color/names';
import { ParameterMappingControl } from '../controls/parameters';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { PluginUIComponent } from '../base';
import { Color } from '../../mol-util/color';
import { ParamMapping } from '../../mol-util/param-mapping';
import { PluginContext } from '../../mol-plugin/context';

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
    target(ctx: PluginContext) { return ctx.canvas3d?.props!; } })({
    values(t, ctx) {
        const renderer = t.renderer;

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
            spin: !!t.trackball.spin,
            camera: t.cameraMode,
            background:  (renderer.backgroundColor === ColorNames.white && t.transparentBackground) 
                ? { name: 'transparent', params: { } }
                : { name: 'opaque', params: { color: renderer.backgroundColor } },
            renderStyle,
            occlusion: t.postprocessing.occlusionEnable,
            outline: t.postprocessing.outlineEnable,
            fog: ctx.canvas3d ? t.cameraFog > 1 : false,
            clipFar: t.cameraClipFar
        };
    },
    update(s, t) {
        t.trackball.spin = s.spin;
        t.cameraMode = s.camera;
        t.transparentBackground = s.background.name === 'transparent';
        t.renderer.backgroundColor = s.background.name === 'transparent' ? ColorNames.white : s.background.params.color;
        switch (s.renderStyle) {
            case 'flat': Object.assign(t.renderer, { lightIntensity: 0, ambientIntensity: 1, roughness: 0.4, metalness: 0 }); break;
            case 'matte':  Object.assign(t.renderer, { lightIntensity: 0.6, ambientIntensity: 0.4, roughness: 1, metalness: 0 }); break;
            case 'glossy':  Object.assign(t.renderer, { lightIntensity: 0.6, ambientIntensity: 0.4, roughness: 0.4, metalness: 0 }); break;
            case 'metallic':  Object.assign(t.renderer, { lightIntensity: 0.6, ambientIntensity: 0.4, roughness: 0.6, metalness: 0.4 }); break;
        }
        t.postprocessing.occlusionEnable = s.occlusion;
        if (s.occlusion) { 
            t.postprocessing.occlusionBias = 0.5;
            t.postprocessing.occlusionRadius = 64;
        }
        t.postprocessing.outlineEnable = s.outline;
        t.cameraFog = s.fog ? 50 : 0;
        t.cameraClipFar = s.clipFar;
    },
    apply(settings, ctx) {
        return PluginCommands.Canvas3D.SetSettings.dispatch(ctx, { settings });
    }
})