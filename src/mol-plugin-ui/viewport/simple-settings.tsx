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
import { ParameterControls } from '../controls/parameters';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { PluginUIComponent } from '../base';
import { Color } from '../../mol-util/color';

const SimpleSettingsParams = {
    spin: Canvas3DParams.trackball.params.spin,
    camera: Canvas3DParams.cameraMode,
    background: PD.MappedStatic('opaque', {
        'transparent': PD.EmptyGroup(),
        'opaque': PD.Group({ color: PD.Color(Color(0xFCFBF9), { description: 'Custom background color' }) }, { isFlat: true })
    }, { description: 'Background of the 3D canvas' }),
    renderStyle: PD.Select('glossy', [['flat', 'Flat'], ['matte', 'Matte'], ['glossy', 'Glossy'], ['metallic', 'Metallic']], { description: 'Style in which the 3D scene is rendered' }),
    occlusion: PD.Boolean(false, { description: 'Darken occluded crevices with the ambient occlusion effect' }),
    outline: PD.Boolean(false, { description: 'Draw outline around 3D objects' }),
    fog: PD.Boolean(false, { description: 'Show fog in the distance' }),
    clipFar: PD.Boolean(true, { description: 'Clip scene in the distance' }),
};

export class SimpleSettingsControl extends PluginUIComponent {
    setSettings = (p: { param: PD.Base<any>, name: keyof typeof SimpleSettingsParams | string, value: any }) => {
        if (p.name === 'spin') {
            if (!this.plugin.canvas3d) return;
            const trackball = this.plugin.canvas3d.props.trackball;
            PluginCommands.Canvas3D.SetSettings.dispatch(this.plugin, { settings: { trackball: { ...trackball, spin: p.value } } });
        } else if (p.name === 'camera') {
            PluginCommands.Canvas3D.SetSettings.dispatch(this.plugin, { settings: { cameraMode: p.value }});
        } else if (p.name === 'background') {
            if (!this.plugin.canvas3d) return;
            const renderer = this.plugin.canvas3d.props.renderer;
            const color: typeof SimpleSettingsParams['background']['defaultValue'] = p.value;
            if (color.name === 'transparent') {
                PluginCommands.Canvas3D.SetSettings.dispatch(this.plugin, { settings: { renderer: { ...renderer, backgroundColor: ColorNames.white }, transparentBackground: true } });
            } else {
                PluginCommands.Canvas3D.SetSettings.dispatch(this.plugin, { settings: { renderer: { ...renderer, backgroundColor: color.params.color }, transparentBackground: false } });
            }
        } else if (p.name === 'renderStyle') {
            if (!this.plugin.canvas3d) return;

            const renderer = this.plugin.canvas3d.props.renderer;
            if (p.value === 'flat') {
                PluginCommands.Canvas3D.SetSettings.dispatch(this.plugin, { settings: {
                    renderer: { ...renderer, lightIntensity: 0, ambientIntensity: 1, roughness: 0.4, metalness: 0 }
                } });
            } else if (p.value === 'matte') {
                PluginCommands.Canvas3D.SetSettings.dispatch(this.plugin, { settings: {
                    renderer: { ...renderer, lightIntensity: 0.6, ambientIntensity: 0.4, roughness: 1, metalness: 0 }
                } });
            } else if (p.value === 'glossy') {
                PluginCommands.Canvas3D.SetSettings.dispatch(this.plugin, { settings: {
                    renderer: { ...renderer, lightIntensity: 0.6, ambientIntensity: 0.4, roughness: 0.4, metalness: 0 }
                } });
            } else if (p.value === 'metallic') {
                PluginCommands.Canvas3D.SetSettings.dispatch(this.plugin, { settings: {
                    renderer: { ...renderer, lightIntensity: 0.6, ambientIntensity: 0.4, roughness: 0.6, metalness: 0.4 }
                } });
            }
        } else if (p.name === 'occlusion') {
            if (!this.plugin.canvas3d) return;
            const postprocessing = this.plugin.canvas3d.props.postprocessing;
            PluginCommands.Canvas3D.SetSettings.dispatch(this.plugin, { settings: {
                postprocessing: { ...postprocessing, occlusionEnable: p.value, occlusionBias: 0.5, occlusionRadius: 64 },
            } });
        } else if (p.name === 'outline') {
            if (!this.plugin.canvas3d) return;
            const postprocessing = this.plugin.canvas3d.props.postprocessing;
            PluginCommands.Canvas3D.SetSettings.dispatch(this.plugin, { settings: {
                postprocessing: { ...postprocessing, outlineEnable: p.value },
            } });
        } else if (p.name === 'fog') {;
            PluginCommands.Canvas3D.SetSettings.dispatch(this.plugin, { settings: {
                cameraFog: p.value ? 50 : 1,
            } });
        } else if (p.name === 'clipFar') {;
            PluginCommands.Canvas3D.SetSettings.dispatch(this.plugin, { settings: {
                cameraClipFar: p.value,
            } });
        }
    }

    get values () {
        const renderer = this.plugin.canvas3d?.props.renderer;

        let renderStyle = 'custom'
        let background: typeof SimpleSettingsParams['background']['defaultValue'] = { name: 'transparent', params: { } }

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

            if (renderer.backgroundColor === ColorNames.white && this.plugin.canvas3d?.props.transparentBackground) {
                background = { name: 'transparent', params: { } }
            } else {
                background = { name: 'opaque', params: { color: renderer.backgroundColor } }
            }
        }

        return {
            spin: !!this.plugin.canvas3d?.props.trackball.spin,
            camera: this.plugin.canvas3d?.props.cameraMode,
            background,
            renderStyle,
            occlusion: this.plugin.canvas3d?.props.postprocessing.occlusionEnable,
            outline: this.plugin.canvas3d?.props.postprocessing.outlineEnable,
            fog: this.plugin.canvas3d ? this.plugin.canvas3d.props.cameraFog > 1 : false,
            clipFar: this.plugin.canvas3d?.props.cameraClipFar
        }
    }

    componentDidMount() {
        this.subscribe(this.plugin.events.canvas3d.settingsUpdated, () => this.forceUpdate());
    }

    render() {
        return <ParameterControls params={SimpleSettingsParams} values={this.values} onChange={this.setSettings} />
    }
}