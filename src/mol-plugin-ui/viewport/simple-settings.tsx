/**
 * Copyright (c) 2019-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { produce } from 'immer';
import { throttleTime } from 'rxjs';
import { Canvas3DParams, Canvas3DProps } from '../../mol-canvas3d/canvas3d';
import { PluginCommands } from '../../mol-plugin/commands';
import { PluginConfig } from '../../mol-plugin/config';
import { StateTransform } from '../../mol-state';
import { Color } from '../../mol-util/color';
import { deepClone } from '../../mol-util/object';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ParamMapping } from '../../mol-util/param-mapping';
import { Mutable } from '../../mol-util/type-helpers';
import { PluginUIComponent } from '../base';
import { PluginUIContext } from '../context';
import { ParameterMappingControl } from '../controls/parameters';
import { ViewportHelpContent } from './help';

export class SimpleSettingsControl extends PluginUIComponent {
    componentDidMount() {
        if (!this.plugin.canvas3d) return;

        this.subscribe(this.plugin.events.canvas3d.settingsUpdated, () => this.forceUpdate());

        this.subscribe(this.plugin.canvas3d!.camera.stateChanged.pipe(throttleTime(500, undefined, { leading: true, trailing: true })), state => {
            if (state.radiusMax !== undefined || state.radius !== undefined) {
                this.forceUpdate();
            }
        });
    }

    render() {
        if (!this.plugin.canvas3d) return null;
        return <>
            <ParameterMappingControl mapping={SimpleSettingsMapping} />
            <ViewportHelpContent />
        </>;
    }
}

const LayoutOptions = {
    'sequence': 'Sequence',
    'log': 'Log',
    'left': 'Left Panel'
};
type LayoutOptions = keyof typeof LayoutOptions

const SimpleSettingsParams = {
    animate: Canvas3DParams.trackball.params.animate,
    camera: Canvas3DParams.camera,
    background: PD.Group({
        color: PD.Color(Color(0xFCFBF9), { label: 'Background', description: 'Custom background color' }),
        transparent: PD.Boolean(false),
        style: Canvas3DParams.postprocessing.params.background,
    }, { pivot: 'color' }),
    lighting: PD.Group({
        occlusion: Canvas3DParams.postprocessing.params.occlusion,
        shadow: Canvas3DParams.postprocessing.params.shadow,
        outline: Canvas3DParams.postprocessing.params.outline,
        fog: Canvas3DParams.cameraFog,
    }, { isFlat: true }),
    clipping: PD.Group<any>({
        ...Canvas3DParams.cameraClipping.params,
    }, { pivot: 'radius' }),
    layout: PD.MultiSelect([] as LayoutOptions[], PD.objectToOptions(LayoutOptions)),
};

type SimpleSettingsParams = typeof SimpleSettingsParams
const SimpleSettingsMapping = ParamMapping({
    params: (ctx: PluginUIContext) => {
        const params = PD.clone(SimpleSettingsParams);
        const controls = ctx.spec.components?.controls;
        if (controls) {
            const options: [LayoutOptions, string][] = [];
            if (controls.top !== 'none') options.push(['sequence', LayoutOptions.sequence]);
            if (controls.bottom !== 'none') options.push(['log', LayoutOptions.log]);
            if (controls.left !== 'none') options.push(['left', LayoutOptions.left]);
            params.layout.options = options;
        }
        const bgStyles = ctx.config.get(PluginConfig.Background.Styles) || [];
        if (bgStyles.length > 0) {
            Object.assign(params.background.params.style, {
                presets: deepClone(bgStyles),
                isFlat: false, // so the presets menu is shown
            });
        }
        return params;
    },
    target(ctx: PluginUIContext) {
        const c = ctx.spec.components?.controls;
        const r = ctx.layout.state.regionState;
        const layout: SimpleSettingsParams['layout']['defaultValue'] = [];
        if (r.top !== 'hidden' && (!c || c.top !== 'none')) layout.push('sequence');
        if (r.bottom !== 'hidden' && (!c || c.bottom !== 'none')) layout.push('log');
        if (r.left !== 'hidden' && (!c || c.left !== 'none')) layout.push('left');
        return { canvas: ctx.canvas3d?.props!, layout };
    }
})({
    values(props, ctx) {
        const { canvas } = props;
        const renderer = canvas.renderer;

        return {
            layout: props.layout,
            animate: canvas.trackball.animate,
            camera: canvas.camera,
            background: {
                color: renderer.backgroundColor,
                transparent: canvas.transparentBackground,
                style: canvas.postprocessing.background,
            },
            lighting: {
                occlusion: canvas.postprocessing.occlusion,
                shadow: canvas.postprocessing.shadow,
                outline: canvas.postprocessing.outline,
                fog: canvas.cameraFog,
            },
            clipping: {
                ...canvas.cameraClipping,
            }
        };
    },
    update(s, props) {
        const canvas = props.canvas as Mutable<Canvas3DProps>;
        canvas.trackball.animate = s.animate;
        canvas.camera = s.camera;
        canvas.transparentBackground = s.background.transparent;
        canvas.renderer.backgroundColor = s.background.color;
        canvas.postprocessing.occlusion = s.lighting.occlusion;
        canvas.postprocessing.shadow = s.lighting.shadow;
        canvas.postprocessing.outline = s.lighting.outline;
        canvas.postprocessing.background = s.background.style;
        canvas.cameraFog = s.lighting.fog;
        canvas.cameraClipping = {
            radius: s.clipping.radius,
            far: s.clipping.far,
            minNear: s.clipping.minNear,
        };

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
});