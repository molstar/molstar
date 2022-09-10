/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { debounceTime } from 'rxjs/operators';
import { PluginStateAnimation } from '../../mol-plugin-state/animation/model';
import { PluginComponent } from '../../mol-plugin-state/component';
import { PluginContext } from '../../mol-plugin/context';
import { Task } from '../../mol-task';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { encodeMp4Animation } from './encoder';

export interface Mp4AnimationInfo {
    width: number,
    height: number
}

export const Mp4AnimationParams = {
    quantization: PD.Numeric(18, { min: 10, max: 51 }, { description: 'Lower is better, but slower.' })
};

export class Mp4Controls extends PluginComponent {
    private currentNames = new Set<string>();
    private animations: PluginStateAnimation[] = [];

    readonly behaviors = {
        animations: this.ev.behavior<PD.Params>({ }),
        current: this.ev.behavior<{ anim: PluginStateAnimation, params: PD.Params, values: any } | undefined>(void 0),
        canApply: this.ev.behavior<PluginStateAnimation.CanApply>({ canApply: false }),
        info: this.ev.behavior<Mp4AnimationInfo>({ width: 0, height: 0 }),
        params: this.ev.behavior<PD.Values<typeof Mp4AnimationParams>>(PD.getDefaultValues(Mp4AnimationParams))
    };

    setCurrent(name?: string) {
        const anim = this.animations.find(a => a.name === name);
        if (!anim) {
            this.behaviors.current.next(anim);
            return;
        }

        const params = anim.params(this.plugin) as PD.Params;
        const values = PD.getDefaultValues(params);

        this.behaviors.current.next({ anim, params, values });
        this.behaviors.canApply.next(anim.canApply?.(this.plugin) ?? { canApply: true });
    }

    setCurrentParams(values: any) {
        this.behaviors.current.next({ ...this.behaviors.current.value!, values });
    }

    get current() {
        return this.behaviors.current.value;
    }

    render() {
        const task = Task.create('Export Animation', async ctx => {
            try {
                const resolution = this.plugin.helpers.viewportScreenshot?.getSizeAndViewport()!;
                const anim = this.current!;
                const movie = await encodeMp4Animation(this.plugin, ctx, {
                    animation: {
                        definition: anim.anim,
                        params: anim.values,
                    },
                    ...resolution,
                    quantizationParameter: this.behaviors.params.value.quantization,
                    pass: this.plugin.helpers.viewportScreenshot?.imagePass!,
                });

                const filename = anim.anim.display.name.toLowerCase().replace(/\s/g, '-').replace(/[^a-z0-9_\-]/g, '');
                return { movie, filename: `${this.plugin.helpers.viewportScreenshot?.getFilename('')}_${filename}.mp4` };
            } catch (e) {
                this.plugin.log.error('Error during animation export');
                throw e;
            }
        });

        return this.plugin.runTask(task, { useOverlay: true });
    }

    private get manager() {
        return this.plugin.managers.animation;
    }

    private syncInfo() {
        const helper = this.plugin.helpers.viewportScreenshot;
        const size = helper?.getSizeAndViewport();
        if (!size) return;

        this.behaviors.info.next({ width: size.viewport.width, height: size.viewport.height });
    }

    private sync() {
        const animations = this.manager.animations.filter(a => a.isExportable);

        const hasAll = animations.every(a => this.currentNames.has(a.name));
        if (hasAll && this.currentNames.size === animations.length) {
            return;
        }

        const params = {
            current: PD.Select(animations[0]?.name,
                animations.map(a => [a.name, a.display.name] as [string, string]),
                { label: 'Animation' })
        };

        const current = this.behaviors.current.value;
        const hasCurrent = !!animations.find(a => a.name === current?.anim.name);

        this.animations = animations;
        if (!hasCurrent) {
            this.setCurrent(animations[0]?.name);
        }
        this.behaviors.animations.next(params);
    }

    private init() {
        if (!this.plugin.canvas3d) return;

        this.subscribe(this.plugin.managers.animation.events.updated.pipe(debounceTime(16)), () => {
            this.sync();
        });

        this.subscribe(this.plugin.canvas3d.resized, () => this.syncInfo());
        this.subscribe(this.plugin.helpers.viewportScreenshot?.events.previewed!, () => this.syncInfo());

        this.subscribe(this.plugin.behaviors.state.isBusy, b => this.updateCanApply(b));
        this.subscribe(this.plugin.managers.snapshot.events.changed, b => this.updateCanApply(b));

        this.sync();
        this.syncInfo();
    }

    private updateCanApply(b?: any) {
        const anim = this.current;
        if (!b && anim) {
            this.behaviors.canApply.next(anim.anim.canApply?.(this.plugin) ?? { canApply: true });
        }
    }

    constructor(private plugin: PluginContext) {
        super();

        this.init();
    }
}