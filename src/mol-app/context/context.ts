/*
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Adapted from LiteMol
 * Copyright (c) 2016 - now David Sehnal, licensed under Apache 2.0, See LICENSE file for more info.
 */

import { UUID } from 'mol-util'
import { PerformanceMonitor } from 'mol-util/performance-monitor';
import { Dispatcher } from '../service/dispatcher'
import { Logger } from '../service/logger'
import { LayoutTarget, LayoutController } from '../controller/layout';
import { ViewportController } from '../controller/visualization/viewport';
import { Stage } from 'mol-view/stage';
import { AnyTransform } from 'mol-view/state/transform';
import { BehaviorSubject } from 'rxjs';
import { AnyEntity } from 'mol-view/state/entity';
import { SequenceViewController } from '../controller/visualization/sequence-view';

export class Settings {
    private settings = new Map<string, any>();

    set(key: string, value: any) {
        this.settings.set(key, value);
    }

    get(key: string) {
        return this.settings.get(key);
    }
}

export class Context {
    id = UUID.create()

    dispatcher = new Dispatcher();
    logger = new Logger(this);
    performance = new PerformanceMonitor();

    stage = new Stage(this);
    viewport = new ViewportController(this);
    layout: LayoutController;
    settings = new Settings();

    // TODO: this is a temporary solution
    components = {
        sequenceView: new SequenceViewController(this)
    };

    currentEntity = new BehaviorSubject(undefined) as BehaviorSubject<AnyEntity | undefined>
    currentTransforms = new BehaviorSubject([] as AnyTransform[])

    createLayout(targets: LayoutTarget[], target: HTMLElement) {
        this.layout = new LayoutController(this, targets, target);
    }

    initStage(canvas: HTMLCanvasElement, container: HTMLDivElement) {
        this.stage.initRenderer(canvas, container)
        return true
    }

    destroy() {
        if (this.stage) {
            this.stage.dispose()
            this.stage = null as any
        }
    }
}