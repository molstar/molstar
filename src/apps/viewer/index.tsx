/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react'
import * as ReactDOM from 'react-dom'

import './index.html'
import 'mol-app/skin/molstar-light.scss'

import { Context } from 'mol-app/context/context';
import { Viewport } from 'mol-app/ui/visualization/viewport'
import { makeEmptyTargets, LayoutRegion } from 'mol-app/controller/layout';
import { Layout } from 'mol-app/ui/layout';
import { LogController } from 'mol-app/controller/misc/log';
import { Log } from 'mol-app/ui/misc/log';
import { JobsController } from 'mol-app/controller/misc/jobs';
import { BackgroundJobs, Overlay } from 'mol-app/ui/misc/jobs';
import { EntityTree } from 'mol-app/ui/entity/tree';
import { EntityTreeController } from 'mol-app/controller/entity/tree';
import { TransformListController } from 'mol-app/controller/transform/list';
import { TransformList } from 'mol-app/ui/transform/list';
import { SequenceView } from 'mol-app/ui/visualization/sequence-view';
import { InteractivityEvents } from 'mol-app/event/basic';
import { MarkerAction } from 'mol-geo/util/marker-data';
import { EveryLoci } from 'mol-model/loci';

const elm = document.getElementById('app')
if (!elm) throw new Error('Can not find element with id "app".')

const ctx = new Context()
const targets = makeEmptyTargets();

targets[LayoutRegion.Main].components.push({
    key: 'molstar-internal-viewport',
    controller: ctx.viewport,
    region: LayoutRegion.Main,
    view: Viewport,
    isStatic: true
});

targets[LayoutRegion.Bottom].components.push({
    key: 'molstar-log',
    controller: new LogController(ctx),
    region: LayoutRegion.Bottom,
    view: Log,
    isStatic: true
});

targets[LayoutRegion.Top].components.push({
    key: 'molstar-sequence-view',
    controller: ctx.components.sequenceView,
    region: LayoutRegion.Top,
    view: SequenceView,
    isStatic: true
});

targets[LayoutRegion.Main].components.push({
    key: 'molstar-background-jobs',
    controller: new JobsController(ctx, 'Background'),
    region: LayoutRegion.Main,
    view: BackgroundJobs,
    isStatic: true
});

targets[LayoutRegion.Root].components.push({
    key: 'molstar-overlay',
    controller: new JobsController(ctx, 'Normal'),
    region: LayoutRegion.Root,
    view: Overlay,
    isStatic: true
});

targets[LayoutRegion.Right].components.push({
    key: 'molstar-transform-list',
    controller: new TransformListController(ctx),
    region: LayoutRegion.Right,
    view: TransformList,
    isStatic: false
});

targets[LayoutRegion.Left].components.push({
    key: 'molstar-entity-tree',
    controller: new EntityTreeController(ctx),
    region: LayoutRegion.Left,
    view: EntityTree,
    isStatic: true
});

ctx.createLayout(targets, elm)
ctx.layout.setState({
    isExpanded: true,
    hideControls: false,
    collapsedControlsLayout: 0
})
// ctx.viewport.setState()

ctx.dispatcher.getStream(InteractivityEvents.HighlightLoci).subscribe(event => {
    ctx.stage.viewer.mark(EveryLoci, MarkerAction.RemoveHighlight)
    if (event && event.data) {
        ctx.stage.viewer.mark(event.data, MarkerAction.Highlight)
    }
})

ReactDOM.render(React.createElement(Layout, { controller: ctx.layout }), elm);
