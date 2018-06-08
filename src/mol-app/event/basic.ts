/*
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Adapted from LiteMol
 * Copyright (c) 2016 - now David Sehnal, licensed under Apache 2.0, See LICENSE file for more info.
 */

import { Event } from './event'
import { Logger } from '../service/logger';
import { Dispatcher } from '../service/dispatcher'
import { LayoutState } from '../controller/layout';
import { ViewportOptions } from '../controller/visualization/viewport';
import { Job } from '../service/job';
import { Element } from 'mol-model/structure'

const Lane = Dispatcher.Lane;

export const LogEvent = Event.create<Logger.Entry>('bs.Log', Lane.Log);

export namespace CommonEvents {
    export const LayoutChanged = Event.create('bs.Common.LayoutChanged', Lane.Slow);
    export const ComponentsChanged = Event.create('bs.Common.ComponentsChanged', Lane.Slow);
}

export namespace JobEvents {
    export const Started = Event.create<Job.Info>('bs.Jobs.Started', Lane.Job);
    export const Completed = Event.create<number>('bs.Jobs.Completed', Lane.Job);
    export const StateUpdated = Event.create<Job.State>('bs.Jobs.StateUpdated', Lane.Busy);
}

export namespace LayoutEvents {
    export const SetState = Event.create<Partial<LayoutState>>('lm.cmd.Layout.SetState', Lane.Slow);
    export const SetViewportOptions = Event.create<ViewportOptions>('bs.cmd.Layout.SetViewportOptions', Lane.Slow);
}

export namespace InteractivityEvents {
    export const HighlightElementLoci = Event.create<Element.Loci | undefined>('bs.Interactivity.HighlightElementLoci', Lane.Slow);
}
