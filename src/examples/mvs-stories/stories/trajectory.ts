/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { MVSData_States } from '../../../extensions/mvs/mvs-data';
import { createMVSBuilder, Root } from '../../../extensions/mvs/tree/mvs/mvs-builder';

const Steps = [
    {
        header: 'Animation Demo',
        key: 'intro',
        description: `### Trajectory
TODO: remove before merging`,
        linger_duration_ms: 2000,
        transition_duration_ms: 500,
        state: (): Root => {
            const builder = createMVSBuilder();

            // builder.download({ url: '/tmp/C2H6.xtc' })
            //     .parse({ format: 'xtc' })
            //     .coordinates({ ref: 'coords' });

            const s = builder.download({ url: '/tmp/26_293_nvt.lammpstrj' }).parse({ format: 'lammpstrj' })
                .modelStructure({ model_index: 2, ref: 'model' });

            s.component({ selector: 'all' })
                .representation({ type: 'ball_and_stick' });

            // const anim = builder.animation();
            // anim.interpolate({
            //     kind: 'scalar',
            //     target_ref: 'model',
            //     property: 'model_index',
            //     start: 0,
            //     end: 7,
            //     duration_ms: 1000,
            // });


            return builder;
        },
    },
];

export function buildStory(): MVSData_States {
    const snapshots = Steps.map((s, i) => {
        const builder = s.state();
        if ((s as any).camera) builder.camera((s as any).camera);

        const description = i > 0 ? `${s.description}\n\n[Go to start](#intro)` : s.description;

        return builder.getSnapshot({
            title: s.header,
            key: s.key,
            description,
            description_format: 'markdown',
            linger_duration_ms: s.linger_duration_ms ?? 500,
            transition_duration_ms: s.transition_duration_ms ?? 1000,
        });
    });

    return {
        kind: 'multiple',
        snapshots,
        metadata: {
            title: 'Animation Showcase',
            version: '1.0',
            timestamp: new Date().toISOString(),
        }
    };
}