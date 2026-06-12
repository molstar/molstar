/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { useEffect, useState } from 'react';
import { createRoot } from 'react-dom/client';
import { useCreatePluginUIViewModel } from '../../extensions/plugin/hooks/use-ui-view-model';
import { useCreatePluginViewModel } from '../../extensions/plugin/hooks/use-view-model';
import { loadAlphaFoldDb, loadPdb, loadMvsState } from '../../extensions/plugin/loaders';
import { PluginCanvas } from '../../extensions/plugin/react';
import { Plugin } from '../../mol-plugin-ui/plugin';
import '../../mol-plugin-ui/skin/light.scss';
import './index.html';
import { createViewerSpec } from '../../apps/viewer/plugin-spec';
import { PluginViewModel } from '../../extensions/plugin/view-model';
import { MVSData } from '../../extensions/mvs';
import { MolViewSpecBehavior } from '../../extensions/mvs/behavior';
import { PluginUIViewModel } from '../../extensions/plugin/ui-view-model';
import { BehaviorSubject } from 'rxjs';
import { useBehavior } from '../../mol-plugin-ui/hooks/use-behavior';

function Root() {
    return <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gridTemplateRows: '1fr 1fr', position: 'absolute', inset: 8, gap: 8 }}>
        <div style={{ gridColumn: 1, gridRow: 1, position: 'relative' }}>
            <div style={{ display: 'flex', flexDirection: 'column', height: '100%' }}>
                <div style={{ padding: 4 }}>Only Canvas</div>
                <div style={{ flexGrow: 1, position: 'relative' }}>
                    <OnlyCanvas />
                </div>
            </div>
        </div>
        <div style={{ gridColumn: 2, gridRow: 1, position: 'relative' }}>
            <div style={{ display: 'flex', flexDirection: 'column', height: '100%' }}>
                <div style={{ padding: 4 }}>Default UI + MolViewSpec</div>
                <div style={{ flexGrow: 1, position: 'relative' }}>
                    <DefaultUI />
                </div>
            </div>
        </div>
        <div style={{ gridColumn: 1, gridRow: 2, position: 'relative' }}>
            <div style={{ display: 'flex', flexDirection: 'column', height: '100%' }}>
                <div style={{ padding: 4 }}>Viewer App Spec with bundled extensions (e.g., pLDDT)</div>
                <div style={{ flexGrow: 1, position: 'relative' }}>
                    <ViewerUI />
                </div>
            </div>
        </div>
        <div style={{ gridColumn: 2, gridRow: 2, position: 'relative' }}>
            <div style={{ display: 'flex', flexDirection: 'column', height: '100%' }}>
                <div style={{ padding: 4 }}>Move UI around reusing a single plugin instance</div>
                <div style={{ flexGrow: 1, position: 'relative' }}>
                    <ReuseModel />
                </div>
            </div>
        </div>
    </div>;
}

function OnlyCanvas() {
    const model = useCreatePluginViewModel();
    useEffect(() => {
        loadPdb(model.plugin, '1tqn');
    }, [model]);

    return <PluginCanvas model={model} />;
}

function DefaultUI() {
    const model = useCreatePluginUIViewModel({
        spec: spec => ({
            ...spec,
            layout: {
                initial: {
                    isExpanded: false,
                    showControls: false,
                }
            },
            behaviors: [
                ...spec.behaviors,
                MolViewSpecBehavior
            ]
        })
    });
    useEffect(() => {
        const builder = MVSData.createBuilder();

        const structure = builder
            .download({ url: 'https://www.ebi.ac.uk/pdbe/entry-files/1cbs.bcif' })
            .parse({ format: 'bcif' })
            .modelStructure({});
        structure
            .component({ selector: 'polymer' })
            .representation({ type: 'cartoon' })
            .color({ color: 'green' });
        structure
            .component({ selector: 'ligand' })
            .representation({ type: 'ball_and_stick' })
            .color({ color: '#cc3399' });

        loadMvsState(model.plugin, builder.getState());
    }, [model]);

    return <Plugin plugin={model.plugin} />;
}

// You can also create your own custom view model by extending PluginUIViewModel or PluginViewModel,
// and use it in useCreatePluginUIViewModel. This allows you to add your own state and methods
// to the view model, and use them in your React components.
class CustomViewerModel extends PluginUIViewModel {
    state = {
        // It is often convenient to keep some state related to the viewer
        // close to it. Use useBehavior hook to access it in the React components.
        loadedId: new BehaviorSubject<string | null>(null),
    };

    loadAlphaFoldDb(id: string) {
        // In React strict mode in dev, the useEffect is called twice, so we
        // check if the same ID is already loaded to avoid reloading it.
        // In a real app, might want to use a more robust solution, e.g., an async queue
        // to serialize the commands.
        if (this.state.loadedId.value === id) return;
        this.state.loadedId.next(id);
        return loadAlphaFoldDb(this.plugin, id);
    }
}

function ViewerUI() {
    const model = useCreatePluginUIViewModel({
        spec: () => createViewerSpec({
            layoutIsExpanded: false,
            layoutShowControls: false,
        }),
        model: spec => new CustomViewerModel({ spec })
    });
    useEffect(() => {
        model.loadAlphaFoldDb('AF-Q9I1F6-F1');
    }, [model]);

    return <>
        <Plugin plugin={model.plugin} />
        <AFIdLabel model={model} />
    </>;
}

function AFIdLabel({ model }: { model: CustomViewerModel }) {
    // useBehavior is used to access the RxJS BehaviorSubject in the model's state.
    // It will re-render the component whenever the value changes.
    const id = useBehavior(model.state.loadedId);
    return <div style={{ position: 'absolute', bottom: 8, left: 8, fontSize: 10, padding: '2px 4px', background: 'rgba(255, 255, 255, 0.8)' }}>
        {id ? `Loaded ${id}` : 'Loading...'}
    </div>;
}

// In simple apps, the model can just be a global variable.
// In more complex apps, it can be stored in a context or a state management solution.
const GlobalModel = new PluginViewModel();

// It doesn't really matter where the initial load is, it can live outside the React tree.
loadPdb(GlobalModel.plugin, '1fdl');

function ReuseModel() {
    const model = GlobalModel;
    const [target, setTarget] = useState<'left' | 'right'>('left');

    const update = async (newTarget: 'left' | 'right', newPdb: string) => {
        // in a real app, this should use an async queue to enque the command
        setTarget(newTarget);
        await model.plugin.clear();
        await loadPdb(model.plugin, newPdb);
    };

    return <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: 8, height: '100%' }}>
        <div style={{ gridColumn: 1, gridRow: 1, position: 'relative', display: 'flex', alignItems: 'center', justifyContent: 'center' }}>
            {target === 'left' && <PluginCanvas model={model} />}
            {target !== 'left' && <button style={{ padding: '2px 4px' }} onClick={() => update('left', '1fdl')}>Show on left (1fdl)</button>}
        </div>
        <div style={{ gridColumn: 2, gridRow: 1, position: 'relative', display: 'flex', alignItems: 'center', justifyContent: 'center' }}>
            {target === 'right' && <PluginCanvas model={model} />}
            {target !== 'right' && <button style={{ padding: '2px 4px' }} onClick={() => update('right', '3hfm')}>Show on right (3hfm)</button>}
        </div>
    </div>;
}

createRoot(document.getElementById('app')!).render(<Root />);