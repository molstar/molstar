/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react';
import { useEffect, useState } from 'react';
import * as ReactDOM from 'react-dom';
import { AlphaOrbitalsExample } from '.';
import { ParameterControls } from '../../mol-plugin-ui/controls/parameters';
import { PluginContextContainer } from '../../mol-plugin-ui/plugin';

export function mountControls(orbitals: AlphaOrbitalsExample, parent: Element) {
    ReactDOM.render(<PluginContextContainer plugin={orbitals.plugin}>
        <Controls orbitals={orbitals} />
    </PluginContextContainer>, parent);
}

function Controls({ orbitals }: { orbitals: AlphaOrbitalsExample }) {
    const params = useBehavior(orbitals.params);
    const values = useBehavior(orbitals.state);

    return <ParameterControls params={params as any} values={values} onChangeValues={(vs: any) => orbitals.state.next(vs)} />;
}


interface Behavior<T> {
    value: T;
    subscribe(f: (v: T) => void): { unsubscribe(): void };
}

export function useBehavior<T>(s: Behavior<T>): T;
// eslint-disable-next-line
export function useBehavior<T>(s: Behavior<T> | undefined): T | undefined;
// eslint-disable-next-line
export function useBehavior<T>(s: Behavior<T> | undefined): T | undefined {
    const [value, setValue] = useState(s?.value);

    useEffect(() => {
        if (!s) return;
        let fst = true;
        const sub = s.subscribe((v) => {
            if (fst) {
                fst = false;
                if (v !== value) setValue(v);
            } else setValue(v);
        });

        return () => {
            sub.unsubscribe();
        };
        // eslint-disable-next-line
    }, [s]);

    return value;
}