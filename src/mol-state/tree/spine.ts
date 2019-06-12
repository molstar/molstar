/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { State } from '../state';
import { StateTransform } from '../transform';
import { StateObject, StateObjectCell } from '../object';

export { StateTreeSpine }

/** The tree spine allows access to ancestor of a node during reconciliation. */
interface StateTreeSpine {
    getAncestorOfType<T extends StateObject.Ctor>(type: T): StateObject.From<T> | undefined;
    getRootOfType<T extends StateObject.Ctor>(type: T): StateObject.From<T> | undefined;
}

namespace StateTreeSpine {
    export class Impl implements StateTreeSpine {
        private current: StateObjectCell | undefined = void 0;
        setCurrent(cell?: StateObjectCell) {
            this.current = cell;
        }

        getAncestorOfType<T extends StateObject.Ctor>(t: T): StateObject.From<T> | undefined {
            if (!this.current) return void 0;
            let cell = this.current;
            while (true) {
                cell = this.cells.get(cell.transform.parent)!;
                if (!cell.obj) return void 0;
                if (cell.obj.type === t.type) return cell.obj as StateObject.From<T>;
                if (cell.transform.ref === StateTransform.RootRef) return void 0;
            }
        }

        getRootOfType<T extends StateObject.Ctor>(t: T): StateObject.From<T> | undefined {
            if (!this.current) return void 0;
            let cell = this.current;
            let ret: StateObjectCell | undefined = void 0;
            while (true) {
                cell = this.cells.get(cell.transform.parent)!;
                if (!cell.obj) return void 0;
                if (cell.obj.type === t.type) {
                    ret = cell;
                }
                if (cell.transform.ref === StateTransform.RootRef) return ret ? ret.obj as StateObject.From<T> : void 0;
            }
        }

        constructor(private cells: State.Cells) {

        }
    }
}