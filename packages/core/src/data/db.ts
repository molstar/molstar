/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Database } from './db/database.js';
import { Table } from './db/table.js';
import { Column } from './db/column.js';
import * as ColumnHelpers from './db/column-helpers.js';

export type DatabaseCollection<T extends Database.Schema> = { [name: string]: Database<T> };

export { Database, Table, Column, ColumnHelpers };
