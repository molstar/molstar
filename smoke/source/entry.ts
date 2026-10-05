import { Color } from '@molstar/core/util/color';
import { Task } from '@molstar/core/task';

export const smoke = { color: Color(0x336699), task: Task.create('source smoke', async () => 'ok') };
