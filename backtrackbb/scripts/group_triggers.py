# -*- coding: utf8 -*-
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import os
from ..mod_setup import configure
from ..bp_types import Trigger, Pick
from ..mod_group_trigs import group_triggers


def main():
    config = configure('group_triggers')

    triggers = []
    for line in open(config.options.trigger_file, 'r'):
        try:
            trigger = Trigger()
            trigger.from_str(line)
            triggers.append(trigger)
        except ValueError:
            try:
                pick = Pick()
                pick.from_str(line)
                triggers[-1].add_pick(pick)
            except ValueError:
                continue

    sorted_trigs = group_triggers(config, triggers)

    trigger_file_out = os.path.basename(config.options.trigger_file)
    trigger_file_out = os.path.splitext(trigger_file_out)[0]
    trigger_file_out = trigger_file_out + '.grouped.dat'

    with open(trigger_file_out, 'w') as f:
        for trigger in sorted_trigs:
            f.write(str(trigger) + '\n')
            picks = sorted(trigger.picks, key=lambda x: x.station)
            for pick in picks:
                f.write(str(pick) + '\n')


if __name__ == '__main__':
    main()
