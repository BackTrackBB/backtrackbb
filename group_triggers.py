#!/usr/bin/env python
# -*- coding: utf8 -*-
import sys
import os
from tatka_modules.parse_config import parse_config
from tatka_modules.bp_types import Trigger, Pick
from tatka_modules.mod_group_trigs import group_triggers

def main():
    progname = sys.argv[0]
    if len(sys.argv) < 3:
        print 'Usage: %s config_file trigger_file' % progname
        sys.exit(1)
    else:
        config_file = sys.argv[1]
        trigger_file = sys.argv[2]

    config = parse_config(config_file)

    triggers = []
    for line in open(trigger_file, 'r'):
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

    trigger_file_out = os.path.basename(trigger_file)
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
