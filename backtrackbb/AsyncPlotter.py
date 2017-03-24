# -*- coding: utf8 -*-
# source: https://gist.github.com/astrofrog/1453933
import multiprocessing as mp
import time


class AsyncPlotter():
    def __init__(self, processes=mp.cpu_count()):
        self.manager = mp.Manager()
        self.nc = self.manager.Value('i', 0)
        self.pids = []
        self.processes = processes

    def async_plotter(self, nc, canvas, filename, processes):
        while nc.value >= processes:
            time.sleep(0.1)
        nc.value += 1
        canvas.print_figure(filename)
        nc.value -= 1

    def save(self, canvas, filename):
        p = mp.Process(target=self.async_plotter,
                       args=(self.nc, canvas, filename, self.processes))
        p.start()
        self.pids.append(p)

    def terminate(self):
        for p in self.pids:
            p.terminate()
            p.join()

    def join(self):
        for p in self.pids:
            p.join()
