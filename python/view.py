#!/usr/bin/env python3

from typing import Tuple
import functools
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


def parse_line(line: str) -> Tuple[float, np.ndarray]:
  data = np.fromstring(line, sep=",")
  t = data[0]
  particles = data[1: -1].reshape((-1, 3))
  return t, particles

def init_animation(ax: plt.Axes, sc):
  ax.set_xlim(-200, 200)
  ax.set_ylim(-200, 200)
  return sc

def next_frame(frame: int, file_, scatter, figure):
  line = file_.readline()
  if line is None or len(line) < 5:
    raise StopIteration

  t, data = parse_line(line)
  figure.suptitle(f"Time: {t:.1f}")

  scatter.set_sizes(np.sqrt(data[:, 0]) * 10)
  scatter.set_offsets(data[:, 1:])
  return scatter


if __name__ == "__main__":
  if not len(sys.argv) > 1:
    raise ValueError('Must provide path to input file as first argument')
  path = sys.argv[1]

  fig, ax = plt.subplots(figsize=(10, 10))
  sc = ax.scatter([], [], s=[])
    
  with open(path) as f:
    update = functools.partial(next_frame, file_=f, scatter=sc, figure=fig)
    init = functools.partial(init_animation, ax=ax, sc=sc)

    animation = FuncAnimation(fig, update, init_func=init, interval=30)
  
    plt.show()