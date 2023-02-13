#!/usr/bin/env python3

from typing import Tuple
import functools
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
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

def next_frame(frame: int, file_, tree_, scatter, figure, ax):
  line = file_.readline()
  if line is None or len(line) < 5:
    raise StopIteration

  t, data = parse_line(line)
  figure.suptitle(f"Time: {t:.1f}")

  # scatter = ax.scatter(data[:, 1], data[:, 2], s=np.sqrt(data[:, 0]) * 10)

  scatter.set_sizes(np.sqrt(data[:, 0]) * 10)
  scatter.set_offsets(data[:, 1:])

  if tree_ is not None:
    line = tree_.readline()
    t2, data = parse_line(line)
    assert(t == t2)
    l = len(ax.patches)
    while(len(ax.patches) > 0):
      ax.patches[0].remove()
    print(l, len(ax.patches))
    for r in data:
      ax.add_patch(patches.Rectangle(r[1:] - r[0], width=r[0] * 2, height=r[0] * 2, fill=False))
  return scatter


if __name__ == "__main__":
  if not len(sys.argv) > 1:
    raise ValueError('Must provide path to input file as first argument')
  path = sys.argv[1]
  if len(sys.argv) > 2:
    tree = sys.argv[2]
  else:
    tree = None

  fig, ax = plt.subplots(figsize=(10, 10))
  sc = ax.scatter([], [], s=[])

  
  if tree is not None:
    t = open(tree)
  else:
    t = None
  with open(path) as f:
    update = functools.partial(next_frame, file_=f, tree_=t, scatter=sc, figure=fig, ax=ax)
    init = functools.partial(init_animation, ax=ax, sc=sc)

    animation = FuncAnimation(fig, update, init_func=init, interval=30)
  
    plt.show()
  if tree is not None:
    t.close()