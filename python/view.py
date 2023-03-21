#!/usr/bin/env python3

from typing import Tuple

import functools
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.animation import FuncAnimation


SIZE_T_BYTES = 8
DATA_DTYPE = np.float64


def get_from_file(file_, ndims) -> Tuple[float, np.ndarray]:
  header = file_.read(SIZE_T_BYTES)
  if header == b'':
    return StopIteration

  bytes_size = np.frombuffer(header, dtype=np.uint64)
  data = file_.read(bytes_size[0])
  data = np.frombuffer(data, DATA_DTYPE)
  t = data[0]
  particles = data[1:].reshape((-1, 1+ndims))
  return t, particles

def init_animation(ax: plt.Axes, sc):
  ax.set_xlim(-250, 250)
  ax.set_ylim(-250, 250)
  return sc

def next_frame(frame: int, file_, tree_, scatter, figure, ax, ndims):
  data = get_from_file(file_, ndims)
  if data == StopIteration:
    return StopIteration
  t, data = data
  figure.suptitle(f"Time: {t:.1f}")

  scatter.set_sizes(np.sqrt(data[:, 0]) * 10)
  scatter.set_offsets(data[:, 1:3])

  if tree_ is not None:
    t2, data = get_from_file(tree_, ndims)
    assert(t == t2)
    ax.patches.clear()
    for r in data:
      ax.add_patch(patches.Rectangle(r[1:3] - r[0], width=r[0] * 2, height=r[0] * 2, fill=False))
  return scatter


if __name__ == "__main__":
  ndims = 3
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
    t = open(tree, "rb")
  else:
    t = None
  with open(path, "rb") as f:
    update = functools.partial(next_frame, file_=f, tree_=t, scatter=sc, figure=fig, ax=ax, ndims=ndims)
    init = functools.partial(init_animation, ax=ax, sc=sc)

    animation = FuncAnimation(fig, update, init_func=init, interval=30)
  
    plt.show()
  if tree is not None:
    t.close()