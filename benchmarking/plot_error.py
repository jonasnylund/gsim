#!/usr/env python3
import os
import sys
from typing import BinaryIO, List, Iterator, Tuple

import matplotlib.pyplot as plt
import numpy as np

NDIMS = 3
SIZE_T_BYTES = 8
DATA_DTYPE = np.float64

SingleState = Tuple[float, np.ndarray]


def get_from_file(file_: BinaryIO, ndims: int) -> Iterator[SingleState]:
  while True:
    header = file_.read(SIZE_T_BYTES)
    if header == b'':
      return StopIteration

    bytes_size = np.frombuffer(header, dtype=np.uint64)
    data = file_.read(bytes_size[0])
    data = np.frombuffer(data, DATA_DTYPE)
    t = data[0]
    particles = data[1:].reshape((-1, 1 + ndims))
    yield t, particles


def read_file(path: str) -> List[SingleState]:
  outputs: List[SingleState] = []
  with open(path, "rb") as f:
    for state in get_from_file(f, NDIMS):
      outputs.append(state)
  return outputs


def get_rmserror(base: str, path: str) -> np.ndarray:
  exact = read_file(os.path.join(base, 'exact/particles.bin'))
  measurement = read_file(os.path.join(path, 'particles.bin'))

  output = []
  for (t0, reference), (t1, measured) in zip(exact, measurement):
    rmse = np.sqrt(np.mean((reference - measured) ** 2))
    output.append((t0, rmse))
  return np.array(output)


if __name__ == '__main__':
  print(sys.argv)
  base_path = sys.argv[1]
  baseline = get_rmserror(base_path, sys.argv[2])

  fig, ax = plt.subplots(figsize=(10, 3))
  for arg in sys.argv[3:]:
    data = get_rmserror(base_path, arg)
    ax.plot(
      data[:, 0],
      data[:, 1] - baseline[:, 1],
      label=os.path.basename(arg).capitalize(),
      )
  ax.legend()
  ax.set_xlabel('Timestep')
  ax.set_ylabel('delta RMSE')

  fig.tight_layout()
  # fig.show()
  plt.show()
