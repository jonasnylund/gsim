import argparse
import os
import pathlib
import time

from gsim.model import py_model  # type: ignore
from viewer import viewer as viewer_lib


def main(args: argparse.Namespace) -> None:
  model = py_model.Model(
    args.dtime,
    args.substeps,
    args.epsilon,
    args.theta,
  )

  viewer = viewer_lib.Viewer()

  model.random_particles(args.num_particles)
  model.initialize()
  start_time = time.time()

  if args.output is not None:
    output_path = os.path.expanduser(args.output)
    if os.path.exists(output_path):
      os.remove(output_path)
  else:
    output_path = None

  while viewer.update():
    print(f'{model.get_time():.1f}', end='\r')
    model.step(1)

    viewer.show(model.get_particles())

    if output_path is not None:
      model.write_particles(os.path.expanduser(output_path))

  model.print_stats()
  print(f'{time.time() - start_time:.3f} s')


if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument(
    '-p',
    '--num_particles',
    type=int,
    default=300,
    help='Number of particles to simulate',
    metavar='NUMBER'
  )
  parser.add_argument(
    '-n',
    '--timesteps',
    type=int,
    default=300,
    help='Number of timesteps to simulate',
    metavar='NUMBER',
  )
  parser.add_argument(
    '--dtime',
    type=float,
    default=0.01,
    help='Time step resolution',
    metavar='NUMBER',
  )
  parser.add_argument(
    '--substeps',
    type=int,
    default=4,
    help='Maximum number of substeps per time step.',
    metavar='INTEGER',
  )
  parser.add_argument(
    '--epsilon',
    type=float,
    default=0.05,
    help='Fudge factor to avoid div/0, giving infinite accelleration',
    metavar='NUMBER',
  )
  parser.add_argument(
    '--theta',
    type=float,
    default=0.3,
    help='Cell width/distance ratio threshold for when to expand the tree during force computations',
    metavar='NUMBER',
  )
  parser.add_argument(
    '-o',
    '--output',
    type=pathlib.Path,
    default=None,
    help='Path to write the output to',
    metavar='PATH',
  )

  main(parser.parse_args())
