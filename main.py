import argparse
import os
import pathlib

from gsim.model import py_model


def main(args: argparse.Namespace) -> None:
  model = py_model.Model(
    args.dtime,
    args.substeps,
    0.01,
    0.3,
  )

  model.random_particles(args.num_particles)
  model.initialize()
  if args.output is not None:
    os.remove(args.output)

  for i in range(args.timesteps):
    print(i, end='\r')
    dt = i + 1 - model.get_time()
    model.step(dt)
    if args.output is not None:
      model.write_particles(os.path.expanduser(args.output))

  model.print_stats()


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
    type=float,
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
    '-o',
    '--output',
    type=pathlib.Path,
    default=None,
    help='Path to write the output to',
    metavar='PATH',
  )

  main(parser.parse_args())