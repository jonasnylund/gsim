#!/usr/env python3

import dataclasses
import os
import sys
import re

import numpy as np
import matplotlib.pyplot as plt

from typing import Any, Dict, List, Tuple, Union


def list_folders(path: str) -> List[str]:
  files = os.listdir(path)
  folders = []
  for file in files:
    if os.path.isdir(os.path.join(path, file)):
      folders.append(os.path.join(path, file))
  return folders


@dataclasses.dataclass
class SingleRun:
  parameters: Dict[str, float]
  timers: Dict[str, float]
  stats: Dict[str, int]
  state: Dict[str, int]


@dataclasses.dataclass
class Scenario:
  iterations: List[SingleRun]

  def parameter(self, name: str) -> float:
    return self.iterations[0].parameters[name]

  def stats(self, name: str) -> int:
    return self.iterations[0].stats[name]
  
  def state(self, name: str) -> int:
    return self.iterations[0].state[name]

  def timers(self, name: str) -> Tuple[float, float]:
    values = []
    for run in self.iterations:
      values.append(run.timers[name])
    return np.mean(values), np.std(values)


def parse_stats(path: str) -> Scenario:
  with open(path) as f:
    lines = f.readlines()
  pattern = re.compile(r"^(?P<name>(\s?[A-Za-z:]+)+)\s+(-\s)?(?P<value>\d+(\.\d+)?)")
  
  dictionaries: Dict[str, Dict[str, Union[float, int]]] = dict()
  iterations: List[SingleRun] = []
  current_id = None
  for line in lines:
    if "--- Final state: ---" in line:
      current_id = "state"
      if dictionaries:
        iterations.append(SingleRun(
          dictionaries["parameters"],
          dictionaries["timers"],
          dictionaries["stats"],
          dictionaries["state"],
          ))
      dictionaries.clear()
    elif "Parameters:" in line:
      current_id = "parameters"
    elif "--- Stats for run: ---" in line:
      current_id = "stats"
    elif "--- Timers: ---" in line:
      current_id = "timers"
    elif current_id is not None:
      current = dictionaries.setdefault(current_id, dict())
      match = pattern.search(line)
      if match is None:
        continue

      key = match["name"].replace(":", "").replace(" ", "_").lower()
      current[key] = float(match["value"])

  if dictionaries:
    iterations.append(SingleRun(
      dictionaries["parameters"],
      dictionaries["timers"],
      dictionaries["stats"],
      dictionaries["state"],
      ))
    
  return Scenario(iterations)


if __name__ == "__main__":
  path = sys.argv[1]
  reference = sys.argv[2] if len(sys.argv) > 2 else None
  FOLDERS = list_folders(path)

  scenarios: Dict[str, Scenario] = dict()
  for folder in FOLDERS:
    name = os.path.basename(folder)
    if "exact" in name:
      continue
    scenarios[name] = parse_stats(os.path.join(folder, "stats.txt"))

  plotdata = dict()
  for name, scenario in scenarios.items():
    mean, stdev = scenario.timers("tree_total")
    print(name, round(mean / 1000, 1), round(stdev / 1000, 1))
    plotdata[name] = (mean / 1000, stdev / 1000)
  
  if reference is not None:
    mean_ref, std_ref = plotdata[reference]
    for name, (mean, stdev) in plotdata.items():
      stdev = std_ref * 1 / mean + stdev * mean_ref / (mean ** 2)
      mean = mean_ref / mean

      plotdata[name] = (mean, stdev)

  plt.bar(
    x = np.arange(len(plotdata)),
    height=[x for x, _ in plotdata.values()],
    tick_label=[str(key) for key in plotdata.keys()],
    yerr=[x for _, x in plotdata.values()],
    capsize=5,
    )

  plt.show()
